/*
 * splinevolume.cpp
 * Compute the value, derivative and hessian of 3d volumes.
 *
 *  Created on: Dec 21, 2019
 *      Author: apedired
 */

#include <mitsuba/render/volume.h>
#include <mitsuba/core/basisspline.h>
#include <mitsuba/core/mstream.h>
#include <mitsuba/core/fresolver.h>
#include <mitsuba/core/properties.h>
#include <mitsuba/core/mmap.h>

#define POWER_ITERATION_STEPS 5
MTS_NAMESPACE_BEGIN

/*!\plugin{splinevolume}{spline-based volume data source that can give both gradients and hessian at a point}
 * \parameters{
 *     \parameter{filename}{\String}{
 *       Specifies the filename of the volume data file to be loaded
 *     }
 *     \parameter{sendData}{\Boolean}{
 *       When this parameter is set to \code{true}, the implementation will
 *       send all volume data to other network render nodes. Otherwise, they
 *       are expected to have access to an identical volume data file that can be
 *       mapped into memory. \default{\code{false}}
 *     }
 *     \parameter{toWorld}{\Transform}{
 *         Optional linear transformation that should be applied to the data
 *     }
 *     \parameter{min, max}{\Point}{
 *         Optional parameter that can be used to re-scale the data so that
 *         it lies in the bounding box between \code{min} and \code{max}.
 *     }
 * }
 *
 * This class implements access to memory-mapped volume data stored on a
 * 3D grid using a simple binary exchange format.
 * The format uses a little endian encoding and is specified as
 * follows:\vspace{3mm}
 *
 * \begin{center}
 * \begin{tabular}{>{\bfseries}p{2cm}p{11cm}}
 * \toprule
 * Position & Content\\
 * \midrule
 * Bytes 1-3&   ASCII Bytes '\code{V}', '\code{O}', and '\code{L}' \\
 * Byte  4&     File format version number (currently 3)\\
 * Bytes 5-8&   Encoding identifier (32-bit integer). The following
 * choices are available:
 * \begin{enumerate}[1.]
 * \item Dense \code{float32}-based representation
 * \item Dense \code{float16}-based representation (\emph{currently not supported by this implementation})
 * \item Dense \code{uint8}-based representation (The range 0..255 will be mapped to 0..1)
 * \item Dense quantized directions. The directions are stored in spherical
 * coordinates with a total storage cost of 16 bit per entry.
 * \end{enumerate}\\
 * Bytes 9-12 &  Number of cells along the X axis (32 bit integer)\\
 * Bytes 13-16 &  Number of cells along the Y axis (32 bit integer)\\
 * Bytes 17-20 &  Number of cells along the Z axis (32 bit integer)\\
 * Bytes 21-24 &  Number of channels (32 bit integer, supported values: 1 or 3)\\
 * Bytes 25-48 &  Axis-aligned bounding box of the data stored in single
 *                precision (order: xmin, ymin, zmin, xmax, ymax, zmax)\\
 * Bytes 49-*  &  Binary data of the volume stored in the specified encoding.
 *                The data are ordered so that the following C-style indexing
 *                operation makes sense after the file has been mapped into memory:\newline
 *                   \ \ \code{data[((zpos*yres + ypos)*xres + xpos)*channels + chan]}\newline
 *                where \code{(xpos, ypos, zpos, chan)} denotes the lookup location.\\
 *
 * \bottomrule
 * \end{tabular}
 * \end{center}
 *
 *
 */
class SplineDataSource : public VolumeDataSource {
public:
	enum EVolumeType {
		EFloat32 = 1,
		EFloat16 = 2,
		EUInt8 = 3,
		EQuantizedDirections = 4
	};

	SplineDataSource(const Properties &props)
		: VolumeDataSource(props) {
		m_volumeToWorld = props.getTransform("toWorld", Transform());
		m_worldToVolume_Rot = Matrix3x3F(m_volumeToWorld.getInverseMatrix().rot3x3());
		m_worldToVolume_RotT = m_worldToVolume_Rot;
		m_worldToVolume_RotT.transpose();
		if (props.hasProperty("min") && props.hasProperty("max")) {
			/* Optionally allow to use an AABB other than
			   the one specified by the grid file */
			m_dataAABB.min = props.getPoint("min");
			m_dataAABB.max = props.getPoint("max");
		}

		/**
		 * When 'sendData' is set to false, only the filename
		 * is transmitted. A following unserialization of the
		 * stream causes the implementation to then look for
		 * the file (which had better exist if unserialization
		 * occurs on a remote machine).
		 */
		 m_sendData = props.getBoolean("sendData", false);


		 loadFromFile(props.getString("filename"));
	}

	SplineDataSource(Stream *stream, InstanceManager *manager)
			: VolumeDataSource(stream, manager) {
		m_volumeToWorld = Transform(stream);
		m_worldToVolume_Rot = Matrix3x3F(stream);
		m_worldToVolume_RotT = Matrix3x3F(stream);
		m_dataAABB = AABB(stream);
		m_sendData = stream->readBool();
		if (m_sendData) {
			m_volumeType = (EVolumeType) stream->readInt();
			m_res = Vector3i(stream);
			m_channels = stream->readInt();
			m_filename = stream->readString();
			size_t volumeSize = getVolumeSize();
			m_data = new uint8_t[volumeSize];
			stream->read(m_data, volumeSize);
		} else {
			fs::path filename = stream->readString();
			loadFromFile(filename);
		}
		configure();
	}

	virtual ~SplineDataSource() {
		if (!m_mmap)
			delete[] m_data;
	}

	size_t getVolumeSize() const {
		size_t nEntries = (size_t) m_res.x
			* (size_t) m_res.y * (size_t) m_res.z;
		switch (m_volumeType) {
			case EFloat32: return 4 * nEntries * m_channels;
			case EFloat16: return 2 * nEntries * m_channels;
			case EUInt8:   return 1 * nEntries * m_channels;
			case EQuantizedDirections:  return 2 * nEntries;
			default:
				Log(EError, "Unknown volume format!");
				return 0;
		}
	}

	void serialize(Stream *stream, InstanceManager *manager) const {
		VolumeDataSource::serialize(stream, manager);

		m_volumeToWorld.serialize(stream);
		m_worldToVolume_Rot.serialize(stream);
		m_worldToVolume_RotT.serialize(stream);
		m_dataAABB.serialize(stream);
		stream->writeBool(m_sendData);

		if (m_sendData) {
			stream->writeInt(m_volumeType);
			m_res.serialize(stream);
			stream->writeInt(m_channels);
			stream->writeString(m_filename.string());
			stream->write(m_data, getVolumeSize());
		} else {
			stream->writeString(m_filename.string());
		}
	}

	void configure() {
		Vector extents(m_dataAABB.getExtents());
		m_worldToVolume = m_volumeToWorld.inverse();
		m_worldToGrid = Transform::scale(Vector(
				(m_res[0] - 1) / extents[0],
				(m_res[1] - 1) / extents[1],
				(m_res[2] - 1) / extents[2])
			) * Transform::translate(-Vector(m_dataAABB.min)) * m_worldToVolume;
//		m_spline.transform(m_worldToVolume);
		m_stepSize = std::numeric_limits<FLOAT>::infinity();
		for (int i=0; i<3; ++i)
			m_stepSize = std::min(m_stepSize, 0.5f * extents[i] / (FLOAT) (m_res[i]-1));
		m_aabb.reset();
		for (int i=0; i<8; ++i)
			m_aabb.expandBy(m_volumeToWorld(m_dataAABB.getCorner(i)));

		/* Precompute cosine and sine lookup tables */
		for (int i=0; i<255; i++) {
			Float angle = (float) i * ((float) M_PI / 255.0f);
			m_cosPhi[i] = std::cos(2.0f * angle);
			m_sinPhi[i] = std::sin(2.0f * angle);
			m_cosTheta[i] = std::cos(angle);
			m_sinTheta[i] = std::sin(angle);
			m_densityMap[i] = i/255.0f;
		}
		m_cosPhi[255] = m_sinPhi[255] = 0;
		m_cosTheta[255] = m_sinTheta[255] = 0;
		m_densityMap[255] = 1.0f;
	}

	void loadFromFile(const fs::path &filename) {
		m_filename = filename;
		fs::path resolved = Thread::getThread()->getFileResolver()->resolve(filename);
		m_mmap = new MemoryMappedFile(resolved);
		ref<MemoryStream> stream = new MemoryStream(m_mmap->getData(), m_mmap->getSize());
		stream->setByteOrder(Stream::ELittleEndian);

		char header[3];
		stream->read(header, 3);
		if (header[0] != 'V' || header[1] != 'O' || header[2] != 'L')
			Log(EError, "Encountered an invalid volume data file "
				"(incorrect header identifier)");
		uint8_t version;
		stream->read(&version, 1);
		if (version != 3)
			Log(EError, "Encountered an invalid volume data file "
				"(incorrect file version)");
		int type = stream->readInt();

		int xres = stream->readInt(),
			yres = stream->readInt(),
			zres = stream->readInt();
		m_res = Vector3i(xres, yres, zres);
		m_channels = stream->readInt();
		std::string format;

		switch (type) {
			case EFloat32:
				if (m_channels != 1 && m_channels != 3)
					Log(EError, "Encountered an unsupported float32 volume data "
						"file (%i channels, only 1 and 3 are supported)",
						m_channels);
				format = "float32";
				break;
			case EFloat16:
				format = "float16";
				Log(EError, "Error: float16 volumes are not yet supported!");
			case EUInt8:
				format = "uint8";
				if (m_channels != 1 && m_channels != 3)
					Log(EError, "Encountered an unsupported uint8 volume data "
						"file (%i channels, only 1 and 3 are supported)", m_channels);
				break;
			case EQuantizedDirections:
				format = "qdir";
				if (m_channels != 3)
					Log(EError, "Encountered an unsupported quantized direction "
							"volume data file (%i channels, only 3 are supported)",
							m_channels);
				break;
			default:
				Log(EError, "Encountered a volume data file of unknown type (type=%i, channels=%i)!", type, m_channels);
		}

		m_volumeType = (EVolumeType) type;

		if (!m_dataAABB.isValid()) {
			FLOAT xmin = stream->readSingle(),
				  ymin = stream->readSingle(),
				  zmin = stream->readSingle();
			FLOAT xmax = stream->readSingle(),
				  ymax = stream->readSingle(),
				  zmax = stream->readSingle();
			m_dataAABB = AABB(Point(xmin, ymin, zmin), Point(xmax, ymax, zmax));
		}

		Log(EDebug, "Mapped \"%s\" into memory: %ix%ix%i (%i channels, format = %s), %s, %s",
			resolved.filename().string().c_str(), m_res.x, m_res.y, m_res.z, m_channels, format.c_str(),
			memString(m_mmap->getSize()).c_str(), m_dataAABB.toString().c_str());
		m_data = (uint8_t *) (((float *) m_mmap->getData()) + 12);


		FLOAT xmin[] = {m_dataAABB.getCorner(0).x, m_dataAABB.getCorner(0).y, m_dataAABB.getCorner(0).z};
		FLOAT xmax[] = {m_dataAABB.getCorner(7).x, m_dataAABB.getCorner(7).y, m_dataAABB.getCorner(7).z};
		int N[] = {xres, yres, zres};
		m_spline.initialize(xmin, xmax,  N);
		m_interpolatableLimits =  AABB(Point(xmin[0], xmin[1], xmin[2]) + Point( 2.0*m_spline.getStride(0)+Epsilon,  2.0*m_spline.getStride(1)+Epsilon,  2.0*m_spline.getStride(2)+Epsilon),
									   Point(xmax[0], xmax[1], xmax[2]) + Point(-2.0*m_spline.getStride(0)-Epsilon, -2.0*m_spline.getStride(1)-Epsilon, -2.0*m_spline.getStride(2)-Epsilon));
		m_maxSDFError = std::sqrt( m_spline.getStride(0)*m_spline.getStride(0) + m_spline.getStride(1)*m_spline.getStride(1) + m_spline.getStride(2)*m_spline.getStride(2));
		const float *floatData = (float *) m_data;

		//convert to FLOAT Data // FIXME? Slow but one time so fine
		FLOAT *FloatData = new FLOAT[m_res.x*m_res.y*m_res.z];
		for(int i=0; i<m_res.x*m_res.y*m_res.z; i++)
			FloatData[i] = (FLOAT)floatData[i];
		m_spline.build(FloatData);
//		for(int z1 = 0; z1 < m_res.z; z1++){
//			for(int y1=0; y1 < m_res.y; y1++){
//				for(int x1=0; x1 < m_res.x; x1++){
//					Point p(x1 * (xmax[0] - xmin[0])/(m_res.x - 1) + xmin[0],
//							y1 * (xmax[1] - xmin[1])/(m_res.y - 1) + xmin[1],
//							z1 * (xmax[2] - xmin[2])/(m_res.z - 1) + xmin[2]
//							);
//					std::cout << value(p) << ", ";
//				}
//				std::cout << std::endl;
//			}
//			std::cout << std::endl;
//		}
//		m_spline.printcoeff3d();

// Display: To validate
//		const float *floatData = (float *) m_data;
//		for(int z1 = 0; z1 < m_res.z; z1++){
//			for(int y1=0; y1 < m_res.y; y1++){
//				for(int x1=0; x1 < m_res.x; x1++)
//					std::cout << FloatData[(z1*m_res.y + y1)*m_res.x + x1] << ", ";
//				std::cout << std::endl;
//			}
//			std::cout << std::endl;
//		}

		delete[] FloatData;
	}

	inline bool insideVolumeLimits(const PointF &pc) const {
		PointF p = m_worldToVolume(pc);
		return (p.x > m_interpolatableLimits.getCorner(0).x && p.x < m_interpolatableLimits.getCorner(7).x &&
				p.y > m_interpolatableLimits.getCorner(0).y && p.y < m_interpolatableLimits.getCorner(7).y &&
				p.z > m_interpolatableLimits.getCorner(0).z && p.z < m_interpolatableLimits.getCorner(7).z );
	}

	inline Float maxSDFError() const {
		return m_maxSDFError;
	}

	inline FLOAT value(const PointF &pc) const {
		PointF p = m_worldToVolume(pc);
		FLOAT x[] = {p.x, p.y, p.z};
//		if(!insideVolume(p))
//			return FLOAT(0.0);
//		return m_spline.value<0, 0, 0>(x);
		return m_spline.value(x);
	}
	inline VectorF gradient(const PointF &pc) const {
		PointF p = m_worldToVolume(pc);
		FLOAT x[] = {p.x, p.y, p.z};
//		if(!insideVolume(p))
//			return VectorF(0.0);
		return (m_worldToVolume_RotT*m_spline.gradient(x));
	}
	inline Matrix3x3F hessian(const PointF &pc) const {
		PointF p = m_worldToVolume(pc);
		FLOAT x[] = {p.x, p.y, p.z};
//		if(!insideVolume(p))
//			return Matrix3x3F(0.0);
		Matrix3x3F H = m_worldToVolume_RotT*m_spline.hessian(x)*m_worldToVolume_Rot;
		return H;
	}

	inline void valueAndGradient(const PointF &pc, FLOAT &f, VectorF &v) const {
		PointF p = m_worldToVolume(pc);
		FLOAT x[] = {p.x, p.y, p.z};
		m_spline.valueAndGradient(x, f, v);

		v = m_worldToVolume_RotT*v;
	}

	inline void gradientAndHessian(const PointF &pc, VectorF &v, Matrix3x3F &M) const {
		PointF p = m_worldToVolume(pc);
		FLOAT x[] = {p.x, p.y, p.z};
		m_spline.gradientAndHessian(x, v, M);
		v = m_worldToVolume_RotT*v;
		M = m_worldToVolume_RotT*M*m_worldToVolume_Rot;
	}

	inline void valueGradientAndHessian(const PointF &pc, FLOAT &f, VectorF &v, Matrix3x3F &M) const {
		PointF p = m_worldToVolume(pc);
		FLOAT x[] = {p.x, p.y, p.z};
		m_spline.valueGradientAndHessian(x, f, v, M);
		v = m_worldToVolume_RotT*v;
		M = m_worldToVolume_RotT*M*m_worldToVolume_Rot;
	}

	inline bool isAcousticRIF() const{
		return false;
	}

	/**
	 * This is needed since Mitsuba might be
	 * compiled with either single/double precision
	 */
	struct float3 {
		float value[3];

		inline float3() { }

		inline float3(float a, float b, float c) {
			value[0] = a; value[1] = b; value[2] = c;
		}

		inline explicit float3(double a, double b, double c) {
			value[0] = (float) a; value[1] = (float) b; value[2] = (float) c;
		}

		inline float3 operator*(FLOAT v) const {
			return float3((float) (value[0]*v), (float) (value[1]*v), (float) (value[2]*v));
		}

		inline float3 operator+(const float3 &f2) const {
			return float3(value[0]+f2.value[0], value[1]+f2.value[1], value[2]+f2.value[2]);
		}

		inline Spectrum toSpectrum() const {
			Spectrum result;
			result.fromLinearRGB(value[0], value[1], value[2]);
			return result;
		}

		inline Vector toVector() const {
			return Vector(value[0], value[1], value[2]);
		}

		float operator[](int i) const {
			return value[i];
		}

		inline Matrix3x3 tensor() const {
			return Matrix3x3(
				value[0]*value[0], value[0]*value[1], value[0]*value[2],
				value[1]*value[0], value[1]*value[1], value[1]*value[2],
				value[2]*value[0], value[2]*value[1], value[2]*value[2]
			);
		}
	};

	Float lookupFloat(const Point &_p) const {
		const Point p = m_worldToGrid.transformAffine(_p);
		const int x1 = math::floorToInt(p.x),
			  y1 = math::floorToInt(p.y),
			  z1 = math::floorToInt(p.z),
			  x2 = x1+1, y2 = y1+1, z2 = z1+1;

		if (x1 < 0 || y1 < 0 || z1 < 0 || x2 >= m_res.x ||
		    y2 >= m_res.y || z2 >= m_res.z)
			return 0;

		const FLOAT fx = p.x - x1, fy = p.y - y1, fz = p.z - z1,
				_fx = 1.0f - fx, _fy = 1.0f - fy, _fz = 1.0f - fz;

		switch (m_volumeType) {
			case EFloat32: {
				const float *floatData = (float *) m_data;
				const FLOAT
					d000 = floatData[(z1*m_res.y + y1)*m_res.x + x1],
					d001 = floatData[(z1*m_res.y + y1)*m_res.x + x2],
					d010 = floatData[(z1*m_res.y + y2)*m_res.x + x1],
					d011 = floatData[(z1*m_res.y + y2)*m_res.x + x2],
					d100 = floatData[(z2*m_res.y + y1)*m_res.x + x1],
					d101 = floatData[(z2*m_res.y + y1)*m_res.x + x2],
					d110 = floatData[(z2*m_res.y + y2)*m_res.x + x1],
					d111 = floatData[(z2*m_res.y + y2)*m_res.x + x2];

				return ((d000*_fx + d001*fx)*_fy +
						(d010*_fx + d011*fx)*fy)*_fz +
					   ((d100*_fx + d101*fx)*_fy +
						(d110*_fx + d111*fx)*fy)*fz;
			}
			case EUInt8: {
				const FLOAT
					d000 = m_densityMap[m_data[(z1*m_res.y + y1)*m_res.x + x1]],
					d001 = m_densityMap[m_data[(z1*m_res.y + y1)*m_res.x + x2]],
					d010 = m_densityMap[m_data[(z1*m_res.y + y2)*m_res.x + x1]],
					d011 = m_densityMap[m_data[(z1*m_res.y + y2)*m_res.x + x2]],
					d100 = m_densityMap[m_data[(z2*m_res.y + y1)*m_res.x + x1]],
					d101 = m_densityMap[m_data[(z2*m_res.y + y1)*m_res.x + x2]],
					d110 = m_densityMap[m_data[(z2*m_res.y + y2)*m_res.x + x1]],
					d111 = m_densityMap[m_data[(z2*m_res.y + y2)*m_res.x + x2]];

				return ((d000*_fx + d001*fx)*_fy +
						(d010*_fx + d011*fx)*fy)*_fz +
					   ((d100*_fx + d101*fx)*_fy +
						(d110*_fx + d111*fx)*fy)*fz;
			}
			default:
				return 0.0f;
		}
	}

	Spectrum lookupSpectrum(const Point &_p) const {
		const Point p = m_worldToGrid.transformAffine(_p);
		const int x1 = math::floorToInt(p.x),
			  y1 = math::floorToInt(p.y),
			  z1 = math::floorToInt(p.z),
			  x2 = x1+1, y2 = y1+1, z2 = z1+1;

		if (x1 < 0 || y1 < 0 || z1 < 0 || x2 >= m_res.x ||
		    y2 >= m_res.y || z2 >= m_res.z)
			return Spectrum(0.0f);

		const FLOAT fx = p.x - x1, fy = p.y - y1, fz = p.z - z1,
				_fx = 1.0f - fx, _fy = 1.0f - fy, _fz = 1.0f - fz;

		switch (m_volumeType) {
			case EFloat32: {
				const float3 *spectrumData = (float3 *) m_data;
				const float3
					&d000 = spectrumData[(z1*m_res.y + y1)*m_res.x + x1],
					&d001 = spectrumData[(z1*m_res.y + y1)*m_res.x + x2],
					&d010 = spectrumData[(z1*m_res.y + y2)*m_res.x + x1],
					&d011 = spectrumData[(z1*m_res.y + y2)*m_res.x + x2],
					&d100 = spectrumData[(z2*m_res.y + y1)*m_res.x + x1],
					&d101 = spectrumData[(z2*m_res.y + y1)*m_res.x + x2],
					&d110 = spectrumData[(z2*m_res.y + y2)*m_res.x + x1],
					&d111 = spectrumData[(z2*m_res.y + y2)*m_res.x + x2];

				return (((d000*_fx + d001*fx)*_fy +
						 (d010*_fx + d011*fx)*fy)*_fz +
						((d100*_fx + d101*fx)*_fy +
						 (d110*_fx + d111*fx)*fy)*fz).toSpectrum();
				}
			case EUInt8: {
				const float3
					d000 = float3(
						m_densityMap[m_data[3*((z1*m_res.y + y1)*m_res.x + x1)+0]],
						m_densityMap[m_data[3*((z1*m_res.y + y1)*m_res.x + x1)+1]],
						m_densityMap[m_data[3*((z1*m_res.y + y1)*m_res.x + x1)+2]]),
					d001 = float3(
						m_densityMap[m_data[3*((z1*m_res.y + y1)*m_res.x + x2)+0]],
						m_densityMap[m_data[3*((z1*m_res.y + y1)*m_res.x + x2)+1]],
						m_densityMap[m_data[3*((z1*m_res.y + y1)*m_res.x + x2)+2]]),
					d010 = float3(
						m_densityMap[m_data[3*((z1*m_res.y + y2)*m_res.x + x1)+0]],
						m_densityMap[m_data[3*((z1*m_res.y + y2)*m_res.x + x1)+1]],
						m_densityMap[m_data[3*((z1*m_res.y + y2)*m_res.x + x1)+2]]),
					d011 = float3(
						m_densityMap[m_data[3*((z1*m_res.y + y2)*m_res.x + x2)+0]],
						m_densityMap[m_data[3*((z1*m_res.y + y2)*m_res.x + x2)+1]],
						m_densityMap[m_data[3*((z1*m_res.y + y2)*m_res.x + x2)+2]]),
					d100 = float3(
						m_densityMap[m_data[3*((z2*m_res.y + y1)*m_res.x + x1)+0]],
						m_densityMap[m_data[3*((z2*m_res.y + y1)*m_res.x + x1)+1]],
						m_densityMap[m_data[3*((z2*m_res.y + y1)*m_res.x + x1)+2]]),
					d101 = float3(
						m_densityMap[m_data[3*((z2*m_res.y + y1)*m_res.x + x2)+0]],
						m_densityMap[m_data[3*((z2*m_res.y + y1)*m_res.x + x2)+1]],
						m_densityMap[m_data[3*((z2*m_res.y + y1)*m_res.x + x2)+2]]),
					d110 = float3(
						m_densityMap[m_data[3*((z2*m_res.y + y2)*m_res.x + x1)+0]],
						m_densityMap[m_data[3*((z2*m_res.y + y2)*m_res.x + x1)+1]],
						m_densityMap[m_data[3*((z2*m_res.y + y2)*m_res.x + x1)+2]]),
					d111 = float3(
						m_densityMap[m_data[3*((z2*m_res.y + y2)*m_res.x + x2)+0]],
						m_densityMap[m_data[3*((z2*m_res.y + y2)*m_res.x + x2)+1]],
						m_densityMap[m_data[3*((z2*m_res.y + y2)*m_res.x + x2)+2]]);

				return (((d000*_fx + d001*fx)*_fy +
						 (d010*_fx + d011*fx)*fy)*_fz +
						((d100*_fx + d101*fx)*_fy +
						 (d110*_fx + d111*fx)*fy)*fz).toSpectrum();

				}
			default: return Spectrum(0.0f);
		}
	}

	Vector lookupVector(const Point &_p) const {
		const Point p = m_worldToGrid.transformAffine(_p);
		const int x1 = math::floorToInt(p.x),
			  y1 = math::floorToInt(p.y),
			  z1 = math::floorToInt(p.z),
			  x2 = x1+1, y2 = y1+1, z2 = z1+1;

		if (x1 < 0 || y1 < 0 || z1 < 0 || x2 >= m_res.x ||
		    y2 >= m_res.y || z2 >= m_res.z)
			return Vector(0.0f);

		const FLOAT fx = p.x - x1, fy = p.y - y1, fz = p.z - z1;
		Vector value;

		#if defined(VINTERP_NEAREST_NEIGHBOR)
			/* Nearest neighbor */
			switch (m_volumeType) {
				case EFloat32: {
					const float3 *vectorData = (float3 *) m_data;
					value = vectorData[
						(((fz < .5) ? z1 : z2) * m_res.y +
						((fy < .5) ? y1 : y2)) * m_res.x +
						((fx < .5) ? x1 : x2)].toVector();
					}
					break;
				case EQuantizedDirections: {
					value = lookupQuantizedDirection(
						(((fz < .5) ? z1 : z2) * m_res.y +
						((fy < .5) ? y1 : y2)) * m_res.x +
						((fx < .5) ? x1 : x2));
					}
					break;
				default:
					return Vector(0.0f);
			}
		#else
			FLOAT _fx = 1.0f - fx, _fy = 1.0f - fy, _fz = 1.0f - fz;

			Matrix3x3 tensor(0.0f);
			switch (m_volumeType) {
				case EFloat32: {
						const float3 *vectorData = (float3 *) m_data;
						for (int k=0; k<8; ++k) {
							uint32_t index = (((k & 4) ? z2 : z1) * m_res.y +
								((k & 2) ? y2 : y1)) * m_res.x + ((k & 1) ? x2 : x1);
							FLOAT factor = ((k & 1) ? fx : _fx) * ((k & 2) ? fy : _fy)
								* ((k & 4) ? fz : _fz);
							Vector d = vectorData[index].toVector();
							tensor(0, 0) += factor * d.x * d.x;
							tensor(0, 1) += factor * d.x * d.y;
							tensor(0, 2) += factor * d.x * d.z;
							tensor(1, 1) += factor * d.y * d.y;
							tensor(1, 2) += factor * d.y * d.z;
							tensor(2, 2) += factor * d.z * d.z;
						}
					}
					break;
				case EQuantizedDirections: {
						for (int k=0; k<8; ++k) {
							uint32_t index = (((k & 4) ? z2 : z1) * m_res.y +
								((k & 2) ? y2 : y1)) * m_res.x + ((k & 1) ? x2 : x1);
							FLOAT factor = ((k & 1) ? fx : _fx) * ((k & 2) ? fy : _fy)
								* ((k & 4) ? fz : _fz);
							Vector d = lookupQuantizedDirection(index);
							tensor(0, 0) += factor * d.x * d.x;
							tensor(0, 1) += factor * d.x * d.y;
							tensor(0, 2) += factor * d.x * d.z;
							tensor(1, 1) += factor * d.y * d.y;
							tensor(1, 2) += factor * d.y * d.z;
							tensor(2, 2) += factor * d.z * d.z;
						}
					}
					break;
				default:
					return Vector(0.0f);
			}

			tensor(1, 0) = tensor(0, 1);
			tensor(2, 0) = tensor(0, 2);
			tensor(2, 1) = tensor(1, 2);

			if (tensor.isZero())
				return Vector(0.0f);

#if 0
			FLOAT lambda[3];
			eig3_noniter(tensor, lambda);
			value = tensor.col(0);
			FLOAT specularity = 1-lambda[1]/lambda[0];
#else
			/* Square the structure tensor for faster convergence */
			tensor *= tensor;

			const FLOAT invSqrt3 = 0.577350269189626f;
			value = Vector(invSqrt3, invSqrt3, invSqrt3);

			/* Determine the dominant eigenvector using
			   a few power iterations */
			for (int i=0; i<POWER_ITERATION_STEPS-1; ++i)
				value = normalize(tensor * value);
			value = tensor * value;
#endif

		#endif

		if (!value.isZero())
			return normalize(m_volumeToWorld(value));
		else
			return Vector(0.0f);
	}

	bool supportsFloatLookups() const { return m_channels == 1; }
	bool supportsSpectrumLookups() const { return m_channels == 3; }
	bool supportsVectorLookups() const { return m_channels == 3; }
	Float getStepSize() const { return m_stepSize; }
	Vector3i getResolution() const {return m_res;}

	Float getMaximumFloatValue() const {
		return 1.0f;
	}

	std::string toString() const {
		std::ostringstream oss;
		oss << "GridVolume[" << endl
			<< "  res = " << m_res.toString() << "," << endl
			<< "  channels = " << m_channels << "," << endl
			<< "  aabb = " << m_dataAABB.toString() << endl
			<< "]";
		return oss.str();
	}

	MTS_DECLARE_CLASS()
protected:
	FINLINE Vector lookupQuantizedDirection(size_t index) const {
		uint8_t theta = m_data[2*index], phi = m_data[2*index+1];
		return Vector(
			m_cosPhi[phi] * m_sinTheta[theta],
			m_sinPhi[phi] * m_sinTheta[theta],
			m_cosTheta[theta]
		);
	}

protected:
	fs::path m_filename;
	uint8_t *m_data;
	basisspline::Spline<3> m_spline;
	bool m_sendData;
	EVolumeType m_volumeType;
	Vector3i m_res;
	int m_channels;
	Transform m_worldToGrid;
	Transform m_worldToVolume;
	Transform m_volumeToWorld;
	Matrix3x3F m_worldToVolume_Rot;
	Matrix3x3F m_worldToVolume_RotT;
	FLOAT m_stepSize;
	AABB m_dataAABB;
	AABB m_interpolatableLimits;
	ref<MemoryMappedFile> m_mmap;
	FLOAT m_cosTheta[256], m_sinTheta[256];
	FLOAT m_cosPhi[256], m_sinPhi[256];
	FLOAT m_densityMap[256];

	// Hack variable to store the maximum distance between two nodes. Useful for the SDF heuristic.
	Float m_maxSDFError;
};

MTS_IMPLEMENT_CLASS_S(SplineDataSource, false, VolumeDataSource);
MTS_EXPORT_PLUGIN(SplineDataSource, "Spline data source");
MTS_NAMESPACE_END
