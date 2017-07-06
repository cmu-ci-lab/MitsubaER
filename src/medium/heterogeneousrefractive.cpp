/*
    author: Adithya Pediredla
*/

#include <mitsuba/core/basisspline.h>
#include <mitsuba/render/volume.h>
#include <mitsuba/render/scene.h>
#include <mitsuba/render/trimesh.h>
#include <mitsuba/core/util.h>
#include <iomanip>

#include "maxexp.h"

#include "WindingNumber/UT_SolidAngle.h"
#include <igl/signed_distance.h>
#include <igl/read_triangle_mesh.h>
#include <igl/parallel_for.h>
#include <igl/readDMAT.h>
#include <igl/readCSV.h>
#include <igl/writeDMAT.h>
//#include "tbb/tbb.h"

#include <Eigen/Core>
#include <cstdlib>

#include <glog/logging.h>
#include <ceres/ceres.h>

using ceres::CostFunction;
using ceres::SizedCostFunction;
using ceres::Problem;
using ceres::Solver;
using ceres::Solve;

Problem problem;

MTS_NAMESPACE_BEGIN

/*!\plugin{heterogeneousrefractive}{heterogeneous refractive medium}
 * \order{1}
 * \parameters{
 *     \parameter{material}{\String}{
 *         Name of a material preset, see
 *         \tblref{medium-coefficients}. \default{\texttt{skin1}}
 *     }
 *     \parameter{sigmaA, sigmaS}{\Spectrum}{
 *         Absorption and scattering
 *         coefficients of the medium in inverse scene units.
 *         These parameters are mutually exclusive with \code{sigmaT} and \code{albedo}
 *         \default{configured based on \code{material}}
 *     }
 *     \parameter{sigmaT, albedo}{\Spectrum}{
 *         Extinction coefficient in inverse scene units
 *         and a (unitless) single-scattering albedo.
 *         These parameters are mutually exclusive with \code{sigmaA} and \code{sigmaS}
 *         \default{configured based on \code{material}}
 *     }
 *     \parameter{\footnotesize{scale}}{\Float}{
 *         Optional scale factor that will be applied to the \code{sigma*} parameters.
 *         It is provided for convenience when accomodating data based on different units,
 *         or to simply tweak the density of the medium. \default{1}
 *     }
 *     \parameter{\Unnamed}{\Phase}{
 *          A nested phase function that describes the directional
 *          scattering properties of the medium. When none is specified,
 *          the renderer will automatically use an instance of
 *          \pluginref{isotropic}.
 *     }
 * }
 *
 * This class implements Eikonal equations along with a flexible homogeneous participating
 * medium with support for arbitrary phase functions and various
 * medium sampling methods.
 * For non-Eikonal part, it provides two different ways of configuring
 * the medium properties. One possibility is to load a material preset
 * using the \code{material} parameter---see \tblref{medium-coefficients}
 * for details. Alternatively, when specifying parameters by hand, they can either
 * be provided using the scattering and absorption coefficients, or
 * by declaring the extinction coefficient and single scattering
 * albedo (whichever is more convenient). Mixing these parameter
 * initialization methods is not allowed.
 * For Eikonal part, it reads the refractive index field, constructs the spline interpolation scheme and uses it to compute the RIF, dRIF, and hRIF.
 *
 * All scattering parameters (named \code{sigma*}) should
 * be provided in inverse scene units. For instance, when a world-space
 * distance of 1 unit corresponds to a meter, the scattering coefficents should
 * have units of inverse meters. For convenience, the \code{scale}
 * parameter can be used to correct the units. For instance, when the scene is
 * in meters and the coefficients are in inverse millimeters, set
 * \code{scale} to \code{1000}.
 *
 * For Eikonal part the RIF is unitless
 *
 * \renderings{
 *    \rendering{A squishy ball rendered with subsurface scattering and
 *    a dielectric BSDF (courtesy of Chanxi Zheng)}{medium_homogeneous_squishy.jpg}
 * }
 *
 * \begin{xml}[caption=Declaration of a forward scattering medium with high albedo]
 * <medium id="myMedium" type="homogeneous">
 *     <spectrum name="sigmaS" value="1"/>
 *     <spectrum name="sigmaA" value="0.05"/>
 *
 *     <phase type="hg">
 *         <float name="g" value="0.7"/>
 *     </phase>
 * </medium>
 * \end{xml}
 *
 * \textbf{Note}: Rendering media that have a spectrally
 * varying extinction coefficient can be tricky, since all
 * commonly used medium sampling methods suffer from high
 * variance in that case. Here, it may often make more sense to render
 * several monochromatic images separately (using only the coefficients for
 * a single channel) and then merge them back into a RGB image. There
 * is a \code{mtsutil} (\secref{mtsutil}) plugin named \code{joinrgb}
 * that will perform this RGB merging process.
 *
 * \begin{table}[h!]
 *     \centering
 *     \vspace{3mm}
 *     {\footnotesize
 *     \begin{tabular}{>{\ttfamily}p{3.8cm}p{.4cm}>{\ttfamily}p{3.8cm}p{.4cm}>{\ttfamily}p{3.8cm}}
 *         \toprule
 *         \rmfamily \small\textbf{Name} &&
 *         \rmfamily \small\textbf{Name} &&
 *         \rmfamily \small\textbf{Name} \\
 *         \cmidrule{1-1} \cmidrule{3-3} \cmidrule{5-5}
 *         Apple && Chicken1 && Chicken2 \\
 *         Cream && Ketchup  && Potato \\
 *         Skimmilk && Skin1 && Skin2 \\
 *         Spectralon && Wholemilk && \\
 *         \cmidrule{1-1} \cmidrule{3-3} \cmidrule{5-5}
 *         Lowfat Milk              &&  Gatorade                &&    White Grapefruit Juice     \\
 *         Reduced Milk             &&  Chardonnay              &&    Shampoo                    \\
 *         Regular Milk             &&  White Zinfandel         &&    Strawberry Shampoo         \\
 *         Espresso                 &&  Merlot                  &&    \mbox{Head \& Shoulders Shampoo}  \\
 *         Mint Mocha Coffee        &&  Budweiser Beer          &&    Lemon Tea Powder           \\
 *         Lowfat Soy Milk          &&  Coors Light Beer        &&    Orange Juice Powder        \\
 *         Regular Soy Milk         &&  Clorox                  &&    Pink Lemonade Powder       \\
 *         Lowfat Chocolate Milk    &&  Apple Juice             &&    Cappuccino Powder          \\
 *         Regular Chocolate Milk   &&  Cranberry Juice         &&    Salt Powder                \\
 *         Coke                     &&  Grape Juice             &&    Sugar Powder               \\
 *         Pepsi                    &&  Ruby Grapefruit Juice   &&    Suisse Mocha               \\
 *         Sprite                   &&                          &&                               \\
 *         \bottomrule
 *     \end{tabular}}
 *     \caption{\label{tbl:medium-coefficients}This
 *          table lists all supported medium material presets. The
 *          top entries are from Jensen et al. \cite{Jensen2001Practical}, and the
 *          bottom ones are from Narasimhan et al. \cite{Narasimhan2006Acquiring}.
 *          They all use units of $\frac{1}{mm}$, so remember to set
 *          \code{scale} appropriately when your scene is not
 *          in units of millimeters.
 *          These material presets can be used with the plugins
 *          \pluginref{homogeneous},\
 *          \pluginref{dipole}, and \
 *          \pluginref{hk}
 *     }
 * \end{table}
 */

template <typename T> FLOAT sgn(T val) {
    return (FLOAT) ((T(0) < val) - (val < T(0)));
}

class HeterogeneousRefractiveMedium;

class DirectConnectionCostFunction : public SizedCostFunction<3, 3>
{
public:
	DirectConnectionCostFunction(const HeterogeneousRefractiveMedium *m, const PointF &p1, const PointF &p2, const bool &isSensorSample, const Matrix3x3F &dpdv0, const Matrix3x3F &dvdv0){
		m_medium = m;
		m_p1 = p1; m_p2 = p2;
		m_dpdv0 = dpdv0; m_dvdv0 = dvdv0;
		m_isSensorSample = isSensorSample;
	}
	virtual bool Evaluate(double const* const* parameters,
	                      double* residuals,
	                      double** jacobians) const;
private:
	const HeterogeneousRefractiveMedium *m_medium;
	PointF m_p1, m_p2;
	Matrix3x3F m_dpdv0, m_dvdv0;
	bool m_isSensorSample;
};

class HeterogeneousRefractiveMedium : public Medium {
public:
	/**
	 * This class supports the following sampling strategies for choosing
	 * a suitable scattering location when sampling the RTE
	 */
	enum ESamplingStrategy {
		EBalance,  /// Exponential distrib.; pick a random channel each time
		ESingle,   /// Exponential distrib.; pick a specified channel
		EManual,   /// Exponential distrib.; manually specify the falloff
		EMaximum   /// Maximum-of-exponential distribution
	};

	HeterogeneousRefractiveMedium(const Properties &props)
			: Medium(props), m_samplingDensity(0.0f), m_maxExpDist(NULL) {

		google::InitGoogleLogging("mitsuba");
		std::string strategy = props.getString("strategy", "balance");

		/* RIF parameters*/
		m_erstepsize = props.getFloat("stepsize", 1e-3);
		m_tol        = props.getFloat("tol2", 1e-6);
		m_rrweight   = props.getFloat("rrweight", 1e-2);
		m_invrrweight= 1/m_rrweight;

		m_precision  = props.getInteger("boundaryprecision", 3);

        m_ceresoptions.check_gradients = props.getBoolean("cerescheckgradients", false);
        m_ceresoptions.gradient_check_relative_precision = props.getFloat("ceresgradient_check_relative_precision", 1e-3);
        m_ceresoptions.max_num_iterations = props.getInteger("ceresmaxiterations", 20); // default

        // These four options are hardcoded. FIXME: Make these programmable?
        m_ceresoptions.minimizer_type = ceres::LINE_SEARCH;
        m_ceresoptions.line_search_direction_type = ceres::BFGS;
        m_ceresoptions.logging_type = ceres::SILENT;
        m_ceresoptions.minimizer_progress_to_stdout = false;

        m_ceresoptions.function_tolerance = props.getFloat("ceresfunctiontolerance", m_tol);
        m_ceresoptions.gradient_tolerance = props.getFloat("ceresgradienttolerance", 0.0);
        m_ceresoptions.parameter_tolerance = props.getFloat("ceresparametertolerance", 0.0);

        m_makeSensorDirectConnections = props.getBoolean("makesensordirectconnections", false);
        m_aggressiveTracing = props.getBoolean("aggressivetracing", false);
		/**
		 * The goal of the medium sampling weight is to be able to
		 * sample medium intarctions according to
		 *    sigma_s(t) * tau(0 <-> t)
		 * as opposed to
		 *    sigma_t(t) * tau(0 <-> t)
		 * See the separate writeup for more details.
		 */
		m_mediumSamplingWeight = props.getFloat("mediumSamplingWeight", -1);
		if (m_mediumSamplingWeight == -1) {
			for (int i=0; i<SPECTRUM_SAMPLES; ++i) {
				/// Record the highest albedo values across channels
				Float albedo = m_sigmaS[i] / m_sigmaT[i];
				if (albedo > m_mediumSamplingWeight && m_sigmaT[i] != 0)
					m_mediumSamplingWeight = albedo;
			}
			if (m_mediumSamplingWeight > 0) {
				/* The medium scatters some light -> place at least half
				   of the samples in it, otherwise we will render lots
				   of spatially varying noise where one pixel has a
				   medium interaction and the neighbors don't */
				m_mediumSamplingWeight = std::max(m_mediumSamplingWeight,
					(Float) 0.5f);
			}
		}

		if (strategy == "balance") {
			m_strategy = EBalance;
		} else if (strategy == "single") {
			m_strategy = ESingle;

			/* By default, choose the lowest-variance channel
			   (the one with the smallest sigma_t, that is) */
			int channel = 0;
			Float smallest = std::numeric_limits<Float>::infinity();
			for (int i=0; i<SPECTRUM_SAMPLES; ++i) {
				if (m_sigmaT[i] < smallest) {
					smallest = m_sigmaT[i];
					channel = i;
				}
			}

			channel = props.getInteger("channel", channel);
			Assert(channel >= 0 && channel < SPECTRUM_SAMPLES);
			m_samplingDensity = m_sigmaT[channel];

			if (props.getBoolean("monochromatic", false)) {
				/* Optionally turn this into a monochromatic medium
				   based on the chosen color channel. This is useful
				   when the whole scene is rendered once per channel
				   and then recombined to create a color image */
				m_sigmaA = Spectrum(m_sigmaA[channel]);
				m_sigmaS = Spectrum(m_sigmaS[channel]);
				m_sigmaT = m_sigmaA + m_sigmaS;
			}
		} else if (strategy == "maximum") {
			m_strategy = EMaximum;
			std::vector<Float> coeffs(SPECTRUM_SAMPLES);
			for (int i=0; i<SPECTRUM_SAMPLES; ++i)
				coeffs[i] = m_sigmaT[i];
			m_maxExpDist = new MaxExpDist(coeffs);
		} else if (strategy == "manual") {
			m_strategy = EManual;
			m_samplingDensity = props.getFloat("samplingDensity");
		} else {
			Log(EError, "Specified an unknown sampling strategy");
		}
	}

	void buildShape(){
		TriMesh *tempTriMesh = dynamic_cast<TriMesh *>(m_shape.get());
		m_isTriMesh = false;
		m_isSphere = false;

		if(tempTriMesh != NULL){
			m_isTriMesh = true;

			// Note: These triangles are already transformed and hence, we do not need to do any more transformations
			int nVertices  = tempTriMesh->getVertexCount();
			int nTriangles = tempTriMesh->getTriangleCount();
			Point *pos = tempTriMesh->getVertexPositions();
			Triangle *triangles = tempTriMesh->getTriangles();

			V.resize(nVertices, 3);
			F.resize(nTriangles, 3);

			for (int i=0; i < nVertices; i++){
				V(i, 0) = pos[i].x; V(i, 1) = pos[i].y; V(i, 2) = pos[i].z;
			}
			for (int i=0; i < nTriangles; i++){
				F(i, 0) = triangles[i].idx[0]; F(i, 1) = triangles[i].idx[1]; F(i, 2) = triangles[i].idx[2];
			}

			U.resize(V.rows());
			for(int i = 0;i<V.rows();i++)
				for(int j = 0;j<3;j++)
					U[i][j] = V(i,j);

			int order = 2;
			m_solid_angle.init(
				F.rows(),
				F.data(),
				V.rows(),
				&U[0],
				order);

		}else{
//			Log(EError, "The shape is a non triangular mesh!");
		}

	}

	HeterogeneousRefractiveMedium(Stream *stream, InstanceManager *manager)
		: Medium(stream, manager), m_maxExpDist(NULL) {
		m_rif = static_cast<VolumeDataSource *>(manager->getInstance(stream));
		m_SDF = static_cast<VolumeDataSource *>(manager->getInstance(stream));
		m_strategy = (ESamplingStrategy) stream->readInt();
		m_samplingDensity = stream->readFloat();
		m_mediumSamplingWeight = stream->readFloat();

		if (m_strategy == EMaximum) {
			std::vector<Float> coeffs(SPECTRUM_SAMPLES);
			for (int i=0; i<SPECTRUM_SAMPLES; ++i)
				coeffs[i] = m_sigmaT[i];
			m_maxExpDist = new MaxExpDist(coeffs);
		}

		configure();
	}

	virtual ~HeterogeneousRefractiveMedium() {
		if (m_maxExpDist)
			delete m_maxExpDist;
	}

	void configure() {
		Medium::configure();
		if (m_rif.get() == NULL)
			Log(EError, "No RIF specified!");
		m_rifAABB = m_rif->getAABB();
		if (m_SDF.get() == NULL)
			Log(EError, "No SDF specified!");
		m_sdfAABB = m_SDF->getAABB();
		if(!m_rif->isAcousticRIF() && m_rifAABB != m_sdfAABB)
			Log(EError, "The bounding boxes of rif and winding number do not match");

		m_albedo = 0;
		for (int i=0; i<SPECTRUM_SAMPLES; ++i) {
			if (m_sigmaT[i] != 0)
				m_albedo = std::max(m_albedo, m_sigmaS[i]/m_sigmaT[i]);
		}
	}

	void serialize(Stream *stream, InstanceManager *manager) const {
		Medium::serialize(stream, manager);
		manager->serialize(stream, m_rif.get());
		manager->serialize(stream, m_SDF.get());
		stream->writeInt(m_strategy);
		stream->writeFloat(m_samplingDensity);
		stream->writeFloat(m_mediumSamplingWeight);
	}

	Spectrum evalTransmittance(const Ray &ray, Sampler *) const {
		Float negLength = ray.mint - ray.maxt;
		Spectrum transmittance;
		for (int i=0; i<SPECTRUM_SAMPLES; ++i)
			transmittance[i] = m_sigmaT[i] != 0
				? math::fastexp(m_sigmaT[i] * negLength) : (Float) 1.0f;
		return transmittance;
	}

	bool sampleDistance(const Ray &ray, MediumSamplingRecord &mRec,
			Sampler *sampler) const {
		FLOAT rand = sampler->next1D(), sampledDistance;
		FLOAT samplingDensity = m_samplingDensity;


//      Code to test the winding number
//		Float minx = m_wnAABB.getCorner(0).x;
//		Float miny = m_wnAABB.getCorner(0).y;
//		Float minz = m_wnAABB.getCorner(0).z;
//		Float maxx = m_wnAABB.getCorner(7).x;
//		Float maxy = m_wnAABB.getCorner(7).y;
//		Float maxz = m_wnAABB.getCorner(7).z;
//		int xres = m_windingnumber->getResolution().x;
//		int yres = m_windingnumber->getResolution().y;
//		int zres = m_windingnumber->getResolution().z;
//		Point p;
//		for (int k=0; k<100; k++){
//			for (int j=0; j<100; j++){
//				for (int i=0; i<100; i++){
//					p.x = i*(maxx- minx)/(xres-1) + minx;
//					p.y = j*(maxy- miny)/(yres-1) + miny;
//					p.z = k*(maxz- minz)/(zres-1) + minz;
//
//					std::cout << m_windingnumber->value(p) << ", ";
//				}
//				std::cout << std::endl;
//			}
//			std::cout << std::endl;
//		}
//		Log(EError, "Crashing for debug");


		if (rand < m_mediumSamplingWeight) {
			rand /= m_mediumSamplingWeight;
			if (m_strategy != EMaximum) {
				/* Choose the sampling density to be used */
				if (m_strategy == EBalance) {
					int channel = std::min((int) (sampler->next1D()
						* SPECTRUM_SAMPLES), SPECTRUM_SAMPLES-1);
					samplingDensity = m_sigmaT[channel];
				}
				sampledDistance = -math::fastlog(1-rand) / samplingDensity;
			} else {
				sampledDistance = m_maxExpDist->sample(1-rand, mRec.pdfSuccess);
			}
		} else {
			/* Don't generate a medium interaction */
			sampledDistance = std::numeric_limits<FLOAT>::infinity();
		}

		bool success = true;
		FLOAT distSurf = 0;
		FLOAT opticalDistance = 0;

		PointF tempP(ray.o);
		VectorF tempV(ray.d);

//		if(!m_windingnumber->insideVolumeLimits(tempP) || fabs(m_windingnumber->value(tempP)) < 1e-2){ // A patch to handle mitsuba's bug. The initial point can itself be outside medium at times.
		if(!m_SDF->insideVolumeLimits(tempP)){ // A patch to handle mitsuba's bug. The initial point can itself be outside medium at times.
			mRec.transmittance = Spectrum(0.0);
			mRec.pdfSuccessRev = mRec.pdfSuccess = 1.0;
			mRec.pdfFailure = 1.0;
			return false;
		}

		FLOAT refStart = m_rif->value(tempP);
		FLOAT refRatioSq = 1.0/(refStart*refStart);

		if(std::isfinite(sampledDistance)){
			tempV *= refStart;
			if(!m_aggressiveTracing){
				success = trace(tempP, tempV, sampledDistance, distSurf, opticalDistance);
			}
			else{
				Float dist_left = sampledDistance;
				Float dist_traced = 0;
				while(dist_left > Epsilon){
					Float sdf = -m_SDF->value(tempP);
					sdf -= m_SDF->maxSDFError();
					if(sdf < Epsilon)
						break;
					Float traceDist = std::min(sdf, dist_left);
					aggressive_trace(tempP, tempV, traceDist, opticalDistance);
					dist_left -= traceDist;
					dist_traced += traceDist;
//					if(!exactInsideShape(tempP))
//						SLog(EError, "aggreessiveTracing strategy Failed");
				}
				success = trace(tempP, tempV, dist_left, distSurf, opticalDistance);
				distSurf += dist_traced;
			}
		}else{
			tempV *= refStart;
			traceTillBoundary(tempP, tempV, distSurf, opticalDistance);
			success = false;
		}

		Float refEnd = m_rif->value(tempP);
		refRatioSq *= refEnd*refEnd;

//		std::cout << "Initial point:" << ray.o.toString() << std::endl;
//		std::cout << "success:" << success << std::endl;
//		std::cout << "sampledDistance:" << sampledDistance << std::endl;
//		std::cout << "final point after trace:" << tempP.toString() << std::endl;


		if(success){
			mRec.t = sampledDistance + ray.mint;
			mRec.opticalLength = opticalDistance;
			mRec.p = Point(tempP);
			mRec.d = Vector(tempV);
			mRec.refRatioSq = refRatioSq;
			mRec.sigmaA = m_sigmaA;
			mRec.sigmaS = m_sigmaS;
			mRec.time = ray.time; // ? Adi: FIXME: This is unclear
			mRec.medium = this;

			/* Fail if there is no forward progress
			   (e.g. due to roundoff errors) */
			if (mRec.p == ray.o)
				success = false;
		}else{
			sampledDistance = distSurf;
			mRec.t = sampledDistance + ray.mint;
			mRec.opticalLength = opticalDistance;
			mRec.p = Point(tempP);
			mRec.d = Vector(tempV);
			mRec.refRatioSq = refRatioSq;
		}

		switch (m_strategy) {
			case EMaximum:
				mRec.pdfFailure = 1-m_maxExpDist->cdf(sampledDistance);
				break;

			case EBalance:
				mRec.pdfFailure = 0;
				mRec.pdfSuccess = 0;
				for (int i=0; i<SPECTRUM_SAMPLES; ++i) {
					FLOAT tmp = math::fastexp(-m_sigmaT[i] * sampledDistance);
					mRec.pdfFailure += tmp;
					mRec.pdfSuccess += m_sigmaT[i] * tmp;
				}
				mRec.pdfFailure /= SPECTRUM_SAMPLES;
				mRec.pdfSuccess /= SPECTRUM_SAMPLES;
				break;

			case ESingle:
			case EManual:
				mRec.pdfFailure = math::fastexp(-samplingDensity * sampledDistance);
				mRec.pdfSuccess = samplingDensity * mRec.pdfFailure;
				break;

			default:
				Log(EError, "Unknown sampling strategy!");
		}

		mRec.transmittance = (m_sigmaT * (-sampledDistance)).exp();
		mRec.pdfSuccessRev = mRec.pdfSuccess = mRec.pdfSuccess * m_mediumSamplingWeight;
		mRec.pdfFailure = m_mediumSamplingWeight * mRec.pdfFailure + (1-m_mediumSamplingWeight);
		mRec.medium = this;
		if (mRec.transmittance.max() < 1e-20)
			mRec.transmittance = Spectrum(0.0f);

		return success;
	}


	void eval(const Ray &ray, const Point &vsp, const Point &vtp, const bool isSensorSample, MediumSamplingRecord &mRec, Sampler *sampler) const {

		FLOAT weight = 1.0;
		VectorF dirToP2;
		VectorF revDirToP1;
		FLOAT opticalDistToP2 = 0.0;
		FLOAT distToP2 = 0.0;

		//Copy to FLOAT variables
		PointF p1(vsp);
		PointF p2(vtp);
		VectorF v(-ray.d); //The ray is defined from vtp to vsp. So, reverse it.

		Float distance = distToP2;
		if(EXPECT_TAKEN(makeDirectConnections(p1, p2, v, isSensorSample, sampler, weight, dirToP2, revDirToP1, opticalDistToP2, distToP2))){
			distance = distToP2;

			switch (m_strategy) {
				case EManual:
				case ESingle: {
						Float temp = math::fastexp(-m_samplingDensity * distance);
						mRec.pdfSuccess = m_samplingDensity * temp;
						mRec.pdfFailure = temp;
					}
					break;

				case EBalance: {
						mRec.pdfSuccess = 0;
						mRec.pdfFailure = 0;
						for (int i=0; i<SPECTRUM_SAMPLES; ++i) {
							Float temp = math::fastexp(-m_sigmaT[i] * distance);
							mRec.pdfSuccess += m_sigmaT[i] * temp;
							mRec.pdfFailure += temp;
						}
						mRec.pdfSuccess /= SPECTRUM_SAMPLES;
						mRec.pdfFailure /= SPECTRUM_SAMPLES;
					}
					break;

				case EMaximum:
					mRec.pdfSuccess = m_maxExpDist->pdf(distance);
					mRec.pdfFailure = 1-m_maxExpDist->cdf(distance);
					break;

				default:
					Log(EError, "Unknown sampling strategy!");
			}

			mRec.transmittance = (m_sigmaT * (-distance)).exp() * weight; // Note that transmittance increased by weight (weight is mostly close to 1 unless russian roulette kicks in)
			mRec.pdfSuccess = mRec.pdfSuccessRev = mRec.pdfSuccess * m_mediumSamplingWeight;
			mRec.pdfFailure = mRec.pdfFailure * m_mediumSamplingWeight + (1-m_mediumSamplingWeight);
		}else{
			// Trick to make failed case return 0.0 spectrum. Make distance can be infinite?
			distance = std::numeric_limits<Float>::infinity();
			mRec.transmittance = Spectrum(0.0f);
			mRec.pdfSuccess = mRec.pdfSuccessRev = 1.0;
			mRec.pdfFailure = 1.0;
		}

		mRec.sigmaA = m_sigmaA;
		mRec.sigmaS = m_sigmaS;
		mRec.time = ray.time;
		mRec.medium = this;
		mRec.distance = distance;
		mRec.opticalLength = opticalDistToP2;
		mRec.drev = Vector(revDirToP1);
		if (mRec.transmittance.max() < 1e-20)
			mRec.transmittance = Spectrum(0.0f);

	}


//	inline VectorF dV(const PointF &p) const{
//		return m_rif->gradient(p);
//	}
//
//	inline void er_step_old(PointF &p, VectorF &v, const FLOAT &stepsize) const{
//		v += VHALF * stepsize * dV(p);
//		p +=         stepsize * v/m_rif->value(p);
//		v += VHALF * stepsize * dV(p);
//	}

	inline void er_step(PointF &p, VectorF &v, const FLOAT &stepsize, FLOAT &opticalDistance) const{
		FLOAT n;
		VectorF G;
		m_rif->valueAndGradient(p, n, G);
		v += VHALF * stepsize * G;
		p +=         stepsize * v/n;
		v += VHALF * stepsize * m_rif->gradient(p);
		opticalDistance += stepsize*n;
	}
	inline void er_step(PointF &p, VectorF &v, const FLOAT &stepsize) const{
		FLOAT n;
		VectorF G;
		m_rif->valueAndGradient(p, n, G);
		v += VHALF * stepsize * G;
		p +=         stepsize * v/n;
		v += VHALF * stepsize * m_rif->gradient(p);
	}

	inline bool trace(PointF &p, VectorF &v, const FLOAT &sampledDistance, FLOAT &distSurf, FLOAT &opticalDistance) const{
	    FLOAT distance = sampledDistance;
	    distSurf = 0;
	    int steps = distance/m_erstepsize;
	    distance  = distance - steps * m_erstepsize;
	    for(int i = 0; i < steps; i++){
	        er_step(p, v, m_erstepsize, opticalDistance);
	        if(!insideShape(p)){
	        	er_step(p, v, -m_erstepsize, opticalDistance);
	        	return false;
	        }
	        distSurf += m_erstepsize;
	    }
	    er_step(p, v, distance, opticalDistance);
        if(!insideShape(p)){
        	er_step(p, v, -distance, opticalDistance);
        	return false;
        }
        distSurf += distance;
        return true;
	}




	//same as trace but does not do inside outside tests and simple goes forward
	inline void aggressive_trace(PointF &p, VectorF &v, const FLOAT &sampledDistance, FLOAT &opticalDistance) const{
	    FLOAT distance = sampledDistance;
	    int steps = distance/m_erstepsize;
	    distance  = distance - steps * m_erstepsize;
	    for(int i = 0; i < steps; i++)
	        er_step(p, v, m_erstepsize, opticalDistance);
	    er_step(p, v, distance, opticalDistance);
	}


	inline bool insideShape(const PointF &p) const{
		return hackForSphere(p);
//		return hackForBox(p);
//		return exactInsideShape(p);
//		return exactInsideShape(p) && (p.x > -.010001); // To take care of inside medium case for US
	}

	inline bool hackForSphere(const PointF &p) const{
//		PointF center(-200, 0, 350);
//		FLOAT radius = 200;
		PointF center(-0.22827, 1.2, 0.152505);
		FLOAT radius = 0.3;
		return ((p - center).lengthSquared() < radius*radius);
	}

	inline bool hackForBox(const PointF &p) const{ // Adi: The obj file is failing to produce accurate WN (for Acquarium Box result). Probably the faces are inverted. FIXME later
		return p.x >= -200 && p.x <= 200 &&
				   p.y >= -200 && p.y <= 200 &&
				   p.z >=   50 && p.z <= 100;
	}

	inline bool approximateInsideShape(const PointF &p) const{
		return std::abs(m_SDF->value(p)) < 0;
	}

	inline bool exactInsideShape(const PointF &p) const{
		double accuracy_scale = 2.0;
		float Pp[3];
		Pp[0] = (float) p.x;
		Pp[1] = (float) p.y;
		Pp[2] = (float) p.z;
		return m_solid_angle.computeSolidAngle(Pp,accuracy_scale)/(4.0*M_PI) > 0.5;
	}

	//FIXME: Make this into bisection search
	inline void traceTillBoundary(PointF &p, VectorF &v, FLOAT &distSurf, FLOAT &opticalDistance) const{

	    distSurf = 0;

	    long int maxsteps = 1e5, i, precision = m_precision; // 1e5 is a large number of steps. Hopefully we will meet the surface before that
	    FLOAT current_stepsize = m_erstepsize;

	    PointF oldp;
	    VectorF oldv;

	    for(i = 0; i < maxsteps; i++){
	    	oldp = p;
	    	oldv = v;
	    	er_step(p, v, current_stepsize, opticalDistance);
	    	if(insideShape(p)){
	    		distSurf += current_stepsize;
	    	}
	    	else{
	    		er_step(p, v, -current_stepsize, opticalDistance);
	    		distSurf -= current_stepsize;
	    		return;
	    		//turns out, this much precision is not necessary with WN interpolation, infact we have to take a step back.
//	    		precision--;
//	    		if(precision < 0)
//	    			break;
//	    		p = oldp;
//	    		v = oldv;
//	    		current_stepsize = current_stepsize / 10;
//	            i  = 0;
//	            maxsteps = 11;
	    	}
	    }
	    if( i == 1e5)
	    	Log(EWarn, "Trace till boundary: Max steps 1e5 are not enough");
	}

//	inline Matrix3x3F d2V(PointF &p, Matrix3x3F &dpdv0) const{
//		return m_rif->hessian(p) * dpdv0;
//	}
//
//	inline Matrix3x3F d2Path(PointF &p, VectorF &v, Matrix3x3F &dpdv0, Matrix3x3F &dvdv0) const{
//		FLOAT invn = 1/m_rif->value(p);
//		return (-invn*invn*Matrix3x3F(v, m_rif->gradient(p))*dpdv0 + invn * dvdv0);
//	}

//	inline void er_derivativestep_old(PointF &p, VectorF &v, Matrix3x3F &dpdv0, Matrix3x3F &dvdv0, const FLOAT &stepsize) const{
//		v     += VHALF * stepsize * dV(p);
//		dvdv0 += VHALF * stepsize * d2V(p, dpdv0);
//		p     +=         stepsize * v/m_rif->value(p);
//		dpdv0 +=         stepsize * d2Path(p, v, dpdv0, dvdv0);
//		v     += VHALF * stepsize * dV(p);
//		dvdv0 += VHALF * stepsize * d2V(p, dpdv0);
//
//	}

	//For optimization. Unwraps all internal code and calls optimized spline interpolants.
	inline void er_derivativestep(PointF &p, VectorF &v, Matrix3x3F &dpdv0, Matrix3x3F &dvdv0, const FLOAT &stepsize) const{
		FLOAT n, invn;
		VectorF G;
		Matrix3x3F H;
		m_rif->valueGradientAndHessian(p, n, G, H);

		v     += VHALF * stepsize * G;
		dvdv0 += VHALF * stepsize * H*dpdv0;
		p     +=         stepsize * v/n;

		m_rif->valueGradientAndHessian(p, n, G, H);
		invn = 1/n;

		dpdv0 +=         stepsize * (-invn*invn*Matrix3x3F(v, G)*dpdv0 + invn * dvdv0);
		v     += VHALF * stepsize * G;
		dvdv0 += VHALF * stepsize * H*dpdv0;
	}

	inline void computefdfBDPT(const VectorF &v_i, const PointF &p1, const PointF &p2, const bool &isSensorSample, Matrix3x3F &dpdv0, Matrix3x3F &dvdv0, VectorF &error, Matrix3x3F &derror) const {

//		if(!m_windingnumber->insideVolumeLimits(p1) || fabs(m_windingnumber->value(p1)) < 1e-2){ // A patch to handle mitsuba's bug. The initial point can itself be outside medium at times.
		if(!m_SDF->insideVolumeLimits(p1)){ // A patch to handle mitsuba's bug. The initial point can itself be outside medium at times.
			error = p1 - p2;
			derror = Matrix3x3F(0.0);
			return;
		}

		//Adithya: Can we detect surely failing case?, if v_i is away from p2? <p2 - p1, v_i> < 0?
		FLOAT currentStepSize = m_erstepsize;
//		int maxSteps = 5*std::floor(dot(pdiff, pdiff)/dot(pdiff, v_i))/currentStepSize; // ADI: Heuristic
		int maxSteps = 1e5; // ADI: Heuristic
		int dec_precision = m_precision;
		long int nBisectionSearches = ceil(dec_precision/log10(2)); // ADI: Make this efficient
		bool leftObject = false;

		PointF p = p1, oldp;
		VectorF v = v_i, oldv;
		Matrix3x3F olddpdv0, olddvdv0;
		bool signOld = std::signbit(dot(p-p2, v)), signNew;

		FLOAT r = m_rif->value(p);
		FLOAT n1 = v_i.length();
		FLOAT n2 = n1*n1;
		FLOAT n3 = n2*n1;
		dvdv0 = (r/n3)*(n2*Matrix3x3F(1.0, 0, 0, 0, 1.0, 0, 0, 0, 1.0) - Matrix3x3F(v, v))*dvdv0;
		v = v/n1 * r;

		for(int i=0; i< maxSteps; i++){
			oldp = p;
			oldv = v;
			olddpdv0 = dpdv0;
			olddvdv0 = dvdv0;
			er_derivativestep(p, v, dpdv0, dvdv0, currentStepSize);
			signNew = std::signbit(dot(p-p2, v));
			// ADI: I believe the if-else-if conditions can be combined and that would be more accurate as well. Test this on MATLAB first!!
			if(signNew != signOld){
	            while(nBisectionSearches > 0){
	                nBisectionSearches--;
	                // load backup
	                p     = oldp    ;
	                v     = oldv    ;
	                dpdv0 = olddpdv0;
	                dvdv0 = olddvdv0;
	                currentStepSize = currentStepSize/2;
	                er_derivativestep(p, v, dpdv0, dvdv0, currentStepSize);
	                signNew = std::signbit(dot(p-p2, v));
	                if(signNew == signOld){
	                	// create better backup
	                    oldp     = p;
	                    oldv     = v;
	                    olddpdv0 = dpdv0;
	                    olddvdv0 = dvdv0;
	                }
	            }
	            break;
			}else if(!insideShape(p)){
	            while(nBisectionSearches > 0){
	                nBisectionSearches--;
	                // load backup
	                p     = oldp    ;
	                v     = oldv    ;
	                dpdv0 = olddpdv0;
	                dvdv0 = olddvdv0;
	                currentStepSize = currentStepSize/2;
	                er_derivativestep(p, v, dpdv0, dvdv0, currentStepSize);
	                if(insideShape(p)){
	                	// create better backup
	                    oldp     = p;
	                    oldv     = v;
	                    olddpdv0 = dpdv0;
	                    olddvdv0 = dvdv0;
	                }
	            }
	            if( (p-p1).lengthSquared() < Epsilon){ //No progress made
	    			error = p1 - p2;
	    			derror = Matrix3x3F(0.0);
	    			return;
	    		}

	            FLOAT nb;
	            VectorF dnb;
	            m_rif->valueAndGradient(p, nb, dnb);
	            VectorF dpdtb = v/nb;
				VectorF N = m_SDF->gradient(p);
				N = normalize(N);
				VectorF dtbdv0 = -dpdv0.preMult(N)/dot(N, dpdtb);

				//boundary hack conditions
				boundaryVelocityDerivative(v, dvdv0, dtbdv0, dnb, N, nb, (FLOAT)1.0);
//				dvdv0 += Matrix3x3F(dnb, dtbdv0);

				FLOAT extra_t = -dot(v, p - p2)/v.lengthSquared();
	            leftObject = true;
	            if(isSensorSample && extra_t < 0){ // for other case, we can go to otherside of extra_t as well,
	            								   // as it improves convergence and intermediate paths can afford to be not super accurate anyway.
	            								   // There is a possibility that it will cause inaccurate paths but they are eliminated in the end
	            	error = p1 - p2; // dError is initialized to zero anyway.
	            	return;
	            }
	            dpdv0 += Matrix3x3F(dpdtb - v, dtbdv0) + extra_t * dvdv0;
				p     += extra_t * v;
				break;
			}else{
				// do nothing
			}
		}
		VectorF dtstardV0;
		VectorF dpdt;
		if(!leftObject){
			VectorF dvdt;
			m_rif->valueAndGradient(p, r, dvdt);
			dpdt = v/r;
			dtstardV0 = -(dpdv0.preMult(v) + dvdv0.preMult(p - p2))/(dot(v, dpdt) + dot(p - p2, dvdt));
		}else{
			dpdt = v;
			dtstardV0 = -(dpdv0.preMult(v) + dvdv0.preMult(p - p2))/dot(v, dpdt);
		}
		Matrix3x3F dQdV = dpdv0 + Matrix3x3F(dpdt, dtstardV0);
		error = p - p2;
		dQdV.transpose();
		derror = dQdV ;
	}

    inline bool computePathLengthsTillClosestP2(const PointF &p1, const PointF &p2, const VectorF &dirToP2, VectorF &revDirToP1, const bool &isSensorSample, FLOAT &opticalDistToP2, FLOAT &distToP2) const {
    	//Adithya: Can we detect surely failing case?, if v_i is away from p2? <p2 - p1, v_i> < 0?
    	distToP2 = 0;
    	opticalDistToP2 = 0;
    	FLOAT currentStepSize = m_erstepsize;
    	int maxSteps = 1e5;

		int dec_precision = m_precision;
		long int nBisectionSearches = ceil(dec_precision/log10(2)); // ADI: Make this efficient

		PointF p = p1, oldp;
		VectorF v = dirToP2, oldv;
		bool signOld = std::signbit(dot(p-p2, v)), signNew;

		for(int i=0; i< maxSteps; i++){
			oldp = p;
			oldv = v;
			er_step(p, v, currentStepSize);
			signNew = std::signbit(dot(p-p2, v));
			if(!insideShape(p)){
				if(!isSensorSample)
					return false;
				else{
					//go close to the boundary
		            while(nBisectionSearches > 0){
		                nBisectionSearches--;
		                // load backup
		                p     = oldp    ;
		                v     = oldv    ;
		                currentStepSize = currentStepSize/2;
		                er_step(p, v, currentStepSize);
		                if(insideShape(p)){
		                    distToP2 += currentStepSize;
		                    opticalDistToP2  += currentStepSize * m_rif->value( VHALF*(p + oldp));
		                	// create better backup
		                    oldp     = p;
		                    oldv     = v;
		                }
		            }
					VectorF N = m_SDF->gradient(p);
					N = normalize(N);

					//Snell's law
					boundaryVelocity(v, N, m_rif->value(p), (FLOAT)1.0);

					FLOAT extra_t = -dot(v, p - p2)/v.lengthSquared();
		            if(extra_t < 0){
		            	return false;
		            }
					p     += extra_t * v;
					opticalDistToP2 += extra_t;// distToP2 does not increment by extra_t as we care only about the geometric distance
				}
				break;
			}

//			if(signNew != signOld || (p-p2).lengthSquared() < m_tol){
			if(signNew != signOld){
	            while(nBisectionSearches > 0){
	                nBisectionSearches--;
	                // load backup
	                p     = oldp    ;
	                v     = oldv    ;
	                currentStepSize = currentStepSize/2;
	                er_step(p, v, currentStepSize);
	                signNew = std::signbit(dot(p-p2, v));
	                if(signNew == signOld){
	                    distToP2 += currentStepSize;
	                    opticalDistToP2  += currentStepSize * m_rif->value( VHALF*(p + oldp));
	                	// create better backup
	                    oldp     = p;
	                    oldv     = v;
	                }
	            }
	            break;
			}else{
                distToP2 += currentStepSize;
                opticalDistToP2  += currentStepSize * m_rif->value( VHALF*(p + oldp));
			}
		}


//		if((p-p2).lengthSquared() > m_tol*100){ // Relax tolerance by (10x)^2 as renormalization can cause minor errors
		if((p-p2).lengthSquared() > m_tol){
			//FAILED CASE
			SLog(EDebug, "Failing after claiming convergence: p: %s; p2: %s; m_tol: %f", p.toString().c_str(), p2.toString().c_str(), m_tol);
			return false;
		}
		revDirToP1 = -normalize(v);
		return true;
    }

    //Hack code to go straight
//    inline void boundaryVelocity(VectorF &v, const VectorF &N, const FLOAT ni, const FLOAT ne) const{
//    }

    inline void boundaryVelocity(VectorF &v, const VectorF &N, const FLOAT ni, const FLOAT ne) const{
    	FLOAT dotp = dot(v, N);
    	FLOAT r = ne/ni;
        r = r*r - 1;
        FLOAT n2 = v.lengthSquared();

        FLOAT sq = (r*n2 + dotp*dotp);
        if(sq < Epsilon){
			// Reflection case
			v = 2*dotp*N - v;
			return;
        }
		sq = std::sqrt(sq);

        v = v - dotp * N + sgn(dotp)*sq*N;
    }

//    inline void boundaryVelocityDerivative(VectorF &v, Matrix3x3F &dvdv0, const VectorF &dtbdv0, const VectorF &dnb, const VectorF &N, const FLOAT ni, const FLOAT ne) const{
//        dvdv0 = (dvdv0 + Matrix3x3F(dnb, dtbdv0));
//    }

    inline void boundaryVelocityDerivative(VectorF &v, Matrix3x3F &dvdv0, const VectorF &dtbdv0, const VectorF &dnb, const VectorF &N, const FLOAT ni, const FLOAT ne) const{
    	FLOAT dotp = dot(v, N);
    	FLOAT r = ne/ni;
        r = r*r - 1;
        FLOAT n2 = v.lengthSquared();

        FLOAT sq = (r*n2 + dotp*dotp);
        if(sq < Epsilon){
			// Reflection case
			v = 2*dotp*N - v;
			dvdv0 = ((FLOAT)2.0*Matrix3x3F(N, N) - Matrix3x3F(1.0, 0, 0, 0, 1.0, 0, 0, 0, 1.0)) * (dvdv0 + Matrix3x3F(dnb, dtbdv0));
			return;
        }
		sq = std::sqrt(sq);

        dvdv0 = (Matrix3x3F(1.0, 0, 0, 0, 1.0, 0, 0, 0, 1.0) - Matrix3x3F(N, N) + sgn(dotp)*Matrix3x3F(N, (r*v + dotp*N)/sq)) * (dvdv0 + Matrix3x3F(dnb, dtbdv0));
        v = v - dotp * N + sgn(dotp)*sq*N;
    }



	inline void uniformSample(const VectorF &in, VectorF &out, Sampler *sampler) const {
		VectorF axisX, axisY;
		coordinateSystem(in, axisX, axisY);
		Vector temp = warp::squareToUniformHemisphere(Point2(sampler->nextFloat(), sampler->nextFloat()));
		out.x = temp.x; out.y = temp.y; out.z = temp.z;
		out = out.x * axisX + out.y * axisY + out.z * in;
	}


	inline bool makeDirectConnections(const PointF &p1, const PointF &p2, const VectorF &d, const bool &isSensorSample, Sampler *sampler, FLOAT &weight, VectorF &dirToP2, VectorF &revDirToP1, FLOAT &opticalDistToP2, FLOAT &distToP2) const {
	    Matrix3x3F dpdv0((FLOAT)0);
	    Matrix3x3F dvdv0((FLOAT)1, 0, 0,
	                    0, 1, 0,
	                    0, 0, 1);

	    VectorF v;
	    VectorF tempSol;
	    int iterations = 1;

//	    if(!m_windingnumber->insideVolumeLimits(p1) || fabs(m_windingnumber->value(p1)) < 1e-2){ // A patch to handle mitsuba's bug. The initial point can itself be outside medium at times.
		if(!m_SDF->insideVolumeLimits(p1)){ // A patch to handle mitsuba's bug. The initial point can itself be outside medium at times.
			return false;
		}

	    FLOAT RIFp = m_rif->value(p1);
	    while(true){
	    	uniformSample(d, v, sampler);
	    	// Hack to match BDPT homogeneous case for validation
//	    	v = p2-p1;
//	    	v = v/v.length();

	    	v *= RIFp;

	    	CostFunction* cost_function = new DirectConnectionCostFunction(this, p1, p2, isSensorSample, dpdv0, dvdv0);
	    	Problem problem;

	    	double x[] = {v.x, v.y, v.z};

	    	problem.AddResidualBlock(cost_function, NULL, x);

	        Solver::Summary summary;
	        Solve(m_ceresoptions, &problem, &summary);

	        if(summary.final_cost < m_tol){
	        	// gradient descent converged. Perform the geometric measurement
	        	if(iterations == 1){
	        		iterations++;
	        		tempSol[0] = x[0];
	        		tempSol[1] = x[1];
	        		tempSol[2] = x[2];
	        		tempSol = normalize(tempSol);
	        	}else
	        		iterations++;
					dirToP2[0] = x[0];
					dirToP2[1] = x[1];
					dirToP2[2] = x[2];
					dirToP2 = normalize(dirToP2);
					if( (tempSol-dirToP2).lengthSquared() < (2*m_tol)){ //using same tolerance here too
						// converged in the sense of Zeltner's paper (or Booth [2007])
						break;
					}
//	            dirToP2[0] = x[0];
//	            dirToP2[1] = x[1];
//	            dirToP2[2] = x[2];
//	            dirToP2 = normalize(dirToP2);
//	            break;
	        }

	        // gradient descent failed, so perform russian roulette
	        if(sampler->nextFloat() < m_rrweight)
	            weight = weight * m_invrrweight;
	        else{
	            dirToP2[0] = x[0];
	            dirToP2[1] = x[1];
	            dirToP2[2] = x[2];
	            dirToP2 = normalize(dirToP2);
	            return false;
	        }
	    }

	    dirToP2 *= RIFp;

	    weight *= (iterations - 1);

	    return computePathLengthsTillClosestP2(p1, p2, dirToP2, revDirToP1, isSensorSample, opticalDistToP2, distToP2); // Can still fail as the ray may go outside the medium
	}

	bool isHomogeneous() const {
		return false;
	}

	bool isheterogeneousrefractive() const {
		return true;
	}

	bool makeSensorDirectConnections() const {return m_makeSensorDirectConnections;};

	Float getRIF(const Point &p) const { return m_rif->value(PointF(p));}

	void addChild(const std::string &name, ConfigurableObject *child) {
		if (child->getClass()->derivesFrom(MTS_CLASS(VolumeDataSource))) {
			VolumeDataSource *volume = static_cast<VolumeDataSource *>(child);

			if (name == "rif") {
				Assert(volume->supportsFloatLookups());
				m_rif = volume;
			} else if (name == "sdf") {
				Assert(volume->supportsFloatLookups());
				m_SDF = volume;
			} else {
				Medium::addChild(name, child);
			}
		} else {
			Medium::addChild(name, child);
		}
	}


	std::string toString() const {
		std::ostringstream oss;
		oss << "HeterogeneousRefractiveMedium[" << endl
			<< "  sigmaA = " << m_sigmaA.toString() << "," << endl
			<< "  sigmaS = " << m_sigmaS.toString() << "," << endl
			<< "  sigmaT = " << m_sigmaT.toString() << "," << endl
			<< "  mediumSamplingWeight = " << m_mediumSamplingWeight << "," << endl
			<< "  samplingDensity = " << m_samplingDensity << "," << endl
			<< "  strategy = ";

		switch (m_strategy) {
			case ESingle: oss << "single," << endl; break;
			case EManual: oss << "manual," << endl; break;
			case EBalance: oss << "balance," << endl; break;
			case EMaximum: oss << "maximum," << endl; break;
		}

		oss << "  phase = " << indent(m_phaseFunction.toString()) << endl
			<< "]";
		return oss.str();
	}

	MTS_DECLARE_CLASS()
private:
	Float m_samplingDensity, m_mediumSamplingWeight;
	ESamplingStrategy m_strategy;
	MaxExpDist *m_maxExpDist;
	Float m_albedo;
	bool m_makeSensorDirectConnections;
	bool m_aggressiveTracing;
protected:
	ref<VolumeDataSource> m_rif;
	ref<VolumeDataSource> m_SDF;
	FLOAT m_erstepsize, m_tol, m_rrweight, m_invrrweight;
	int m_precision;
	AABB m_rifAABB;
	AABB m_sdfAABB;
	Solver::Options m_ceresoptions;

	// For Exact inside outside tests (and winding number calculation)
	bool m_isTriMesh;
	bool m_isSphere;
	Eigen::Matrix<float, Eigen::Dynamic,3,Eigen::RowMajor> V;
	std::vector<HDK_Sample::UT_Vector3T<float> > U;
	Eigen::Matrix<int,Eigen::Dynamic,3,Eigen::RowMajor> F;
	HDK_Sample::UT_SolidAngle<float, float> m_solid_angle; // double for UT_SolidAngle is not defined

	// For Exact signed distance function calculation

};

bool DirectConnectionCostFunction::Evaluate(double const* const* parameters,
                      double* residuals,
                      double** jacobians) const{
    VectorF v_i;
    PointF p1 = m_p1;
    PointF p2 = m_p2;
    Matrix3x3F dpdv0 = m_dpdv0;
    Matrix3x3F dvdv0 = m_dvdv0;

    VectorF error(0.0);
    Matrix3x3F derror(0.0);
    bool isSensorSample = m_isSensorSample;

    v_i.x = parameters[0][0];
    v_i.y = parameters[0][1];
    v_i.z = parameters[0][2];
    m_medium->computefdfBDPT(v_i, p1, p2, isSensorSample, dpdv0, dvdv0, error, derror);

    residuals[0] = error.x;
    residuals[1] = error.y;
    residuals[2] = error.z;
    if (jacobians != NULL && jacobians[0] != NULL){
        jacobians[0][0] = derror.m[0][0];
        jacobians[0][1] = derror.m[1][0];
        jacobians[0][2] = derror.m[2][0];
        jacobians[0][3] = derror.m[0][1];
        jacobians[0][4] = derror.m[1][1];
        jacobians[0][5] = derror.m[2][1];
        jacobians[0][6] = derror.m[0][2];
        jacobians[0][7] = derror.m[1][2];
        jacobians[0][8] = derror.m[2][2];
    }
	return true;
}

MTS_IMPLEMENT_CLASS_S(HeterogeneousRefractiveMedium, false, Medium)
MTS_EXPORT_PLUGIN(HeterogeneousRefractiveMedium, "Heterogeneous Refractive Medium");
MTS_NAMESPACE_END
