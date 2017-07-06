/*
    This file is part of Mitsuba, a physically based rendering system.

    Copyright (c) 2007-2014 by Wenzel Jakob and others.

    Mitsuba is free software; you can redistribute it and/or modify
    it under the terms of the GNU General Public License Version 3
    as published by the Free Software Foundation.

    Mitsuba is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with this program. If not, see <http://www.gnu.org/licenses/>.
*/

#include <mitsuba/render/bsdf.h>
#include <mitsuba/hw/basicshader.h>
#include <mitsuba/render/medium.h>
#include "ior.h"

MTS_NAMESPACE_BEGIN

class HSmoothDielectric : public BSDF {
public:
	HSmoothDielectric(const Properties &props) : BSDF(props) {
		m_specularReflectance = new ConstantSpectrumTexture(
			props.getSpectrum("specularReflectance", Spectrum(1.0f)));
		m_specularTransmittance = new ConstantSpectrumTexture(
			props.getSpectrum("specularTransmittance", Spectrum(1.0f)));
	}

	HSmoothDielectric(Stream *stream, InstanceManager *manager)
			: BSDF(stream, manager) {
		m_specularReflectance = static_cast<Texture *>(manager->getInstance(stream));
		m_specularTransmittance = static_cast<Texture *>(manager->getInstance(stream));
		configure();
	}

	bool isheterogeneousbsdf() const {return true;}
	bool ishdielectric() const {return true;}

	void serialize(Stream *stream, InstanceManager *manager) const {
		BSDF::serialize(stream, manager);

		manager->serialize(stream, m_specularReflectance.get());
		manager->serialize(stream, m_specularTransmittance.get());
	}

	void configure() {
		/* Verify the input parameters and fix them if necessary */
		m_specularReflectance = ensureEnergyConservation(
			m_specularReflectance, "specularReflectance", 1.0f);
		m_specularTransmittance = ensureEnergyConservation(
			m_specularTransmittance, "specularTransmittance", 1.0f);

		m_components.clear();
		m_components.push_back(EDeltaReflection | EFrontSide | EBackSide
			| (m_specularReflectance->isConstant() ? 0 : ESpatiallyVarying));
		m_components.push_back(EDeltaTransmission | EFrontSide | EBackSide | ENonSymmetric
			| (m_specularTransmittance->isConstant() ? 0 : ESpatiallyVarying));

		m_usesRayDifferentials =
			m_specularReflectance->usesRayDifferentials() ||
			m_specularTransmittance->usesRayDifferentials();

		BSDF::configure();
	}

	void addChild(const std::string &name, ConfigurableObject *child) {
		if (child->getClass()->derivesFrom(MTS_CLASS(Texture))) {
			if (name == "specularReflectance")
				m_specularReflectance = static_cast<Texture *>(child);
			else if (name == "specularTransmittance")
				m_specularTransmittance = static_cast<Texture *>(child);
			else
				BSDF::addChild(name, child);
		} else {
			BSDF::addChild(name, child);
		}
	}

	/// Reflection in local coordinates
	inline Vector reflect(const Vector &wi) const {
		return Vector(-wi.x, -wi.y, wi.z);
	}

	Spectrum eval(const BSDFSamplingRecord &bRec, EMeasure measure) const {
		NotImplementedError("eval");
		return Spectrum(0.0);
	}

	Float pdf(const BSDFSamplingRecord &bRec, EMeasure measure) const{
		NotImplementedError("pdf");
		return Float(0.0);
	}

	Spectrum sample(BSDFSamplingRecord &bRec, Float &pdf, const Point2 &sample) const{
		NotImplementedError("sample");
		return Spectrum(0.0);
	}

	Spectrum sample(BSDFSamplingRecord &bRec, const Point2 &sample) const {
		NotImplementedError("sample");
		return Spectrum(0.0);
	}

	Float getEta() const {
		NotImplementedError("getETA");
		return 1.0f;
	}

	inline void getEtaInvEta(const Point &p, Float &eta, Float &invEta) const{
		eta = m_shape->getInteriorMedium()->getRIF(p);
		invEta = 1/eta;
	}

	/// Refraction in local coordinates
	inline Vector refract(const Vector &wi, const Point &p, Float cosThetaT) const {
		Float eta, invEta;
		getEtaInvEta(p, eta, invEta);
		Float scale = -(cosThetaT < 0 ? invEta : eta);
		return Vector(scale*wi.x, scale*wi.y, cosThetaT);
	}

	Spectrum eval(const BSDFSamplingRecord &bRec, const Point &p, EMeasure measure) const {
		bool sampleReflection   = (bRec.typeMask & EDeltaReflection)
				&& (bRec.component == -1 || bRec.component == 0) && measure == EDiscrete;
		bool sampleTransmission = (bRec.typeMask & EDeltaTransmission)
				&& (bRec.component == -1 || bRec.component == 1) && measure == EDiscrete;

		Float eta, invEta;
		getEtaInvEta(p, eta, invEta);

		Float cosThetaT;
		Float F = fresnelDielectricExt(Frame::cosTheta(bRec.wi), cosThetaT, eta);

		if (Frame::cosTheta(bRec.wi) * Frame::cosTheta(bRec.wo) >= 0) {
			if (!sampleReflection || std::abs(dot(reflect(bRec.wi), bRec.wo)-1) > DeltaEpsilon)
				return Spectrum(0.0f);

			return m_specularReflectance->eval(bRec.its) * F;
		} else {
			if (!sampleTransmission || std::abs(dot(refract(bRec.wi, p, cosThetaT), bRec.wo)-1) > DeltaEpsilon)
				return Spectrum(0.0f);

			/* Radiance must be scaled to account for the solid angle compression
			   that occurs when crossing the interface. */
			Float factor = (bRec.mode == ERadiance)
				? (cosThetaT < 0 ? invEta : eta) : 1.0f;

			return m_specularTransmittance->eval(bRec.its)  * factor * factor * (1 - F);
		}
	}

	Float pdf(const BSDFSamplingRecord &bRec, const Point &p, EMeasure measure) const {
		bool sampleReflection   = (bRec.typeMask & EDeltaReflection)
				&& (bRec.component == -1 || bRec.component == 0) && measure == EDiscrete;
		bool sampleTransmission = (bRec.typeMask & EDeltaTransmission)
				&& (bRec.component == -1 || bRec.component == 1) && measure == EDiscrete;

		Float eta, invEta;
		getEtaInvEta(p, eta, invEta);

		Float cosThetaT;
		Float F = fresnelDielectricExt(Frame::cosTheta(bRec.wi), cosThetaT, eta);

		if (Frame::cosTheta(bRec.wi) * Frame::cosTheta(bRec.wo) >= 0) {
			if (!sampleReflection || std::abs(dot(reflect(bRec.wi), bRec.wo)-1) > DeltaEpsilon)
				return 0.0f;

			return sampleTransmission ? F : 1.0f;
		} else {
			if (!sampleTransmission || std::abs(dot(refract(bRec.wi, p, cosThetaT), bRec.wo)-1) > DeltaEpsilon)
				return 0.0f;

			return sampleReflection ? 1-F : 1.0f;
		}
	}

	Spectrum sample(BSDFSamplingRecord &bRec, Float &pdf, const Point &p, const Point2 &sample) const {
		bool sampleReflection   = (bRec.typeMask & EDeltaReflection)
				&& (bRec.component == -1 || bRec.component == 0);
		bool sampleTransmission = (bRec.typeMask & EDeltaTransmission)
				&& (bRec.component == -1 || bRec.component == 1);

		Float eta, invEta;
		getEtaInvEta(p, eta, invEta);

		Float cosThetaT;
		Float F = fresnelDielectricExt(Frame::cosTheta(bRec.wi), cosThetaT, eta);

		if (sampleTransmission && sampleReflection) {
			if (sample.x <= F) {
				bRec.sampledComponent = 0;
				bRec.sampledType = EDeltaReflection;
				bRec.wo = reflect(bRec.wi);
				bRec.eta = 1.0f;
				pdf = F;

				return m_specularReflectance->eval(bRec.its);
			} else {
				bRec.sampledComponent = 1;
				bRec.sampledType = EDeltaTransmission;
				bRec.wo = refract(bRec.wi, p, cosThetaT);
				bRec.eta = cosThetaT < 0 ? eta : invEta;
				pdf = 1-F;

				/* Radiance must be scaled to account for the solid angle compression
				   that occurs when crossing the interface. */
				Float factor = (bRec.mode == ERadiance)
					? (cosThetaT < 0 ? invEta : eta) : 1.0f;

				return m_specularTransmittance->eval(bRec.its) * (factor * factor);
			}
		} else if (sampleReflection) {
			bRec.sampledComponent = 0;
			bRec.sampledType = EDeltaReflection;
			bRec.wo = reflect(bRec.wi);
			bRec.eta = 1.0f;
			pdf = 1.0f;

			return m_specularReflectance->eval(bRec.its) * F;
		} else if (sampleTransmission) {
			bRec.sampledComponent = 1;
			bRec.sampledType = EDeltaTransmission;
			bRec.wo = refract(bRec.wi, p, cosThetaT);
			bRec.eta = cosThetaT < 0 ? eta : invEta;
			pdf = 1.0f;

			/* Radiance must be scaled to account for the solid angle compression
			   that occurs when crossing the interface. */
			Float factor = (bRec.mode == ERadiance)
				? (cosThetaT < 0 ? invEta : eta) : 1.0f;

			return m_specularTransmittance->eval(bRec.its) * (factor * factor * (1-F));
		}

		return Spectrum(0.0f);
	}

	Spectrum sample(BSDFSamplingRecord &bRec, const Point &p, const Point2 &sample) const {
		bool sampleReflection   = (bRec.typeMask & EDeltaReflection)
				&& (bRec.component == -1 || bRec.component == 0);
		bool sampleTransmission = (bRec.typeMask & EDeltaTransmission)
				&& (bRec.component == -1 || bRec.component == 1);

		Float eta, invEta;
		getEtaInvEta(p, eta, invEta);

		Float cosThetaT;
		Float F = fresnelDielectricExt(Frame::cosTheta(bRec.wi), cosThetaT, eta);

		if (sampleTransmission && sampleReflection) {
			if (sample.x <= F) {
				bRec.sampledComponent = 0;
				bRec.sampledType = EDeltaReflection;
				bRec.wo = reflect(bRec.wi);
				bRec.eta = 1.0f;

				return m_specularReflectance->eval(bRec.its);
			} else {
				bRec.sampledComponent = 1;
				bRec.sampledType = EDeltaTransmission;
				bRec.wo = refract(bRec.wi, p, cosThetaT);
				bRec.eta = cosThetaT < 0 ? eta : invEta;

				/* Radiance must be scaled to account for the solid angle compression
				   that occurs when crossing the interface. */
				Float factor = (bRec.mode == ERadiance)
					? (cosThetaT < 0 ? invEta : eta) : 1.0f;

				return m_specularTransmittance->eval(bRec.its) * (factor * factor);
			}
		} else if (sampleReflection) {
			bRec.sampledComponent = 0;
			bRec.sampledType = EDeltaReflection;
			bRec.wo = reflect(bRec.wi);
			bRec.eta = 1.0f;

			return m_specularReflectance->eval(bRec.its) * F;
		} else if (sampleTransmission) {
			bRec.sampledComponent = 1;
			bRec.sampledType = EDeltaTransmission;
			bRec.wo = refract(bRec.wi, p, cosThetaT);
			bRec.eta = cosThetaT < 0 ? eta : invEta;

			/* Radiance must be scaled to account for the solid angle compression
			   that occurs when crossing the interface. */
			Float factor = (bRec.mode == ERadiance)
				? (cosThetaT < 0 ? invEta : eta) : 1.0f;

			return m_specularTransmittance->eval(bRec.its) * (factor * factor * (1-F));
		}

		return Spectrum(0.0f);
	}

	Float getRoughness(const Intersection &its, int component) const {
		return 0.0f;
	}

	std::string toString() const {
		std::ostringstream oss;
		oss << "SmoothDielectric[" << endl
			<< "  id = \"" << getID() << "\"," << endl
			<< "  specularReflectance = " << indent(m_specularReflectance->toString()) << "," << endl
			<< "  specularTransmittance = " << indent(m_specularTransmittance->toString()) << endl
			<< "]";
		return oss.str();
	}

	Shader *createShader(Renderer *renderer) const;

	MTS_DECLARE_CLASS()
private:
	ref<Texture> m_specularTransmittance;
	ref<Texture> m_specularReflectance;
};


class HSmoothDielectricShader : public Shader {
public:
	HSmoothDielectricShader(Renderer *renderer) :
		Shader(renderer, EBSDFShader) {
		m_flags = ETransparent;
	}

	Float getAlpha() const {
		return 0.3f;
	}

	void generateCode(std::ostringstream &oss,
			const std::string &evalName,
			const std::vector<std::string> &depNames) const {
		oss << "vec3 " << evalName << "(vec2 uv, vec3 wi, vec3 wo) {" << endl
			<< "    if (cosTheta(wi) < 0.0 || cosTheta(wo) < 0.0)" << endl
			<< "    	return vec3(0.0);" << endl
			<< "    return vec3(inv_pi * cosTheta(wo));" << endl
			<< "}" << endl
			<< endl
			<< "vec3 " << evalName << "_diffuse(vec2 uv, vec3 wi, vec3 wo) {" << endl
			<< "    return " << evalName << "(uv, wi, wo);" << endl
			<< "}" << endl;
	}


	MTS_DECLARE_CLASS()
};

Shader *HSmoothDielectric::createShader(Renderer *renderer) const {
	return new HSmoothDielectricShader(renderer);
}

MTS_IMPLEMENT_CLASS(HSmoothDielectricShader, false, Shader)
MTS_IMPLEMENT_CLASS_S(HSmoothDielectric, false, BSDF)
MTS_EXPORT_PLUGIN(HSmoothDielectric, "Heterogeneous Smooth Dielectric");
MTS_NAMESPACE_END
