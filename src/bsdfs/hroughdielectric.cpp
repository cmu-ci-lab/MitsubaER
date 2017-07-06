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
#include <mitsuba/render/sampler.h>
#include <mitsuba/hw/basicshader.h>
#include <mitsuba/render/medium.h>
#include "microfacet.h"
#include "ior.h"

MTS_NAMESPACE_BEGIN

class HRoughDielectric : public BSDF {
public:
	HRoughDielectric(const Properties &props) : BSDF(props) {
		m_specularReflectance = new ConstantSpectrumTexture(
			props.getSpectrum("specularReflectance", Spectrum(1.0f)));
		m_specularTransmittance = new ConstantSpectrumTexture(
			props.getSpectrum("specularTransmittance", Spectrum(1.0f)));

		MicrofacetDistribution distr(props);
		m_type = distr.getType();
		m_sampleVisible = distr.getSampleVisible();

		m_alphaU = new ConstantFloatTexture(distr.getAlphaU());
		if (distr.getAlphaU() == distr.getAlphaV())
			m_alphaV = m_alphaU;
		else
			m_alphaV = new ConstantFloatTexture(distr.getAlphaV());
	}

	HRoughDielectric(Stream *stream, InstanceManager *manager)
	 : BSDF(stream, manager) {
		m_type = (MicrofacetDistribution::EType) stream->readUInt();
		m_sampleVisible = stream->readBool();
		m_alphaU = static_cast<Texture *>(manager->getInstance(stream));
		m_alphaV = static_cast<Texture *>(manager->getInstance(stream));
		m_specularReflectance = static_cast<Texture *>(manager->getInstance(stream));
		m_specularTransmittance = static_cast<Texture *>(manager->getInstance(stream));

		configure();
	}

	bool isheterogeneousbsdf() const {return true;}
	bool ishroughdielectric() const {return true;}


	void serialize(Stream *stream, InstanceManager *manager) const {
		BSDF::serialize(stream, manager);

		stream->writeUInt((uint32_t) m_type);
		stream->writeBool(m_sampleVisible);
		manager->serialize(stream, m_alphaU.get());
		manager->serialize(stream, m_alphaV.get());
		manager->serialize(stream, m_specularReflectance.get());
		manager->serialize(stream, m_specularTransmittance.get());
	}

	void configure() {
		unsigned int extraFlags = 0;
		if (m_alphaU != m_alphaV)
			extraFlags |= EAnisotropic;

		if (!m_alphaU->isConstant() || !m_alphaV->isConstant())
			extraFlags |= ESpatiallyVarying;

		m_components.clear();
		m_components.push_back(EGlossyReflection | EFrontSide
			| EBackSide | EUsesSampler | extraFlags
			| (m_specularReflectance->isConstant() ? 0 : ESpatiallyVarying));
		m_components.push_back(EGlossyTransmission | EFrontSide
			| EBackSide | EUsesSampler | ENonSymmetric | extraFlags
			| (m_specularTransmittance->isConstant() ? 0 : ESpatiallyVarying));

		/* Verify the input parameters and fix them if necessary */
		m_specularReflectance = ensureEnergyConservation(
			m_specularReflectance, "specularReflectance", 1.0f);
		m_specularTransmittance = ensureEnergyConservation(
			m_specularTransmittance, "specularTransmittance", 1.0f);

		m_usesRayDifferentials =
			m_alphaU->usesRayDifferentials() ||
			m_alphaV->usesRayDifferentials() ||
			m_specularReflectance->usesRayDifferentials() ||
			m_specularTransmittance->usesRayDifferentials();

		BSDF::configure();
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


	Spectrum eval(const BSDFSamplingRecord &bRec, const Point &p, EMeasure measure) const {
		Float p_eta, p_invEta;
		getEtaInvEta(p, p_eta, p_invEta);

		if (measure != ESolidAngle || Frame::cosTheta(bRec.wi) == 0)
			return Spectrum(0.0f);

		/* Determine the type of interaction */
		bool reflect = Frame::cosTheta(bRec.wi)
			* Frame::cosTheta(bRec.wo) > 0;

		Vector H;
		if (reflect) {
			/* Stop if this component was not requested */
			if ((bRec.component != -1 && bRec.component != 0)
				|| !(bRec.typeMask & EGlossyReflection))
				return Spectrum(0.0f);

			/* Calculate the reflection half-vector */
			H = normalize(bRec.wo+bRec.wi);
		} else {
			/* Stop if this component was not requested */
			if ((bRec.component != -1 && bRec.component != 1)
				|| !(bRec.typeMask & EGlossyTransmission))
				return Spectrum(0.0f);

			/* Calculate the transmission half-vector */
			Float eta = Frame::cosTheta(bRec.wi) > 0
				? p_eta : p_invEta;

			H = normalize(bRec.wi + bRec.wo*eta);
		}

		/* Ensure that the half-vector points into the
		   same hemisphere as the macrosurface normal */
		H *= math::signum(Frame::cosTheta(H));

		/* Construct the microfacet distribution matching the
		   roughness values at the current surface position. */
		MicrofacetDistribution distr(
			m_type,
			m_alphaU->eval(bRec.its).average(),
			m_alphaV->eval(bRec.its).average(),
			m_sampleVisible
		);

		/* Evaluate the microfacet normal distribution */
		const Float D = distr.eval(H);
		if (D == 0)
			return Spectrum(0.0f);

		/* Fresnel factor */
		const Float F = fresnelDielectricExt(dot(bRec.wi, H), p_eta);

		/* Smith's shadow-masking function */
		const Float G = distr.G(bRec.wi, bRec.wo, H);

		if (reflect) {
			/* Calculate the total amount of reflection */
			Float value = F * D * G /
				(4.0f * std::abs(Frame::cosTheta(bRec.wi)));

			return m_specularReflectance->eval(bRec.its) * value;
		} else {
			Float eta = Frame::cosTheta(bRec.wi) > 0.0f ? p_eta : p_invEta;

			/* Calculate the total amount of transmission */
			Float sqrtDenom = dot(bRec.wi, H) + eta * dot(bRec.wo, H);
			Float value = ((1 - F) * D * G * eta * eta
				* dot(bRec.wi, H) * dot(bRec.wo, H)) /
				(Frame::cosTheta(bRec.wi) * sqrtDenom * sqrtDenom);

			/* Missing term in the original paper: account for the solid angle
			   compression when tracing radiance -- this is necessary for
			   bidirectional methods */
			Float factor = (bRec.mode == ERadiance)
				? (Frame::cosTheta(bRec.wi) > 0 ? p_invEta : p_eta) : 1.0f;

			return m_specularTransmittance->eval(bRec.its)
				* std::abs(value * factor * factor);
		}
	}

	Float pdf(const BSDFSamplingRecord &bRec, const Point &p, EMeasure measure) const {
		if (measure != ESolidAngle)
			return 0.0f;

		Float p_eta, p_invEta;
		getEtaInvEta(p, p_eta, p_invEta);

		/* Determine the type of interaction */
		bool hasReflection   = ((bRec.component == -1 || bRec.component == 0)
							  && (bRec.typeMask & EGlossyReflection)),
		     hasTransmission = ((bRec.component == -1 || bRec.component == 1)
							  && (bRec.typeMask & EGlossyTransmission)),
		     reflect         = Frame::cosTheta(bRec.wi)
				             * Frame::cosTheta(bRec.wo) > 0;

		Vector H;
		Float dwh_dwo;

		if (reflect) {
			/* Zero probability if this component was not requested */
			if ((bRec.component != -1 && bRec.component != 0)
				|| !(bRec.typeMask & EGlossyReflection))
				return 0.0f;

			/* Calculate the reflection half-vector */
			H = normalize(bRec.wo+bRec.wi);

			/* Jacobian of the half-direction mapping */
			dwh_dwo = 1.0f / (4.0f * dot(bRec.wo, H));
		} else {
			/* Zero probability if this component was not requested */
			if ((bRec.component != -1 && bRec.component != 1)
				|| !(bRec.typeMask & EGlossyTransmission))
				return 0.0f;

			/* Calculate the transmission half-vector */
			Float eta = Frame::cosTheta(bRec.wi) > 0
				? p_eta : p_invEta;

			H = normalize(bRec.wi + bRec.wo*eta);

			/* Jacobian of the half-direction mapping */
			Float sqrtDenom = dot(bRec.wi, H) + eta * dot(bRec.wo, H);
			dwh_dwo = (eta*eta * dot(bRec.wo, H)) / (sqrtDenom*sqrtDenom);
		}

		/* Ensure that the half-vector points into the
		   same hemisphere as the macrosurface normal */
		H *= math::signum(Frame::cosTheta(H));

		/* Construct the microfacet distribution matching the
		   roughness values at the current surface position. */
		MicrofacetDistribution sampleDistr(
			m_type,
			m_alphaU->eval(bRec.its).average(),
			m_alphaV->eval(bRec.its).average(),
			m_sampleVisible
		);

		/* Trick by Walter et al.: slightly scale the roughness values to
		   reduce importance sampling weights. Not needed for the
		   Heitz and D'Eon sampling technique. */
		if (!m_sampleVisible)
			sampleDistr.scaleAlpha(1.2f - 0.2f * std::sqrt(
				std::abs(Frame::cosTheta(bRec.wi))));

		/* Evaluate the microfacet model sampling density function */
		Float prob = sampleDistr.pdf(math::signum(Frame::cosTheta(bRec.wi)) * bRec.wi, H);

		if (hasTransmission && hasReflection) {
			Float F = fresnelDielectricExt(dot(bRec.wi, H), p_eta);
			prob *= reflect ? F : (1-F);
		}

		return std::abs(prob * dwh_dwo);
	}

	Spectrum sample(BSDFSamplingRecord &bRec, const Point &p, const Point2 &_sample) const {
		Point2 sample(_sample);

		bool hasReflection = ((bRec.component == -1 || bRec.component == 0)
							  && (bRec.typeMask & EGlossyReflection)),
		     hasTransmission = ((bRec.component == -1 || bRec.component == 1)
							  && (bRec.typeMask & EGlossyTransmission)),
		     sampleReflection = hasReflection;

		if (!hasReflection && !hasTransmission)
			return Spectrum(0.0f);

		Float p_eta, p_invEta;
		getEtaInvEta(p, p_eta, p_invEta);

		/* Construct the microfacet distribution matching the
		   roughness values at the current surface position. */
		MicrofacetDistribution distr(
			m_type,
			m_alphaU->eval(bRec.its).average(),
			m_alphaV->eval(bRec.its).average(),
			m_sampleVisible
		);

		/* Trick by Walter et al.: slightly scale the roughness values to
		   reduce importance sampling weights. Not needed for the
		   Heitz and D'Eon sampling technique. */
		MicrofacetDistribution sampleDistr(distr);
		if (!m_sampleVisible)
			sampleDistr.scaleAlpha(1.2f - 0.2f * std::sqrt(
				std::abs(Frame::cosTheta(bRec.wi))));

		/* Sample M, the microfacet normal */
		Float microfacetPDF;
		const Normal m = sampleDistr.sample(math::signum(Frame::cosTheta(bRec.wi)) * bRec.wi, sample, microfacetPDF);
		if (microfacetPDF == 0)
			return Spectrum(0.0f);

		Float cosThetaT;
		Float F = fresnelDielectricExt(dot(bRec.wi, m), cosThetaT, p_eta);
		Spectrum weight(1.0f);

		if (hasReflection && hasTransmission) {
			if (bRec.sampler->next1D() > F)
				sampleReflection = false;
		} else {
			weight = Spectrum(hasReflection ? F : (1-F));
		}

		if (sampleReflection) {
			/* Perfect specular reflection based on the microfacet normal */
			bRec.wo = reflect(bRec.wi, m);
			bRec.eta = 1.0f;
			bRec.sampledComponent = 0;
			bRec.sampledType = EGlossyReflection;

			/* Side check */
			if (Frame::cosTheta(bRec.wi) * Frame::cosTheta(bRec.wo) <= 0)
				return Spectrum(0.0f);

			weight *= m_specularReflectance->eval(bRec.its);
		} else {
			if (cosThetaT == 0)
				return Spectrum(0.0f);

			/* Perfect specular transmission based on the microfacet normal */
			bRec.wo = refract(bRec.wi, m, p_eta, cosThetaT);
			bRec.eta = cosThetaT < 0 ? p_eta : p_invEta;
			bRec.sampledComponent = 1;
			bRec.sampledType = EGlossyTransmission;

			/* Side check */
			if (Frame::cosTheta(bRec.wi) * Frame::cosTheta(bRec.wo) >= 0)
				return Spectrum(0.0f);

			/* Radiance must be scaled to account for the solid angle compression
			   that occurs when crossing the interface. */
			Float factor = (bRec.mode == ERadiance)
				? (cosThetaT < 0 ? p_invEta : p_eta) : 1.0f;

			weight *= m_specularTransmittance->eval(bRec.its) * (factor * factor);
		}

		if (m_sampleVisible)
			weight *= distr.smithG1(bRec.wo, m);
		else
			weight *= std::abs(distr.eval(m) * distr.G(bRec.wi, bRec.wo, m)
				* dot(bRec.wi, m) / (microfacetPDF * Frame::cosTheta(bRec.wi)));

		return weight;
	}

	Spectrum sample(BSDFSamplingRecord &bRec, Float &pdf, const Point &p, const Point2 &_sample) const {
		Point2 sample(_sample);

		bool hasReflection = ((bRec.component == -1 || bRec.component == 0)
							  && (bRec.typeMask & EGlossyReflection)),
		     hasTransmission = ((bRec.component == -1 || bRec.component == 1)
							  && (bRec.typeMask & EGlossyTransmission)),
		     sampleReflection = hasReflection;

		if (!hasReflection && !hasTransmission)
			return Spectrum(0.0f);

		Float p_eta, p_invEta;
		getEtaInvEta(p, p_eta, p_invEta);

		/* Construct the microfacet distribution matching the
		   roughness values at the current surface position. */
		MicrofacetDistribution distr(
			m_type,
			m_alphaU->eval(bRec.its).average(),
			m_alphaV->eval(bRec.its).average(),
			m_sampleVisible
		);

		/* Trick by Walter et al.: slightly scale the roughness values to
		   reduce importance sampling weights. Not needed for the
		   Heitz and D'Eon sampling technique. */
		MicrofacetDistribution sampleDistr(distr);
		if (!m_sampleVisible)
			sampleDistr.scaleAlpha(1.2f - 0.2f * std::sqrt(
				std::abs(Frame::cosTheta(bRec.wi))));

		/* Sample M, the microfacet normal */
		Float microfacetPDF;
		const Normal m = sampleDistr.sample(math::signum(Frame::cosTheta(bRec.wi)) * bRec.wi, sample, microfacetPDF);
		if (microfacetPDF == 0)
			return Spectrum(0.0f);
		pdf = microfacetPDF;

		Float cosThetaT;
		Float F = fresnelDielectricExt(dot(bRec.wi, m), cosThetaT, p_eta);
		Spectrum weight(1.0f);

		if (hasReflection && hasTransmission) {
			if (bRec.sampler->next1D() > F) {
				sampleReflection = false;
				pdf *= 1-F;
			} else {
				pdf *= F;
			}
		} else {
			weight *= hasReflection ? F : (1-F);
		}

		Float dwh_dwo;
		if (sampleReflection) {
			/* Perfect specular reflection based on the microfacet normal */
			bRec.wo = reflect(bRec.wi, m);
			bRec.eta = 1.0f;
			bRec.sampledComponent = 0;
			bRec.sampledType = EGlossyReflection;

			/* Side check */
			if (Frame::cosTheta(bRec.wi) * Frame::cosTheta(bRec.wo) <= 0)
				return Spectrum(0.0f);

			weight *= m_specularReflectance->eval(bRec.its);

			/* Jacobian of the half-direction mapping */
			dwh_dwo = 1.0f / (4.0f * dot(bRec.wo, m));
		} else {
			if (cosThetaT == 0)
				return Spectrum(0.0f);

			/* Perfect specular transmission based on the microfacet normal */
			bRec.wo = refract(bRec.wi, m, p_eta, cosThetaT);
			bRec.eta = cosThetaT < 0 ? p_eta : p_invEta;
			bRec.sampledComponent = 1;
			bRec.sampledType = EGlossyTransmission;

			/* Side check */
			if (Frame::cosTheta(bRec.wi) * Frame::cosTheta(bRec.wo) >= 0)
				return Spectrum(0.0f);

			/* Radiance must be scaled to account for the solid angle compression
			   that occurs when crossing the interface. */
			Float factor = (bRec.mode == ERadiance)
				? (cosThetaT < 0 ? p_invEta : p_eta) : 1.0f;

			weight *= m_specularTransmittance->eval(bRec.its) * (factor * factor);

			/* Jacobian of the half-direction mapping */
			Float sqrtDenom = dot(bRec.wi, m) + bRec.eta * dot(bRec.wo, m);
			dwh_dwo = (bRec.eta*bRec.eta * dot(bRec.wo, m)) / (sqrtDenom*sqrtDenom);
		}

		if (m_sampleVisible)
			weight *= distr.smithG1(bRec.wo, m);
		else
			weight *= std::abs(distr.eval(m) * distr.G(bRec.wi, bRec.wo, m)
				* dot(bRec.wi, m) / (microfacetPDF * Frame::cosTheta(bRec.wi)));

		pdf *= std::abs(dwh_dwo);

		return weight;
	}

	void addChild(const std::string &name, ConfigurableObject *child) {
		if (child->getClass()->derivesFrom(MTS_CLASS(Texture))) {
			if (name == "alpha")
				m_alphaU = m_alphaV = static_cast<Texture *>(child);
			else if (name == "alphaU")
				m_alphaU = static_cast<Texture *>(child);
			else if (name == "alphaV")
				m_alphaV = static_cast<Texture *>(child);
			else if (name == "specularReflectance")
				m_specularReflectance = static_cast<Texture *>(child);
			else if (name == "specularTransmittance")
				m_specularTransmittance = static_cast<Texture *>(child);
			else
				BSDF::addChild(name, child);
		} else {
			BSDF::addChild(name, child);
		}
	}

	inline void getEtaInvEta(const Point &p, Float &eta, Float &invEta) const{
		eta = m_shape->getInteriorMedium()->getRIF(p);
		invEta = 1/eta;
	}


	Float getRoughness(const Intersection &its, int component) const {
		return 0.5f * (m_alphaU->eval(its).average()
			+ m_alphaV->eval(its).average());
	}

	std::string toString() const {
		std::ostringstream oss;
		oss << "HRoughDielectric[" << endl
			<< "  id = \"" << getID() << "\"," << endl
			<< "  distribution = " << MicrofacetDistribution::distributionName(m_type) << "," << endl
			<< "  sampleVisible = " << m_sampleVisible << "," << endl
			<< "  alphaU = " << indent(m_alphaU->toString()) << "," << endl
			<< "  alphaV = " << indent(m_alphaV->toString()) << "," << endl
			<< "  specularReflectance = " << indent(m_specularReflectance->toString()) << "," << endl
			<< "  specularTransmittance = " << indent(m_specularTransmittance->toString()) << endl
			<< "]";
		return oss.str();
	}

	Shader *createShader(Renderer *renderer) const;

	MTS_DECLARE_CLASS()
private:
	MicrofacetDistribution::EType m_type;
	ref<Texture> m_specularTransmittance;
	ref<Texture> m_specularReflectance;
	ref<Texture> m_alphaU, m_alphaV;
	bool m_sampleVisible;
};

/* Fake glass shader -- it is really hopeless to visualize
   this material in the VPL renderer, so let's try to do at least
   something that suggests the presence of a transparent boundary */
class HRoughDielectricShader : public Shader {
public:
	HRoughDielectricShader(Renderer *renderer) :
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

Shader *HRoughDielectric::createShader(Renderer *renderer) const {
	return new HRoughDielectricShader(renderer);
}

MTS_IMPLEMENT_CLASS(HRoughDielectricShader, false, Shader)
MTS_IMPLEMENT_CLASS_S(HRoughDielectric, false, BSDF)
MTS_EXPORT_PLUGIN(HRoughDielectric, "Heterogeneous Rough dielectric BSDF");
MTS_NAMESPACE_END
