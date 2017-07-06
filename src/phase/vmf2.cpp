/*
    This file is part of Mitsuba, a physically based rendering system.

    Copyright (c) 2007-2011 by Wenzel Jakob and others.

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

#include <mitsuba/render/phase.h>
#include <mitsuba/render/sampler.h>
#include <mitsuba/core/frame.h>
#include <boost/math/special_functions.hpp>
#include "phaseutils.h"

MTS_NAMESPACE_BEGIN

/**
 * Mixture of two forward von Mises-Fisher Phase Function
 */
class vMF2PhaseFunction : public PhaseFunction {
public:
	vMF2PhaseFunction(const Properties &props)
		: PhaseFunction(props) {
		m_kapa = props.getFloat("kappa");
		m_kapa2 = props.getFloat("kappa2");
		m_weight = props.getFloat("weight", 1.0f);
	}

	vMF2PhaseFunction(Stream *stream, InstanceManager *manager)
		: PhaseFunction(stream, manager) {
		m_kapa = stream->readDouble();
		m_kapa2 = stream->readDouble();
		m_weight = stream->readDouble();
	}

	virtual ~vMF2PhaseFunction() { }


	void serialize(Stream *stream, InstanceManager *manager) const {
		PhaseFunction::serialize(stream, manager);
		stream->writeDouble(m_kapa);
		stream->writeDouble(m_kapa2);
		stream->writeDouble(m_weight);
	}

	void configure() {
		norm = vmfNorm(m_kapa);
		norm2 = vmfNorm(m_kapa2);
		Log(EInfo, toString().c_str());
	}

	Float sample(PhaseFunctionSamplingRecord &pRec,
			Sampler *sampler) const {
		Point2 sample(sampler->next2D());

		if ( sampler->next1D() < m_weight ) {
			// sample first vMF
			Vector dir;
			Float theta = 2.0f*M_PI*sampler->next1D();
            Float xi = sampler->next1D();
			dir.z = static_cast<Float>(vmfInverseCdf(xi, m_kapa));
			Float v = std::sqrt(std::max(1.0f - dir.z*dir.z, 0.0f));
			dir.x = std::cos(theta)*v;
			dir.y = std::sin(theta)*v;
			pRec.wo = Frame(-pRec.wi).toWorld(dir);
		} else {
			// sample second vMF
			Vector dir;
			Float theta = 2.0f*M_PI*sampler->next1D();
            Float xi = sampler->next1D();
			dir.z = static_cast<Float>(vmfInverseCdf(xi, m_kapa2));
			Float v = std::sqrt(std::max(1.0f - dir.z*dir.z, 0.0f));
			dir.x = std::cos(theta)*v;
			dir.y = std::sin(theta)*v;
			pRec.wo = Frame(-pRec.wi).toWorld(dir);
		}
		return 1.0f;
	}

	Float sample(PhaseFunctionSamplingRecord &pRec,
			Float &pdf, Sampler *sampler) const {
		vMF2PhaseFunction::sample(pRec, sampler);
		pdf = vMF2PhaseFunction::eval(pRec);
		return 1.0f;
	}

	Float eval(const PhaseFunctionSamplingRecord &pRec) const {
		double val = dot(-pRec.wi, pRec.wo);
		return static_cast<Float>(m_weight * vmfPdf(val, m_kapa, norm) \
				+ (1.0 - m_weight) * vmfPdf(val, m_kapa2, norm2));
	}

	std::string toString() const {
		std::ostringstream oss;
		oss << "vMF2PhaseFunction[" << endl
			<< "   kappa = " << m_kapa << endl
			<< "   kappa2 = " << m_kapa2 << endl
			<< "   weight = " << m_weight << endl
			<< "]";
		return oss.str();
	}

	MTS_DECLARE_CLASS()

private:
	double m_kapa, m_weight, m_kapa2;
	double norm, norm2;
};


MTS_IMPLEMENT_CLASS_S(vMF2PhaseFunction, false, PhaseFunction)
MTS_EXPORT_PLUGIN(vMF2PhaseFunction, "Mixture of two von Mises-Fisher phase function");
MTS_NAMESPACE_END
