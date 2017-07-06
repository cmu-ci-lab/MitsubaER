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
#include <mitsuba/core/warp.h>
#include <boost/math/special_functions.hpp>
#include "phaseutils.h"

MTS_NAMESPACE_BEGIN

/**
 * von Mises-Fisher Phase Function
 */
class vMFPhaseFunction : public PhaseFunction {
public:
	vMFPhaseFunction(const Properties &props) 
		: PhaseFunction(props) {
		m_kapa = props.getFloat("kappa");
		m_weight = props.getFloat("weight", 1.0f);
	}

	vMFPhaseFunction(Stream *stream, InstanceManager *manager) 
		: PhaseFunction(stream, manager) {
		m_kapa = stream->readDouble();
		m_weight = stream->readDouble();
	}

	virtual ~vMFPhaseFunction() { }


	void serialize(Stream *stream, InstanceManager *manager) const {
		PhaseFunction::serialize(stream, manager);
		stream->writeDouble(m_kapa);
		stream->writeDouble(m_weight);
	}

	void configure() {
		norm = vmfNorm(m_kapa);
		Log(EInfo, toString().c_str());
	}

	Float sample(PhaseFunctionSamplingRecord &pRec,
			Sampler *sampler) const {
		Point2 sample(sampler->next2D());
		
		if ( sampler->next1D() < m_weight ) {
			// sample front
			Vector dir;
			Float theta = 2.0f*M_PI*sampler->next1D();
            Float xi = sampler->next1D();
			dir.z = static_cast<Float>(vmfInverseCdf(xi, m_kapa));
			Float v = std::sqrt(std::max(1.0f - dir.z*dir.z, 0.0f));
			dir.x = std::cos(theta)*v;
			dir.y = std::sin(theta)*v;
			pRec.wo = Frame(-pRec.wi).toWorld(dir);
		} else {
			// sample back (uniform)
			pRec.wo = warp::squareToUniformSphere(sampler->next2D());
		}
		return 1.0f;
	}

	Float sample(PhaseFunctionSamplingRecord &pRec,
			Float &pdf, Sampler *sampler) const {
		vMFPhaseFunction::sample(pRec, sampler);
		pdf = vMFPhaseFunction::eval(pRec);
		return 1.0f;
	}

	Float eval(const PhaseFunctionSamplingRecord &pRec) const {
		double val = dot(-pRec.wi, pRec.wo);
		return m_weight * static_cast<Float>(vmfPdf(val, m_kapa, norm)) \
				+ (1.0 - m_weight) / (4.0 * M_PI);
	}

	std::string toString() const {
		std::ostringstream oss;
		oss << "vMFPhaseFunction[" << endl
			<< "   kappa = " << m_kapa << endl
			<< "   weight = " << m_weight << endl
			<< "]";
		return oss.str();
	}

	MTS_DECLARE_CLASS()

private:
	double m_kapa, m_weight;
	double norm;
};


MTS_IMPLEMENT_CLASS_S(vMFPhaseFunction, false, PhaseFunction)
MTS_EXPORT_PLUGIN(vMFPhaseFunction, "von Mises-Fisher phase function");
MTS_NAMESPACE_END
