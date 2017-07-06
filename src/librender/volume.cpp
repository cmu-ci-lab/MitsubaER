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

#include <mitsuba/render/volume.h>

#if defined(__LINUX__)
#include <sys/mman.h>
#include <fcntl.h>
#endif

MTS_NAMESPACE_BEGIN

VolumeDataSource::VolumeDataSource(Stream *stream, InstanceManager *manager) :
	ConfigurableObject(stream, manager) {
	m_aabb = AABB(stream);
}

VolumeDataSource::VolumeDataSource(const Properties &props) : ConfigurableObject(props) { }

VolumeDataSource::~VolumeDataSource() { }

void VolumeDataSource::serialize(Stream *stream, InstanceManager *manager) const {
	ConfigurableObject::serialize(stream, manager);
	m_aabb.serialize(stream);
}

bool VolumeDataSource::insideVolumeLimits(const PointF &p) const{
	Log(EError, "'%s': does not implement insideVolumeLimits()!", getClass()->getName().c_str());
	return false;
}

Vector3i VolumeDataSource::getResolution() const{
	Log(EError, "'%s': does not implement getResolution()!", getClass()->getName().c_str());
	return Vector3i(0);
}

Float VolumeDataSource::maxSDFError() const{
	Log(EError, "'%s': does not implement maxSDFError()!", getClass()->getName().c_str());
	return Float(0);
}

FLOAT VolumeDataSource::value(const PointF &p) const{
	Log(EError, "'%s': does not implement value()!", getClass()->getName().c_str());
	return 0.0;
}

VectorF VolumeDataSource::gradient(const PointF &p) const{
	Log(EError, "'%s': does not implement gradient()!", getClass()->getName().c_str());
	return VectorF(0.0);
}

Matrix3x3F VolumeDataSource::hessian(const PointF &p) const{
	Log(EError, "'%s': does not implement hessian()!", getClass()->getName().c_str());
	return Matrix3x3F(0.0);
}

void VolumeDataSource::valueAndGradient(const PointF &p, FLOAT &f, VectorF &v) const{
	Log(EError, "'%s': does not implement valueAndGradient()!", getClass()->getName().c_str());
}
void VolumeDataSource::gradientAndHessian(const PointF &p, VectorF &v, Matrix3x3F &M) const{
	Log(EError, "'%s': does not implement gradientAndHessian()!", getClass()->getName().c_str());
}
void VolumeDataSource::valueGradientAndHessian(const PointF &p, FLOAT &f, VectorF &v, Matrix3x3F &M) const{
	Log(EError, "'%s': does not implement valueGradientAndHessian()!", getClass()->getName().c_str());
}

bool VolumeDataSource::isAcousticRIF() const{
	Log(EError, "'%s': does not implement isAcousticRIF()!", getClass()->getName().c_str());
	return false;
}

Float VolumeDataSource::lookupFloat(const Point &p) const {
	Log(EError, "'%s': does not implement lookupFloat()!", getClass()->getName().c_str());
	return 0;
}

Spectrum VolumeDataSource::lookupSpectrum(const Point &p) const {
	Log(EError, "'%s': does not implement lookupSpectrum()!", getClass()->getName().c_str());
	return Spectrum(0.0f);
}

Vector VolumeDataSource::lookupVector(const Point &p) const {
	Log(EError, "'%s': does not implement lookupVector()!", getClass()->getName().c_str());
	return Vector();
}

bool VolumeDataSource::supportsFloatLookups() const {
	return false;
}

bool VolumeDataSource::supportsSpectrumLookups() const {
	return false;
}

bool VolumeDataSource::supportsVectorLookups() const {
	return false;
}

MTS_IMPLEMENT_CLASS(VolumeDataSource, true, ConfigurableObject)
MTS_NAMESPACE_END

