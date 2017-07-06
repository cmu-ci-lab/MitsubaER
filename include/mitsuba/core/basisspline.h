/*
 * Spline.h
 * Interpolate 1D, 2D, 3D with basis splines
 *
 *  Created on: Dec 16, 2019
 *      Author: apedired
 */

#pragma once
#if !defined(__BASISSPLINE_H_)
#define __BASISSPLINE_H_

#include <iostream>

#include <mitsuba/core/transform.h>
#include <mitsuba/core/platform.h>
#include <mitsuba/core/constants.h>
#include <mitsuba/core/matrix.h>
#include <math.h>

MTS_NAMESPACE_BEGIN

namespace basisspline{

inline int modulo(int a, int b) {
    int r = a % b;
    return (r < 0) ? r+b : r;
}

template <typename T> int sgn(T val) {
    return (T(0) < val) - (val < T(0));
}

namespace detail{

template<int D>
inline FLOAT kernel(FLOAT x);

template<>
inline FLOAT kernel<0>(FLOAT x){
	x = std::abs(x);
	if(x > 2)
		return (FLOAT)(0.0);
	if(x > 1)
		return (ONE_SIXTH*(2-x)*(2-x)*(2-x));
	return (TWO_THIRD - x*x + VHALF*x*x*x);
//   	FLOAT y = 0;
//	if(x > 2 || x < -2)
//		return y;
//	if(x > -2){
//		y += (x+2)*(x+2)*(x+2);
//		if(x > -1){
//			y += -4*(x+1)*(x+1)*(x+1);
//			if(x > 0){
//				y += 6*x*x*x;
//				if(x > 1){
//					y += -4*(x-1)*(x-1)*(x-1);
//				}
//			}
//		}
//	}
//	return ONE_SIXTH*y;
}
template<>
inline FLOAT kernel<1>(FLOAT x){
	int s = sgn(x);
	x = std::abs(x);
	if(x > 2)
		return (FLOAT)(0.0);
	if(x > 1)
		return s*(-VHALF*(2-x)*(2-x));
	return s*((1.5*x - 2)*x);
//   	FLOAT y = 0;
//	if(x > 2 || x < -2)
//		return y;
//	if(x > -2){
//		y += (x+2)*(x+2);
//		if(x > -1){
//			y += -4*(x+1)*(x+1);
//			if(x > 0){
//				y += 6*x*x;
//				if(x > 1){
//					y += -4*(x-1)*(x-1);
//				}
//			}
//		}
//	}
//	return VHALF*y;
}
template<>
inline FLOAT kernel<2>(FLOAT x){
	x = std::abs(x);
	if(x > 2)
		return (FLOAT)(0.0);
	if(x > 1)
		return 2-x;
	return 3*x-2;
//   	FLOAT y = 0;
//	if(x > 2 || x < -2)
//		return y;
//	if(x > -2){
//		y += (x+2);
//		if(x > -1){
//			y += -4*(x+1);
//			if(x > 0){
//				y += 6*x;
//				if(x > 1){
//					y += -4*(x-1);
//				}
//			}
//		}
//	}
//	return y;
}

}

template<int DIM>
class Spline{
public:

	Spline(){ // empty constructor, which waits for a proper initialization. HACK: FIXME
	}
    inline void initialize(const FLOAT xmin[DIM], const FLOAT xmax[DIM],  const int N[DIM]){ // to compensate for above HACK. FIXME;
        uint datasize = 1;
        for(int i=0; i < DIM; i++){
            this->xmin[i] = xmin[i];
            this->xmax[i] = xmax[i];
            this->N[i] = N[i];
            xres[i] = (N[i]-1)/(xmax[i] - xmin[i]);
            dxres[i] = xres[i];
            dxres2[i] = dxres[i] * dxres[i];
            datasize= datasize*N[i];
        }
        coeff = new FLOAT[datasize];
        built = false; // Note: currently not testing when the value is asked as that would make code slow
        z1 = -2 + std::sqrt(3);
    }

    Spline(const FLOAT xmin[DIM], const FLOAT xmax[DIM],  const int N[DIM]){
        uint datasize = 1;
        for(int i=0; i < DIM; i++){
            this->xmin[i] = xmin[i];
            this->xmax[i] = xmax[i];
            this->N[i] = N[i];
            xres[i] = (N[i]-1)/(xmax[i] - xmin[i]);
            datasize= datasize*N[i];
        }
        coeff = new FLOAT[datasize];
        built = false; // Note: currently not testing when the value is asked as that would make code slow
        z1 = -2 + std::sqrt(3);
        
    }
    ~Spline(){
        delete[] coeff;
    }

    inline void printcoeff3d() const{
    	for(int k=0; k<N[2]; k++){
        	for(int j=0; j<N[1]; j++){
            	for(int i=0; i<N[0]; i++){
                	std::cout << coeff[i + j*N[0] + k*N[0]*N[1]] << ", ";
            	}
            	std::cout << std::endl;
        	}
        	std::cout << std::endl;
    	}
    }


    inline void build(const FLOAT data[]){
        if(DIM == 1){
            build1d(data, 0, 1, N[0], coeff);
            built = true;
            return;
        }
        if(DIM == 2){
            build2d(data);
            built = true;
            return;
        }
        if(DIM == 3){
            build3d(data);
            built = true;
            return;
        }
	    std::cerr << "Error: Current implementation of spline interpolation works only upto 3 dimensions \n";
        exit (EXIT_FAILURE);
    }

    template<int DX>
    inline FLOAT value(const FLOAT x[]) const{
        if(DIM != 1){
	        std::cerr << "Error: 1D-Value called for " << DIM << "dimensions \n";
            exit (EXIT_FAILURE);
        }

        FLOAT y[DIM]; // make a duplicate
        y[0] = x[0];
        return value1d<DX>(y);
    }

    template<int DX, int DY>
    inline FLOAT value(const FLOAT x[]) const{
        if(DIM != 2){
	        std::cerr << "Error: 2D-Value called for " << DIM << "dimensions \n";
            exit (EXIT_FAILURE);
        }

        FLOAT y[DIM]; // make a duplicate
        y[0] = x[0]; y[1] = x[1];
        return value2d<DX, DY>(y);
    }

    template<int DX, int DY, int DZ>
    inline FLOAT value(const FLOAT x[]) const{
        if(DIM != 3){
	        std::cerr << "Error: 3D-Value called for " << DIM << "dimensions \n";
            exit (EXIT_FAILURE);
        }

        FLOAT y[DIM]; // make a duplicate
        y[0] = x[0]; y[1] = x[1]; y[2] = x[2]; 
        return value3d<DX, DY, DZ>(y);
    }

    // special functions that work only inside ERCRDR. comment outside
    inline Matrix3x3F hessian2d(const FLOAT y[]) const{
        if(DIM != 2){
	        std::cerr << "Error: 2D-Value called for " << DIM << "dimensions \n";
            exit (EXIT_FAILURE);
        }

        FLOAT x[DIM]; // make a duplicate
        x[0] = y[0]; x[1] = y[1];

        FLOAT Hxx = 0;
        FLOAT Hxy = 0;
        FLOAT Hyy = 0;

        convertToX(x);

        int wrap_index1, wrap_index2;
        for(int index1 = ceil(x[0]-2); index1 <= floor(x[0]+2); index1++){
            wrap_index1 = modulo(index1, 2*N[0]-2);
            if(wrap_index1 >= N[0]){
                wrap_index1 = 2*N[0] - 2 - wrap_index1;
            }
            for(int index2 = ceil(x[1]-2); index2 <= floor(x[1]+2); index2++){
                wrap_index2 = modulo(index2, 2*N[1]-2);
                if(wrap_index2 >= N[1]){
                    wrap_index2 = 2*N[1] - 2 - wrap_index2;
                }
				Hxx += coeff[wrap_index1+wrap_index2*N[0]] * kernel<2>(x[0]-index1) * kernel<0>(x[1]-index2);
				Hxy += coeff[wrap_index1+wrap_index2*N[0]] * kernel<1>(x[0]-index1) * kernel<1>(x[1]-index2);
				Hyy += coeff[wrap_index1+wrap_index2*N[0]] * kernel<0>(x[0]-index1) * kernel<2>(x[1]-index2);
            }
        }

        Hxx *= dxres2[0];
        Hyy *= dxres2[1];
        Hxy *= dxres[0]*dxres[1];

        return Matrix3x3F(0, 0,   0,
        				 0, Hyy, Hxy,
    					 0, Hxy, Hxx);
    }

    inline VectorF gradient2d(const FLOAT y[]) const{
        if(DIM != 2){
	        std::cerr << "Error: 2D-Value called for " << DIM << "dimensions \n";
            exit (EXIT_FAILURE);
        }

        FLOAT x[DIM]; // make a duplicate
        x[0] = y[0]; x[1] = y[1];
        VectorF v(0);

        convertToX(x);

        int wrap_index1, wrap_index2;
        for(int index1 = ceil(x[0]-2); index1 <= floor(x[0]+2); index1++){
            wrap_index1 = modulo(index1, 2*N[0]-2);
            if(wrap_index1 >= N[0]){
                wrap_index1 = 2*N[0] - 2 - wrap_index1;
            }
            for(int index2 = ceil(x[1]-2); index2 <= floor(x[1]+2); index2++){
                wrap_index2 = modulo(index2, 2*N[1]-2);
                if(wrap_index2 >= N[1]){
                    wrap_index2 = 2*N[1] - 2 - wrap_index2;
                }
				v.y += coeff[wrap_index1+wrap_index2*N[0]] * kernel<0>(x[0]-index1) * kernel<1>(x[1]-index2);
				v.z += coeff[wrap_index1+wrap_index2*N[0]] * kernel<1>(x[0]-index1) * kernel<0>(x[1]-index2);
            }
        }
		v.y *= dxres[1];
		v.z *= dxres[0];
		return v;
    }


    inline FLOAT value(const FLOAT x1[]) const{
        FLOAT x[DIM]; // make a duplicate
        x[0] = x1[0]; x[1] = x1[1]; x[2] = x1[2];
        convertToX(x);
        FLOAT v = 0;
        for(int index1 = ceil(x[0]-2); index1 <= floor(x[0]+2); index1++){
            for(int index2 = ceil(x[1]-2); index2 <= floor(x[1]+2); index2++){
                for(int index3 = ceil(x[2]-2); index3 <= floor(x[2]+2); index3++){
                    v += coeff[index1+index2*N[0]+index3*N[0]*N[1]] * kernel<0>(x[0]-index1) * kernel<0>(x[1]-index2) * kernel<0>(x[2]-index3);
                }
            }
        }
        return v;
    }


    inline VectorF gradient(FLOAT x[]) const{
        if(DIM != 3){
	        std::cerr << "Error: 3D-Value called for " << DIM << "dimensions \n";
            exit (EXIT_FAILURE);
        }

        VectorF v(0.0);

        convertToX(x);

        FLOAT precomputeK0x;
        FLOAT precomputeK0y;
        FLOAT precomputeK0z;
        FLOAT precomputeCoeff;

//        int wrap_index1, wrap_index2, wrap_index3; // hack as usage never warrants this case
        for(int index1 = ceil(x[0]-2); index1 <= floor(x[0]+2); index1++){
//            wrap_index1 = modulo(index1, 2*N[0]-2);
//            if(wrap_index1 >= N[0]){
//                wrap_index1 = 2*N[0] - 2 - wrap_index1;
//            }
            for(int index2 = ceil(x[1]-2); index2 <= floor(x[1]+2); index2++){
//                wrap_index2 = modulo(index2, 2*N[1]-2);
//                if(wrap_index2 >= N[1]){
//                    wrap_index2 = 2*N[1] - 2 - wrap_index2;
//                }
                for(int index3 = ceil(x[2]-2); index3 <= floor(x[2]+2); index3++){
//                    wrap_index3 = modulo(index3, 2*N[2]-2);
//                    if(wrap_index3 >= N[2]){
//                        wrap_index3 = 2*N[2] - 2 - wrap_index3;
//                    }
//                	precomputeCoeff = coeff[wrap_index1+wrap_index2*N[0]+wrap_index3*N[0]*N[1]];
                    precomputeCoeff = coeff[index1+index2*N[0]+index3*N[0]*N[1]];
                    precomputeK0x   = kernel<0>(x[0]-index1);
                    precomputeK0y   = kernel<0>(x[1]-index2);
                    precomputeK0z   = kernel<0>(x[2]-index3);
                    v.x += precomputeCoeff * kernel<1>(x[0]-index1) * precomputeK0y * precomputeK0z;
                    v.y += precomputeCoeff * precomputeK0x * kernel<1>(x[1]-index2) * precomputeK0z;
                    v.z += precomputeCoeff * precomputeK0x * precomputeK0y * kernel<1>(x[2]-index3);
                }
            }
        }
		v.x *= dxres[0];
		v.y *= dxres[1];
		v.z *= dxres[2];
		return v;
    }


    inline Matrix3x3F hessian(FLOAT x[]) const{

        FLOAT Hxx = 0;
        FLOAT Hyy = 0;
        FLOAT Hzz = 0;
        FLOAT Hxy = 0;
        FLOAT Hyz = 0;
        FLOAT Hzx = 0;

        convertToX(x);

        FLOAT precomputeK0x;
        FLOAT precomputeK0y;
        FLOAT precomputeK0z;
        FLOAT precomputeK1x;
        FLOAT precomputeK1y;
        FLOAT precomputeK1z;
        FLOAT precomputeCoeff;


//        int wrap_index1, wrap_index2, wrap_index3; // hack as usage never warrants this case
        for(int index1 = ceil(x[0]-2); index1 <= floor(x[0]+2); index1++){
//            wrap_index1 = modulo(index1, 2*N[0]-2);
//            if(wrap_index1 >= N[0]){
//                wrap_index1 = 2*N[0] - 2 - wrap_index1;
//            }
            for(int index2 = ceil(x[1]-2); index2 <= floor(x[1]+2); index2++){
//                wrap_index2 = modulo(index2, 2*N[1]-2);
//                if(wrap_index2 >= N[1]){
//                    wrap_index2 = 2*N[1] - 2 - wrap_index2;
//                }
                for(int index3 = ceil(x[2]-2); index3 <= floor(x[2]+2); index3++){
//                    wrap_index3 = modulo(index3, 2*N[2]-2);
//                    if(wrap_index3 >= N[2]){
//                        wrap_index3 = 2*N[2] - 2 - wrap_index3;
//                    }

//                    precomputeCoeff = coeff[wrap_index1+wrap_index2*N[0]+wrap_index3*N[0]*N[1]];
                    precomputeCoeff = coeff[index1+index2*N[0]+index3*N[0]*N[1]];
                    precomputeK0x   = kernel<0>(x[0]-index1);
                    precomputeK0y   = kernel<0>(x[1]-index2);
                    precomputeK0z   = kernel<0>(x[2]-index3);

                    precomputeK1x   = kernel<1>(x[0]-index1);
                    precomputeK1y   = kernel<1>(x[1]-index2);
                    precomputeK1z   = kernel<1>(x[2]-index3);

					Hxx += precomputeCoeff * kernel<2>(x[0]-index1) * precomputeK0y * precomputeK0z;
					Hyy += precomputeCoeff * precomputeK0x * kernel<2>(x[1]-index2) * precomputeK0z;
					Hzz += precomputeCoeff * precomputeK0x * precomputeK0y * kernel<2>(x[2]-index3);

					Hxy += precomputeCoeff * precomputeK1x * precomputeK1y * precomputeK0z;
					Hyz += precomputeCoeff * precomputeK0x * precomputeK1y * precomputeK1z;
					Hzx += precomputeCoeff * precomputeK1x * precomputeK0y * precomputeK1z;
                }
            }
        }

        Hxx *= dxres2[0];
        Hyy *= dxres2[1];
        Hzz *= dxres2[2];

        Hxy *= dxres[0]*dxres[1];
        Hyz *= dxres[1]*dxres[2];
        Hzx *= dxres[2]*dxres[0];

        return Matrix3x3F(Hxx, Hxy, Hzx,
        				 Hxy, Hyy, Hyz,
						 Hzx, Hyz, Hzz);
    }

	inline void valueAndGradient(FLOAT x[], FLOAT &f, VectorF &v) const{

        v.x = v.y = v.z = 0;
        f = 0;

        convertToX(x);

        FLOAT precomputeK0x;
        FLOAT precomputeK0y;
        FLOAT precomputeK0z;
        FLOAT precomputeCoeff;


        for(int index1 = ceil(x[0]-2); index1 <= floor(x[0]+2); index1++){
            for(int index2 = ceil(x[1]-2); index2 <= floor(x[1]+2); index2++){
                for(int index3 = ceil(x[2]-2); index3 <= floor(x[2]+2); index3++){
                    precomputeCoeff = coeff[index1+index2*N[0]+index3*N[0]*N[1]];
                    precomputeK0x   = kernel<0>(x[0]-index1);
                    precomputeK0y   = kernel<0>(x[1]-index2);
                    precomputeK0z   = kernel<0>(x[2]-index3);

                    f   += precomputeCoeff * precomputeK0x * precomputeK0y * precomputeK0z;

                    v.x += precomputeCoeff * kernel<1>(x[0]-index1) * precomputeK0y * precomputeK0z;
                    v.y += precomputeCoeff * precomputeK0x * kernel<1>(x[1]-index2) * precomputeK0z;
                    v.z += precomputeCoeff * precomputeK0x * precomputeK0y * kernel<1>(x[2]-index3);
                }
            }
        }

		v.x *= dxres[0];
		v.y *= dxres[1];
		v.z *= dxres[2];
    }

	inline void gradientAndHessian(FLOAT x[], VectorF &v, Matrix3x3F &H) const{

        v.x = v.y = v.z = 0;

        FLOAT Hxx = 0;
        FLOAT Hyy = 0;
        FLOAT Hzz = 0;
        FLOAT Hxy = 0;
        FLOAT Hyz = 0;
        FLOAT Hzx = 0;

        convertToX(x);

        FLOAT precomputeK0x;
        FLOAT precomputeK0y;
        FLOAT precomputeK0z;
        FLOAT precomputeK1x;
        FLOAT precomputeK1y;
        FLOAT precomputeK1z;
        FLOAT precomputeCoeff;


        for(int index1 = ceil(x[0]-2); index1 <= floor(x[0]+2); index1++){
            for(int index2 = ceil(x[1]-2); index2 <= floor(x[1]+2); index2++){
                for(int index3 = ceil(x[2]-2); index3 <= floor(x[2]+2); index3++){
                    precomputeCoeff = coeff[index1+index2*N[0]+index3*N[0]*N[1]];
                    precomputeK0x   = kernel<0>(x[0]-index1);
                    precomputeK0y   = kernel<0>(x[1]-index2);
                    precomputeK0z   = kernel<0>(x[2]-index3);

                    precomputeK1x   = kernel<1>(x[0]-index1);
                    precomputeK1y   = kernel<1>(x[1]-index2);
                    precomputeK1z   = kernel<1>(x[2]-index3);

					Hxx += precomputeCoeff * kernel<2>(x[0]-index1) * precomputeK0y * precomputeK0z;
					Hyy += precomputeCoeff * precomputeK0x * kernel<2>(x[1]-index2) * precomputeK0z;
					Hzz += precomputeCoeff * precomputeK0x * precomputeK0y * kernel<2>(x[2]-index3);

					Hxy += precomputeCoeff * precomputeK1x * precomputeK1y * precomputeK0z;
					Hyz += precomputeCoeff * precomputeK0x * precomputeK1y * precomputeK1z;
					Hzx += precomputeCoeff * precomputeK1x * precomputeK0y * precomputeK1z;

                    v.x += precomputeCoeff * precomputeK1x * precomputeK0y * precomputeK0z;
                    v.y += precomputeCoeff * precomputeK0x * precomputeK1y * precomputeK0z;
                    v.z += precomputeCoeff * precomputeK0x * precomputeK0y * precomputeK1z;
                }
            }
        }

		v.x *= dxres[0];
		v.y *= dxres[1];
		v.z *= dxres[2];

        Hxx *= dxres2[0];
        Hyy *= dxres2[1];
        Hzz *= dxres2[2];

        Hxy *= dxres[0]*dxres[1];
        Hyz *= dxres[1]*dxres[2];
        Hzx *= dxres[2]*dxres[0];

        H = Matrix3x3F(Hxx, Hxy, Hzx,
        			   Hxy, Hyy, Hyz,
					   Hzx, Hyz, Hzz);
    }

	inline void valueGradientAndHessian(FLOAT x[], FLOAT &f, VectorF &v, Matrix3x3F &H) const{

		v.x = v.y = v.z = 0;
        f = 0;

        FLOAT Hxx = 0;
        FLOAT Hyy = 0;
        FLOAT Hzz = 0;
        FLOAT Hxy = 0;
        FLOAT Hyz = 0;
        FLOAT Hzx = 0;

        convertToX(x);

        FLOAT precomputeK0x;
        FLOAT precomputeK0y;
        FLOAT precomputeK0z;
        FLOAT precomputeK1x;
        FLOAT precomputeK1y;
        FLOAT precomputeK1z;
        FLOAT precomputeCoeff;


        for(int index1 = ceil(x[0]-2); index1 <= floor(x[0]+2); index1++){
            for(int index2 = ceil(x[1]-2); index2 <= floor(x[1]+2); index2++){
                for(int index3 = ceil(x[2]-2); index3 <= floor(x[2]+2); index3++){
                    precomputeCoeff = coeff[index1+index2*N[0]+index3*N[0]*N[1]];
                    precomputeK0x   = kernel<0>(x[0]-index1);
                    precomputeK0y   = kernel<0>(x[1]-index2);
                    precomputeK0z   = kernel<0>(x[2]-index3);

                    precomputeK1x   = kernel<1>(x[0]-index1);
                    precomputeK1y   = kernel<1>(x[1]-index2);
                    precomputeK1z   = kernel<1>(x[2]-index3);

                    f += precomputeCoeff * precomputeK0x * precomputeK0y * precomputeK0z;

					Hxx += precomputeCoeff * kernel<2>(x[0]-index1) * precomputeK0y * precomputeK0z;
					Hyy += precomputeCoeff * precomputeK0x * kernel<2>(x[1]-index2) * precomputeK0z;
					Hzz += precomputeCoeff * precomputeK0x * precomputeK0y * kernel<2>(x[2]-index3);

					Hxy += precomputeCoeff * precomputeK1x * precomputeK1y * precomputeK0z;
					Hyz += precomputeCoeff * precomputeK0x * precomputeK1y * precomputeK1z;
					Hzx += precomputeCoeff * precomputeK1x * precomputeK0y * precomputeK1z;

                    v.x += precomputeCoeff * precomputeK1x * precomputeK0y * precomputeK0z;
                    v.y += precomputeCoeff * precomputeK0x * precomputeK1y * precomputeK0z;
                    v.z += precomputeCoeff * precomputeK0x * precomputeK0y * precomputeK1z;
                }
            }
        }

		v.x *= dxres[0];
		v.y *= dxres[1];
		v.z *= dxres[2];

        Hxx *= dxres2[0];
        Hyy *= dxres2[1];
        Hzz *= dxres2[2];

        Hxy *= dxres[0]*dxres[1];
        Hyz *= dxres[1]*dxres[2];
        Hzx *= dxres[2]*dxres[0];

        H = Matrix3x3F(Hxx, Hxy, Hzx,
        			   Hxy, Hyy, Hyz,
					   Hzx, Hyz, Hzz);
    }

//	inline void transform(const Transform &worldToVolume){
//        if(DIM != 3){
//	        std::cerr << "Error: transform is defined only for " << DIM << "dimensions \n";
//            exit (EXIT_FAILURE);
//        }
//        PointF p(xres[0], xres[1], xres[2]);
//        p = worldToVolume(p);
////        volumeToWorld(p);
//        for(int i=0; i<DIM; i++){
//        	dxres[i] = p[i];
//        	dxres2[i] = dxres[i]*dxres[i];
//        }
//	}

	inline FLOAT getStride(int dim) const{
		return 1.0/xres[dim];
	}

    /*    template<int DX, int DY, int DZ>
    inline FLOAT value(const FLOAT x[]) const{
        FLOAT y[DIM]; // make a duplicate
        for(int i=0; i<DIM; i++){
            y[i] = x[i];
        }
        if(DIM == 1)
            return value1d<DX>(y);
        if(DIM == 2)
            return value2d<DX, DY>(y);
        if(DIM == 3)
            return value3d<DX, DY, DZ>(y);
	    std::cerr << "Error: Current implementation of spline interpolation works only upto 3 dimensions \n";
        exit (EXIT_FAILURE);
    }
*/
private:
//public:
    FLOAT xmin[DIM];
    FLOAT xmax[DIM];
    FLOAT xres[DIM];
    FLOAT dxres[DIM];
    FLOAT dxres2[DIM];
    FLOAT *coeff;
    int N[DIM];
    bool  built;
    FLOAT z1;


    inline void convertToX(FLOAT x[DIM]) const{
        for(int i=0; i<DIM; i++)
            x[i] = (x[i] - xmin[i]) * xres[i] ;
    }


    template<int D>
    inline FLOAT kernel(FLOAT x) const{
    	return detail::kernel<D>(x);
    }

//
//    template<int D>
//    inline FLOAT kernel(FLOAT x) const{
//    	FLOAT y = 0;
//	    if(D == 0){
//	        if(x > 2 || x < -2)
//	            return y;
//	        if(x > -2){
//	            y += (x+2)*(x+2)*(x+2);
//	            if(x > -1){
//	                y += -4*(x+1)*(x+1)*(x+1);
//	                if(x > 0){
//	                    y += 6*x*x*x;
//	                    if(x > 1){
//	                        y += -4*(x-1)*(x-1)*(x-1);
//	                    }
//	                }
//	            }
//	        }
//	        return y/6;
//	    }else if(D == 1){
//	        if(x > 2 || x < -2)
//	            return y;
//	        if(x > -2){
//	            y += (x+2)*(x+2);
//	            if(x > -1){
//	                y += -4*(x+1)*(x+1);
//	                if(x > 0){
//	                    y += 6*x*x;
//	                    if(x > 1){
//	                        y += -4*(x-1)*(x-1);
//	                    }
//	                }
//	            }
//	        }
//	        return y/2;
//	    }else if(D == 2){
//	        if(x > 2 || x < -2)
//	            return y;
//	        if(x > -2){
//	            y += (x+2);
//	            if(x > -1){
//	                y += -4*(x+1);
//	                if(x > 0){
//	                    y += 6*x;
//	                    if(x > 1){
//	                        y += -4*(x-1);
//	                    }
//	                }
//	            }
//	        }
//	        return y;
//	    }else{
//	        std::cerr << "Error: More than double derivative is not defined \n";
//	    }
//    }

    template<int DX>
    inline FLOAT value1d(FLOAT x[]) const{
        convertToX(x);
        FLOAT v = 0;
        int wrap_index;
        for(int index = ceil(x[0]-2); index <= floor(x[0]+2); index++){
            wrap_index = modulo(index, 2*N[0]-2);
            if(wrap_index >= N[0]){
                wrap_index = 2*N[0] - 2 - wrap_index;
            }
            v += coeff[wrap_index] * kernel<DX>(x[0]-index);
        }
        if(DX == 1)
            v *= dxres[0];
        if(DX == 2)
            v *= dxres2[0];
        return v;
    }

    template<int DX, int DY>
    inline FLOAT value2d(FLOAT x[]) const{
        convertToX(x);
        FLOAT v = 0;
        int wrap_index1, wrap_index2;
        for(int index1 = ceil(x[0]-2); index1 <= floor(x[0]+2); index1++){
            wrap_index1 = modulo(index1, 2*N[0]-2);
            if(wrap_index1 >= N[0]){
                wrap_index1 = 2*N[0] - 2 - wrap_index1;
            }
            for(int index2 = ceil(x[1]-2); index2 <= floor(x[1]+2); index2++){
                wrap_index2 = modulo(index2, 2*N[1]-2);
                if(wrap_index2 >= N[1]){
                    wrap_index2 = 2*N[1] - 2 - wrap_index2;
                }
                v += coeff[wrap_index1+wrap_index2*N[0]] * kernel<DX>(x[0]-index1) * kernel<DY>(x[1]-index2);
            }
        }
        if(DX == 1)
            v *= dxres[0];
        if(DX == 2)
            v *= dxres2[0];
        if(DY == 1)
            v *= dxres[1];
        if(DY == 2)
            v *= dxres2[1];
        return v;
    }

    template<int DX, int DY, int DZ>
    inline FLOAT value3d(FLOAT x[]) const{
        convertToX(x);
        FLOAT v = 0;
//        int wrap_index1, wrap_index2, wrap_index3;
        for(int index1 = ceil(x[0]-2); index1 <= floor(x[0]+2); index1++){
//            wrap_index1 = modulo(index1, 2*N[0]-2);
//            if(wrap_index1 >= N[0]){
//                wrap_index1 = 2*N[0] - 2 - wrap_index1;
//            }
            for(int index2 = ceil(x[1]-2); index2 <= floor(x[1]+2); index2++){
//                wrap_index2 = modulo(index2, 2*N[1]-2);
//                if(wrap_index2 >= N[1]){
//                    wrap_index2 = 2*N[1] - 2 - wrap_index2;
//                }
                for(int index3 = ceil(x[2]-2); index3 <= floor(x[2]+2); index3++){
//                    wrap_index3 = modulo(index3, 2*N[2]-2);
//                    if(wrap_index3 >= N[2]){
//                        wrap_index3 = 2*N[2] - 2 - wrap_index3;
//                    }
//                    v += coeff[wrap_index1+wrap_index2*N[0]+wrap_index3*N[0]*N[1]] * kernel<DX>(x[0]-index1) * kernel<DY>(x[1]-index2) * kernel<DZ>(x[2]-index3);
                    v += coeff[index1+index2*N[0]+index3*N[0]*N[1]] * kernel<DX>(x[0]-index1) * kernel<DY>(x[1]-index2) * kernel<DZ>(x[2]-index3);
                }
            }
        }
        if(DX == 1)
            v *= dxres[0];
        if(DX == 2)
            v *= dxres2[0];
        if(DY == 1)
            v *= dxres[1];
        if(DY == 2)
            v *= dxres2[1];
        if(DZ == 1)
            v *= dxres[2];
        if(DZ == 2)
            v *= dxres2[2];
        return v;
    }

    // build 1d is called from build2d and build3d so it is written in a generic style without affecting the contents of the class
    inline void build1d(const FLOAT data[], const int &offset, const int &stride, const int &size, FLOAT out[]) const{
        FLOAT *cp = new FLOAT[size];
        FLOAT *cn = new FLOAT[size];
        cp[0] = 0;

        for(int i=0; i<size; i++)
            cp[0] += data[offset+i*stride] * pow(z1, i);

        for(int i=size-2; i>0; i--)
            cp[0] += data[offset+i*stride] * pow(z1, 2*size-2-i);

        cp[0] /= (1-pow(z1, 2*size-2));

        for(int i=1; i<size; i++)
            cp[i] = data[offset + i*stride] + z1 * cp[i-1];


        cn[size-1] = z1/(z1*z1-1) * (cp[size-1] + z1*cp[size-2]); 

        for(int i=size-2; i>=0; i--)
            cn[i] = z1*(cn[i+1] - cp[i]);

        for(int i=0; i<size; i++)
            out[i] = 6*cn[i];


        delete[] cp;
        delete[] cn;
    }

    inline void build2d(const FLOAT data[]){
        const int size = std::max(N[0], N[1]);
        FLOAT *temp = new FLOAT[size];

        // Build along rows
        for(int i=0; i<N[0]; i++){
            build1d(data, i, N[0], N[1], temp);
            for(int j=0; j<N[1]; j++){
                coeff[j*N[0]+i] = temp[j]; 
            }
        }

        // Build along columns
        for(int i=0; i<N[1]; i++){
            build1d(coeff, i*N[0], 1, N[0], temp);
            for(int j=0; j<N[0]; j++){
                coeff[i*N[0]+j] = temp[j]; 
            }
        }

        delete[] temp;
    }

    inline void build3d(const FLOAT data[]){
        FLOAT *temp = new FLOAT[std::max(std::max(N[0], N[1]), N[2])];

        for(int k=0; k<N[2]; k++)
            for(int i=0; i<N[0]; i++){
                build1d(data, k*N[0]*N[1] + i, N[0], N[1], temp);
                for(int t=0; t < N[1]; t++)
                    coeff[i+t*N[0]+k*N[0]*N[1]] = temp[t];
            }

        for(int k=0; k<N[2]; k++)
            for(int j=0; j<N[1]; j++){
                build1d(coeff, k*N[0]*N[1] + j*N[0], 1, N[0], temp);
                for(int t=0; t < N[0]; t++)
                    coeff[t+j*N[0]+k*N[0]*N[1]] = temp[t];
            }

        for(int i=0; i<N[0]; i++)
            for(int j=0; j<N[1]; j++){
                build1d(coeff, j*N[0] + i, N[0]*N[1], N[2], temp);
                for(int t=0; t < N[2]; t++)
                    coeff[i+j*N[0]+t*N[0]*N[1]] = temp[t];
            }

        delete[] temp;
    }
};

}	/* namespace spline */

MTS_NAMESPACE_END

#endif /* SPLINE_H_ */


