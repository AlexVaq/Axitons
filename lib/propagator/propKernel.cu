#include "enum-vars.h"
#include "cudaErrors.h"

#define	TwoPi	(2.*M_PI)

template<typename Float>
static __device__ __forceinline__ Float modPi (const Float x, const Float OneOvPi, const Float TwoPiZ)
{
	const Float tmp = x*OneOvPi;

	if (tmp >=  1.)
		return (x-TwoPiZ);

	if (tmp <  -1.)
		return (x+TwoPiZ);

	return x;
}

template<typename Float, const int nNeig>
static __device__ __forceinline__ constexpr Float C(int nIdx) {

	switch (nNeig) {

		default:
		case 1: {
			return	((Float) 1.0);
			break;
		}

		case 2: {
			constexpr Float C[2] = { 4.0/3.0, -1.0/12.0 };
			return	C[nIdx];
			break;
		}

		case 3: {
			constexpr Float C[3] = { 1.5,     -3.0/20.0,  1.0/90.0 };
			return	C[nIdx];
			break;
		}

		case 4: {
			constexpr Float C[4] = { 8.0/5.0, -1.0/5.0,   8.0/315.0, -1.0/560.0 };
			return	C[nIdx];
			break;
		}
	}
}


template<typename Float, const int nNeig, const bool wMod>
static __device__ __forceinline__ void	propagateCoreGpu(const uint idx, const Float * __restrict__ field, Float * __restrict__ dev, Float * __restrict__ misc, const Float gm, const Float zQ,
							 const Float iz, const Float dzc, const Float dzd, const Float ood2, const int Lx, const Float zP, const Float tPz)
{
	Float mel = 0.0, pPc, a, f0n;
	f0n = field[idx];

	if (idx != 0) {
		pPc = 1.0/((Float) idx);
		if (idx < Lx - nNeig){
			#pragma unroll
			for (int nIdx=1; nIdx<=nNeig; nIdx++)
			{
				auto rIdx  = __sad (idx, nIdx, 0);
				// mel       += (field[idx+nIdx]*(1.0 + nIdx*pPc) + field[rIdx]*((Float) rIdx)*pPc - 2.0*f0n)*C<Float,nNeig>(nIdx-1);
				mel       += (field[idx+nIdx]*(1.0 + nIdx*pPc) + field[rIdx]*(1.0 - nIdx*pPc) - 2.0*f0n)*C<Float,nNeig>(nIdx-1);
			}
		} else {
			// mel += (field[idx-2]-f0n)*(1.0 + pPc) - 2.0*(field[idx-1]-f0n)*(1 + 2.0*pPc);
			mel +=  -(field[idx-3]-f0n)*(1.0 + 2.0*pPc/3.0)  + (field[idx-2]-f0n)*(4.0 + 3.0*pPc) - (field[idx-1]-f0n)*(5.0 + 6.0*pPc);
		}
	} else {
		#pragma unroll
		for (int nIdx=1; nIdx<=nNeig; nIdx++) {
			mel += (field[nIdx] - f0n)*2.0*C<Float,nNeig>(nIdx-1);
		}
	}

	// if (idx != 0) { // FIXME SLOW AS HELL
	// 	pPc = 1.0/((Float) idx);
	// 	#pragma unroll
	// 	for (int nIdx=1; nIdx<=nNeig; nIdx++)
	// 	{
	// 		auto rIdx  = __sad (idx, nIdx, 0);
	// 		mel       += (field[idx+nIdx]*(1.0 + nIdx*pPc) + field[rIdx]*((Float) rIdx)*pPc - 2.0*f0n)*C<Float,nNeig>(nIdx-1);
	// 	}
	// } else {
	// 	#pragma unroll
	// 	for (int nIdx=1; nIdx<=nNeig; nIdx++) {
	// 		mel += (field[nIdx] - f0n)*2.0*C<Float,nNeig>(nIdx-1);
	// 	}
	// }


/*
	if (idx > nNeig) { // FIXME SLOW AS HELL
		pPc = 1.0/cIdx;
		#pragma unroll
		for (int nIdx=1; nIdx<nNeig; nIdx++)
			mel += (field[idx+nIdx]*(((Float) nIdx) + pPc) + field[idx-nIdx]*(((Float) nIdx) - pPc))*C[nIdx] - f0n*C[nIdx];
	} else {
		for (int cIdx=0;     i<idx;   i++)
			mel += (field[cIdx+idx]*(((Float) nIdx) + pPc) + field[cIdx-idx]*(((Float) nIdx) - pPc))*C[cIdx] - f0n*C[cIdx];

		mel += (field[idx*2]*(((Float) idx) + pPc) + f0n*(((Float) idx) - pPc - 1))*C[idx];

		for (int cIdx=idx+1; i<nNeig; i++)
			mel += (field[idx+cIdx]*(((Float) cIdx) + pPc) + field[idx-cIdx]*(((Float) cIdx) - pPc))*C[cIdx] - f0n*C[cIdx];
	}
*/
	a = mel*ood2 - zQ*sin(f0n*iz);

	mel	 = dev[idx];

	if ((idx > Lx*0.8) && (gm > 0.0) ) {	// FIXME No hardcode this 0.8
		// variable measured with respect to axion mass?
		Float nGm = gm*(idx-0.8*Lx)/((Float) idx)*5;
		a         = (a - nGm*mel)/(1. + 0.5*nGm*dzc);
	}

	mel	+= a*dzc;
	dev[idx] = mel;
	mel	*= dzd;
	f0n	+= mel;

	misc[idx] = f0n; //modPi(tmp, zP, tPz);
}

template<typename Float, const int nNeig, const bool wMod>
__global__ void	propagateKernel(const Float * __restrict__ field, Float * __restrict__ dev, Float * __restrict__ misc, const Float gamma, const Float zQ, const Float dzc, const Float dzd,
				const Float ood2, const Float iz, const int Lx, const Float zP=0, const Float tPz=0)
{
	uint idx = threadIdx.x + blockDim.x*(blockIdx.x + gridDim.x*blockIdx.y);
	//uint idx = Vo + (threadIdx.x + blockDim.x*blockIdx.x) + Sf*(threadIdx.y + blockDim.y*blockIdx.y);

	// if	(idx >= Lx - nNeig)
	// 	return;
	if	(idx < Lx)
		propagateCoreGpu<Float,nNeig,wMod>(idx, field, dev, misc, gamma, zQ, iz, dzc, dzd, ood2, Lx, zP, tPz);
}

void	propGpu(const void * __restrict__ field, void * __restrict__ dev, void * __restrict__ misc, const double z, const double dz, const double c, const double d, const double ood2,
		const double aMass2, const int Lx, FieldPrecision precision, const int nNeig, const double gamma, const int xBlock, const int yBlock, const int zBlock)
{
	#define	BLSIZE 256
	dim3 gridSize((Lx+BLSIZE-1)/BLSIZE,1,1);
	dim3 blockSize(BLSIZE,1,1);

	if (precision == DoublePrecision)
	{
		const double dzc   = dz*c;
		const double dzd   = dz*d;
		const double zQ    = aMass2*z*z*z;//axionmass2((double) zR, nQcd, zthres, zrestore)*zR*zR*zR;
		const double iZ    = 1./z;
		const double gdzc2 = gamma*sqrt(aMass2)*z;

		switch (nNeig) {
			default:
			case 1:
			propagateKernel<double,1,false><<<gridSize,blockSize,0>>>((const double *) field, (double *) dev, (double *) misc, gdzc2, zQ, dzc, dzd, ood2, iZ, Lx);
			break;

			case 2:
			propagateKernel<double,2,false><<<gridSize,blockSize,0>>>((const double *) field, (double *) dev, (double *) misc, gdzc2, zQ, dzc, dzd, ood2, iZ, Lx);
			break;

			case 3:
			propagateKernel<double,3,false><<<gridSize,blockSize,0>>>((const double *) field, (double *) dev, (double *) misc, gdzc2, zQ, dzc, dzd, ood2, iZ, Lx);
			break;

			case 4:
			propagateKernel<double,4,false><<<gridSize,blockSize,0>>>((const double *) field, (double *) dev, (double *) misc, gdzc2, zQ, dzc, dzd, ood2, iZ, Lx);
			break;
		}
	}
	else if (precision == SinglePrecision)
	{
		const float dzc   = dz*c;
		const float dzd   = dz*d;
		const float zQ    = (float) (aMass2*z*z*z);//axionmass2((double) zR, nQcd, zthres, zrestore)*zR*zR*zR;
		const float iZ    = 1./z;
		const float gdzc2 = (float) gamma*sqrt(aMass2)*z;

		switch (nNeig) {
			default:
			case 1:
			propagateKernel<float, 1,false><<<gridSize,blockSize,0>>>((const float *) field, (float *) dev, (float *) misc, gdzc2, zQ, dzc, dzd, (float) ood2, iZ, Lx);
			break;

			case 2:
			propagateKernel<float, 2,false><<<gridSize,blockSize,0>>>((const float *) field, (float *) dev, (float *) misc, gdzc2, zQ, dzc, dzd, (float) ood2, iZ, Lx);
			break;

			case 3:
			propagateKernel<float, 3,false><<<gridSize,blockSize,0>>>((const float *) field, (float *) dev, (float *) misc, gdzc2, zQ, dzc, dzd, (float) ood2, iZ, Lx);
			break;

			case 4:
			propagateKernel<float, 4,false><<<gridSize,blockSize,0>>>((const float *) field, (float *) dev, (float *) misc, gdzc2, zQ, dzc, dzd, (float) ood2, iZ, Lx);
			break;
		}
	}
}
/*
void	propModGpu(const void * __restrict__ field, void * __restrict__ dev, void * __restrict__ misc, const double z, const double dz, const double c, const double d, const double ood2,
		   const double aMass2, const int Lx, FieldPrecision precision, const int xBlock, const int yBlock, const int zBlock)
{
	dim3 gridSize((Lx+BLSIZE-1)/BLSIZE,1,1);
	dim3 blockSize(BLSIZE,1,1);
	//dim3 gridSize((Sf+xBlock-1)/xBlock,(Lz2+yBlock-1)/yBlock,1);
	//dim3 blockSize(xBlock,yBlock,1);

	if (precision == DoublePrecision)
	{
		const double dzc  = dz*c;
		const double dzd  = dz*d;
		const double zQ   = aMass2*z*z*z;//xionmass2((double) zR, nQcd, zthres, zrestore)*zR*zR*zR;
		const double iZ   = 1./z;
		const double tPz  = 2.*M_PI*z;
		propagateKernel<double,true><<<gridSize,blockSize,0>>>((const double*) field, (double*) dev, (double*) misc, zQ, dzc, dzd, ood2, iZ, Lx, M_1_PI*iZ, tPz);
	}
	else if (precision == SinglePrecision)
	{
		const float dzc = dz*c;
		const float dzd = dz*d;
		const float zQ = (float) (aMass2*z*z*z);//axionmass2((double) zR, nQcd, zthres, zrestore)*zR*zR*zR;
		const float iZ   = 1./z;
		const float tPz  = 2.*M_PI*z;
		propagateKernel<float, true><<<gridSize,blockSize,0>>>((const float *) field, (float *) dev, (float *) misc, zQ, dzc, dzd, ood2, iZ, Lx, M_1_PI*iZ, tPz);
	}
}
*/
void	propGpu(const void * __restrict__ field, void * __restrict__ dev, void * __restrict__ misc, const double z, const double dz, const double c, const double d, const double ood2,
		const double aMass2, const int Lx, FieldPrecision precision, const int nNeig, const double gamma, const int xBlock, const int yBlock, const int zBlock, const FieldType wMod)
{
	switch (wMod) {

		case	FieldCompact:
			printf ("Compact propagator not implemented\n");
			//propModGpu(field, dev, misc, z, dz, c, d, ood2, aMass2, Lx, precision, gamma, xBlock, yBlock, zBlock);
			break;

		case	FieldNonCompact:
			propGpu	  (field, dev, misc, z, dz, c, d, ood2, aMass2, Lx, precision, nNeig, gamma, xBlock, yBlock, zBlock);
			break;
	}

	CudaCheckError();

	return;
}
