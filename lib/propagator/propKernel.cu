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

template<typename Float, const bool wMod>
static __device__ __forceinline__ void	propagateCoreGpu(const uint idx, const Float * __restrict__ field, Float * __restrict__ dev, Float * __restrict__ misc, const Float zQ,
							 const Float iz, const Float dzc, const Float dzd, const Float ood2, const int Lx, const Float zP, const Float tPz)
{
	Float mel, pFc, a, tmp;

	if (idx != 0) { // FIXME SLOW AS HELL
		tmp = field[idx];

		//if (wMod) {
		//	dev = modPi(mPr - tmp, zP, tPz);
		//} else
			pFc = 1 + ((Float) (2*idx+1))/(idx*idx);
			mel = pFc*field[idx+2] - (1 + pFc)*field[idx+1] + tmp;
	} else
		mel = 0.;

	a = mel*ood2 - zQ*sin(tmp*iz);

	mel	 = dev[idx];
	mel	+= a*dzc;
	dev[idx] = mel;
	mel	*= dzd;
	tmp	+= mel;

	misc[idx] = tmp; //modPi(tmp, zP, tPz);
}

template<typename Float, const bool wMod>
__global__ void	propagateKernel(const Float * __restrict__ field, Float * __restrict__ dev, Float * __restrict__ misc, const Float zQ, const Float dzc, const Float dzd,
				const Float ood2, const Float iz, const int Lx, const Float zP=0, const Float tPz=0)
{
	uint idx = threadIdx.x + blockDim.x*(blockIdx.x + gridDim.x*blockIdx.y);
	//uint idx = Vo + (threadIdx.x + blockDim.x*blockIdx.x) + Sf*(threadIdx.y + blockDim.y*blockIdx.y);

	if	(idx >= Lx-2)
		return;

	propagateCoreGpu<Float,wMod>(idx, field, dev, misc, zQ, iz, dzc, dzd, ood2, Lx, zP, tPz);
}

void	propGpu(const void * __restrict__ field, void * __restrict__ dev, void * __restrict__ misc, const double z, const double dz, const double c, const double d, const double ood2,
		const double aMass2, const int Lx, FieldPrecision precision, const int xBlock, const int yBlock, const int zBlock)
{
	#define	BLSIZE 256
	dim3 gridSize((Lx+BLSIZE-1)/BLSIZE,1,1);
	dim3 blockSize(BLSIZE,1,1);

	if (precision == DoublePrecision)
	{
		const double dzc  = dz*c;
		const double dzd  = dz*d;
		const double zQ   = aMass2*z*z*z;//axionmass2((double) zR, nQcd, zthres, zrestore)*zR*zR*zR;
		const double iZ   = 1./z;
		propagateKernel<double,false><<<gridSize,blockSize,0>>>((const double *) field, (double *) dev, (double *) misc, zQ, dzc, dzd, ood2, iZ, Lx);
	}
	else if (precision == SinglePrecision)
	{
		const float dzc = dz*c;
		const float dzd = dz*d;
		const float zQ = (float) (aMass2*z*z*z);//axionmass2((double) zR, nQcd, zthres, zrestore)*zR*zR*zR;
		const float iZ   = 1./z;
		propagateKernel<float, false><<<gridSize,blockSize,0>>>((const float *) field, (float *) dev, (float *) misc, zQ, dzc, dzd, (float) ood2, iZ, Lx);
	}
}

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

void	propGpu(const void * __restrict__ field, void * __restrict__ dev, void * __restrict__ misc, const double z, const double dz, const double c, const double d, const double ood2,
		const double aMass2, const int Lx, FieldPrecision precision, const int xBlock, const int yBlock, const int zBlock, const FieldType wMod)
{
	switch (wMod) {
	
		case	FieldCompact:
			propModGpu(field, dev, misc, z, dz, c, d, ood2, aMass2, Lx, precision, xBlock, yBlock, zBlock);
			break;

		case	FieldNonCompact:
			propGpu	  (field, dev, misc, z, dz, c, d, ood2, aMass2, Lx, precision, xBlock, yBlock, zBlock);
			break;
	}

	CudaCheckError();

	return;
}
