#ifndef	PropagatorClassGuard
#define	PropagatorClassGuard

#include <cmath>
#include <string>
#include <cstring>
#include <functional>
#include "propagator/propBase.h"
#include "fields/fields.h"
#include "cosmos/cosmos.h"
#include "utils/utils.h"

#include <cuda.h>
#include <cuda_runtime.h>
#include <cuda_device_runtime_api.h>
#include "cudaErrors.h"
#include "propagator/propKernel.h"

#ifdef	__GNUG__
#pragma GCC optimize ("unroll-loops")
#endif

template<const int nStages, const bool lastStage>
class	PropClass : public PropBase {
  protected:

    Cosmos		*background;
    Axiton		*field;
    const size_t	Lx;
    const double	ood2;

    double	c[nStages + (lastStage == true ? 1 : 0)];
    double	d[nStages];

    const int	nNeig;

    template<FieldExpansion fExp>
    inline void	propCore	(const double);

  public:

    inline		 PropClass	(Cosmos *background, Axiton *field, int nNeig);
    inline		~PropClass	() override {};

    inline void	setCoeff	(const double * __restrict__ nC, const double * __restrict__ nD) {
        for(int i=0; i<nStages; i++) {
            c[i] = nC[i];
            d[i] = nD[i];
        }
        if (lastStage) {
            c[nStages] = nC[nStages];
        }
    }


    inline void	propagate	(const double)	override;

    inline double	cFlops		()		override;
    inline double	cBytes		()		override;
};

template<const int nStages, const bool lastStage>
PropClass<nStages, lastStage>::PropClass(Cosmos *background, Axiton *field, int nNeig) : background(background), field(field), nNeig(nNeig), Lx(background->CosmosLatt()),
    ood2(1./(background->Delta()*background->Delta())) {

    /*	Default block size gives just one block	*/

    xBest = xBlock = 256;
    yBest = yBlock = 1;
    zBest = zBlock = 1;
}

template<const int nStages, const bool lastStage>
template<FieldExpansion fExp>
void	PropClass<nStages, lastStage>::propCore	(const double dz) {

    auto wMod = field->Field();

#pragma unroll
    for (int s = 0; s<nStages; s+=2) {
        const double	c1 = c[s], c2 = c[s+1], d1 = d[s], d2 = d[s+1];

        auto R   = field->R<fExp>();
        auto maa = background->AxionMassSq(R);
        propGpu(field->fieldGpu(), field->devGpu(), field->miscGpu(), R, dz, c1, d1, ood2, maa, Lx, field->Precision(), nNeig, background->Gamma(), xBlock, yBlock, zBlock, wMod);
        cudaDeviceSynchronize();        // This is not strictly necessary, but simplifies things a lot

        field->zUpdate(dz*d1);
        R   = field->R<fExp>();
        maa = background->AxionMassSq(R);
        propGpu(field->miscGpu(), field->devGpu(), field->fieldGpu(), R, dz, c2, d2, ood2, maa, Lx, field->Precision(), nNeig, background->Gamma(), xBlock, yBlock, zBlock, wMod);
        cudaDeviceSynchronize();        // This is not strictly necessary, but simplifies things a lot
        field->zUpdate(dz*d2);
    }

    if (lastStage) {
        const double c0  = c[nStages];
        auto         R   = field->R<fExp>();
        auto         maa = background->AxionMassSq(R);
        propGpu(field->fieldGpu(), field->devGpu(), field->fieldGpu(), R, dz, c0, 0., ood2, maa, Lx, field->Precision(), nNeig, background->Gamma(), xBlock, yBlock, zBlock, wMod);
        cudaDeviceSynchronize();        // This is not strictly necessary, but simplifies things a lot
    }
}

template<const int nStages, const bool lastStage>
void	PropClass<nStages, lastStage>::propagate	(const double dz) {

    switch (background->Expansion()) {
    case	Minkowski:
        propCore<Minkowski> (dz);
        break;

    case	Radiation:
        propCore<Radiation> (dz);
        break;
    }

    field->setStatus (FieldBaseDev, FieldGpu);
}

template<const int nStages, const bool lastStage>
double	PropClass<nStages, lastStage>::cFlops	() {

    switch (field->Field()) {
    case FieldCompact:
        return	(1e-9 * ((double) field->Size()) * (21.0 * ((double) nStages) + (lastStage ? 19.0 : 0.)));
        break;

    case FieldNonCompact:
        return	(1e-9 * ((double) field->Size()) * (19.0 * ((double) nStages) + (lastStage ? 17.0 : 0.)));
        break;
    }

    return	0.;
}

template<const int nStages, const bool lastStage>
double	PropClass<nStages, lastStage>::cBytes	() {

    return	(1e-9 * ((double) (field->Size()*field->Precision())) * (6.0 * ((double) nStages) + (lastStage ? 4.0 : 0.)));
}

#endif
