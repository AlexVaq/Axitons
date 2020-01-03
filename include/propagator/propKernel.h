#ifndef	PropagatorGpu_
	#define	PropagatorGpu

	void	propGpu	(const void * __restrict__ field, void * __restrict__ dev, void * __restrict__ misc, const double R, const double dz, const double c, const double d,
			 const double ood2, const double aMass2, const int Lx, FieldPrecision precision, const int xBlock, const int yBlock, const int zBlock, const FieldType wMod);
#endif
