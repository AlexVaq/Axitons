#ifndef	AxionMassGuard
#define	AxionMassGuard

#include <cmath>

class	AxionMass {
  private:

    double	zThRes;
    double	zRestore;
    double	indi3;
    double	nQcd;

  public:

    AxionMass(double zThRes, double zRestore, double indi3, double nQcd) : zThRes(zThRes), zRestore(zRestore), indi3(indi3), nQcd(nQcd) {}

    double	operator()(double RNow)
    {
        double aMass;

        if ((RNow > zThRes) &&  (zThRes < zRestore))
        {
            aMass = indi3*indi3*pow(zThRes, nQcd);
            if (RNow > zRestore)
                aMass *= pow(RNow/zRestore, nQcd);
        }
        else
            aMass = indi3*indi3*pow(RNow, nQcd);

        return aMass;
    }
};
#endif
