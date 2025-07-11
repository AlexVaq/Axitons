
#ifndef	CosmosGuard
#define	CosmosGuard

#include <memory>
#include "enum-vars.h"
#include "utils/iParms.h"
#include "cosmos/axionMass.h"

class	Cosmos {
  private:

    int	   nSize;
    double lSize;
    double delta;
    double gamma;

    FieldExpansion fExp;
    std::unique_ptr<AxionMass> aMass;

    iParms			   myParms;

  public:

    Cosmos(iParms myParms) : nSize(myParms.nSize), lSize(myParms.lSize),
                             delta(lSize/((double) nSize)),gamma(myParms.gamma),
                             fExp(myParms.fExp),   myParms(myParms)
    {
        aMass = std::make_unique<AxionMass> (myParms.zThRes, myParms.zRestore,
                                             myParms.indi3, myParms.nQcd);
    }

    FieldExpansion Expansion  () const   { return fExp;  }
    int		   CosmosLatt () const   { return nSize; }
    double	   CosmosSize () const   { return lSize; }
    double	   Delta      () const   { return delta; }
    double	   Gamma      () const   { return gamma; }
    double	   AxionMassSq(double R) { return (*aMass)(R); }
    iParms&	   InitParms  ()         { return myParms; }
};
#endif
