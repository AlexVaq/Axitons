
#ifndef	CosmosGuard
	#define	CosmosGuard

	#include <memory>
	#include "enum-vars.h"
	#include "utils/iParms.h"
	#include "cosmos/axionMass.h"

	class	Cosmos
	{
		private:

		int			   nSize;
		double			   lSize;
		double			   delta;
		FieldExpansion		   fExp;
		std::unique_ptr<AxionMass> aMass;

		public:

				Cosmos(iParms myParms) : nSize(myParms.nSize), lSize(myParms.lSize), delta(lSize/nSize), fExp(myParms.fExp) {

	 		aMass = std::make_unique<AxionMass> (myParms.zThRes, myParms.zRestore, myParms.indi3, myParms.nQcd);
		}

		FieldExpansion	Expansion  ()	      { return fExp;        }
		int		CosmosLatt ()	      { return nSize;       }
		double		CosmosSize ()	      { return lSize;       }
		double		Delta      ()	      { return delta;       }
		double		AxionMassSq(double R) { return (*aMass)(R); }
	};
#endif
