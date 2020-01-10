#ifndef	GeneratorGuard
	#define	GeneratorGuard

	#include "enum-vars.h"
	#include "fields/fields.h"
	#include "cosmos/cosmos.h"
	#include <gen/spline.h>

	class	Generator {

		private:

		InitialCond	icType;
		Axiton			*field;
		Cosmos			*cosmo;

		tk::spline sm, sv;

		public:

			Generator (InitialCond icType, Axiton *field, Cosmos *cosmo) : icType(icType), field(field), cosmo(cosmo) {}

		void	Construct (double parm1, int parm2, double zInit);
		void  SplineSetup();
		void 	fillGen();
	};
#endif
