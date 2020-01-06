#ifndef	GeneratorGuard
	#define	GeneratorGuard

	#include "enum-vars.h"
	#include "fields/fields.h"

	class	Generator {

		private:

		InitialCond	icType;
		Axiton		*field;

		public:

			Generator (InitialCond icType, Axiton *field) : icType(icType), field(field) {}

		void	Construct (double parm1, int parm2, double zInit);
	};
#endif
