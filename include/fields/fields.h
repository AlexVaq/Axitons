#ifndef	axitonClassGuard
	#define	axitonClassGuard

	#include "enum-vars.h"
	#include "utils/tunable.h"

	class	Axiton : public Tunable {

		private:

			void		*hField;
			void		*dField;
			void		*hDev;
			void		*dDev;
			void		*hMisc;
			void		*dMisc;

			FieldPrecision	prec;
			int		length;
			int		bytes;

			double		zVar;

			FieldType	fMod;

		public:

					 Axiton      (int len, FieldPrecision pr);
					~Axiton      ();

			void		*fieldCpu    ()         { return hField; }
			const void	*fieldCpu    () const   { return hField; }
			void		*fieldGpu    ()         { return dField; }
			const void	*fieldGpu    () const   { return dField; }
			void		*devCpu      ()         { return hDev;   }
			const void	*devCpu      () const   { return hDev;   }
			void		*devGpu      ()         { return dDev;   }
			const void	*devGpu      () const   { return dDev;   }
			void		*miscCpu     ()         { return hMisc;  }
			const void	*miscCpu     () const   { return hMisc;  }
			void		*miscGpu     ()         { return dMisc;  }
			const void	*miscGpu     () const   { return dMisc;  }

			int		Size         ()         { return length; }
			int		Data         ()         { return bytes;  }
			FieldPrecision	Precision    ()         { return prec;   }
			FieldType	Field        () const   { return fMod;   }

			double		z	     ()         { return zVar;   }
			const double	z	     () const   { return zVar;   }
			void		zUpdate	     (double z) { zVar += z;     }

			template<FieldExpansion fExp>
			double		R	     ()         { switch(fExp) { case Minkowski: { return 1.0; break; } case Radiation: { return zVar; break; } } }
			template<FieldExpansion fExp>
			const double	R	     () const   { switch(fExp) { case Minkowski: { return 1.0; break; } case Radiation: { return zVar; break; } } }

			void		transferField(FieldIndex fIdx, MemDirection mIdx);
	};
#endif
