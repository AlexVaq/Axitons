#ifndef	PropagatorBaseGuard
	#define	PropagatorBaseGuard

	#include "enum-vars.h"
	#include "utils/utils.h"
	#include "fields/fields.h"

	class	PropBase : public Tunable
	{
		protected:

		std::string	baseName;

		public:

			 PropBase() {};
		virtual	~PropBase() {};

		inline void	setBaseName	(const char *bName) { baseName.assign(bName); }
		inline void	getBaseName	() 		    { name = baseName; }

		virtual void	propagate	(const double) = 0;

		virtual double	cFlops		() = 0;
		virtual double	cBytes		() = 0;
	};
#endif
