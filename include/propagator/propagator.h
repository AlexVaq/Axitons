#ifndef	PropagatorGuard
	#define	PropagatorGuard

	#include "enum-vars.h"
	#include "cosmos/cosmos.h"
	#include "fields/fields.h"

	void	initPropagator	(PropagatorType pType, Cosmos *bck, Axiton *field, int nNeig);
	void	propagate	(const double dz);
#endif
