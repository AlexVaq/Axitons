#ifndef	InitialParametersGuard
	#define	InitialParametersGuard

	#include"enum-vars.h"

	struct	iParms
	{
		int		procArgs;

		int		nSize;
		FieldPrecision	fPrec;
		PropagatorType	pType;
		FieldType	fMod;
		double		wDz;

		InitialCond	cType;
		int		parm1;
		double		parm2;

		FieldExpansion	fExp;
		double		lSize;
		double		indi3;
		double		nQcd;
		double		zThRes;
		double		zRestore;
		double		zInit;
		double		zEnd;

		int		nSteps;
		int		dump;
		int		fIndex;
		double		wTime;
	};
#endif
