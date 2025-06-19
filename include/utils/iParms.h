#ifndef	InitialParametersGuard
#define	InitialParametersGuard

#include <string>
#include"enum-vars.h"

struct	iParms {
    int		procArgs;

    int		nSize;
    FieldPrecision	fPrec;
    PropagatorType	pType;
    FieldType	fMod;
    double		wDz;

    InitialCond	cType;
    double		parm1;
    int		parm2;

    FieldExpansion	fExp;
    double		lSize;
    double		indi3;
    double		nQcd;
    double		zThRes;
    double		zRestore;
    double		zInit;
    double		zEnd;

    int		nNeig;
    int		nSteps;
    int		dump;
    int		fIndex;
    double		wTime;
    double		gamma;

    std::string	outputName;
    std::string	outputDir;
    std::string	wisdomDir;
};
#endif
