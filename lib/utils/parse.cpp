#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <cmath>
#include <limits>
#include <sys/stat.h>
#include <vector>

#include "enum-vars.h"
#include "utils/initGpu.h"
#include "utils/iParms.h"
#include "utils/memAlloc.h"

#include "cosmos/cosmos.h"

char outName[128] = "axiton\0";
char outDir[1024] = "out/m/\0";
char wisDir[1024] = "./\0";

bool mCreateOut = false;

void	createOutput() {
	struct stat tStat;

	if (mCreateOut == false)
		return;

	if (stat("out", &tStat) != 0) {
		auto  dErr = mkdir("out", S_IRWXU | S_IRWXG | S_IROTH | S_IXOTH);
		if (dErr == -1) {
			printf("Error: can't write on filesystem\n");
			exit(1);
		}
	}

	if (stat("out/m", &tStat) != 0) {
		auto  dErr = mkdir("out/m", S_IRWXU | S_IRWXG | S_IROTH | S_IXOTH);
		if (dErr == -1) {
			printf("Error: can't write on filesystem\n");
			exit(1);
		}
	}
}

void	PrintUsage(char *name)
{
	printf("\nUsage: %s [Options]\n\n", name);

	printf("\nOptions:\n");

	printf("\nSize of the grid:\n");
	printf("  --size  [int]                 Number of lattice points along x and y (Lx). Local size is Lx^2 x Lz (default 128).\n");
	printf("                                Splitting occurs in the z-dimension, so the total lattice is Lx^2 x (zgrid * Lz).\n");

	printf("\nSimulation parameters:\n");
	printf("  --prec  double/single         Precision of the axion field simulation (default single)\n");
	printf("  --prop  leap/rkn4/om2/om4     Numerical propagator to be used for molecular dynamics (default, use rkn4).\n");
	printf("  --steps [int]                 Number of steps of the simulation (default 500).\n");
	printf("  --nn    [int]                 Number of neighbours in the laplacian (1-4, default 1).\n");
	printf("  --wDz   [float]               Adaptive time step dz = wDz/frequency [l/raxion3D].\n");

	printf("\nPhysical parameters:\n");
	printf("  --zi    [float]               Initial value of the conformal time z (default 0.5).\n");
	printf("  --zf    [float]               Final value of the conformal time z (default 1.0).\n");
	printf("  --Rc    [float]               Critical value of scale factor R, (mass_A^2 = constant for R>Rc) (default 1.e5).\n");
	printf("  --lsize [float]               Physical size of the system (default 4.0).\n");
	printf("  --qcd   [float]               Exponent of topological susceptibility (default 7).\n");
	printf("  --ind3  [float]               Factor multiplying axion mass^2 (default, 1).\n");
	printf("                                Setting 0.0 turns on massless Axion mode.\n");
	printf("  --gm    [float]               Damping for the field at large r (last 10\% of the points, default 1.0).\n");


	printf("\nInitial conditions:\n");
	printf("  --ctype iflat/sinc2           Initial configuration, either with smoothing or with FFT and a maximum momentum\n");
	printf("\n");
	printf("  --kmax  [int]                 Maximum momentum squared for the generation of the configuration with --ctype kmax/tkachev (default 2)\n");
	printf("  --kcr   [float]               kritical kappa (default 1.0).\n");
	printf("\n");
	printf("  --sIter [int]                 Number of smoothing steps for the generation of the configuration with --ctype smooth (default 40)\n");
	printf("  --alpha [float]               alpha parameter for the smoothing (default 0.143).\n");
	printf("\n");
	printf("  --index [idx]                 Loads HDF5 file at out/dump as initial conditions (default, don't load).\n");

	printf("\nOutput:\n");
	printf("--name  [filename]              Uses filename to name the output files in out/dump, instead of the default \"axion\"\n");
	printf("--dump  [int]                   frequency of the output (default 100).\n");
	printf("--wTime [float]                 Simulates during approx. [float] hours and then writes the configuration to disk.\n");

	printf("\nLogging:\n");
	printf("--help                          Prints this message.\n");

	return;
}

iParms	parseArgs (int argc, char *argv[])
{
	bool	passed;
	int	procArgs = 0;

	iParms	defaultParms;

	// Defaults
	defaultParms.nSize    = 16384;
	defaultParms.fExp     = Radiation;
	defaultParms.fMod     = FieldNonCompact;
	defaultParms.wDz      = 0.8;
	defaultParms.gamma    = 1.0;

	defaultParms.cType    = IcFlat;
	defaultParms.parm1    = 1.0;
	defaultParms.parm2    = 64;

	defaultParms.lSize    = 4.0;
	defaultParms.indi3    = 1.0;
	defaultParms.nQcd     = 7.0;
	defaultParms.zThRes   = 16.0;
	defaultParms.zRestore = 1.0e+10;
	defaultParms.zInit    = 0.5;
	defaultParms.zEnd     = 0.6;

	defaultParms.pType    = PropagatorRKN4;
	defaultParms.fPrec    = SinglePrecision;
	defaultParms.nNeig    = 1;
	defaultParms.nSteps   = 1000;
	defaultParms.dump     = 100;
	defaultParms.fIndex   = -1;
	defaultParms.wTime    = std::numeric_limits<std::size_t>::max();

	for (int i=1; i<argc; i++)
	{
		passed = false;

		if (!strcmp(argv[i], "--help"))
		{
			PrintUsage(argv[0]);
			exit(0);
		}

		if (!strcmp(argv[i], "--wTime"))
		{
			if (i+1 == argc)
			{
				printf("Error: I need a value for the walltime.\n");
				exit(1);
			}

			double wTime = atof(argv[i+1]);

			if (wTime < 0.)
			{
				printf("Error: Walltime must be larger than or equal to 0.\n");
				exit(1);
			}

			defaultParms.wTime = wTime*3600000000;	// Walltime is processed in microseconds, but the expected precision is much worse

			i++;
			procArgs++;
			passed = true;
			goto endFor;
		}

		if (!strcmp(argv[i], "--size"))
		{
			if (i+1 == argc)
			{
				printf("Error: I need a size.\n");
				exit(1);
			}

			int sizeN = atoi(argv[i+1]);

			if (sizeN < 2)
			{
				printf("Error: Size must be larger than 2.\n");
				exit(1);
			}

			defaultParms.nSize = sizeN;

			i++;
			procArgs++;
			passed = true;
			goto endFor;
		}

		if (!strcmp(argv[i], "--gm"))
		{
			if (i+1 == argc)
			{
				printf("Error: I need a value for gamma.\n");
				exit(1);
			}

			double gamma = atof(argv[i+1]);

			if (gamma < 0.)
			{
				printf("Error: Gamma must be larger than or equal to 0.\n");
				exit(1);
			}

			defaultParms.gamma = gamma;

			i++;
			procArgs++;
			passed = true;
			goto endFor;
		}

		if (!strcmp(argv[i], "--kcr"))
		{
			if (i+1 == argc)
			{
				printf("Error: I need a value for the critical kappa.\n");
				exit(1);
			}

			double kCrit = atof(argv[i+1]);

			if (kCrit < 0.)
			{
				printf("Error: Critical kappa must be larger than or equal to 0.\n");
				exit(1);
			}

			defaultParms.parm1 = kCrit;

			i++;
			procArgs++;
			passed = true;
			goto endFor;
		}

		if (!strcmp(argv[i], "--alpha"))
		{
			if (i+1 == argc)
			{
				printf("Error: I need a value for alpha.\n");
				exit(1);
			}

			double alpha = atof(argv[i+1]);

			if ((alpha < 0.) || (alpha > 1.))
			{
				printf("Error: Alpha parameter must belong to the [0,1] interval.\n");
				exit(1);
			}

			defaultParms.parm1 = alpha;

			i++;
			procArgs++;
			passed = true;
			goto endFor;
		}

		if (!strcmp(argv[i], "--zi"))
		{
			if (i+1 == argc)
			{
				printf("Error: I need a value for the initial conformal time.\n");
				exit(1);
			}

			double zInit = atof(argv[i+1]);

			if (zInit < 0.)
			{
				printf("Error: Initial conformal time must be larger than 0.\n");
				exit(1);
			}

			defaultParms.zInit = zInit;

			i++;
			procArgs++;
			passed = true;
			goto endFor;
		}

		if (!strcmp(argv[i], "--logi"))
		{
			if (i+1 == argc)
			{
				printf("Error: I need a value for the initial conformal time.\n");
				exit(1);
			}

			double logi = atof(argv[i+1]);

			i++;
			procArgs++;
			passed = true;
			goto endFor;
		}

		if (!strcmp(argv[i], "--zf"))
		{
			if (i+1 == argc)
			{
				printf("Error: I need a value for the Final conformal time.\n");
				exit(1);
			}

			double zFinl = atof(argv[i+1]);

			if (zFinl < 0.)
			{
				printf("Error: Final conformal time must be larger than 0.\n");
				exit(1);
			}

			defaultParms.zEnd = zFinl;

			i++;
			procArgs++;
			passed = true;
			goto endFor;
		}

		if (!strcmp(argv[i], "--lsize"))
		{
			if (i+1 == argc)
			{
				printf("Error: I need a value for the physical size of the universe.\n");
				exit(1);
			}

			double sizeL = atof(argv[i+1]);

			if (sizeL <= 0.)
			{
				printf("Error: Physical size must be greater than zero.\n");
				exit(1);
			}

			defaultParms.lSize = sizeL;

			i++;
			procArgs++;
			passed = true;
			goto endFor;
		}

		if (!strcmp(argv[i], "--ind3"))
		{
			if (i+1 == argc)
			{
				printf("Error: I need a value for axion mass ind3.\n");
				exit(1);
			}

			double indi3 = atof(argv[i+1]);

			if (indi3 < 0.)
			{
				printf("Error: Indi3 must be greater than or equal to zero.\n");
				exit(1);
			}

			defaultParms.indi3 = indi3;

			i++;
			procArgs++;
			passed = true;
			goto endFor;
		}


		if (!strcmp(argv[i], "--wDz"))
		{
			if (i+1 == argc)
			{
				printf("Error: I need a value for the adaptive time step.\n");
				exit(1);
			}

			double wDz = atof(argv[i+1]);

			if (wDz < 0.)
			{
				printf("Error: backwards propagation?\n");
				exit(1);
			}

			defaultParms.wDz = wDz;

			i++;
			procArgs++;
			passed = true;
			goto endFor;
		}


		if (!strcmp(argv[i], "--qcd"))
		{
			if (i+1 == argc)
			{
				printf("Error: I need an exponent for the susceptibility nQcd!.\n");
				exit(1);
			}

			double nQcd = atof(argv[i+1]);

			if (nQcd < 0)
			{
				printf("Error: The exponent of the top. susceptibility nQcd must be equal or greater than 0.\n");
				exit(1);
			}

			defaultParms.nQcd = nQcd;

			i++;
			procArgs++;
			passed = true;
			goto endFor;
		}

		if (!strcmp(argv[i], "--steps"))
		{
			if (i+1 == argc)
			{
				printf("Error: I need a number of steps.\n");
				exit(1);
			}

			int nSteps = atoi(argv[i+1]);

			if (nSteps < 0)
			{
				printf("Error: Number of steps must be > 0.\n");
				exit(1);
			}

			defaultParms.nSteps = nSteps;

			i++;
			procArgs++;
			passed = true;
			goto endFor;
		}

		if (!strcmp(argv[i], "--dump"))
		{
			if (i+1 == argc)
			{
				printf("Error: I need a print rate.\n");
				exit(1);
			}

			int dump = atoi(argv[i+1]);

			if (dump < 0)
			{
				printf("Error: Print rate must be equal or greater than zero.\n");
				exit(1);
			}

			defaultParms.dump = dump;

			i++;
			procArgs++;
			passed = true;
			goto endFor;
		}

		if (!strcmp(argv[i], "--name"))
		{
			if (i+1 == argc)
			{
				printf("Error: I need a name for the files.\n");
				exit(1);
			}

			if (strlen(argv[i+1]) > 96)
			{
				printf("Error: name too long, keep it under 96 characters\n");
				exit(1);
			}

			strcpy (outName, argv[i+1]);

			i++;
			procArgs++;
			passed = true;
			goto endFor;
		}

		if (!strcmp(argv[i], "--kmax"))
		{
			if (i+1 == argc)
			{
				printf("Error: I need an integer value for the maximum momentum.\n");
				exit(1);
			}

			int kMax = atoi(argv[i+1]);

			if (kMax < 0)
			{
				printf("Error: the maximum momentum must be equal or greater than zero.\n");
				exit(1);
			}

			defaultParms.parm2 = kMax;

			i++;
			procArgs++;
			passed = true;
			goto endFor;
		}

		if (!strcmp(argv[i], "--sIter"))
		{
			if (i+1 == argc)
			{
				printf("Error: I need a number of iterations for the smoothing.\n");
				exit(1);
			}

			int iter = atoi(argv[i+1]);

			if (iter < 0)
			{
				printf("Error: Number of iterations must be equal or greater than zero.\n");
				exit(1);
			}

			defaultParms.parm2 = iter;

			i++;
			procArgs++;
			passed = true;
			goto endFor;
		}

		if (!strcmp(argv[i], "--index"))
		{
			if (i+1 == argc)
			{
				printf("Error: I need an index for the file.\n");
				exit(1);
			}

			int fIndex = atoi(argv[i+1]);

			if (fIndex < 0)
			{
				printf("Error: Filename index must be equal or greater than zero.\n");
				exit(1);
			}

			defaultParms.fIndex = fIndex;
			defaultParms.cType  = IcRead;

			i++;
			procArgs++;
			passed = true;
			goto endFor;
		}

		if (!strcmp(argv[i], "--nn"))
		{
			if (i+1 == argc)
			{
				printf("Error: I need a number of neighbours.\n");
				exit(1);
			}

			int nNeig = atoi(argv[i+1]);

			if (nNeig < 1 || nNeig > 4)
			{
				printf("Error: The number of neighbours must be comprised between 1 and 4.\n");
				exit(1);
			}

			defaultParms.nNeig = nNeig;

			i++;
			procArgs++;
			passed = true;
			goto endFor;
		}

		if (!strcmp(argv[i], "--ctype"))
		{
			if (i+1 == argc)
			{
				printf("Error: I need a value for the configuration type (flat/r2/hyper/sinc2).\n");
				exit(1);
			}

			InitialCond	cType = IcNone;

			if (!strcmp(argv[i+1], "flat"))
				cType = IcFlat;
			else if (!strcmp(argv[i+1], "r2"))
				cType = IcOvr2;
			else if (!strcmp(argv[i+1], "hyper"))
				cType = IcOvX;
			else if (!strcmp(argv[i+1], "sinc2"))
				cType = IcSinc2;
			else
			{
				printf("Error: Unrecognized configuration type %s\n", argv[i+1]);
				exit(1);
			}

			if (defaultParms.cType == IcRead)
				printf("Conflicting types (reading configuration from disk and setting initial conditions). Overriding with data from disk.\n");
			else
				defaultParms.cType = cType;

			i++;
			procArgs++;
			passed = true;
			goto endFor;
		}

		if (!strcmp(argv[i], "--prec"))
		{
			FieldPrecision sPrec;

			if (i+1 == argc)
			{
				printf("Error: I need a value for the precision (double/single/mixed).\n");
				exit(1);
			}

			if (!strcmp(argv[i+1], "double"))
			{
				sPrec = DoublePrecision;
			}
			else if (!strcmp(argv[i+1], "single"))
			{
				sPrec = SinglePrecision;
			}
			else
			{
				printf("Error: Unrecognized precision %s\n", argv[i+1]);
				exit(1);
			}

			defaultParms.fPrec = sPrec;

			i++;
			procArgs++;
			passed = true;
			goto endFor;
		}

		if (!strcmp(argv[i], "--prop"))
		{
			if (i+1 == argc)
			{
				printf("Error: I need a propagator class (leap/leap4/rkn4/om2/om4).\n");
				exit(1);
			}

			PropagatorType	pType;

			if (!strcmp(argv[i+1], "leap"))
			{
				pType = PropagatorLeapFrog2;
			}
			else if (!strcmp(argv[i+1], "leap4"))
			{
				pType = PropagatorLeapFrog4;
			}
			else if (!strcmp(argv[i+1], "rkn4"))
			{
				pType = PropagatorRKN4;
			}
			else if (!strcmp(argv[i+1], "om2"))
			{
				pType = PropagatorOmelyan2;
			}
			else if (!strcmp(argv[i+1], "om4"))
			{
				pType = PropagatorOmelyan4;
			}
			else
			{
				printf("Error: unrecognized propagator %s\n", argv[i+1]);
				exit(1);
			}

			defaultParms.pType = pType;

			i++;
			procArgs++;
			passed = true;
			goto endFor;
		}

		if (!strcmp(argv[i], "--cax"))
		{
			defaultParms.fMod = FieldCompact;

			procArgs++;
			passed = true;
			goto endFor;
		}

		endFor:

		if (!passed)
		{
			PrintUsage(argv[0]);
			printf("\n\nUnrecognized option %s\n", argv[i]);
			exit(1);
		}

	}

	/*	Set the output directory, according to an environmental variable	*/

	if (const char *outPath = std::getenv("AXITONS_OUTPUT")) {
		if (strlen(outPath) < 1022) {
			struct stat tStat;
			if (stat(outPath, &tStat) == 0 && S_ISDIR(tStat.st_mode)) {
				strcpy(outDir, outPath);
			} else {
				printf("Path %s doesn't exist, using default\n", outPath);
				mCreateOut = true;
			}
		}
	} else {
		mCreateOut = true;
	}

	/*	Set the directory where the FFTW wisdom is/will be stored		*/

	if (const char *wisPath = std::getenv("AXITONS_WISDOM")) {
		if (strlen(wisPath) < 1022) {
			struct stat tStat;
			if (stat(wisPath, &tStat) == 0 && S_ISDIR(tStat.st_mode)) {
				strcpy(wisDir, wisPath);
			} else {
				printf("Path %s doesn't exist, using default\n", wisPath);
			}
		}
	}

	defaultParms.procArgs = procArgs;

	return	defaultParms;
}

Cosmos	initAxitons(int argc, char *argv[]) {

	iParms myParms = parseArgs(argc, argv);

	initGpu();

	createOutput();

	printf ("Output folder set to %s\n", outDir);
	printf ("FFTW wisdom folder set to %s\n", wisDir);

	myParms.outputDir  = outDir;
	myParms.outputName = outName;
	myParms.wisdomDir  = wisDir;

	return	Cosmos(myParms);
}

void	endAxitons() {
	printMemStats();

	return;
}
