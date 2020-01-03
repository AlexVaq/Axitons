#include <cmath>
#include <cstring>
#include <chrono>

#include <complex>
#include <vector>

#include "propagator/allProp.h"
#include "utils/utils.h"
#include "fields/fields.h"

#include <iostream>

#define	ScaleFactor 1.5

using namespace std;

int	main (int argc, char *argv[])
{
	Cosmos myCosmos = initAxitons(argc, argv);

	//--------------------------------------------------
	//       READING INITIAL CONDITIONS
	//--------------------------------------------------

	Field *axiton;
	char fileName[256];

	std::cout << zInit << std::endl;

	if ((fIndex == -1) && (myCosmos.ICData().cType == CONF_NONE)) {
		LogOut("Error: Neither initial conditions nor configuration to be loaded selected. Empty field.\n");
	} else {
		if (fIndex == -1)
			//This generates initial conditions
			axiton = new Axiton (&myCosmos, sizeN, sPrec, zInit, fTypeP, lType);
		else
		{
			//This reads from an Axion.$fIndex file
			readConf(&myCosmos, &axiton, fIndex);
			if (axiton == nullptr)
			{
				printf ("Error reading HDF5 file\n");
				exit (0);
			}
		}
	}

	//--------------------------------------------------
	//          SETTING BASE PARAMETERS
	//--------------------------------------------------

	double delta = axiton->Delta();
	double dz;

	if (nSteps == 0)
		dz = 0.;
	else
		dz = (zFinl - zInit)/((double) nSteps);

	LogOut("--------------------------------------------------\n");
	LogOut("           INITIAL CONDITIONS                     \n\n");

	LogOut("Length =  %2.5f\n", myCosmos.PhysSize());
	LogOut("N      =  %ld\n",   axion->Length());
	LogOut("Nz     =  %ld\n",   axion->Depth());
	LogOut("zGrid  =  %ld\n",   zGrid);
	LogOut("dx     =  %2.5f\n", delta);
	LogOut("dz     =  %2.5f\n", dz);
	LogOut("LL     =  %2.5f\n", myCosmos.Lambda());
	LogOut("--------------------------------------------------\n");

	//--------------------------------------------------
	//   THE TIME ITERATION LOOP
	//--------------------------------------------------

	LogOut("--------------------------------------------------\n");
	LogOut("           STARTING COMPUTATION                   \n");
	LogOut("--------------------------------------------------\n");

	std::chrono::high_resolution_clock::time_point start, current, old;
	std::chrono::milliseconds elapsed;

	int counter = 0;
	int index = 0;

	if (fIndex == -1)
	{
		LogOut ("Dumping configuration %05d ...", index);
		writeConf(axion, index);
		LogOut ("Done!\n");
	}
	else
		index = fIndex;

	if (dump > nSteps)
		dump = nSteps;

	int nLoops;

	if (dump == 0)
		nLoops = 0;
	else
		nLoops = (int)(nSteps/dump);

	LogOut ("Start redshift loop\n");

	start = std::chrono::high_resolution_clock::now();
	old = start;

	initPropagator (pType, axiton);

	LogOut("--------------------------------------------------\n");
	LogOut("            TUNING PROPAGATOR                     \n");
	LogOut("--------------------------------------------------\n");

	tunePropagator (axiton);

	for (int zloop = 0; zloop < nLoops; zloop++)
	{
		//--------------------------------------------------
		// THE TIME ITERATION SUB-LOOP
		//--------------------------------------------------

		index++;

		for (int zsubloop = 0; zsubloop < dump; zsubloop++)
			propagate (dz);

/*
		createMeas  (axion, index);
		writeMapHdf5(axion);
		writeEDens  (axion, MAP_ALL);
		writeString (axion, strDen);
		writeEnergy (axion, eRes);
		writePoint  (axion);
		cDensityMap (axion);
		destroyMeas ();
*/
	} // zloop

	current = std::chrono::high_resolution_clock::now();
	elapsed = std::chrono::duration_cast<std::chrono::milliseconds>(current - start);

	LogOut("\n PROGRAMM FINISHED\n");

	LogOut ("Transferring configuration to host\n");
	axion->transferCpu(FIELD_MV);

	writeConf(axion, index);

	LogOut("z_final = %f\n", *axion->zV());
	LogOut("#_steps = %i\n", counter);
	LogOut("#_prints = %i\n", index);
	LogOut("Total time: %2.3f s\n", elapsed.count()*1.e-3);

	delete axiton;

	endAxiton();

	return 0;
}
