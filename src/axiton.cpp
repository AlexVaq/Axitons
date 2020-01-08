#include <cmath>
#include <cstring>
#include <chrono>
#include <iostream>
#include <vector>

#include "axiton.h"
#include <cuda.h>
#include <cuda_runtime.h>
#include "cudaErrors.h"


using namespace std;

int	main (int argc, char *argv[])
{
	Cosmos myCosmos = initAxitons(argc, argv);
	iParms &myParms = myCosmos.InitParms();

	Axiton		*axiton;

	axiton = new Axiton (myParms.nSize, myParms.fPrec);

	Hdf5ReadWriter	IOHandler(myParms);
	Generator	AxitonFactory(myParms.cType, axiton);

	auto R   = axiton->R<Radiation>();
	printf (" - z = %f R = %f o2 = %e\n",  axiton->z(), R, 1.0/myCosmos.Delta());
	AxitonFactory.Construct(myParms.parm1*myParms.zInit, myParms.parm2, myParms.zInit);

	IOHandler.writeConf(&myCosmos, axiton);
	IOHandler.nextFile ();

	printf ("Transferring configuration to device\n"); fflush(stdout);
	axiton->transferField(FieldBase | FieldDev, HostToDevice);
	CudaCheckError();

	double delta = myCosmos.Delta();
	double dz;

	if (myParms.wDz == 0.0){
		if (myParms.nSteps == 0)
			dz = 0.;
		else
			dz = (myParms.zEnd - myParms.zInit)/((double) myParms.nSteps);
		printf ("User specified time step (dz = %.2e)\n",dz); fflush(stdout);
	} else {
		dz = myParms.wDz*delta/sqrt(12.);
		printf ("Dynamical time step (dz = %.2e)\n",dz); fflush(stdout);
	}

	std::chrono::high_resolution_clock::time_point start, current, old;
	std::chrono::milliseconds elapsed;

	int counter = 0;
	auto dump = myParms.dump;

	if (dump > myParms.nSteps)
		dump = myParms.nSteps;

	int nLoops;

	if (dump == 0)
		nLoops = 0;
	else
		nLoops = (int)(myParms.nSteps/dump);

	printf ("Start redshift loop\n"); fflush(stdout);

	start = std::chrono::high_resolution_clock::now();
	old = start;

	printf ("Initial z = %f\n",  axiton->z());
	initPropagator (myParms.pType, &myCosmos, axiton, myParms.nNeig);

	bool endLoop = false;
	for (int zloop = 0; zloop < nLoops; zloop++)
	{
		if (endLoop){
			propagate (dz);
			break;
		}

		for (int zsubloop = 0; zsubloop < dump; zsubloop++) {
			propagate (dz);
		}

		axiton->transferField(FieldBase | FieldDev, DeviceToHost);
		IOHandler.nextFile();
		IOHandler.writeConf(&myCosmos, axiton);
		auto R   = axiton->R<Radiation>();
		auto maa = myCosmos.AxionMassSq(R);
		printf (" - z = %f R = %f mA = %e o2 = %e\n",  axiton->z(), R, sqrt(maa),1.0/myCosmos.Delta());

		if (axiton->z() + dz > myParms.zEnd)
			{
				dz = myParms.zEnd - axiton->z();
				endLoop = true;
			}
	}

	current = std::chrono::high_resolution_clock::now();
	elapsed = std::chrono::duration_cast<std::chrono::milliseconds>(current - start);

	axiton->transferField(FieldBase | FieldDev, DeviceToHost);
	IOHandler.nextFile();
	IOHandler.writeConf(&myCosmos, axiton);

	printf ("Final z  = %f\n", axiton->z());
	printf ("#_steps  = %i\n", nLoops*dump);
	printf ("#_prints = %i\n", IOHandler.currentIndex());
	printf ("Total time: %2.3f s\n", elapsed.count()*1.e-3);

	delete axiton;

	endAxitons();

	return 0;
}
