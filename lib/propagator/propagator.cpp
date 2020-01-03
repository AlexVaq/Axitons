#include <cstdio>
#include <cstdlib>
#include <memory>
#include <chrono>
#include <string>
#include "fields/fields.h"
#include "enum-vars.h"
#include "propagator/propClass.h"
#include "utils/utils.h"

std::unique_ptr<PropBase> prop;

class	PropLeap : public PropClass<2, true> {

	public:
		PropLeap(Cosmos *bck, Axiton *field) :
		PropClass<2, true>(bck, field) {
		//	Set up Leapfrog parameters

		double nC[3] = { 0.5, 0.5, 0.0 };
		double nD[2] = { 1.0, 0.0 };

		this->setCoeff(nC, nD);
	}
};

class	PropMLeap : public PropClass<4, true> {

	public:
		PropMLeap(Cosmos *bck, Axiton *field) :
		PropClass<4, true>(bck, field) {

		//	Set up Leapfrog parameters

		double nC[5] = { 0.125, 0.25, 0.25, 0.25, 0.125 };
		double nD[4] = { 0.25,  0.25, 0.25, 0.25 };

		this->setCoeff(nC, nD);
	}
};

class	PropOmelyan2 : public PropClass<2, true> {

	public:
		PropOmelyan2(Cosmos *bck, Axiton *field) :
		PropClass<2, true>(bck, field) {
		constexpr double chi = +0.19318332750378360;

		//	Set up Omelyan parameters for BABAB

		double nC[3] = { chi, 1.-2.*chi, chi };
		double nD[2] = { 0.5, 0.5 };

		this->setCoeff(nC, nD);
	}
};

class	PropOmelyan4 : public PropClass<4, true> {

	public:
		PropOmelyan4(Cosmos *bck, Axiton *field) :
		PropClass<4, true>(bck, field) {
		constexpr double xi  = +0.16449865155757600;
		constexpr double lb  = -0.02094333910398989;
		constexpr double chi = +1.23569265113891700;

		//	Set up Omelyan parameters for BABABABAB

		double nC[5] = { xi, chi, 1.-2.*(xi+chi), chi, xi };
		double nD[4] = { 0.5*(1.-2.*lb), lb, lb, 0.5*(1.-2.*lb) };

		this->setCoeff(nC, nD);
	}
};

class	PropRKN4 : public PropClass<4, false> {

	public:
		PropRKN4(Cosmos *bck, Axiton *field) :
		PropClass<4, false>(bck, field) {
		//	Set up RKN parameters for BABABABA

		const double nC[4] = { +0.1344961992774310892, -0.2248198030794208058, +0.7563200005156682911, +0.3340036032863214255 };
		const double nD[4] = { +0.5153528374311229364, -0.085782019412973646,  +0.4415830236164665242, +0.1288461583653841854 };

		this->setCoeff(nC, nD);
	}
};


void	initPropagator	(PropagatorType pType, Cosmos *bck, Axiton *field) {

	bool wasTuned = false;

	unsigned int xBlock, yBlock, zBlock;

	if (prop != nullptr)
		if (prop->IsTuned()) {
			wasTuned = true;
			xBlock = prop->TunedBlockX();
			yBlock = prop->TunedBlockY();
			zBlock = prop->TunedBlockZ();
		}

	switch (pType) {
		case PropagatorOmelyan2:
		{
			prop = std::make_unique<PropOmelyan2>	(bck, field);
			break;
		}

		case PropagatorOmelyan4:
		{
			prop = std::make_unique<PropOmelyan4>	(bck, field);
			break;
		}

		case PropagatorLeapFrog2:
		{
			prop = std::make_unique<PropLeap>	(bck, field);
			break;
		}

		case PropagatorLeapFrog4:
		{
			prop = std::make_unique<PropMLeap>	(bck, field);
			break;
		}

		case PropagatorRKN4:
		{
			prop = std::make_unique<PropRKN4>	(bck, field);
			break;
		}
	}

	if (wasTuned) {
		prop->SetBlockX(xBlock);
		prop->SetBlockY(yBlock);
		prop->SetBlockZ(zBlock);
		prop->UpdateBestBlock();
	}

	//printf	("Propagator %ssuccessfully initialized", prop->Name().c_str());
	printf	("Propagator successfully initialized\n");
}

void	propagate	(const double dz)
{
	(prop->propagate)(dz);

	return;
}

#if 0
void	tunePropagator (Scalar *field) {
	// Hash CPU model so we don't mix different cache files
	LogMsg(VERB_NORMAL,"\n");
 	LogMsg(VERB_NORMAL,"[tp] Tune propagator!\n");

	int  myRank   = commRank();
	bool newFile  = false, found = false;

	if (prop == nullptr) {
		LogError("Error: propagator not initialized, can't be tuned.");
		return;
	}

	Profiler &prof = getProfiler(PROF_TUNER);

	std::chrono::high_resolution_clock::time_point start, end;
	size_t bestTime, lastTime, cTime;

	LogMsg (VERB_HIGH, "[tp] Started tuner");
	prof.start();

	if (field->Device() == DEV_CPU)
		prop->InitBlockSize(field->Length(), field->Depth(), field->DataSize(), field->DataAlign());
	else
		prop->InitBlockSize(field->Length(), field->Depth(), field->DataSize(), field->DataAlign(), true);

	/*	Check for a cache file	*/

	if (myRank == 0) {
		FILE *cacheFile;
		char tuneName[2048];
		// if (pType == PROP_BASE)
		// 	sprintf (tuneName, "%s/tuneCache.dat", wisDir);
		// else if (pType == PROP_NNEIG)
		// 	sprintf (tuneName, "%s/tuneCache%d.dat", wisDir,field->getNg());
		if (pType & PROP_BASE)
			sprintf (tuneName, "%s/tuneCache.dat", wisDir);
		else if (pType & PROP_NNEIG){
			LogMsg(VERB_HIGH,"%01d",field->getNg());
			sprintf (tuneName, "%s/tuneCache%01d.dat", wisDir,field->getNg());
			LogMsg(VERB_HIGH,"%s",tuneName);
			}
		if ((cacheFile = fopen(tuneName, "r")) == nullptr) {
LogMsg(VERB_HIGH,"[tp] new cache!!");
LogMsg (VERB_NORMAL, "Missing tuning cache file %s, will create a new one", tuneName);
			newFile = true;
		} else {
			int	     rMpi, rThreads;
			size_t       rLx, rLz, Nghost;
			unsigned int rBx, rBy, rBz, fType, myField ;
			if      (field->Field() == FIELD_SAXION)
				myField = 0;
			else if (field->Field() == FIELD_AXION)
				myField = 1;
			else if (field->Field() == FIELD_NAXION)
				myField = 2;
			else if (field->Field() == FIELD_PAXION)
				myField = 3;

			char	     mDev[8];

			std::string tDev(field->Device() == DEV_GPU ? "Gpu" : "Cpu");
LogMsg(VERB_HIGH,"[tp] Reading cache file %s",tuneName);
			do {
				fscanf (cacheFile, "%s %d %d %lu %lu %u %u %u %u %lu\n", reinterpret_cast<char*>(&mDev), &rMpi, &rThreads, &rLx, &rLz, &fType, &rBx, &rBy, &rBz, &Nghost);
				std::string fDev(mDev);
LogMsg(VERB_HIGH,"[tp] Read: MPI %d, threads %d, Ng %d, Lx,Lz (%d,%d) rBx,y,z (%d,%d,%d)  ",rMpi, rThreads, Nghost, rLx, rLz, rBx, rBy, rBz);
				if (rMpi == commSize() && rThreads == omp_get_max_threads() && rLx == field->Length() && rLz == field->Depth() && fType == myField && fDev == tDev && Nghost == field->getNg()) {
					if ((field->Device() == DEV_CPU && (rBx <= prop->BlockX() && rBy <= field->Surf()/prop->BlockX() && rBz <= field->Depth())) ||
					    (field->Device() == DEV_GPU	&& (rBx <= prop->MaxBlockX() && rBy <= prop->MaxBlockY() && rBz <= prop->MaxBlockZ()))) {
						found = true;
LogMsg(VERB_HIGH,"[tp] X!!");
						prop->SetBlockX(rBx);
LogMsg(VERB_HIGH,"[tp] Y!!");
						prop->SetBlockY(rBy);
LogMsg(VERB_HIGH,"[tp] Z!!");
						prop->SetBlockZ(rBz);
LogMsg(VERB_HIGH,"[tp] update best block!!");
						prop->UpdateBestBlock();
					}
				}
			}	while(!feof(cacheFile) && !found);

			fclose (cacheFile);
LogMsg(VERB_HIGH,"[tp] cache file closed!!");
		}
	}
LogMsg(VERB_HIGH,"[tp] BCAST!");

	MPI_Bcast (&found, sizeof(found), MPI_BYTE, 0, MPI_COMM_WORLD);

	commSync();

	// If a cache file was found, we broadcast the best block and exit
	if (found) {
LogMsg(VERB_HIGH,"[tp] optimum found!");
		unsigned int block[3];

		if (myRank == 0) {
			block[0] = prop->TunedBlockX();
			block[1] = prop->TunedBlockY();
			block[2] = prop->TunedBlockZ();
		}

		MPI_Bcast (&block, sizeof(int)*3, MPI_BYTE, 0, MPI_COMM_WORLD);
		commSync();

		if (myRank != 0) {
			prop->SetBlockX(block[0]);
			prop->SetBlockY(block[1]);
			prop->SetBlockZ(block[2]);
			prop->UpdateBestBlock();
		}

LogMsg (VERB_NORMAL, "Tuned values read from cache file. Best block %u x %u x %u", prop->TunedBlockX(), prop->TunedBlockY(), prop->TunedBlockZ());
LogMsg (VERB_HIGH,   "Chosen block %u x %u x %u\n", prop->BlockX(), prop->BlockY(), prop->BlockZ());
		prop->Tune();
		prof.stop();
		prof.add(prop->Name(), 0., 0.);
		return;
	}

	// Otherwise we start tuning

LogMsg (VERB_HIGH,   "[tp] Start tuning ... ");

	start = std::chrono::high_resolution_clock::now();
	propagate(field, 0.);
	end   = std::chrono::high_resolution_clock::now();

	cTime = std::chrono::duration_cast<std::chrono::nanoseconds>(end-start).count();

	// If there is an error in GPU propagation, we set the time to an absurd value
	#ifdef USE_GPU
	if (field->Device() == DEV_GPU) {
		auto gErr = cudaGetLastError();

		if (gErr != cudaSuccess)
			cTime = std::numeric_limits<std::size_t>::max();
	}
	#endif

	MPI_Allreduce(&cTime, &bestTime, 1, MPI_UNSIGNED_LONG_LONG, MPI_MAX, MPI_COMM_WORLD);

	if (field->Device() == DEV_GPU && cTime == std::numeric_limits<std::size_t>::max())
		LogMsg (VERB_HIGH, "Block %u x %u x %u gave an error and couldn't run on the GPU", prop->BlockX(), prop->BlockY(), prop->BlockZ());
	else
		LogMsg (VERB_HIGH, "Block %u x %u x %u done in %lu ns", prop->BlockX(), prop->BlockY(), prop->BlockZ(), bestTime);

	prop->AdvanceBlockSize();

	while (!prop->IsTuned()) {

		start = std::chrono::high_resolution_clock::now();
		propagate(field, 0.);
		end   = std::chrono::high_resolution_clock::now();

		cTime = std::chrono::duration_cast<std::chrono::nanoseconds>(end-start).count();

		#ifdef USE_GPU
		if (field->Device() == DEV_GPU) {
			auto gErr = cudaGetLastError();

			if (gErr != cudaSuccess)
				cTime = std::numeric_limits<std::size_t>::max();
		}
		#endif

    MPI_Allreduce(&cTime, &lastTime, 1, MPI_UNSIGNED_LONG_LONG, MPI_MAX, MPI_COMM_WORLD);

		if (field->Device() == DEV_GPU && cTime == std::numeric_limits<std::size_t>::max())
			LogMsg (VERB_HIGH, "Block %u x %u x %u gave an error and couldn't run on the GPU", prop->BlockX(), prop->BlockY(), prop->BlockZ());
		else{
			LogMsg (VERB_HIGH, "Test Block %u x %u x %u done in %lu ns", prop->BlockX(), prop->BlockY(), prop->BlockZ(), lastTime);
			LogMsg (VERB_HIGH, "Best Block %u x %u x %u done in %lu ns", prop->TunedBlockX(), prop->TunedBlockY(), prop->TunedBlockZ(), bestTime);
		}

		if (lastTime < bestTime) {
			bestTime = lastTime;
			prop->UpdateBestBlock();
			LogMsg (VERB_HIGH, "Best block updated");
		}

		prop->AdvanceBlockSize();
	}

	prop->getBaseName();

	char loli[2048];

	switch (field->Field()) {
		case FIELD_SAXION:
			if (pType & (PROP_BASE | PROP_NNEIG)){
				sprintf (loli, "N %01d Ng %01d Saxion", Nng, field->getNg());
				prop->appendName(loli);
			} else {
				prop->appendName("Saxion");
			}
			break;

		case FIELD_AXION:
		if (pType & (PROP_BASE | PROP_NNEIG)){
				sprintf (loli, "N %01d Ng %01d Axion", Nng, field->getNg());
				prop->appendName(loli);
			} else {
				prop->appendName("Axion");
			}
			break;

		case FIELD_AXION_MOD:
			if (pType & (PROP_BASE | PROP_NNEIG)){
				sprintf (loli, "N %01d Ng %01d Axion Mod", Nng, field->getNg());
				prop->appendName(loli);
			} else {
				prop->appendName("Axion Mod");
			}
			break;

		case FIELD_NAXION:
			if (pType & (PROP_BASE | PROP_NNEIG)){
					sprintf (loli, "N %01d Ng %01d Naxion", Nng, field->getNg());
					prop->appendName(loli);
				} else {
					prop->appendName("Naxion");
				}
				break;

		case FIELD_PAXION:
			if (pType & (PROP_BASE | PROP_NNEIG)){
					sprintf (loli, "N %01d Ng %01d Paxion", Nng, field->getNg());
					prop->appendName(loli);
				} else {
					prop->appendName("Paxion");
				}
				break;
		default:
			LogError ("Error: invalid field type");
			prof.stop();
			return;
	}

	Profiler &propProf = getProfiler(PROF_PROP);
	propProf.reset(prop->Name());

	prop->SetBestBlock();
	LogMsg (VERB_NORMAL, "Propagator tuned! Best block %u x %u x %u in %lu ns", prop->TunedBlockX(), prop->TunedBlockY(), prop->TunedBlockZ(), bestTime);

	/*	Write cache file if necessary, block of rank 0 prevails		*/

	if (myRank == 0) {
		FILE *cacheFile;
		char tuneName[2048];
		// sprintf (tuneName, "%s/tuneCache.dat", wisDir);
		if (pType & PROP_BASE)
			sprintf (tuneName, "%s/tuneCache.dat", wisDir);
		else if (pType & PROP_NNEIG){
			LogMsg(VERB_HIGH,"%01d",field->getNg());
			sprintf (tuneName, "%s/tuneCache%01d.dat", wisDir,field->getNg());
			LogMsg(VERB_HIGH,"%s",tuneName);
			}

		// We distinguish between opening and appending a new line
		if (!newFile) {
			if ((cacheFile = fopen(tuneName, "a")) == nullptr) {
				LogError ("Error: can't open cache file, can't save tuning results");
				commSync();
				prof.stop();
				prof.add(prop->Name(), 0., 0.);
			}
		} else {
			if ((cacheFile = fopen(tuneName, "w")) == nullptr) {
				LogError ("Error: can't create cache file, can't save tuning results");
				commSync();
				prof.stop();
				prof.add(prop->Name(), 0., 0.);
			}
		}

		unsigned int fType ;
		if      (field->Field() == FIELD_SAXION)
			fType = 0;
		else if (field->Field() == FIELD_AXION)
			fType = 1;
		else if (field->Field() == FIELD_NAXION)
			fType = 2;
		else if (field->Field() == FIELD_PAXION)
			fType = 3;

		std::string myDev(field->Device() == DEV_GPU ? "Gpu" : "Cpu");
		fprintf (cacheFile, "%s %d %d %lu %lu %u %u %u %u %lu\n", myDev.c_str(), commSize(), omp_get_max_threads(), field->Length(), field->Depth(),
									fType, prop->TunedBlockX(), prop->TunedBlockY(), prop->TunedBlockZ(), field->getNg());
		fclose  (cacheFile);
	}
LogMsg (VERB_NORMAL, "\n");

	commSync();
	prof.stop();
	prof.add(prop->Name(), 0., 0.);
}
#endif
