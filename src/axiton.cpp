#include <cmath>
#include <cstring>
#include <chrono>
#include <iostream>
#include <vector>

#include "axiton.h"
#include <cuda.h>
#include <cuda_runtime.h>
#include "cudaErrors.h"

#include <fftw3.h>

void    printsample  (FILE *fichero, Axiton *axiton, size_t *idxprint, int npoints, double mA, double delta);
void    loadsample(int N, int n, double *fa, double *fra, Axiton *axiton, size_t *idxprint, int Nidxprint, double delta);
void    printfftsample(int N, FILE *file_samp, int Nidxprint, double *in, fftw_complex *out, fftw_plan p);
void    addfftsample(int N, FILE *file_samp, int Nidxprint, double *fa, double *fra, double *in, fftw_complex *out, fftw_plan p, double *Ak, double *Amk);
void    saveps(int N, FILE *file_samp, int Nidxprint, double *Ak, double *Amk);
double  getR(Axiton *axiton, iParms &myParms);

using namespace std;

int	main (int argc, char *argv[])
{
	Cosmos myCosmos = initAxitons(argc, argv);
	iParms &myParms = myCosmos.InitParms();

	Axiton		*axiton;

	axiton = new Axiton (myParms.nSize, myParms.fPrec);

	Hdf5ReadWriter	IOHandler(myParms);
	Generator	AxitonFactory(myParms.cType, axiton, &myCosmos);
	

	/* Extra output */
	FILE *file_samp ;
	char out3Name[2048];
	sprintf (out3Name, "sample.txt");
	file_samp = fopen(out3Name,"w+");

	int Nidxprint = 4;
	size_t *idxprint = (size_t *) malloc(Nidxprint*sizeof(size_t));

	/* FFTs hardcoded length -openMP? */
	size_t N = 1024;
	// for the FFT
	double *in  = (double *) fftw_malloc(N*sizeof(double));
	fftw_complex *out = (fftw_complex *) fftw_malloc((N/2+1)*2*sizeof(fftw_complex));
	// for the field
	double *fa   = (double *) fftw_malloc(N*(Nidxprint)*sizeof(double));
	// for the gradient
	double *fra  = (double *) fftw_malloc(N*(Nidxprint)*sizeof(double));
	// to accumulate Ak, A-k
	double *Ak  = (double *) fftw_malloc((N/2+1)*(Nidxprint)*sizeof(double));
	double *Amk = (double *) fftw_malloc((N/2+1)*(Nidxprint)*sizeof(double));
	
	fftw_plan p = fftw_plan_dft_r2c_1d(N, in, out, FFTW_ESTIMATE);
	FILE *file_fft ;
	char out4Name[2048];
	sprintf (out4Name, "samplefft.txt");
	file_fft = fopen(out4Name,"w+");
	memset(Ak, 0, (N/2+1)*(Nidxprint)*sizeof(double));
	memset(Amk, 0, (N/2+1)*(Nidxprint)*sizeof(double));

	FILE *file_info;
	char out5Name[2048];
	sprintf (out5Name,"info.txt");
	file_info = fopen(out5Name,"w+");

	double R;
	switch (myParms.fExp){
		case Radiation:
	       	R = axiton->R<Radiation>();
		break;
		case Minkowski:
		R = axiton->R<Minkowski>();
		break;
	}
	//auto R   = axiton->R<myParms.fExp>();
	printf (" - z = %f R = %f o2 = %e\n",  axiton->z(), R, 1.0/myCosmos.Delta());
	AxitonFactory.Construct(myParms.parm1*myParms.zInit, myParms.parm2, myParms.zInit);

	printf("Current file %d",IOHandler.currentIndex());
	IOHandler.writeConf(&myCosmos, axiton);

	printf ("Transferring configuration to device\n"); fflush(stdout);
	axiton->transferField(FieldBase | FieldDev, HostToDevice);
	CudaCheckError();

	double delta = myCosmos.Delta();
	double dz;
	double Ri;
	switch (myParms.fExp){
		case Radiation:
			Ri = axiton->R<Radiation>();
			break;
		case Minkowski:
		    	Ri = axiton->R<Minkowski>();
			break;
	}
        auto maai = myCosmos.AxionMassSq(Ri);

	if (myParms.wDz == 0.0){
		if (myParms.nSteps == 0)
			dz = 0.;
		else
			dz = (myParms.zEnd - myParms.zInit)/((double) myParms.nSteps);
		printf ("User specified time step (dz = %.2e)\n",dz); fflush(stdout);
	} else {
		dz = myParms.wDz/sqrt(maai*Ri*Ri + 12./(delta*delta));
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

	
	size_t rim = (size_t) (3.0/sqrt(maai)/delta);
	
	idxprint[0] = 2;
	idxprint[1] = rim;
	idxprint[2] = axiton->Size()/50;
	idxprint[3] = axiton->Size()/2;
	
	for (int i = 0; i < Nidxprint;i++)
	{
		printf("idxprint %d %d \n",i,idxprint[i]);
		fprintf(file_info,"%d ",idxprint[i]);
	}
	if (0)
		printsample(file_samp, axiton, idxprint, Nidxprint,sqrt(maai),delta);

	
	int n = 0;
	bool endLoop = false;
	for (int zloop = 0; zloop < nLoops; zloop++)
	{
		if (endLoop){
			propagate (dz);
			break;
		}

		for (int zsubloop = 0; zsubloop < dump; zsubloop++) {
			propagate (dz);
			// to print in sample
			if(0) 
			{
				auto maa = myCosmos.AxionMassSq(R);
				double R = getR(axiton,myParms);
				idxprint[1] =  (size_t) (3.0/sqrt(maa)/(delta*R));
				//printf("%d %d idxprintrim %d\n",zloop,zsubloop,idxprint[1]);
				printsample(file_samp, axiton, idxprint, Nidxprint,sqrt(maa), delta);
			}
			if (0)
			{
				//printf("n=%d",n);
				loadsample(N, n, fa, fra, axiton, idxprint, Nidxprint, delta);
				n++; 
			
				if (n%N == 0){
					printf("outputfft n %d N %d...", n, N);
					//printfftsample(N, file_fft, Nidxprint, in, out, p);
					addfftsample(N, file_fft, Nidxprint, fa, fra, in, out, p, Ak, Amk);
					n = 0;
					printf("done!\n");
			}
			}
		}
		//printf("Current file %d",IOHandler.currentIndex());
		axiton->transferField(FieldBase | FieldDev, DeviceToHost);
		IOHandler.nextFile();
		IOHandler.writeConf(&myCosmos, axiton);
		double R = getR(axiton,myParms);
		//switch (myParms.fExp){
		//	case Radiation:
		//	R = axiton->R<Radiation>();
		//	break;
		//	case Minkowski:
		//	R = axiton->R<Minkowski>();
		//	break;
		//}
		//auto R   = axiton->R<Radiation>();
		auto maa = myCosmos.AxionMassSq(R);
		printf (" - z = %f (dz %.2e) R = %f mA = %e o2 = %e\n",  axiton->z(), dz, R, sqrt(maa),1.0/myCosmos.Delta());

		if (myParms.wDz != 0.0)
		{
			// printf ("!!!! MASS ALARM !!!! \n");
			dz = myParms.wDz/sqrt(maa*R*R +12/delta/delta);
		}


		if (axiton->z() + dz > myParms.zEnd)
			{
				dz = myParms.zEnd - axiton->z();
				endLoop = true;
			}
	}

	saveps(N, file_fft, Nidxprint, Ak, Amk);
	current = std::chrono::high_resolution_clock::now();
	elapsed = std::chrono::duration_cast<std::chrono::milliseconds>(current - start);

	axiton->transferField(FieldBase | FieldDev, DeviceToHost);
	IOHandler.nextFile();
	IOHandler.writeConf(&myCosmos, axiton);
	
	free(idxprint);
	fftw_destroy_plan(p);
    	fftw_free(in); 
	fftw_free(out);
	//fftw_free(fa); 
	//fftw_free(fra); 
	//fftw_free(Ak); 
	//fftw_free(Amk); 


	printf ("Final z  = %f\n", axiton->z());
	printf ("#_steps  = %i\n", nLoops*dump);
	printf ("#_prints = %i\n", IOHandler.currentIndex());
	printf ("Total time: %2.3f s\n", elapsed.count()*1.e-3);

	fclose(file_samp);
	fclose(file_fft);
	fclose(file_info);
	delete axiton;

	endAxitons();

	return 0;
}

void printsample(FILE *file_samp, Axiton *axiton, size_t *idxprint, int Nidxprint, double mA, double delta)
{
	double buff[3];
	fprintf(file_samp,"%f %f ",axiton->z(), mA);

	for (int i=0;i<Nidxprint;i++)
	{
		//printf("cosa %d %d\n",i,idxprint[i]);
		cudaMemcpy(&(buff[0]),&(static_cast<double*>(axiton->fieldGpu())[idxprint[i]-1]),3*sizeof(double),cudaMemcpyDeviceToHost);
		//cudaMemcpy(&(buff[1]),&(static_cast<double*>(axiton->fieldGpu())[idxprint[i]-1]),3*sizeof(double),cudaMemcpyDeviceToHost);
		// derivative of a\times r required
		fprintf(file_samp,"%.10e %.10e ",buff[1],((idxprint[i]+1)*buff[2]-(idxprint[i]-1)*buff[0])/(2.));
	}
	fprintf(file_samp," \n");
	fflush(file_samp);
}

void loadsample(int N, int n, double *fa, double *fra, Axiton *axiton, size_t *idxprint, int Nidxprint, double delta)
{
	double buff[3];
	for (int i=0;i<Nidxprint;i++)
	{
		cudaMemcpy(&(fa[i*N+n]),&(static_cast<double*>(axiton->fieldGpu())[idxprint[i]]),sizeof(double),cudaMemcpyDeviceToHost);
		cudaMemcpy(&(buff[0]),&(static_cast<double*>(axiton->devGpu())[idxprint[i]-1]),3*sizeof(double),cudaMemcpyDeviceToHost);
		// mises the 1/delta!
		fra[i*N+n] = (buff[2]-buff[0])/(2.*delta);
	}
}

void printfftsample(int N, FILE *file_samp, int Nidxprint, double *fa, double *fra, double *in, fftw_complex *out, fftw_plan p)
{
	for (int i=0;i<Nidxprint;i++)
	{
		// load field
		memcpy(in,fa,N*sizeof(double));
		fftw_execute(p);
		// print
		for (int j = 0;j<N/2+1;j++)
			fprintf(file_samp,"%e %e ",out[j][0],out[j][1]);
		fprintf(file_samp,"\n");
		// load derivative
		memcpy(in,fra,N*sizeof(double));
		fftw_execute(p);
		for (int j = 0;j<N/2+1;j++)
			fprintf(file_samp,"%e %e ",out[j][0],out[j][1]);
		fprintf(file_samp,"\n");
	}
}

void addfftsample(int N, FILE *file_samp, int Nidxprint, double *fa, double *fra, double *in, fftw_complex *out, fftw_plan p, double *Ak, double *Amk)
{
	size_t Nh = N/2+1;
	for (int i=0;i<Nidxprint;i++)
	{
		// load derivative field
		memcpy(in,fra,N*sizeof(double));
		fftw_execute(p);
		memcpy(&(out[(N/2+1)]),out,(N/2+1)*sizeof(fftw_complex));
		
		// load field 
		memcpy(in,fa,N*sizeof(double));
		fftw_execute(p);
		
		// calculate Ak Amk
		// Ak = [ikr fft(a) + fft(ra)]
		// Amk = [ikr fft(a) - fft(ra)]
		// call fft(a)  = a + ib
		//      fft(rs) = c + id
		for (int j = 0;j<N/2+1;j++){
			Ak[i*Nh+j] += out[j][0]*out[j][0] + out[j][1]*out[j][1];
			Amk[i*Nh+j] += out[Nh+j][0]*out[Nh+j][0] + out[Nh+j][1]*out[Nh+j][1];
		}
	}

}


void saveps(int N, FILE *file_samp, int Nidxprint, double *Ak, double *Amk)
{
	for (int i=0;i<Nidxprint*2;i++)
	{
		for (int j = 0;j<N/2+1;j++)
			fprintf(file_samp,"%e ",Ak[i*(N/2+1)+j]);
		fprintf(file_samp,"\n");
	for (int j = 0;j<N/2+1;j++)
			fprintf(file_samp,"%e ",Amk[i*(N/2+1)+j]);
		fprintf(file_samp,"\n");
	}
}



double getR(Axiton *axiton, iParms &myParms)
{
	switch (myParms.fExp){
	case Radiation:
	return axiton->R<Radiation>();
	break;
	case Minkowski:
	return axiton->R<Minkowski>();
	break;
	}
}

