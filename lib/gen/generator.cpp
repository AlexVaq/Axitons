#include<cstdlib>
#include<cstring>

#include <cmath>
#include <functional>

#include "enum-vars.h"
#include "gen/generator.h"
#include "fields/fields.h"

template<typename Float>
void	fillArray(void	*array, size_t size, std::function<Float(size_t)> func) {

	Float *data = static_cast<Float*>(array);

	#pragma omp parallel for schedule(static)
	for (size_t i=0; i<size; i++)
		data[i] = func(i);
}

void	fillFlat (Axiton *field, double value) {

	switch (field->Precision()) {
		case SinglePrecision:

		// fillArray<float>  (field->fieldCpu(), field->Size(), [&] (size_t x) -> float  { return (x == 0) ? (float) value : ((float)  value)/(((float)  x)*((float)  x)); } );
		fillArray<float>  (field->fieldCpu(), field->Size(), [&] (size_t x) -> float  { return (float) value; });
		fillArray<float>  (field->devCpu  (), field->Size(), [&] (size_t x) -> float  { return 0.f; });
		fillArray<float>  (field->miscCpu (), field->Size(), [&] (size_t x) -> float  { return 0.f; });
		break;

		case DoublePrecision:

		fillArray<double> (field->fieldCpu(), field->Size(), [&] (size_t x) -> double { return value; } );
		fillArray<double> (field->devCpu  (), field->Size(), [&] (size_t x) -> double { return 0.; });
		fillArray<double> (field->miscCpu (), field->Size(), [&] (size_t x) -> double { return 0.; });
		break;
	}
}

void	fillOvr2 (Axiton *field, double value) {

	switch (field->Precision()) {
		case SinglePrecision:

		fillArray<float>  (field->fieldCpu(), field->Size(), [&] (size_t x) -> float  { return (x == 0) ? (float) value : ((float)  value)/(((float)  x)*((float)  x)); } );
		fillArray<float>  (field->devCpu  (), field->Size(), [&] (size_t x) -> float  { return 0.f; });
		fillArray<float>  (field->miscCpu (), field->Size(), [&] (size_t x) -> float  { return 0.f; });
		break;

		case DoublePrecision:

		fillArray<double> (field->fieldCpu(), field->Size(), [&] (size_t x) -> double { return (x == 0) ? value : value/(((double)  x)*((double)  x)); } );
		fillArray<double> (field->devCpu  (), field->Size(), [&] (size_t x) -> double { return 0.; });
		fillArray<double> (field->miscCpu (), field->Size(), [&] (size_t x) -> double { return 0.; });
		break;
	}
}

void	fillSinc2(Axiton *field, double value, int coef) {

	double	cf = ((double) coef)/((double) field->Size());

	switch (field->Precision()) {
		case SinglePrecision:

		fillArray<float>  (field->fieldCpu(), field->Size(), [&] (size_t x) -> float  { return (x == 0) ? (float) value : ((float)  value)*(sin(cf*x)*sin(cf*x)/(cf*cf*x*x)); } );
		fillArray<float>  (field->devCpu  (), field->Size(), [&] (size_t x) -> float  { return 0.f; });
		fillArray<float>  (field->miscCpu (), field->Size(), [&] (size_t x) -> float  { return 0.f; });
		break;

		case DoublePrecision:

		fillArray<double> (field->fieldCpu(), field->Size(), [&] (size_t x) -> double { return (x == 0) ? value: value*(sin(cf*x)*sin(cf*x)/(cf*cf*x*x)); } );
		fillArray<double> (field->devCpu  (), field->Size(), [&] (size_t x) -> double { return 0.; });
		fillArray<double> (field->miscCpu (), field->Size(), [&] (size_t x) -> double { return 0.; });
		break;
	}
}

void	fillOvX(Axiton *field, double value, int coef) {

	double	cf = ((double) coef)/((double) field->Size());

	switch (field->Precision()) {
		case SinglePrecision:

		fillArray<float>  (field->fieldCpu(), field->Size(), [&] (size_t x) -> float  { return ((float)  value)/(1.0 + ((float) cf)*x); } );
		fillArray<float>  (field->devCpu  (), field->Size(), [&] (size_t x) -> float  { return 0.f; });
		fillArray<float>  (field->miscCpu (), field->Size(), [&] (size_t x) -> float  { return 0.f; });
		break;

		case DoublePrecision:

		fillArray<double> (field->fieldCpu(), field->Size(), [&] (size_t x) -> double { return value/(1.0 + cf*x); } );
		fillArray<double> (field->devCpu  (), field->Size(), [&] (size_t x) -> double { return 0.; });
		fillArray<double> (field->miscCpu (), field->Size(), [&] (size_t x) -> double { return 0.; });
		break;
	}
}

void	Generator::fillGen() {


	switch (field->Precision()) {
		case SinglePrecision:

		fillArray<float>  (field->fieldCpu(), field->Size(), [&] (size_t x) -> float  { return ((float) sm(x*cosmo->Delta())); } );
		fillArray<float>  (field->devCpu  (), field->Size(), [&] (size_t x) -> float  { return ((float) sv(x*cosmo->Delta())); } );
		fillArray<float>  (field->miscCpu (), field->Size(), [&] (size_t x) -> float  { return 0.f; });
		break;

		case DoublePrecision:

		fillArray<double> (field->fieldCpu(), field->Size(), [&] (size_t x) -> double { return sm(x*cosmo->Delta()); } );
		fillArray<double> (field->devCpu  (), field->Size(), [&] (size_t x) -> double { return sv(x*cosmo->Delta()); });
		fillArray<double> (field->miscCpu (), field->Size(), [&] (size_t x) -> double { return 0.; });
		break;
	}
}

void	Generator::Construct (double parm1, int parm2, double zInit) {

	switch (icType) {

		case IcFlat:
		fillFlat (field, parm1);
		break;

		case IcOvr2:
		fillOvr2 (field, parm1);
		break;

		case IcOvX:
		fillOvX  (field, parm1, parm2);
		break;

		case IcSinc2:
		fillSinc2(field, parm1, parm2);
		break;

		case IcGen:
			SplineSetup();
			fillGen();
		break;

		case IcNone:
		default:
		break;
	}

	field->setStatus    (FieldBaseDev, FieldCpu);
	field->transferField(FieldBaseDev, HostToDevice);
	field->zUpdate	    (zInit);
}





//##############################################################################

void Generator::SplineSetup()
{

  /*Read Cosmology*/

  // char cosName[2048];
  // if (const char *cosPath = std::getenv("AXITONS_DIR")) {
  //   if (strlen(cosPath) < 1022) {
  //     struct stat tStat;
  //     if (stat(cosPath, &tStat) == 0 && S_ISDIR(tStat.st_mode)) {
  //       strcpy(cosName, cosPath);
  //     } else {
  //       printf("Path %s doesn't exist, using default\n", cosPath);
  //     }
  //   }
  // }
  // sprintf(cosName, "%s%s", cosName,  "./ics.txt");
  std::vector<double>	xV, mV, vV;
  FILE *cFile = nullptr;
  if (((cFile  = fopen("./ics.txt", "r")) == nullptr)){
    return ;
  }
  else
  {
    double x, m, v;

    char buffer[200];
    fgets(buffer, 200, cFile); // reads header
    fgets(buffer, 200, cFile); //reads 2nd line space

    fscanf (cFile ,"%lf %lf %lf", &x, &m, &v);
    while(!feof(cFile)){
      xV.push_back(x);
      mV.push_back(m);
      vV.push_back(v);
    	fscanf (cFile ,"%lf %lf %lf", &x, &m, &v);
    }
  }

	// for (int i; i < xV.size();i++){
  //   xV[i] *= xV.back();
	// 	xV[i] *= xV.back();
  // }

  sm.set_points(xV,mV);
	sv.set_points(xV,vV);
	printf("user defined ICs\n");
}
