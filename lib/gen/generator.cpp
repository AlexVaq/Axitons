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

		fillArray<float>  (field->fieldCpu(), field->Size(), [&] (size_t x) -> float { return (float) value; } );
		fillArray<float>  (field->devCpu  (), field->Size(), [&] (size_t x) -> float { return 0.f; });
		fillArray<float>  (field->miscCpu (), field->Size(), [&] (size_t x) -> float { return 0.f; });
		break;

		case DoublePrecision:

		fillArray<double> (field->fieldCpu(), field->Size(), [&] (size_t x) -> double { return value; } );
		fillArray<double> (field->devCpu  (), field->Size(), [&] (size_t x) -> double { return 0.; });
		fillArray<double> (field->miscCpu (), field->Size(), [&] (size_t x) -> double { return 0.; });
		break;
	}
}

void	fillSinc2(Axiton *field, double value, int coef) {

	switch (field->Precision()) {
		case SinglePrecision:

		fillArray<float>  (field->fieldCpu(), field->Size(), [&] (size_t x) -> float { return (x == 0) ? (float) value : (float) value*(sin(coef*x)/(coef*x)); } );
		fillArray<float>  (field->devCpu  (), field->Size(), [&] (size_t x) -> float { return 0.f; });
		fillArray<float>  (field->miscCpu (), field->Size(), [&] (size_t x) -> float { return 0.f; });
		break;

		case DoublePrecision:

		fillArray<double> (field->fieldCpu(), field->Size(), [&] (size_t x) -> double { (x == 0) ? value : value*(sin(coef*x)/(coef*x)); } );
		fillArray<double> (field->devCpu  (), field->Size(), [&] (size_t x) -> double { return 0.; });
		fillArray<double> (field->miscCpu (), field->Size(), [&] (size_t x) -> double { return 0.; });
		break;
	}
}

void	Generator::Construct (double parm1, int parm2) {

	switch (icType) {

		case IcFlat:
		fillFlat (field, parm1);
		break;

		case IcSinc2:
		fillSinc2(field, parm1, parm2);
		break;

		case IcNone:
		default:
		break;
	}
}
