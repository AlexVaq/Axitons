#include "fields/fields.h"
#include "utils/utils.h"

#include <stdio.h>

#include <cuda.h>
#include <cuda_runtime.h>
#include "cudaErrors.h"


	Axiton::Axiton (int len, FieldPrecision pr) : length(len), prec(pr), bytes(len * pr) {

	printf("Constructing fields...\n"); fflush(stdout);

	trackAlloc (&hField, bytes);
	trackAlloc (&hDev,   bytes);
	trackAlloc (&hMisc,  bytes);

	deviceAlloc(&dField, bytes);
	deviceAlloc(&dDev,   bytes);
	deviceAlloc(&dMisc,  bytes);

	printf("Field constructed\n"); fflush(stdout);
}

	Axiton::~Axiton () {

	trackFree (dMisc);
	trackFree (dDev);
	trackFree (dField);

	trackFree (hMisc);
	trackFree (hDev);
	trackFree (hField);
}


void	Axiton::transferField(FieldIndex fIdx, MemDirection mIdx)
{
	if (fIdx & FieldBase) {
		if	(mIdx == DeviceToHost) 
			cudaMemcpy(hField, dField,  bytes, cudaMemcpyDeviceToHost);
		else if (mIdx == HostToDevice)
			cudaMemcpy(dField, hField,  bytes, cudaMemcpyHostToDevice);
	}

	if (fIdx & FieldDev) {
		if	(mIdx == DeviceToHost) 
			cudaMemcpy(hDev,   dDev,    bytes, cudaMemcpyDeviceToHost);
		else if (mIdx == HostToDevice)
			cudaMemcpy(dDev,   hDev,    bytes, cudaMemcpyHostToDevice);
	}

	if (fIdx & FieldMisc) {
		if	(mIdx == DeviceToHost) 
			cudaMemcpy(hMisc,  dMisc,   bytes, cudaMemcpyDeviceToHost);
		else if (mIdx == HostToDevice)
			cudaMemcpy(dMisc,  hMisc,   bytes, cudaMemcpyHostToDevice);
	}
}

