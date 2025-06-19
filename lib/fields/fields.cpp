#include "fields/fields.h"
#include "utils/utils.h"

#include <stdio.h>

#include <cuda.h>
#include <cuda_runtime.h>
#include "cudaErrors.h"


Axiton::Axiton (int len, FieldPrecision pr, FieldType fType) : length(len), prec(pr), bytes(len * pr), fMod(fType) {

    printf("Constructing fields...\n");
    fflush(stdout);

    trackAlloc (&hField, bytes);
    trackAlloc (&hDev,   bytes);
    trackAlloc (&hMisc,  bytes);

    deviceAlloc(&dField, bytes);
    deviceAlloc(&dDev,   bytes);
    deviceAlloc(&dMisc,  bytes);

    fStatus = FieldUndefined;
    dStatus = FieldUndefined;
    zVar	= 0.0;

    printf("Field constructed\n");
    fflush(stdout);
}

Axiton::~Axiton () {

    trackFree (dMisc);
    trackFree (dDev);
    trackFree (dField);

    trackFree (hMisc);
    trackFree (hDev);
    trackFree (hField);
}


void	Axiton::transferField(FieldIndex fIdx, MemDirection mIdx) {
    if (fIdx & FieldBase) {
        if	  ((mIdx == DeviceToHost) && (fStatus & FieldGpu)) {
            cudaMemcpy(hField, dField,  bytes, cudaMemcpyDeviceToHost);
            fStatus |= FieldCpu;
        } else if ((mIdx == HostToDevice) && (fStatus & FieldCpu)) {
            cudaMemcpy(dField, hField,  bytes, cudaMemcpyHostToDevice);
            fStatus |= FieldGpu;
        }
    }

    if (fIdx & FieldDev) {
        if	  ((mIdx == DeviceToHost) && (dStatus & FieldGpu)) {
            cudaMemcpy(hDev,   dDev,    bytes, cudaMemcpyDeviceToHost);
            dStatus |= FieldCpu;
        } else if ((mIdx == HostToDevice) && (dStatus & FieldCpu)) {
            cudaMemcpy(dDev,   hDev,    bytes, cudaMemcpyHostToDevice);
            dStatus |= FieldGpu;
        }
    }

    if (fIdx & FieldMisc) {
        if	  (mIdx == DeviceToHost) {
            cudaMemcpy(hMisc,  dMisc,   bytes, cudaMemcpyDeviceToHost);
        } else if (mIdx == HostToDevice) {
            cudaMemcpy(dMisc,  hMisc,   bytes, cudaMemcpyHostToDevice);
        }
    }
}

