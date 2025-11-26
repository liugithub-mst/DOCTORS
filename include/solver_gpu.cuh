
#ifndef SOLVER_GPU_CUH
#define SOLVER_GPU_CUH
#include "macros.h"

SOL_T* gsSolverIsoGPU(struct solver_metadata metadata, const std::vector<double>* uFlux_vector);
SOL_T* gsSolverHarmonicGPU(struct solver_metadata metadata, const std::vector<double>* uFlux_vector);
float* getufluxArray(const std::vector<double>* uflux_vector, int num_energies, int num_voxels);
RAY_T* raytraceIsoGPU(struct solver_metadata metadata);
struct solver_metadata *transferMetadataToDevice(struct solver_metadata metadata);

#endif
