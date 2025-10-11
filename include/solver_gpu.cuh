
#ifndef SOLVER_GPU_CUH
#define SOLVER_GPU_CUH

float *gsSolverIsoGPU(struct solver_metadata metadata, const std::vector<double>* uFlux_vector);
float *gsSolverHarmonicGPU(struct solver_metadata metadata, const std::vector<double>* uFlux_vector);
float *getufluxArray(const std::vector<double>* uflux_vector, int num_energies, int num_voxels);
float *raytraceIsoGPU(struct solver_metadata metadata);
struct solver_metadata *transferMetadataToDevice(struct solver_metadata metadata);

#endif
