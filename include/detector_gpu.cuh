
#ifndef DETECTOR_GPU_CUH
#define DETECTOR_GPU_CUH

// TODO: Use separate data structures for footprint (FP) and non-FP method
struct rayDetInt{
  // detector-index that ray hits (incident detector-pixel)
  int idx_det;
  // distance between point of intersection and detector pixel-center
  float dist;

  // displacement \delta between voxel-center projection (mid-point of footprint) and (incident) detector-pixel
  // this is NOT absolute valued
  float delta_x;
  float delta_z;
  // voxel-footprint width
  float W_x;
  float W_z;
  // range of detector-pixels indices that overlap with voxel-footprint
  int jx_min, jx_max;
  int jz_min, jz_max;

};

// Direct transmission flux delivered to the detector
// version based on voxel-detector footprint (FP)
float *getPrimaryImg_GPU(struct detector_metadata metadata, float* uflux);
// REMOVE LATER: non-footprint (non-FP) version
float *getPrimaryImg_nonFP_GPU(struct detector_metadata metadata, float* uflux); 

// Compute first-order and Multi-order (2nd and higher order) flux delivered to the detector
float **getScatterImgIso_GPU(float* uFlux, float* cFlux, struct solver_metadata metadata_s, struct detector_metadata metadata_d);
// Compute first-order and Multi-order (2nd and higher order) flux delivered to the detector
float **getScatterImgAnIso_GPU(float* uFlux, float* cFlux, struct solver_metadata metadata_s, struct detector_metadata metadata_d);

#endif
