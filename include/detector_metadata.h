
#ifndef DETECTOR_METADATA_H
#define DETECTOR_METADATA_H

#define NV_VOXEL_DENSITY_THRESHOLD 0.005 // threshold in g/cm^3 for being considered a non-vaccum voxel

struct detector_metadata{

  // Detector-array parameters
  float det_lowerX, det_lowerZ;
  float det_upperX, det_upperZ;
  float det_Dx, det_Dz;
  float det_cenY;
  float det_mag;
  int det_Nx, det_Nz, num_det;

  // Mesh-grid parameters
  int mesh_Nx, mesh_Ny, mesh_Nz, num_voxels;
  float mesh_Dx, mesh_Dy, mesh_Dz;

  // Source position parameters
  float source_cenX, source_cenY, source_cenZ;

  // List of non-void voxels in the mesh
  int *nv_voxels;
  int num_nv_voxels;

  // number of energies (needed since ray-trace is energy-wise)
  int num_energies;

  // spectrum (needed for fringe detector-trace case where no voxels project onto given detector pixel)
  float *sIntensity;

};

#endif
