
#ifndef SOLVER_METADATA_H
#define SOLVER_METADATA_H

struct solver_metadata
{
  int num_materials;
  int num_angles;
  int num_voxels;
  int num_energies;
  int pn;
  float voxel_vol;
  float voxel_area_yz, voxel_area_xz, voxel_area_xy;
  int Nx, Ny, Nz;
  float sx, sy, sz;                      // XL: source position
  float *sIntensity;                     // XL: source strength
  float *XNodes, *YNodes, *ZNodes;       // XL: voxel position
  float sphi, stheta;                    // XL: source fan or cone angle (degrees)
  float isox, isoy, isoz;                // XL: rotation center position
  int stype;                             // XL: source type

  // int *zoneId; //size: #voxels
  int *zoneIdMapped; //size: #voxels

  // We can compute Ayz, Axz, Axy nd quad_d can on the fly.
  float *mu, *zi, *eta; //size: #angle-directions for scatter
  float *quad_wt; // size: #angle-directions for scatter

  /*
  float *Ayz, *Axz, *Axy;  // size: #angle-directions for scatter
  int *quad_d;    // size: #angle-directions for scatter
  */

  float *atomDensity; // size: #voxels
  float *sigma_s; //size: #materials * #energy-groups * #energy-groups * #legendre-coefficients
  float *sigma_total; //size:  #materials * #energy-groups
};

struct sweep_pattern
{
  int d; //0-7
  int num_stages; // #sub-sweeps
  int num_voxels; // #voxels in entire 3D sweep
  int *num_voxels_per_stage;
  int *starting_index_per_stage;
  int *voxel_list;
};

#endif
