
#ifndef SWEEP_H
#define SWEEP_H

class sweep
{
public:
  int dix; // can be +1 or -1
  int diy; // can be +1 or -1
  int diz; // can be +1 or -1
  int Nx;  // num_voxels in x direction
  int Ny;  // num_voxels in y direction
  int Nz;  // num_voxels in z direction
  int xjmp; // Ny*Nz
  int yjmp; // Nz
  std::vector< std::vector<int> > pattern; // sweep pattern, for each stage of the sweep provides a list of voxels to update, #lists = #stages

  sweep(int d, int nx, int ny, int nz); // d has to be range of 0-7
  void compute_pattern();

  // small routines
  int flatten(int ix, int iy, int iz);
  std::vector<int> getNeighbors(int ind);
  bool num_not_present_in_list(std::vector<int> num_list, int num);
};

// VS Added
std::vector<sweep> compute_sweep_list(int nx, int ny, int nz); // get all 8 possible voxel-sweeping pattern

#endif
