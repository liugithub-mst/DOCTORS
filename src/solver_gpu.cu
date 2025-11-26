#include <iostream>
#include <stdio.h>
#include <stdlib.h>
#include <vector>
#include <cuda_runtime.h>

#include "solver_metadata.h"
#include "sweep.h"
#include "solver_gpu.cuh"

#include "macros.h"
#include "device_utils.cuh"
#include "legendre_gpu.cuh"


/*********** BEGIN: ROUTINES FOR ISOTROPIC LBTE SOLVER *******************/

// The input metadata is on the CPU (host).
struct solver_metadata *transferMetadataToDevice(struct solver_metadata metadata)
{
  struct solver_metadata s_host, *s_device;

  // Scalar parameters stay on the host
  s_host.num_voxels = metadata.num_voxels;
  s_host.num_angles = metadata.num_angles;
  s_host.num_energies = metadata.num_energies;
  s_host.num_materials = metadata.num_materials;
  s_host.pn = metadata.pn;
  s_host.voxel_vol = metadata.voxel_vol;
  s_host.voxel_area_yz = metadata.voxel_area_yz;
  s_host.voxel_area_xz = metadata.voxel_area_xz;
  s_host.voxel_area_xy = metadata.voxel_area_xy;

  s_host.Nx = metadata.Nx;
  s_host.Ny = metadata.Ny;
  s_host.Nz = metadata.Nz;

  s_host.sx = metadata.sx;  // source x position
  s_host.sy = metadata.sy;  // source y position
  s_host.sz = metadata.sz;  // source z position


  // Dynamic arrays are allocated memory on the device, but their pointers stay on the host
  // Quadrature size: #angles
  cudaMalloc((void **)&s_host.mu, sizeof(float)*s_host.num_angles);
  cudaMemcpy(s_host.mu, metadata.mu, sizeof(float)*s_host.num_angles, cudaMemcpyHostToDevice);

  cudaMalloc((void **)&s_host.zi, sizeof(float)*s_host.num_angles);
  cudaMemcpy(s_host.zi, metadata.zi, sizeof(float)*s_host.num_angles, cudaMemcpyHostToDevice);

  cudaMalloc((void **)&s_host.eta, sizeof(float)*s_host.num_angles);
  cudaMemcpy(s_host.eta, metadata.eta, sizeof(float)*s_host.num_angles, cudaMemcpyHostToDevice);

  cudaMalloc((void **)&s_host.quad_wt, sizeof(float)*s_host.num_angles);
  cudaMemcpy(s_host.quad_wt, metadata.quad_wt, sizeof(float)*s_host.num_angles, cudaMemcpyHostToDevice);

  cudaMalloc((void**)&s_host.XNodes, sizeof(float) * (s_host.Nx + 1));  
  cudaMemcpy(s_host.XNodes, metadata.XNodes, sizeof(float) * (s_host.Nx + 1), cudaMemcpyHostToDevice);

  cudaMalloc((void**)&s_host.YNodes, sizeof(float) * (s_host.Ny + 1));  
  cudaMemcpy(s_host.YNodes, metadata.YNodes, sizeof(float) * (s_host.Ny + 1), cudaMemcpyHostToDevice);

  cudaMalloc((void**)&s_host.ZNodes, sizeof(float) * (s_host.Nz + 1));  
  cudaMemcpy(s_host.ZNodes, metadata.ZNodes, sizeof(float) * (s_host.Nz + 1), cudaMemcpyHostToDevice);

  cudaMalloc((void**)&s_host.sIntensity, sizeof(float) * s_host.num_energies);  
  cudaMemcpy(s_host.sIntensity, metadata.sIntensity, sizeof(float) * s_host.num_energies, cudaMemcpyHostToDevice);


  // Density and (mapped) zone-Id, size: #voxels
  cudaMalloc((void **)&s_host.atomDensity, sizeof(float)*s_host.num_voxels);
  cudaMemcpy(s_host.atomDensity, metadata.atomDensity, sizeof(float)*s_host.num_voxels, cudaMemcpyHostToDevice);

  cudaMalloc((void **)&s_host.zoneIdMapped, sizeof(int)*s_host.num_voxels);
  cudaMemcpy(s_host.zoneIdMapped, metadata.zoneIdMapped, sizeof(int)*s_host.num_voxels, cudaMemcpyHostToDevice);

  // Differential scatter cross-section: size: #materials * #(src)energies * #(sink)energies
  cudaMalloc((void **)&s_host.sigma_s, sizeof(float)*(s_host.num_materials * s_host.num_energies * s_host.num_energies * (s_host.pn+1)));
  cudaMemcpy(s_host.sigma_s, metadata.sigma_s, sizeof(float)*(s_host.num_materials * s_host.num_energies * s_host.num_energies * (s_host.pn+1)), cudaMemcpyHostToDevice);

  // Total cross-section, size: #materials * #energies
  cudaMalloc((void **)&s_host.sigma_total, sizeof(float)*(s_host.num_materials * s_host.num_energies));
  cudaMemcpy(s_host.sigma_total, metadata.sigma_total, sizeof(float)*(s_host.num_materials * s_host.num_energies), cudaMemcpyHostToDevice);

  // Copy strcuture to device. Now bot the scalar parameters, pointers and the arrays are on device,
  cudaMalloc((void**)&s_device, sizeof(struct solver_metadata));
  cudaMemcpy(s_device, &s_host, sizeof(struct solver_metadata), cudaMemcpyHostToDevice);

  return s_device;
}


struct sweep_pattern *transferSweepPatternToDevice(struct sweep_pattern sweep_struct)
{
  struct sweep_pattern s_host, *s_device;

  s_host.d = sweep_struct.d;
  s_host.num_stages = sweep_struct.num_stages;
  s_host.num_voxels = sweep_struct.num_voxels;

  // Dynamic arrays are allocated memory on the device, but their pointers stay on the host
  cudaMalloc((void **)&s_host.num_voxels_per_stage, sizeof(int)*s_host.num_stages);
  cudaMemcpy(s_host.num_voxels_per_stage, sweep_struct.num_voxels_per_stage, sizeof(int)*s_host.num_stages, cudaMemcpyHostToDevice);

  cudaMalloc((void **)&s_host.starting_index_per_stage, sizeof(int)*s_host.num_stages);
  cudaMemcpy(s_host.starting_index_per_stage, sweep_struct.starting_index_per_stage, sizeof(int)*s_host.num_stages, cudaMemcpyHostToDevice);

  cudaMalloc((void **)&s_host.voxel_list, sizeof(int)*s_host.num_voxels);
  cudaMemcpy(s_host.voxel_list, sweep_struct.voxel_list, sizeof(int)*s_host.num_voxels, cudaMemcpyHostToDevice);

  // Copy strcuture to device. Now bot the scalar parameters, pointers and the arrays are on device,
  cudaMalloc((void**)&s_device, sizeof(struct sweep_pattern));
  cudaMemcpy(s_device, &s_host, sizeof(struct sweep_pattern), cudaMemcpyHostToDevice);

  return s_device;
}


struct sweep_pattern convertSweepClassToStruct(sweep sweep_obj)
{
    struct sweep_pattern sweep_struct;
    int i, j, istart, num_voxels;
    num_voxels = sweep_obj.Nx*sweep_obj.Ny*sweep_obj.Nz;


    sweep_struct.d = 4*(sweep_obj.dix+1)/2 + 2*(sweep_obj.diy+1)/2 + (sweep_obj.diz+1)/2; 
    sweep_struct.num_stages = sweep_obj.pattern.size();
    sweep_struct.num_voxels = num_voxels;
    sweep_struct.num_voxels_per_stage = new int[sweep_struct.num_stages];
    sweep_struct.starting_index_per_stage = new int[sweep_struct.num_stages];
    sweep_struct.voxel_list = new int[num_voxels]; //each voxel is covered only once

    sweep_struct.num_voxels_per_stage[0] = sweep_obj.pattern[0].size(); // this will be = 1
    sweep_struct.starting_index_per_stage[0] = 0;

    for(i=1; i<sweep_struct.num_stages; i++)
    {
      sweep_struct.num_voxels_per_stage[i] = sweep_obj.pattern[i].size();
      sweep_struct.starting_index_per_stage[i] = sweep_struct.starting_index_per_stage[i-1] + sweep_struct.num_voxels_per_stage[i-1];
    }

    for(i=0; i<sweep_struct.num_stages; i++)
    {
      istart = sweep_struct.starting_index_per_stage[i];
      for(j=0; j<sweep_struct.num_voxels_per_stage[i]; j++)
        sweep_struct.voxel_list[istart + j] = sweep_obj.pattern[i][j];
      if(istart >= num_voxels)
      {
        std::cout<<"Error in convertSweepClassToStruct(): Index exceeded number of voxels"<<std::endl;
        exit(-1);
      }
    }

    return sweep_struct;
}


float *getufluxArray(const std::vector<double>* uflux_vector, int num_energies, int num_voxels)
{
  float *uflux_array = new float[num_energies * num_voxels];

  for(int ie=0; ie<num_energies; ie++)
  for(int i=0;  i<num_voxels; i++)
    uflux_array[ie*num_voxels+i] = (float) ((*uflux_vector)[ie*num_voxels + i]);

  return uflux_array;
}

__global__ void compute_totalSource_Iso(float *totalSource, float *uFlux, float *cFlux, struct solver_metadata *metadata, int ie)
{
  int ir = blockIdx.x * blockDim.x + threadIdx.x;

  int im;
  int im_step   = metadata->num_energies * metadata->num_energies * (metadata->pn+1);
  int src_step  = metadata->num_energies * (metadata->pn+1);
  int sink_step = metadata->pn+1;
  int il=0;

  if(ir < metadata->num_voxels)
  {
    // initialize
    totalSource[ir]=0;

    // materil mapped ID
    im = metadata->zoneIdMapped[ir];

    // uncollidex-flux source term (downscatter + inscatter)
    for(int iie = 0; iie <= ie; iie++)
      totalSource[ir] += uFlux[iie*metadata->num_voxels + ir] * metadata->sigma_s[im * im_step + iie * src_step + ie * sink_step + il] * metadata->atomDensity[ir] * metadata->voxel_vol;

    // collided-flux (downscatter only)
    for(int iie = 0; iie < ie; iie++)
      totalSource[ir] += cFlux[iie*metadata->num_voxels + ir] * metadata->sigma_s[im * im_step + iie * src_step + ie * sink_step + il] * metadata->atomDensity[ir] * metadata->voxel_vol;
  }
}

//compute scatter within same energy group.
__global__ void compute_inscatter_Iso(float *inscatter, float *cFlux, struct solver_metadata *metadata, int ie)
{
  int ir = blockIdx.x * blockDim.x + threadIdx.x;

  int im;
  int im_step   = metadata->num_energies * metadata->num_energies * (metadata->pn+1);
  int src_step  = metadata->num_energies * (metadata->pn+1);
  int sink_step = metadata->pn+1;
  int il=0;

  if(ir < metadata->num_voxels)
  {
    im = metadata->zoneIdMapped[ir];
    inscatter[ir] = cFlux[ie*metadata->num_voxels + ir] *  metadata->sigma_s[im * im_step + ie * src_step + ie * sink_step + il]  *  metadata->atomDensity[ir] * metadata->voxel_vol;
  }
}

__global__ void compute_totalCrossSection(SOL_T *totalCrossSection, struct solver_metadata *metadata, int ie)
{
    int ir = blockIdx.x * blockDim.x + threadIdx.x;
    int im;

    if(ir < metadata->num_voxels)
    {
      im = metadata->zoneIdMapped[ir];
      double sigma = metadata->sigma_total[im * metadata->num_energies + ie];
      double density = metadata->atomDensity[ir];

      // Basic sanity checks
      if (sigma < 0.0 || isnan(sigma)) sigma = 0.0;
      if (density < 0.0 || isnan(density)) density = 0.0;

      totalCrossSection[ir] = sigma * density * metadata->voxel_vol;

      if (totalCrossSection[ir] < 1e-8f) 
      {
         printf("Warning: totalXS[%d] = %f (sigma = %f, density = %f)\n", ir, totalCrossSection[ir], sigma, density);
      }
    }
    __syncthreads(); 
}

// Sum the flux across differnt angles, each computed by a different "block" of threads.
__global__ void sumFluxAcrossAngles(SOL_T *nextFlux, SOL_T *nextFlux_block, struct solver_metadata *metadata)
{
    int ir = blockIdx.x * blockDim.x + threadIdx.x;

    if(ir < metadata->num_voxels)
    {
      for(int iang=0; iang < metadata->num_angles; iang++)
        nextFlux[ir] += nextFlux_block[iang * metadata->num_voxels + ir];
    }
    __syncthreads(); 
}

// Find difference between previous and current iteration LBTE solution
/*  Debug 
__global__ void updateDiff(SOL_T *norm_diff, SOL_T *diff, SOL_T *nextFlux, SOL_T *cFlux, struct solver_metadata *metadata, int ie)
{
  int ir = blockIdx.x * blockDim.x + threadIdx.x;
  
  if(ir < metadata->num_voxels)
  {
    diff[ir] = fabsf(nextFlux[ir]-cFlux[ie * metadata->num_voxels + ir]);

    if(nextFlux[ir] != 0)
      norm_diff[ir] = diff[ir]/fabsf(nextFlux[ir]);
    else
      norm_diff[ir] = 0; // to take care of divide by zero issues
  }
  __syncthreads();
}
*/
__global__ void updateDiff(SOL_T *norm_diff, SOL_T *diff, SOL_T *nextFlux, SOL_T *cFlux, struct solver_metadata *metadata, int ie)
{ 
  int ir = blockIdx.x * blockDim.x + threadIdx.x;
  
  if(ir < metadata->num_voxels)
  {
    diff[ir] = fabs(nextFlux[ir]-cFlux[ie * metadata->num_voxels + ir]);

    if(fabs(nextFlux[ir]) > TINY_T)
      norm_diff[ir] = diff[ir]/fabs(nextFlux[ir]);
    else
      norm_diff[ir] = 0.0; // to take care of divide by zero issues
  }
  __syncthreads();
}

// copies flux summed across angles from a local buffer (energy-independent) to a global buffer (energy-dependent)
__global__ void updatecFluxAtEnergy(SOL_T *cFlux, SOL_T *nextFlux, struct solver_metadata *metadata, int ie)
{
  int ir = blockIdx.x * blockDim.x + threadIdx.x;

  if(ir < metadata->num_voxels)
    cFlux[ie * metadata->num_voxels + ir] = nextFlux[ir];

  __syncthreads();
}


// single lBTE iteration at given energy
__global__ void LBTE_Iteration_Iso(float *nextFlux_block, float *totalSource, float *inscatter, float *totalCrossSection,
                                   float *outboundFluxX, float *outboundFluxY, float *outboundFluxZ,
                                   struct sweep_pattern **sweep_list, struct solver_metadata *metadata, int ie)
{
  // each block solves LBTE for a specific scatter-angle in the quadrature
  int iang = blockIdx.x;
  // within a block
  int Nthreads = blockDim.x;
  int tid = threadIdx.x;

  // Get parameters for specific angle within quadrature
  float Ayz = 2 * fabsf(metadata->mu[iang]) * metadata->voxel_area_yz;
  float Axz = 2 * fabsf(metadata->zi[iang]) * metadata->voxel_area_xz;
  float Axy = 2 * fabsf(metadata->eta[iang])* metadata->voxel_area_xy;

  // cartesian directions: (+/-)x, (+/-)y, (+/-) z.
  int dix = (metadata->mu[iang] < 0) ? -1 : 1;
  int diy = (metadata->zi[iang] < 0) ? -1 : 1;
  int diz = (metadata->eta[iang]< 0) ? -1 : 1;
  int pattern_ind  = 4*(dix+1)/2 + 2*(diy+1)/2 + (diz+1)/2;    // Range: 0->7

  int iang_offset = iang * metadata->Nx * metadata->Ny * metadata->Nz;
  int xjmp = metadata->Ny * metadata->Nz;
  int yjmp = metadata->Nz;

  float influxX, influxY, influxZ;
  float outx, outy, outz;
  float numer, denom, angFlux;

  int voxel_ind_local, voxel_ind_global, vbegin, vend;
  int ir, ix, iy, iz, stage, cflag;


  // Begin voxel-sweep
  // Thread synchronization within block needed for this for-loop
  for(stage=0; stage < sweep_list[pattern_ind]->num_stages; stage++)
  {
    vbegin = (tid*sweep_list[pattern_ind]->num_voxels_per_stage[stage])/Nthreads;
    vend   = ((tid+1)*sweep_list[pattern_ind]->num_voxels_per_stage[stage])/Nthreads;

    // NO thread synchronization needed for this loop
    for(voxel_ind_local=vbegin; voxel_ind_local < vend; voxel_ind_local++)
    {
          // index within the full (global) voxel-list that goes across stages
          voxel_ind_global = sweep_list[pattern_ind]->starting_index_per_stage[stage] + voxel_ind_local;

          // Voxel 3D coordinates (flattened into 1 coordinate)
          ir = sweep_list[pattern_ind]->voxel_list[voxel_ind_global];
          ix = ir/xjmp;
          iy = (ir%xjmp)/yjmp;
          iz = (ir%xjmp)%yjmp;

          // Extract neighbors in a 3-point neighborhood
          if(dix > 0)
          {
              if(ix == 0)                                               // If this is a boundary cell
                  influxX = 0;                                       // then the in-flux is zero
              else                                                      // otherwise
                  influxX = outboundFluxX[iang_offset + (ix-1)*xjmp + iy*yjmp + iz];  // the in-flux is the out-flux from the previous cell
          }
          else
          {
              if(ix == (metadata->Nx-1))
                  influxX = 0;
              else
                  influxX = outboundFluxX[iang_offset + (ix+1)*xjmp + iy*yjmp + iz];
          }

          if(diy > 0)
          {
              if(iy == 0)
                  influxY = 0;
              else
                  influxY = outboundFluxY[iang_offset + ix*xjmp + (iy-1)*yjmp + iz];
          }
          else
          {
              if(iy == (metadata->Ny-1))
                  influxY = 0;
              else
                  influxY = outboundFluxY[iang_offset + ix*xjmp + (iy+1)*yjmp + iz];
          }

          if(diz > 0)
          {
              if(iz == 0)
                  influxZ = 0;
              else
                  influxZ = outboundFluxZ[iang_offset + ix*xjmp + iy*yjmp + iz-1];
          }
          else
          {
              if(iz == (metadata->Nz-1))
                  influxZ = 0;
              else
                  influxZ = outboundFluxZ[iang_offset + ix*xjmp + iy*yjmp + iz+1];
          }

          numer = totalSource[ir] + inscatter[ir] + Ayz*influxX + Axz*influxY + Axy*influxZ;

          denom = totalCrossSection[ir] + Ayz + Axz + Axy;

          // update voxel
          angFlux = numer/denom;

          // value at voxel-boundaries (used as incoming flux for next stage of sweep)
          outx = 2*angFlux - influxX;
          outy = 2*angFlux - influxY;
          outz = 2*angFlux - influxZ;

          cflag = 0;
          if(outx < 0)
          {
              cflag = 1;
              outx = 0;
          }
          if(outy < 0)
          {
              cflag = 1;
              outy = 0;
          }
          if(outz < 0)
          {
              cflag = 1;
              outz = 0;
          }
          if(cflag)
          {
              angFlux = (influxX + influxY + influxZ + outx + outy + outz)/6.0;
          }

          outboundFluxX[iang_offset + ir] = outx;
          outboundFluxY[iang_offset + ir] = outy;
          outboundFluxZ[iang_offset + ir] = outz;
          nextFlux_block[iang_offset + ir] = metadata->quad_wt[iang] * angFlux;

    } 
    // Synhcronize threads in the block before next stage of the voxel-sweep
    __syncthreads();
  }

}
/* XL: comment out 
float *gsSolverIsoGPU(struct solver_metadata metadata, const std::vector<double>* uFlux_vector)
{
  struct solver_metadata *metadata_device;
  struct sweep_pattern sweep_temp, *sweep_device[8], **sweep_list_device;
  int ie, deviceNum=0; // gpu-id

  // device pointers - allocate on device. NO copy from host-->device needed
  float *outboundFluxX, *outboundFluxY, *outboundFluxZ; //size: #angles * #voxels. (Optional) re-set to 0 prior to every LBTE iteration
  float *nextFlux_block; //size: #angles * #voxels . Should be re-set to 0 prior to every LBTE iteration

  float *inscatter;    // size: #voxels. This needs to be recomputed at every energy level.
  float *totalSource;  // size: #voxels. compute on device. This needs to be recomputed (after initialization to 0) at every energy level.
  float *nextFlux;     // size: #voxels. accumulated flux over all angles.
                       // reset to 0 prior to every LBTE iteration.
  float *totalCrossSection; // size:#voxels. Meeds to be compute once every energy-lvel.


  float *uFlux_device; // size: #energies * #voxels. Just copy once from host-->device (before the LBTE solver starts).
  float *cFlux_device; // size: #energies * #voxels. Just intialize once to 0 (before LBTE solver starts)
                       // After every LBTE iteration update this based on nextFlux.

  float *diff_device, *norm_diff_device;  // after every LBTE iteration copy this value from device-->host

  // host pointers
  float *cFlux_host = new float[metadata.num_energies * metadata.num_voxels];
  float *uFlux_host; // use getufluxArray to set this

  // Convergence parameters
  int maxIterations = 25;
  float epsilon = 0.001;
  // variables for tracking convergenxe
  int i, it=0;
  float maxDiff=0, totDiff=0;
  float *diff, *norm_diff;

  // Precompute sweep patterns
  std::cout<<"Precomputing 3-D voxel-sweep patterns for different directions ..."<<std::endl;
  std::vector<sweep> sweep_list = compute_sweep_list(metadata.Nx, metadata.Ny, metadata.Nz);


  if (cudaSetDevice(deviceNum) != cudaSuccess) {
      fprintf(stderr, "Error initializing device %d\n", deviceNum);
      exit(-1);
  }
  else
      printf("Using gpu_device %d\n", deviceNum);

  //----metadata transfer--
  std::cout<<"Transferring metadata to device ..."<<std::endl;
  metadata_device = transferMetadataToDevice(metadata);

  //---sweep pattern transfer ---
  std::cout<<"Transferring sweep patterns to device ..."<<std::endl;
  // sweep_device[d] is a pointer that resides on the host but points to a memory location on the device (GPU).
  for(int d=0; d<8; d++)
  {
    sweep_temp = convertSweepClassToStruct(sweep_list[d]); // class-->structure
    sweep_device[d] = transferSweepPatternToDevice(sweep_temp); // transfer data of structure to device
  }
  // List of all 8  sweep-patterns represneted by a **sweep_list_device, where sweep_list_device[d], d=0...8-1, is a pointer on the device ...
  // that represents a particular sweep-pattern.
  cudaMalloc((void **)&sweep_list_device, 8*sizeof(struct sweep_pattern*));
  cudaMemcpy(sweep_list_device, &sweep_device[0], 8*sizeof(struct sweep_pattern*), cudaMemcpyHostToDevice);

  //--allocate arrays on device--
  // initialization / transfer may need to be done repeatedly during the LBTE solver.
  std::cout<<"Allocating vectors on device for the LBTE solver ..."<<std::endl;
  cudaMalloc((void **)&outboundFluxX, sizeof(float)*metadata.num_voxels*metadata.num_angles);
  cudaMalloc((void **)&outboundFluxY, sizeof(float)*metadata.num_voxels*metadata.num_angles);
  cudaMalloc((void **)&outboundFluxZ, sizeof(float)*metadata.num_voxels*metadata.num_angles);
  cudaMalloc((void **)&nextFlux_block, sizeof(float)*metadata.num_voxels*metadata.num_angles);

  cudaMalloc((void **)&inscatter, sizeof(float)*metadata.num_voxels);
  cudaMalloc((void **)&totalSource, sizeof(float)*metadata.num_voxels);
  cudaMalloc((void **)&nextFlux, sizeof(float)*metadata.num_voxels);
  cudaMalloc((void **)&totalCrossSection, sizeof(float)*metadata.num_voxels);

  cudaMalloc((void **)&uFlux_device, sizeof(float)*metadata.num_voxels*metadata.num_energies);
  cudaMalloc((void **)&cFlux_device, sizeof(float)*metadata.num_voxels*metadata.num_energies);

  // set (u/c)fluxes on device before LBTE solver begins
  uFlux_host = getufluxArray(uFlux_vector, metadata.num_energies, metadata.num_voxels);
  cudaMemcpy(uFlux_device, uFlux_host, sizeof(float)*metadata.num_voxels*metadata.num_energies, cudaMemcpyHostToDevice);
  cudaMemset(cFlux_device, 0, sizeof(float)*metadata.num_voxels*metadata.num_energies);

  // convgerence tracking
  diff = new float[metadata.num_voxels];
  norm_diff = new float[metadata.num_voxels];
  cudaMalloc((void **)&diff_device, sizeof(float)*metadata.num_voxels);
  cudaMalloc((void **)&norm_diff_device, sizeof(float)*metadata.num_voxels);

  // Grid and Block dimensions
  // for precomputing source-term, downscatter and in-scatter
  dim3 block_precomp(BLOCK_SIZE_PRECOMP);
  dim3 grid_precomp(metadata.num_voxels/BLOCK_SIZE_PRECOMP+1);
  // for LBTE iteration
  dim3 block_LBTE(BLOCK_SIZE_LBTE);
  dim3 grid_LBTE(metadata.num_angles);

  bool notTrackedHighestEnergy = true;
  int ira, highestEnergy = 0;
  while(notTrackedHighestEnergy)
  {
      for(ira = 0; ira < metadata.num_voxels; ira++)
      {
          if((*uFlux_vector)[(unsigned int) (highestEnergy*metadata.num_voxels + ira)] > 0)
            break;
      }
      if(ira==metadata.num_voxels)
      {
          std::cout << "No external source or downscatter, skipping energy group " << highestEnergy << std::endl;
          highestEnergy++;
      }
      else
          notTrackedHighestEnergy = false;

      if(highestEnergy >= metadata.num_energies)
      {
          std::cout << "Zero flux everywhere from the raytracer" << std::endl;
          return NULL;
      }
  }

  // LBTE solver
  for(ie=highestEnergy; ie<metadata.num_energies; ie++)
  {
    std::cout<<"### Energy group "<<ie<<" ###"<<std::endl;
    // LBTE iteration count for given energy
    it=0;
    maxDiff=1; // just set so that it enters the loop

    // Precompute source-term from uncollided-flux and downscatter from collided flux
    compute_totalSource_Iso<<< grid_precomp, block_precomp  >>>(totalSource, uFlux_device, cFlux_device, metadata_device, ie);
    // Precompute tota scatter cross-section
    compute_totalCrossSection<<< grid_precomp, block_precomp >>>(totalCrossSection, metadata_device, ie);

    // Single LBTE iteration
    while((it<maxIterations) && (maxDiff > epsilon))
    {
      // Prior to every LBTE iteration (before sweep of angular-direction and voxel-grid) do the following re-initialization.
      // initialization required every iteration
      cudaMemset(nextFlux_block, 0, sizeof(float)*metadata.num_voxels*metadata.num_angles);
      cudaMemset(nextFlux, 0, sizeof(float)*metadata.num_voxels);

      // this initialization should be done once before LBTE solver. re-initialization before every LBTE iteration is optional.
      cudaMemset(outboundFluxX, 0, sizeof(float)*metadata.num_voxels*metadata.num_angles);
      cudaMemset(outboundFluxY, 0, sizeof(float)*metadata.num_voxels*metadata.num_angles);
      cudaMemset(outboundFluxZ, 0, sizeof(float)*metadata.num_voxels*metadata.num_angles);

      // Pre-compute in-scatter from collided-flux
      compute_inscatter_Iso<<< grid_precomp, block_precomp >>>(inscatter, cFlux_device, metadata_device, ie);

      // LBTE - angular-sweep and voxel-sweep
      LBTE_Iteration_Iso<<< grid_LBTE, block_LBTE >>>(nextFlux_block, totalSource, inscatter, totalCrossSection,  outboundFluxX, outboundFluxY, outboundFluxZ, sweep_list_device, metadata_device, ie);

      // At end of LBTE iteration  do ..
      // 1) nextFlux <-- sum nextFlux_block over all angles
      // 2) find the difference between previous and current iteration LBTE solution (to check for convergence)
      // 3) cFlux_device[energy, :] <-- nextFlux
      sumFluxAcrossAngles<<< grid_precomp, block_precomp >>>(nextFlux, nextFlux_block, metadata_device);
      updateDiff<<< grid_precomp, block_precomp >>>(norm_diff_device, diff_device, nextFlux, cFlux_device, metadata_device, ie);
      updatecFluxAtEnergy<<< grid_precomp, block_precomp >>>(cFlux_device, nextFlux, metadata_device, ie);

      // Track convergence here
      cudaMemcpy(diff, diff_device, sizeof(float)*metadata.num_voxels, cudaMemcpyDeviceToHost);
      cudaMemcpy(norm_diff, norm_diff_device, sizeof(float)*metadata.num_voxels, cudaMemcpyDeviceToHost);
      totDiff=0;
      maxDiff=0;
      for(i = 0; i < metadata.num_voxels; i++)
      {
          maxDiff = std::max(maxDiff, norm_diff[i]);
          totDiff += diff[i];
      }

      it++;

      std::cout<<"Iteration "<<it<<": Relative error = "<<maxDiff<<std::endl;
    } // END: while-loop corresponding to LBTE iterations at given energy-group

    if(!(it < maxIterations))
        std::cout << "Max iterations hit" << std::endl;
    else if(!(maxDiff > epsilon))
        std::cout << "Converged on relative error" << std::endl;
    else
        std::cout << "Converged on precision bound" << std::endl;

    // For updating Progress
    std::cout << ((ie+1) * 100)/metadata.num_energies << "%% finished" << std::endl;

  } //END: for(ie=0; ...)

  // Copy LBTE solution across all energy levels
  cudaMemcpy(cFlux_host, cFlux_device, sizeof(float)*metadata.num_energies*metadata.num_voxels, cudaMemcpyDeviceToHost);

  // De-allocate memory
  cudaFree(outboundFluxX);
  cudaFree(outboundFluxY);
  cudaFree(outboundFluxZ);
  cudaFree(nextFlux_block);
  cudaFree(cFlux_device);
  cudaFree(uFlux_device);

  cudaFree(nextFlux);
  cudaFree(inscatter);
  cudaFree(totalSource);
  cudaFree(totalCrossSection);

  free(uFlux_host);

  return cFlux_host;

}
*/
/*********** END: ROUTINES FOR ISOTROPIC LBTE SOLVER *******************/

/*********** BEGIN: ROUTINES FOR ANISOTROPIC LBTE SOLVER ***************/

// Function to convert uncollided flux to flux moments and generate in-moments
__global__ void update_umomentsKernel(SOL_T* uFlux, struct solver_metadata* metadata, SOL_T* umoments, SOL_T* cmoments, 
                                      SOL_T* inmoments, int numMoments, int highestE, int currentE)
{
    unsigned int xIndx = blockIdx.x;
    unsigned int yIndx = blockIdx.y;
    unsigned int zIndx = blockIdx.z;
    //unsigned int ilIndx = threadIdx.x;

    unsigned int ir0 = xIndx * metadata->Ny * metadata->Nz + yIndx * metadata->Nz + zIndx;

    // Assume uniform mesh
    float dx = metadata->XNodes[1] - metadata->XNodes[0];
    float dy = metadata->YNodes[1] - metadata->YNodes[0];
    float dz = metadata->ZNodes[1] - metadata->ZNodes[0];

    // Real position (x, y, z) of a voxel
    float x = metadata->XNodes[xIndx] + dx / 2;
    float y = metadata->YNodes[yIndx] + dy / 2;
    float z = metadata->ZNodes[zIndx] + dz / 2;

    SOL_T deltaX = x - metadata->sx;
    SOL_T deltaY = y - metadata->sy;
    SOL_T deltaZ = z - metadata->sz;

    SOL_T srcToCellDist = sqrt(deltaX * deltaX + deltaY * deltaY + deltaZ * deltaZ);
    //if (srcToCellDist < 1e-8) srcToCellDist = 1e-8;  //XL: Debug
    
    SOL_T xcos = deltaX / srcToCellDist;     // sin(theta)sin(phi)
    SOL_T ycos = deltaY / srcToCellDist;     // cos(theta) chose Y axis as the polar axis as the X-ray is shooting down along Y axis
    SOL_T zcos = deltaZ / srcToCellDist;     // sin(theta)cos(phi)
    
    // Clamp xcos and other cosines so it could not exceed [-1, 1]
    //xcos = fmin(fmax(xcos, -1.0f), 1.0f);
    //ycos = fmin(fmax(ycos, -1.0f), 1.0f);
    //zcos = fmin(fmax(zcos, -1.0f), 1.0f);

    SOL_T phi = atan2(xcos, zcos);
    SOL_T theta = acos(ycos);                // chose Y axis as the polar axis as the X-ray is shooting down along Y axis

    for (unsigned int il = 0; il <= metadata->pn; il++)
    {
        // Regular uncollided flux moment, if m = 0
        umoments[currentE * metadata->num_voxels * numMoments + ir0 * numMoments + il * il] = uFlux[currentE * metadata->num_voxels + ir0] * sqrt(2.0*il + 1.0) * computAssocLegendre(il, 0, ycos);

        //Generate regular in-moments
        for (unsigned int iep = highestE; iep <= currentE; iep++)
        {
            SOL_T total_zeroMoments = umoments[iep * metadata->num_voxels * numMoments + ir0 * numMoments + il * il] +
                                      cmoments[iep * metadata->num_voxels * numMoments + ir0 * numMoments + il * il];

            int im = metadata->zoneIdMapped[ir0];
            int im_step = metadata->num_energies * metadata->num_energies * (metadata->pn + 1);
            int src_step = metadata->num_energies * (metadata->pn + 1);
            int sink_step = metadata->pn + 1;
            SOL_T coeff = metadata->sigma_s[im * im_step + iep * src_step + currentE * sink_step + il] * metadata->atomDensity[ir0] * metadata->voxel_vol;
            inmoments[ir0 * numMoments + il * il] += coeff * total_zeroMoments;  // zero in-moments
        }

        for (unsigned int imm = 1; imm <= il; imm++)
        {
            umoments[currentE * metadata->num_voxels * numMoments + ir0 * numMoments + il * il + 2 * imm - 1] =
                uFlux[currentE * metadata->num_voxels + ir0] * normConst(il, imm) * computAssocLegendre(il, imm, ycos) * sin(imm * phi);  
            umoments[currentE * metadata->num_voxels * numMoments + ir0 * numMoments + il * il + 2 * imm] =
                uFlux[currentE * metadata->num_voxels + ir0] * normConst(il, imm) * computAssocLegendre(il, imm, ycos) * cos(imm * phi);  

            // Generate cosine and sine in-moments
            for (unsigned int iep = highestE; iep <= currentE; iep++)
            {
                SOL_T total_sinMoments = umoments[iep * metadata->num_voxels * numMoments + ir0 * numMoments + il * il + 2 * imm - 1] +
                                         cmoments[iep * metadata->num_voxels * numMoments + ir0 * numMoments + il * il + 2 * imm - 1];
                SOL_T total_cosMoments = umoments[iep * metadata->num_voxels * numMoments + ir0 * numMoments + il * il + 2 * imm] +
                                         cmoments[iep * metadata->num_voxels * numMoments + ir0 * numMoments + il * il + 2 * imm];

                int im = metadata->zoneIdMapped[ir0];
                int im_step = metadata->num_energies * metadata->num_energies * (metadata->pn + 1);
                int src_step = metadata->num_energies * (metadata->pn + 1);
                int sink_step = metadata->pn + 1;
                SOL_T coeff = metadata->sigma_s[im * im_step + iep * src_step + currentE * sink_step + il] * metadata->atomDensity[ir0] * metadata->voxel_vol;
                inmoments[ir0 * numMoments + il * il + 2 * imm - 1] += coeff * total_sinMoments;  // sine in-moments
                inmoments[ir0 * numMoments + il * il + 2 * imm] += coeff * total_cosMoments;      // cosine in-moments
            }
        }

    }
}

// Update the collided flux moment across differnt angles, each computed by a different "block" of threads.
__global__ void updateFluxMoment(SOL_T* nextMoment, SOL_T* nextFlux_block, struct solver_metadata* metadata, int numMoments)
{
    int ir = blockIdx.x * blockDim.x + threadIdx.x;

    if (ir < metadata->num_voxels)
    {
        for (int iang = 0; iang < metadata->num_angles; iang++)
        {
            SOL_T costheta = metadata->eta[iang];  
            SOL_T phi = atan2(metadata->mu[iang], metadata->zi[iang]);
            for (unsigned int il = 0; il <= metadata->pn; il++)
            {
                nextMoment[ir * numMoments + il * il] += nextFlux_block[iang * metadata->num_voxels + ir] * computAssocLegendre(il, 0, costheta);
                for (unsigned int im = 1; im <= il; im++)
                {
                    nextMoment[ir * numMoments + il * il + 2 * im - 1] += nextFlux_block[iang * metadata->num_voxels + ir] *
                                                                          computAssocLegendre(il, im, costheta) * sin(im * phi);
                    nextMoment[ir * numMoments + il * il + 2 * im] += nextFlux_block[iang * metadata->num_voxels + ir] *
                                                                      computAssocLegendre(il, im, costheta) * cos(im * phi);
                }
            }

        }

    }
    
}

// Copies updated collided flux moment from a local buffer (energy-independent) to a global buffer (energy-dependent)
__global__ void updatecMomentAtEnergy(SOL_T* cMoment, SOL_T* nextMoment, struct solver_metadata* metadata, int ie, int numMoments)
{
    unsigned int xIndx = blockIdx.x;
    unsigned int yIndx = blockIdx.y;
    unsigned int zIndx = blockIdx.z;
    unsigned int ilIndx = threadIdx.x;

    unsigned int ir = xIndx * metadata->Ny * metadata->Nz * numMoments + yIndx * metadata->Nz * numMoments + zIndx * numMoments + ilIndx;

    if (ir < metadata->num_voxels * numMoments)
        cMoment[ie * metadata->num_voxels * numMoments + ir] = nextMoment[ir];
    else
        printf("updatecMomentAtEnergy error: ir = %d\n", ir);
}

// Single LBTE iteration at given energy using Spherical Harmonic method
__global__ void LBTE_Iteration_Harmonic(SOL_T* nextFlux_block, SOL_T* totalCrossSection, SOL_T* cmoments, SOL_T* inmoments,
                                        SOL_T* outboundFluxX, SOL_T* outboundFluxY, SOL_T* outboundFluxZ, int numMoments,
                                        struct sweep_pattern** sweep_list, struct solver_metadata* metadata, int ie)
{
    // each block solves LBTE for a specific scatter-angle in the quadrature
    int iang = blockIdx.x;
    // within a block
    int Nthreads = blockDim.x;
    int tid = threadIdx.x;

    // Get parameters for specific angle within quadrature
    SOL_T Ayz = 2 * fabs(metadata->mu[iang]) * metadata->voxel_area_yz;
    SOL_T Axz = 2 * fabs(metadata->zi[iang]) * metadata->voxel_area_xz;
    SOL_T Axy = 2 * fabs(metadata->eta[iang]) * metadata->voxel_area_xy;

    // Find the angle (theta, phi)
    SOL_T theta = acos(metadata->eta[iang]);  // Choose Y axis as the polar axis
    SOL_T phi = atan2(metadata->mu[iang], metadata->zi[iang]); 

    // cartesian directions: (+/-)x, (+/-)y, (+/-) z.
    int dix = (metadata->mu[iang] < 0) ? -1 : 1;
    int diy = (metadata->eta[iang] < 0) ? -1 : 1;
    int diz = (metadata->zi[iang] < 0) ? -1 : 1;
    int pattern_ind = 4 * (dix + 1) / 2 + 2 * (diy + 1) / 2 + (diz + 1) / 2;    // Range: 0->7

    int iang_offset = iang * metadata->Nx * metadata->Ny * metadata->Nz;
    int xjmp = metadata->Ny * metadata->Nz;
    int yjmp = metadata->Nz;

    SOL_T influxX, influxY, influxZ;
    SOL_T outx, outy, outz;
    SOL_T numer, denom, angFlux;

    int voxel_ind_local, voxel_ind_global, vbegin, vend;
    int ir, ix, iy, iz, stage, cflag;


    // Begin voxel-sweep
    // Thread synchronization within block needed for this for-loop
    for (stage = 0; stage < sweep_list[pattern_ind]->num_stages; stage++)
    {
        vbegin = (tid * sweep_list[pattern_ind]->num_voxels_per_stage[stage]) / Nthreads;
        vend = ((tid + 1) * sweep_list[pattern_ind]->num_voxels_per_stage[stage]) / Nthreads;

        // NO thread synchronization needed for this loop
        for (voxel_ind_local = vbegin; voxel_ind_local < vend; voxel_ind_local++)
        {
            // index within the full (global) voxel-list that goes across stages
            voxel_ind_global = sweep_list[pattern_ind]->starting_index_per_stage[stage] + voxel_ind_local;

            // Voxel 3D coordinates (flattened into 1 coordinate)
            ir = sweep_list[pattern_ind]->voxel_list[voxel_ind_global];
            ix = ir / xjmp;
            iy = (ir % xjmp) / yjmp;
            iz = (ir % xjmp) % yjmp;

            // Extract neighbors in a 3-point neighborhood
            if (dix > 0)
            {
                if (ix == 0)                                               // If this is a boundary cell
                    influxX = 0;                                       // then the in-flux is zero
                else                                                      // otherwise
                    influxX = outboundFluxX[iang_offset + (ix - 1) * xjmp + iy * yjmp + iz];  // the in-flux is the out-flux from the previous cell
            }
            else
            {
                if (ix == (metadata->Nx - 1))
                    influxX = 0;
                else
                    influxX = outboundFluxX[iang_offset + (ix + 1) * xjmp + iy * yjmp + iz];
            }

            if (diy > 0)
            {
                if (iy == 0)
                    influxY = 0;
                else
                    influxY = outboundFluxY[iang_offset + ix * xjmp + (iy - 1) * yjmp + iz];
            }
            else
            {
                if (iy == (metadata->Ny - 1))
                    influxY = 0;
                else
                    influxY = outboundFluxY[iang_offset + ix * xjmp + (iy + 1) * yjmp + iz];
            }

            if (diz > 0)
            {
                if (iz == 0)
                    influxZ = 0;
                else
                    influxZ = outboundFluxZ[iang_offset + ix * xjmp + iy * yjmp + iz - 1];
            }
            else
            {
                if (iz == (metadata->Nz - 1))
                    influxZ = 0;
                else
                    influxZ = outboundFluxZ[iang_offset + ix * xjmp + iy * yjmp + iz + 1];
            }

            // Calcualte down-scattering source
            SOL_T totalSource = 0.0;
            SOL_T coeff = 0.0;
           
            for (unsigned int il = 0; il <= metadata->pn; il++)
            {
                int im = metadata->zoneIdMapped[ir];
                int im_step = metadata->num_energies * metadata->num_energies * (metadata->pn + 1);
                int src_step = metadata->num_energies * (metadata->pn + 1);
                int sink_step = metadata->pn + 1;

                coeff = metadata->sigma_s[im * im_step + ie * src_step + ie * sink_step + il] * metadata->atomDensity[ir] * metadata->voxel_vol;
                
                // downscatter from the uncollided and collided flux including group inscatter
                totalSource += (2.0 * il + 1.0) * computAssocLegendre(il, 0, cos(theta)) * (coeff * cmoments[ie * metadata->num_voxels * numMoments + ir * numMoments + il * il]
                                                                                       + inmoments[ir * numMoments + il * il]);

                
                 SOL_T n_const = 1.0; //XL: Debug chnage 2.0 * (2.0 * il + 1.0); // Normalization constant
                for (unsigned int imm = 1; imm <= il; imm++)
                {
                    totalSource +=  n_const * doubleFactorial(il-im) / doubleFactorial(il+im) * computAssocLegendre(il, imm, cos(theta)) * sin(imm * phi) *
                                   (coeff * cmoments[ie * metadata->num_voxels * numMoments + ir * numMoments + il * il + 2 * imm -1]
                                   + inmoments[ir * numMoments + il * il + 2 * imm -1]);
                    totalSource +=  n_const * doubleFactorial(il-im) / doubleFactorial(il+im) * computAssocLegendre(il, imm, cos(theta)) * cos(imm * phi) *
                                   (coeff * cmoments[ie * metadata->num_voxels * numMoments + ir * numMoments + il * il + 2 * imm]
                                   + inmoments[ir * numMoments + il * il + 2 * imm]);

                }

            }

            numer = totalSource + Ayz * influxX + Axz * influxY + Axy * influxZ;

            //XL: Debug
            denom = fmax(totalCrossSection[ir] + Ayz + Axz + Axy, 1e-30);

            // update voxel
            angFlux = numer / denom;
            if(isnan(angFlux) || isinf(angFlux))
            {
                 printf("NaN detected at voxel %d angle %d\n", ir, iang);
            }

            // value at voxel-boundaries (used as incoming flux for next stage of sweep)
            outx = 2 * angFlux - influxX;
            outy = 2 * angFlux - influxY;
            outz = 2 * angFlux - influxZ;

            cflag = 0;
            if (outx < 0)
            {
                cflag = 1;
                outx = 0;
            }
            if (outy < 0)
            {
                cflag = 1;
                outy = 0;
            }
            if (outz < 0)
            {
                cflag = 1;
                outz = 0;
            }
            if (cflag)
            {
                angFlux = (influxX + influxY + influxZ + outx + outy + outz) / 6.0;
            }

            outboundFluxX[iang_offset + ir] = outx;
            outboundFluxY[iang_offset + ir] = outy;
            outboundFluxZ[iang_offset + ir] = outz;
            nextFlux_block[iang_offset + ir] = metadata->quad_wt[iang] * angFlux;

        } // for(voxel_ind_local=vbegin ...)


        // Synhcronize threads in the block before next stage of the voxel-sweep
        __syncthreads();
    }// END: for for(int stage=0; ...)

}


// Anisotropic Solver 
SOL_T* gsSolverHarmonicGPU(struct solver_metadata metadata, const std::vector<SOL_T>* uFlux_vector)
{
    // Do some input checks
    if (metadata.pn > 8)
    {
        std::cout << " Pn Check failed, Pn = " << metadata.pn << std::endl;
        std::cout << " Maximum Pn = 8 " << std::endl;
        return NULL;
    }

    struct solver_metadata* metadata_device;
    struct sweep_pattern sweep_temp, * sweep_device[8], ** sweep_list_device;
    int ie, deviceNum = 0; // gpu-id

    // The number of flux moments used in Spherial Harmonic expansion
    const unsigned int momentCount = (metadata.pn + 1) * (metadata.pn + 1);
    std::cout << "MomentCount = " << momentCount << std::endl;

    // device pointers - allocate on device. NO copy from host-->device needed
    SOL_T* outboundFluxX, *outboundFluxY, *outboundFluxZ; //size: #angles * #voxels. (Optional) re-set to 0 prior to every LBTE iteration
    SOL_T* nextFlux_block; //size: #angles * #voxels . Should be re-set to 0 prior to every LBTE iteration

    SOL_T* nextFlux;     // size: #voxels. accumulated flux over all angles.
                         // reset to 0 prior to every LBTE iteration.
    SOL_T* nextMoment;   // size: #voxels * #moments. update moments over all angles.
                         // reset to 0 prior to every LBTE iteration.
    SOL_T* totalCrossSection; // size:#voxels. Meeds to be compute once every energy-lvel.

    SOL_T* uFlux_device; // size: #energies * #voxels. Just copy once from host-->device (before the LBTE solver starts).
    SOL_T* cFlux_device; // size: #energies * #voxels. Just intialize once to 0 (before LBTE solver starts)
                         // After every LBTE iteration update this based on nextFlux.

    SOL_T* umoments_device; // size: #energies * #voxels * momentCount. Need to compute on device (before the LBTE solver starts).
    SOL_T* cmoments_device; // size: #energies * #voxels * momentCount. Just intialize once to 0 (before LBTE solver starts)
    SOL_T* inmoments_device; // size: #voxels * momentCount. Sum of moments over energy groups. Just intialize once to 0 (before LBTE solver starts)

    SOL_T* diff_device, *norm_diff_device;  // after every LBTE iteration copy this value from device-->host
    

    // host pointers
    SOL_T* cFlux_host = new SOL_T[metadata.num_energies * metadata.num_voxels]{};
    SOL_T* uFlux_host; // use getufluxArray to set this
    SOL_T* cmoments_host = new SOL_T[metadata.num_energies * metadata.num_voxels * momentCount]{}; // collided flux moments
   

    // Convergence parameters
    int maxIterations = 25;
    float epsilon = 0.001;

    // variables for tracking convergenxe
    int i, it = 0;
    SOL_T maxDiff = 0.0, totDiff = 0.0, totDiffPre = 0.0, totDiff2 = 0.0, rmsDiff = 0.0, avgDiff = 0.0;
    SOL_T *diff, *norm_diff;

    // Precompute sweep patterns
    std::cout << "Precomputing 3-D voxel-sweep patterns for different directions ..." << std::endl;
    std::vector<sweep> sweep_list = compute_sweep_list(metadata.Nx, metadata.Ny, metadata.Nz);


    if (cudaSetDevice(deviceNum) != cudaSuccess) {
        fprintf(stderr, "Error initializing device %d\n", deviceNum);
        exit(-1);
    }
    else
        printf("Using gpu_device %d\n", deviceNum);

    // Total memory requested on Device
    SOL_T req_mem_device = 0.0;

    //----metadata transfer--
    std::cout << "Transferring metadata to device ..." << std::endl;
    metadata_device = transferMetadataToDevice(metadata);
    req_mem_device  = sizeof(float) * metadata.num_angles * 4 + sizeof(float) * (metadata.Nx + 1) * 3 + sizeof(float) * metadata.num_energies
                    + sizeof(float) * metadata.num_voxels * 2 + sizeof(float) * metadata.num_materials * metadata.num_energies * metadata.num_energies * (metadata.pn + 1)
                    + sizeof(float) * metadata.num_materials * metadata.num_energies + sizeof(struct solver_metadata);

    //---sweep pattern transfer ---
    std::cout << "Transferring sweep patterns to device ..." << std::endl;
    // sweep_device[d] is a pointer that resides on the host but points to a memory location on the device (GPU).
    for (int d = 0; d < 8; d++)
    {
        sweep_temp = convertSweepClassToStruct(sweep_list[d]); // class-->structure
        sweep_device[d] = transferSweepPatternToDevice(sweep_temp); // transfer data of structure to device
    }
    // List of all 8  sweep-patterns represneted by a **sweep_list_device, where sweep_list_device[d], d=0...8-1, is a pointer on the device ...
    // that represents a particular sweep-pattern.
    cudaMalloc((void**)&sweep_list_device, 8 * sizeof(struct sweep_pattern*));
    cudaMemcpy(sweep_list_device, &sweep_device[0], 8 * sizeof(struct sweep_pattern*), cudaMemcpyHostToDevice);
    req_mem_device += 8 * sizeof(struct sweep_pattern*);

    //--allocate arrays on device--
    // initialization / transfer may need to be done repeatedly during the LBTE solver.
    std::cout << "Allocating vectors on device for the LBTE solver ..." << std::endl;
    cudaMalloc((void**)&outboundFluxX, sizeof(SOL_T) * metadata.num_voxels * metadata.num_angles);
    cudaMalloc((void**)&outboundFluxY, sizeof(SOL_T) * metadata.num_voxels * metadata.num_angles);
    cudaMalloc((void**)&outboundFluxZ, sizeof(SOL_T) * metadata.num_voxels * metadata.num_angles);
    cudaMalloc((void**)&nextFlux_block, sizeof(SOL_T) * metadata.num_voxels * metadata.num_angles);
    req_mem_device += 4. * sizeof(SOL_T) * metadata.num_voxels * metadata.num_angles;

    cudaMalloc((void**)&nextFlux, sizeof(SOL_T) * metadata.num_voxels);
    cudaMalloc((void**)&totalCrossSection, sizeof(SOL_T) * metadata.num_voxels);
    req_mem_device += 2. * sizeof(SOL_T) * metadata.num_voxels;

    cudaMalloc((void**)&uFlux_device, sizeof(SOL_T) * metadata.num_voxels * metadata.num_energies);
    cudaMalloc((void**)&cFlux_device, sizeof(SOL_T) * metadata.num_voxels * metadata.num_energies);
    req_mem_device += 2. * sizeof(SOL_T) * metadata.num_voxels * metadata.num_energies;

    cudaMalloc((void**)&cmoments_device, sizeof(SOL_T)* metadata.num_energies* metadata.num_voxels* momentCount);  // collided moments on device
    cudaMalloc((void**)&umoments_device, sizeof(SOL_T)* metadata.num_energies* metadata.num_voxels* momentCount);  // uncollided moments on device
    cudaMalloc((void**)&inmoments_device, sizeof(SOL_T)* metadata.num_voxels* momentCount);                        // sum of c/uc moments over energy
    cudaMalloc((void**)&nextMoment, sizeof(SOL_T) * metadata.num_voxels * momentCount);                            // temp moments in LBTE inner iteration
    req_mem_device += 2. * sizeof(SOL_T) * metadata.num_voxels * metadata.num_energies * momentCount;
    req_mem_device += 2. * sizeof(SOL_T) * metadata.num_voxels * momentCount;
    
    // set (u/c)fluxes on device before LBTE solver begins
    //uFlux_host = getufluxArray(uFlux_vector, metadata.num_energies, metadata.num_voxels);  // Convert double vector to float array
    cudaMemcpy(uFlux_device, uFlux_vector->data(), sizeof(SOL_T) * metadata.num_voxels * metadata.num_energies, cudaMemcpyHostToDevice);
    cudaMemset(cFlux_device, 0, sizeof(SOL_T) * metadata.num_voxels * metadata.num_energies);

    // set (u/c)moments on device before LBTE solver begins
    cudaMemset(cmoments_device, 0, sizeof(SOL_T) * metadata.num_voxels * metadata.num_energies * momentCount);
    cudaMemset(umoments_device, 0, sizeof(SOL_T) * metadata.num_voxels * metadata.num_energies * momentCount);
    cudaMemset(inmoments_device, 0, sizeof(SOL_T) * metadata.num_voxels * momentCount);

    // convgerence tracking make it double
    diff = new SOL_T[metadata.num_voxels];
    norm_diff = new SOL_T[metadata.num_voxels];
    cudaMalloc((void**)&diff_device, sizeof(SOL_T) * metadata.num_voxels);
    cudaMalloc((void**)&norm_diff_device, sizeof(SOL_T) * metadata.num_voxels);
    req_mem_device += 2 * sizeof(SOL_T) * metadata.num_voxels;

    // Check total requested memory on Device
    std::cout << "Total requested memory on GPU is: " << req_mem_device / 1024. / 1024. / 1024. << " GB\n" << std::endl;

    // Grid and Block dimensions
    // for precomputing source-term, downscatter and in-scatter
    dim3 block_precomp(BLOCK_SIZE_PRECOMP);
    dim3 grid_precomp(metadata.num_voxels / BLOCK_SIZE_PRECOMP + 1);

    // for LBTE iteration
    dim3 block_LBTE(BLOCK_SIZE_LBTE);
    dim3 grid_LBTE(metadata.num_angles);

    // Grid and Block dimensions for generating umoments
    dim3 block_umoment(1); // XL: current assume a fixed source
    dim3 grid_umoment(metadata.Nx, metadata.Ny, metadata.Nz);
    std::cout << "uGrid: " << grid_umoment.x << "x" << grid_umoment.y << "x" << grid_umoment.z << ", uBlock: "
        << block_umoment.x << "x" << block_umoment.y << std::endl;

    // Grid and Block dimensions for updating cmoments
    dim3 block_cmoment(momentCount);
    dim3 grid_cmoment(metadata.Nx, metadata.Ny, metadata.Nz);
    std::cout << "cGrid: " << grid_cmoment.x << "x" << grid_cmoment.y << "x" << grid_cmoment.z << ", cBlock: "
        << block_cmoment.x << "x" << block_cmoment.y << std::endl;

    bool notTrackedHighestEnergy = true;
    int ira, highestEnergy = 0;
    while (notTrackedHighestEnergy)
    {
        for (ira = 0; ira < metadata.num_voxels; ira++)
        {
            if ((*uFlux_vector)[(unsigned int)(highestEnergy * metadata.num_voxels + ira)] > 0)
                break;
        }
        if (ira == metadata.num_voxels)
        {
            std::cout << "No external source or downscatter, skipping energy group " << highestEnergy << std::endl;
            highestEnergy++;
        }
        else
            notTrackedHighestEnergy = false;

        if (highestEnergy >= metadata.num_energies)
        {
            std::cout << "Zero flux everywhere from the raytracer" << std::endl;
            return NULL;
        }
    }

    // LBTE solver
    //for (ie = highestEnergy; ie < 1; ie++)  //XL: debug
    for (ie = highestEnergy; ie < metadata.num_energies; ie++)
    {
        std::cout << "### Energy group " << ie << " ###" << std::endl;
        // LBTE iteration count for given energy
        it = 0;
        maxDiff = 1.0; // just set so that it enters the loop
        rmsDiff = 1.0;
        totDiff = 1.0e10;
        totDiffPre = 1.0e11;

        // Precompute total scatter cross-section
        compute_totalCrossSection <<< grid_precomp, block_precomp >>> (totalCrossSection, metadata_device, ie);

        // Initialize the inmoments to zero for every energy group
        cudaMemset(inmoments_device, 0, sizeof(SOL_T) * metadata.num_voxels * momentCount);

        // Update uncollided flux moment for energy group ie and generate in-moments
        update_umomentsKernel <<< grid_umoment, block_umoment >>> (uFlux_device, metadata_device, umoments_device, cmoments_device,
                                                                  inmoments_device, momentCount, highestEnergy, ie);                                                                 
        cudaDeviceSynchronize();

        // Single LBTE iteration
        while ((it < maxIterations) && (maxDiff > epsilon))
       //while ((it < maxIterations) && (maxDiff > epsilon) && (totDiff/totDiffPre < 1.0))  //XL: Debug use rms for convergence
        {
            // Prior to every LBTE iteration (before sweep of angular-direction and voxel-grid) do the following re-initialization.
            // initialization required every iteration
            cudaMemset(nextFlux_block, 0, sizeof(SOL_T) * metadata.num_voxels * metadata.num_angles);
            cudaMemset(nextFlux, 0, sizeof(SOL_T) * metadata.num_voxels);
            cudaMemset(nextMoment, 0, sizeof(SOL_T) * metadata.num_voxels * momentCount);

            // this initialization should be done once before LBTE solver. re-initialization before every LBTE iteration is optional.
            cudaMemset(outboundFluxX, 0, sizeof(SOL_T) * metadata.num_voxels * metadata.num_angles);
            cudaMemset(outboundFluxY, 0, sizeof(SOL_T) * metadata.num_voxels * metadata.num_angles);
            cudaMemset(outboundFluxZ, 0, sizeof(SOL_T) * metadata.num_voxels * metadata.num_angles);

            // Pre-compute in-scatter from collided-flux
            //compute_inscatter_Iso << < grid_precomp, block_precomp >> > (inscatter, cFlux_device, metadata_device, ie);

            // LBTE - angular-sweep and voxel-sweep
            LBTE_Iteration_Harmonic <<< grid_LBTE, block_LBTE >>> (nextFlux_block, totalCrossSection, cmoments_device, inmoments_device,
                                                                   outboundFluxX, outboundFluxY, outboundFluxZ, momentCount,
                                                                   sweep_list_device, metadata_device, ie);
            cudaDeviceSynchronize();
            // At end of LBTE iteration  do ..
            // 1) nextFlux <-- sum nextFlux_block over all angles
            // 2) nextMoment  <-- Update collided flux moment using angular flux over all angles
            // 3) find the difference between previous and current iteration LBTE scalar flux solution (to check for convergence)
            // 4) cFlux_device[energy, :] <-- nextFlux
            // 5) cmoment_device[energy, :]  <-- nextMoment
            sumFluxAcrossAngles <<< grid_precomp, block_precomp >>> (nextFlux, nextFlux_block, metadata_device);
            updateFluxMoment <<< grid_precomp, block_precomp >>> (nextMoment, nextFlux_block, metadata_device, momentCount);
            updateDiff <<< grid_precomp, block_precomp >>> (norm_diff_device, diff_device, nextFlux, cFlux_device, metadata_device, ie);
            updatecFluxAtEnergy <<< grid_precomp, block_precomp >>> (cFlux_device, nextFlux, metadata_device, ie);
            updatecMomentAtEnergy <<< grid_cmoment, block_cmoment >>> (cmoments_device, nextMoment, metadata_device, ie, momentCount);

            // Track convergence here
            cudaMemcpy(diff, diff_device, sizeof(SOL_T) * metadata.num_voxels, cudaMemcpyDeviceToHost);
            cudaMemcpy(norm_diff, norm_diff_device, sizeof(SOL_T) * metadata.num_voxels, cudaMemcpyDeviceToHost);
            
            totDiffPre = totDiff;
            totDiff = 0.0;
            maxDiff = -1.0e35f;
            totDiff2 = 0.0;
            for (i = 0; i < metadata.num_voxels; i++)
            {
                maxDiff = std::max(maxDiff, norm_diff[i]);
                totDiff += norm_diff[i];
                totDiff2 += norm_diff[i] * norm_diff[i]; // Sum of squares
            }
            avgDiff = totDiff / metadata.num_voxels;
            rmsDiff = sqrt(totDiff2 / metadata.num_voxels);
            
            // increment iteration count
            it++;

            std::cout<<"Iteration "<< it << ": Max relative error = " << maxDiff << " avgDiff = " << avgDiff << std::endl;
    
        } // END: while-loop corresponding to LBTE iterations at given energy-group

        if (!(it < maxIterations))
            std::cout << "Max iterations hit" << std::endl;
        else if (!(maxDiff > epsilon))
            std::cout << "Converged on relative error" << std::endl;
        else
            std::cout << "Converged on precision bound" << std::endl;

        // For updating Progress
        std::cout << ((ie + 1) * 100) / metadata.num_energies << "%% finished" << std::endl;

    } //END: for(ie=0; ...)

    // Copy LBTE solution across all energy levels
    cudaMemcpy(cFlux_host, cFlux_device, sizeof(SOL_T) * metadata.num_energies * metadata.num_voxels, cudaMemcpyDeviceToHost);
    //cudaMemcpy(cmoments_host, cmoments_device, sizeof(SOL_T) * metadata.num_energies * metadata.num_voxels * momentCount, cudaMemcpyDeviceToHost);

    // De-allocate memory
    cudaFree(outboundFluxX);
    cudaFree(outboundFluxY);
    cudaFree(outboundFluxZ);
    cudaFree(nextFlux_block);
    cudaFree(cFlux_device);
    cudaFree(uFlux_device);

    cudaFree(nextFlux);
    cudaFree(nextMoment);
    cudaFree(totalCrossSection);

    cudaFree(cmoments_device);
    cudaFree(umoments_device);
    cudaFree(inmoments_device);

    cudaFree(norm_diff_device);
    cudaFree(diff_device);

    free(uFlux_host);

    return cFlux_host;
    //return cmoments_host;
}

__global__ void raytraceKernel(struct solver_metadata* metadata, RAY_T* uflux)
{
    unsigned int xIndxStart = blockIdx.x;
    unsigned int yIndxStart = blockIdx.y;
    unsigned int zIndxStart = blockIdx.z;
    unsigned int sIndxStart = threadIdx.x; 

    const unsigned short DIRECTION_X = 1;
    const unsigned short DIRECTION_Y = 2;
    const unsigned short DIRECTION_Z = 3;

    RAY_T tiny = 1.0E-35f;
    RAY_T huge = 1.0E35f;

    unsigned int ir0 = xIndxStart * metadata->Ny * metadata->Nz + yIndxStart * metadata->Nz + zIndxStart;

    RAY_T* meanFreePaths = new RAY_T[metadata->num_energies];

    // Assume uniform mesh
    float dx = metadata->XNodes[1] - metadata->XNodes[0]; 
    float dy = metadata->YNodes[1] - metadata->YNodes[0];
    float dz = metadata->ZNodes[1] - metadata->ZNodes[0];
    
    // Real position (x, y, z) of a voxel
    float x = metadata->XNodes[xIndxStart] + dx / 2;
    float y = metadata->YNodes[yIndxStart] + dy / 2;
    float z = metadata->ZNodes[zIndxStart] + dz / 2;

    // Start raytracing through the geometry
    int xIndx = xIndxStart;
    int yIndx = yIndxStart;
    int zIndx = zIndxStart;
    int sIndx = sIndxStart;

    RAY_T srcToCellX = metadata->sx - x;
    RAY_T srcToCellY = metadata->sy - y;
    RAY_T srcToCellZ = metadata->sz - z;

    RAY_T srcToCellDist = sqrt(srcToCellX * srcToCellX + srcToCellY * srcToCellY + srcToCellZ * srcToCellZ);
    RAY_T srcToPtDist = 0.0; // Used later in the while loop

    RAY_T xcos = srcToCellX / srcToCellDist;
    RAY_T ycos = srcToCellY / srcToCellDist;
    RAY_T zcos = srcToCellZ / srcToCellDist;

    int xBoundIndx = (xcos > 0 ? xIndx + 1 : xIndx);
    int yBoundIndx = (ycos > 0 ? yIndx + 1 : yIndx);
    int zBoundIndx = (zcos > 0 ? zIndx + 1 : zIndx);

    // Clear MFP array to zeros
    for (unsigned int i = 0; i < metadata->num_energies; i++)
    {
        meanFreePaths[i] = 0.0f;
    }

    bool exhaustedRay = false;
    int iter = 0;

    while (!exhaustedRay)
    {
        unsigned int ir = xIndx * metadata->Ny * metadata->Nz + yIndx * metadata->Nz + zIndx;

        // Determine the distance to cell boundaries
        RAY_T tx = (fabs(xcos) < tiny ? huge : (metadata->XNodes[xBoundIndx] - x) / xcos); // Distance traveled when next cell is entered
        RAY_T ty = (fabs(ycos) < tiny ? huge : (metadata->YNodes[yBoundIndx] - y) / ycos);
        RAY_T tz = (fabs(zcos) < tiny ? huge : (metadata->ZNodes[zBoundIndx] - z) / zcos);

        // Determine the shortest distance traveled [cm] before any surface is crossed
        RAY_T tmin;
        unsigned short dirHitFirst;
        
        if (tx < ty && tx < tz)
        {
            tmin = tx;
            dirHitFirst = DIRECTION_X;
        }
        else if (ty < tz)
        {
            tmin = ty;
            dirHitFirst = DIRECTION_Y;
        }
        else
        {
            tmin = tz;
            dirHitFirst = DIRECTION_Z;
        }

        // Calculate distance from cell to source
        srcToCellX = metadata->sx - x;
        srcToCellY = metadata->sy - y;
        srcToCellZ = metadata->sz - z;
        srcToPtDist = sqrt(srcToCellX * srcToCellX + srcToCellY * srcToCellY + srcToCellZ * srcToCellZ);
        if (srcToPtDist < tmin) // Determine if cell is source cell
        {
            tmin = srcToPtDist;
            exhaustedRay = true;
            dirHitFirst = 0;  // assign a null direction
        }

        // Update mfp array
        unsigned int zid = metadata->zoneIdMapped[ir];
        for (unsigned int ie = 0; ie < metadata->num_energies; ie++)
        {
            //                   [cm] * [b]                                                      * [atom/b-cm]
            meanFreePaths[ie] += tmin * metadata->sigma_total[zid * metadata->num_energies + ie] * metadata->atomDensity[ir];
        }

        // Update cell indices and positions
        if (dirHitFirst == DIRECTION_X) // x direction
        {
            x = metadata->XNodes[xBoundIndx];
            y += tmin * ycos;
            z += tmin * zcos;
            if (xcos > 0)
            {
                xIndx++;
                xBoundIndx = xIndx + 1;
            }
            else
            {
                xIndx--;
                xBoundIndx = xIndx;
            }
            if (xIndx < 0 || xIndx >(metadata->Nx - 1))
                exhaustedRay = true;
        }
        else if (dirHitFirst == DIRECTION_Y) // y direction
        {
            x += tmin * xcos;
            y = metadata->YNodes[yBoundIndx];
            z += tmin * zcos;
            if (ycos > 0)
            {
                yIndx++;
                yBoundIndx = yIndx + 1;
            }
            else
            {
                yIndx--;
                yBoundIndx = yIndx;
            }
            if (yIndx < 0 || yIndx >(metadata->Ny - 1))
                exhaustedRay = true;
        }
        else if (dirHitFirst == DIRECTION_Z) // z direction
        {
            x += tmin * xcos;
            y += tmin * ycos;
            z = metadata->ZNodes[zBoundIndx];
            if (zcos > 0)
            {
                zIndx++;
                zBoundIndx = zIndx + 1;
            }
            else
            {
                zIndx--;
                zBoundIndx = zIndx;
            }
            if (zIndx < 0 || zIndx >(metadata->Nz - 1))
                exhaustedRay = true;
        }

        iter++;

    } // End of while loop

    for (unsigned int ie = 0; ie < metadata->num_energies; ie++)
    {
        RAY_T flx = metadata->sIntensity[ie] * exp(-meanFreePaths[ie]) / (4 * 3.1415926 * srcToCellDist * srcToCellDist);
        uflux[ie * metadata->num_voxels + ir0] += static_cast<SOL_T>(flx);
    }

    delete[] meanFreePaths;
}

RAY_T* raytraceIsoGPU(struct solver_metadata metadata)
{
    struct solver_metadata* metadata_device;

    int deviceNum = 0; // gpu-id

    // device pointers - allocate on device. NO copy from host-->device needed
    RAY_T* uFlux_device; // size: #energies * #voxels. Just copy once from host-->device (before the LBTE solver starts).

    // host pointers
    RAY_T* uFlux_host = new RAY_T[metadata.num_energies * metadata.num_voxels];

    if (cudaSetDevice(deviceNum) != cudaSuccess) {
        fprintf(stderr, "Error initializing device %d\n", deviceNum);
        exit(-1);
    }
    else
        printf("Using gpu_device %d\n", deviceNum);

    //----metadata transfer--
    std::cout << "Transferring metadata to device ..." << std::endl;
    metadata_device = transferMetadataToDevice(metadata);

    //--allocate arrays on device--
    std::cout << "Allocating vectors on device for the Ray Tracing ..." << std::endl;
    cudaMalloc((void**)&uFlux_device, sizeof(RAY_T) * metadata.num_voxels * metadata.num_energies);
    cudaMemset(uFlux_device, 0, sizeof(RAY_T) * metadata.num_voxels * metadata.num_energies);

    // Grid and Block dimensions
    dim3 block_raytrace(1); // Current assume a fixed isotropic source
    dim3 grid_raytrace(metadata.Nx, metadata.Ny, metadata.Nz);
    std::cout << "Grid: " << grid_raytrace.x << "x" << grid_raytrace.y << "x" << grid_raytrace.z << ", Block: "
        << block_raytrace.x << "x" << block_raytrace.y << std::endl;

    // Launch raytace kernel
    raytraceKernel <<<grid_raytrace, block_raytrace >>>(metadata_device, uFlux_device);
    
    // Copy raytrace solution across all energy levels
    cudaMemcpy(uFlux_host, uFlux_device, sizeof(RAY_T) * metadata.num_energies * metadata.num_voxels, cudaMemcpyDeviceToHost);

    // De-allocate memory
    cudaFree(uFlux_device);
    
    cudaFree(metadata_device);

    // write collided-flux out to file extaneral to this routine
    return uFlux_host;
}