

#include <iostream>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "quadrature.h"
#include "mesh.h"
#include "xsection.h"
#include "sourceparams.h"

#include "solver_metadata.h"
#include "solver_metadata_extract.h"

void solver_metadata_extract(struct solver_metadata *metadata, const Quadrature* quad, const Mesh* mesh,
                             const XSection* xs, const SourceParams *srcPar, int pn)
{
  std::vector<int> zoneId_list;
  int im_step, src_step, sink_step;
  int ir, im, ie_src, ie_sink, il;
  unsigned int k;

  int num_angles = (int)quad->angleCount();
  float *quad_wt = new float[num_angles];
  float *mu = new float[num_angles];
  float *zi = new float[num_angles];
  float *eta = new float[num_angles];

  float* XNodes = new float[mesh->xNodeCt];  
  for (int inode = 0; inode < (int)mesh->xNodeCt; inode++)
  {
      XNodes[inode] = mesh->xNodes[inode];
  }
  float* YNodes = new float[mesh->yNodeCt];  
  for (int inode = 0; inode < (int)mesh->yNodeCt; inode++)
  {
      YNodes[inode] = mesh->yNodes[inode];
  }
  float* ZNodes = new float[mesh->zNodeCt];  
  for (int inode = 0; inode < (int)mesh->zNodeCt; inode++)
  {
      ZNodes[inode] = mesh->zNodes[inode];
  }


  for(int qIndx = 0; qIndx < num_angles; qIndx++)
  {
      quad_wt[qIndx] = quad->wt[qIndx];
      mu[qIndx] = quad->mu[qIndx];
      zi[qIndx] = quad->zi[qIndx];
      eta[qIndx] = quad->eta[qIndx];
  }


  int num_voxels = (int)mesh->voxelCount();
  float voxel_vol  = mesh->vol[0]; // uniform grid assumption
  float voxel_area_yz = mesh->dy[0] * mesh->dz[0];
  float voxel_area_xz = mesh->dx[0] * mesh->dz[0];
  float voxel_area_xy = mesh->dx[0] * mesh->dy[0];


  float *atomDensity = new float[num_voxels];
  int *zoneIdMapped = new int[num_voxels];

  for(ir=0; ir<num_voxels; ir++)
    atomDensity[ir] = mesh->atomDensity[ir];

  // initialize.
  zoneId_list.push_back(mesh->zoneId[0]); // list of unique material-IDs.
                                          // if zoneId is unqiue, insert it into list. The index into the list (here 0) is the mapped index for this zoneId (here mesh-zoneId[0]).
  zoneIdMapped[0] = 0; // voxel-wise mapped index. maps voxel index --> material ID.
  // loop over all remaining voxels
  for(ir=1; ir<num_voxels; ir++)
  {
    // check if zoneId of this voxel is present in the list, if not append it
    // should be a short list, so we can use a simple search routine
    for(k=0; k<zoneId_list.size(); k++)
    {
      if(zoneId_list[k] == mesh->zoneId[ir])
        break;
    }
    // no match found
    if(k==zoneId_list.size())
      zoneId_list.push_back(mesh->zoneId[ir]);

    // assigned a mapped id for this voxel
    zoneIdMapped[ir] = k;
  }

  int num_materials = (int)zoneId_list.size();
  int num_energies = (int)xs->groupCount();
  float *sigma_s = new float[num_materials * num_energies * num_energies * (pn+1)];
  float *sigma_total = new float[num_materials * num_energies];

  im_step = num_energies * num_energies * (pn+1);
  src_step = num_energies * (pn+1);
  sink_step = pn+1;

  // Representation of 2-D arrays: GPUs prefer flattened 1-D arrays to a an array of 2-D arrays (unequal lengths)
  for(im=0; im<num_materials; im++)
  for(ie_src=0; ie_src<num_energies; ie_src++)
  {
    // energies are in decending order
    // No downscatter when sink-energy is greater than source-energy.
    memset(&sigma_s[im * im_step + ie_src * src_step],  0, sizeof(float)*(ie_src * sink_step));
    // Downscatter since sink energy is less than or equal to source-energy.
    for(ie_sink=ie_src; ie_sink<num_energies; ie_sink++)
    {
      for(il=0; il<=pn; il++)
        sigma_s[im * im_step + ie_src * src_step + ie_sink * sink_step + il] = xs->scatXs2d(zoneId_list[im], ie_src, ie_sink, il);
    }

    sigma_total[im * num_energies + ie_src] = xs->totXs1d(zoneId_list[im], ie_src);
  }

  // Source intensity
  float* sIntensity = new float[num_energies];
  for (int ieg = 0; ieg < num_energies; ieg++)
  {
      sIntensity[ieg] = srcPar->spectraIntensity[ieg];
  }

  // Set structure
  metadata->num_materials = num_materials;
  metadata->num_angles = num_angles;
  metadata->num_voxels = num_voxels;
  metadata->num_energies = num_energies;
  metadata->voxel_vol = voxel_vol;
  metadata->voxel_area_yz = voxel_area_yz ;
  metadata->voxel_area_xz = voxel_area_xz ;
  metadata->voxel_area_xy = voxel_area_xy ;
  metadata->pn = pn;
  metadata->Nx = (int)mesh->xElemCt;
  metadata->Ny = (int)mesh->yElemCt;
  metadata->Nz = (int)mesh->zElemCt;
  metadata->sx = srcPar->sourceX;  
  metadata->sy = srcPar->sourceY;  
  metadata->sz = srcPar->sourceZ;  
  metadata->sIntensity = sIntensity; 
  metadata->XNodes = XNodes;       
  metadata->YNodes = YNodes;       
  metadata->ZNodes = ZNodes;       
  metadata->sphi = srcPar->sourcePhi; 
  metadata->stheta = srcPar->sourceTheta; 
  metadata->stype = srcPar->sourceType;  
  metadata->isox = srcPar->isoX;         
  metadata->isoy = srcPar->isoY;         
  metadata->isoz = srcPar->isoZ;        

  printf("#materials = %d, #angles=%d, #voxels=%d and #energies=%d \n", metadata->num_materials, metadata->num_angles, metadata->num_voxels, metadata->num_energies);
  printf("(Nx, Ny, Nz) = (%d, %d, %d) \n", metadata->Nx, metadata->Ny, metadata->Nz);

  metadata->zoneIdMapped = zoneIdMapped;
  metadata->mu = mu;
  metadata->zi = zi;
  metadata->eta = eta;
  metadata->quad_wt = quad_wt;
  metadata->atomDensity = atomDensity;
  metadata->sigma_s = sigma_s;
  metadata->sigma_total = sigma_total;

}
