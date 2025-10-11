
#include <stdlib.h>
#include <iostream>
#include <vector>
#include "sweep.h"

sweep::sweep(int d, int nx, int ny, int nz)
{
  Nx=nx;
  Ny=ny;
  Nz=nz;

  if((d<0) || (d>7))
  {
    std::cout<<"Error in sweep::sweep(): argument d must be in the range 0-->7"<<std::endl;
    exit(-1);
  }

  // Map: d=0-->(-1,-1,-1) and d=7-->(+1,+1,+1) based on a 2*2*2 cube
  dix = 2*(d/4)-1;
  diy = 2*((d%4)/2)-1;
  diz = 2*((d%4)%2)-1;

  // derived
  xjmp = Ny*Nz;
  yjmp = Nz;

  // Compute sweep pattern
  compute_pattern();
}

int sweep::flatten(int ix, int iy, int iz)
{
  return (ix*xjmp + iy*yjmp + iz);
}

std::vector<int> sweep::getNeighbors(int ind)
{
  std::vector<int> neighbors;

  int ix = ind/xjmp;
  int iy = (ind%xjmp)/yjmp;
  int iz = (ind%xjmp)%yjmp;

  if(((ix+dix)>0) && ((ix+dix)<Nx))
    neighbors.push_back(flatten(ix+dix,iy,iz));

  if(((iy+diy)>0) && ((iy+diy)<Ny))
    neighbors.push_back(flatten(ix,iy+diy,iz));

  if(((iz+diz)>0) && ((iz+diz)<Nz))
    neighbors.push_back(flatten(ix,iy,iz+diz));

  return neighbors;
}

bool sweep::num_not_present_in_list(std::vector<int> num_list, int num)
{
  bool flag = true;
  for(unsigned int i=0; i < num_list.size(); i++)
  {
    if(num_list[i] == num)
    {
      flag = false;
      break;
    }
  }
  return flag;
}


void sweep::compute_pattern(void)
{
  int ixStart = (dix > 0) ? 0:(Nx-1);
  int iyStart = (diy > 0) ? 0:(Ny-1);
  int izStart = (diz > 0) ? 0:(Nz-1);

  std::vector<int> pt_list; // each pint in the list is "flattened", i.e. represented by single int rather than 3 coordinates
  std::vector<int> pt_list_new;
  std::vector<int> neighbors;

  // initialize: 1st stage of sweep has only 1 voxel
  pt_list.push_back(flatten(ixStart, iyStart, izStart));

  //std::cout<<"---Begin computation of sweep pattern--"<<std::endl;

  while(pt_list.size()!=0)
  {
    //std::cout<<"Sweep size = "<<pt_list.size()<<std::endl;

    // Insert voxels for a given stage into pattern
    pattern.push_back(pt_list);

    // grow list of voxels to sweep based on 3-pt neighborhoood
    for(unsigned int j=0; j < pt_list.size(); j++)
    {
      neighbors = getNeighbors(pt_list[j]);
      for(unsigned int k=0; k < neighbors.size(); k++)
      {
        if(num_not_present_in_list(pt_list_new, neighbors[k]))
          pt_list_new.push_back(neighbors[k]);
      }
    }
    pt_list = pt_list_new;
    pt_list_new.clear();
  }
}


std::vector<sweep> compute_sweep_list(int nx, int ny, int nz)
{
  std::vector<sweep> sweep_list;
  for(int d=0; d<8; d++)
  {
    sweep sweep_obj(d, nx, ny, nz);
    sweep_list.push_back(sweep_obj);
  }
  return sweep_list;
}
