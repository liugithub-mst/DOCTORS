#ifndef MESH_H_
#define MESH_H_

#include <vector>

#include "quadrature.h"

// VS: Uniform spatial grid, here x is slowest changing in (x,y,z)

class Mesh
{
public:
    Mesh();

    int material;

    unsigned int xElemCt;       // Number of voxels in x dimension
    unsigned int yElemCt;       // Number of voxels in y dimension
    unsigned int zElemCt;       // Number of voxels in z dimension

    unsigned int xNodeCt;       // Number of nodes in x dimension ~xElemCt+1
    unsigned int yNodeCt;       // Number of nodes in y dimension
    unsigned int zNodeCt;       // Number of nodes in z dimension

    std::vector<float> xNodes;  // Coordinates between mesh elements in x
    std::vector<float> yNodes;  // Coordinates between mesh elements in y
    std::vector<float> zNodes;  // Coordinates between mesh elements in z

    std::vector<float> dx;      // mesh segment size in x
    std::vector<float> dy;      // mesh segment size in y
    std::vector<float> dz;      // mesh segment size in z

    std::vector<float> Axy;     // Replaces DA (Already has the x2 factored in
    std::vector<float> Ayz;     // Replaces DB
    std::vector<float> Axz;     // Replaces DC

    // Always initialized
    std::vector<int> zoneId;    // related to material identity
    std::vector<float> vol;     // volume of each voxel

    // Only initialized when reading CT data
    std::vector<unsigned short> ct;  // 16 bit unsigned values 0...65535

    // density in [g/cm^3]
    std::vector<float> density;

    // atom density in [atom/b-cm]
    std::vector<float> atomDensity;

    inline unsigned long long int voxelCount() const
                                     { return xElemCt * yElemCt * zElemCt; }

    inline unsigned long int xjmp() const { return yElemCt * zElemCt; }
    inline unsigned int yjmp() const { return zElemCt; }

    inline float getMaxX() const { return m_maxX; }
    inline float getMaxY() const { return m_maxY; }
    inline float getMaxZ() const { return m_maxZ; }

    void calcAreas(const Quadrature *quad, const int eGroups);
    void initCtVariables();

    int getFlatIndex(unsigned int xindx, unsigned int yindx, unsigned int zindx) const;
    int getZoneIdAt(unsigned int xindx, unsigned int yindx, unsigned int zindx) const; //VS: included unsigned
    float getPhysicalDensityAt(unsigned int xindx, unsigned int yindx, unsigned int zindx) const; // VS: included unsigned
    float getAtomDensityAt(unsigned int xindx, unsigned int yindx, unsigned int zindx) const; //VS: included unsigned

    void uniform(const int xelems, const int yelems, const int zelems,
                 const float xLen, const float yLen, const float zLen);

    // Added below utilities to set the physical properties directly such as zoneId (material Id) and atomic-density
    // Needed for python interface to this C++ backend
    // Assume ::uniform() or equivalent has been run prior to below
    // set voxel-wise material label
    inline void set_zoneId(int *material_id)
    {
        unsigned long long int voxel_count =  voxelCount();
        zoneId.resize(voxel_count, 0);
        for(unsigned long long int i=0; i<voxel_count; i++)
            zoneId[i] = material_id[i];
    }
    // set voxel-wise atomic density (scaled down by 1E+24)
    inline void set_atomDensity(float *atomDensity_scaled)
    {
        unsigned long long int voxel_count =  voxelCount();
        atomDensity.resize(voxel_count, 0);
        for(unsigned long long int i=0; i<voxel_count; i++)
            atomDensity[i] = atomDensity_scaled[i];
    }
    // added just for sake of completeness, not needed. set voxel-wise density.
    inline void set_density(float *massDensity)
    {
        unsigned long long int voxel_count =  voxelCount();
        density.resize(voxel_count, 0);
        for(unsigned long long int i=0; i<voxel_count; i++)
            density[i] = massDensity[i];
    }

private:
    float m_maxX;  // Physical length in x
    float m_maxY;  // Physical length in y
    float m_maxZ;  // Physical length in z

};
#endif // MESH_H_
