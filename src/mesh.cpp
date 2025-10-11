#define _USE_MATH_DEFINES
#include <cmath>
#include <cassert>

#include <iostream>
#include <cassert>

#include "mesh.h"

Mesh::Mesh() : material(0)
{

}

void Mesh::uniform(const int xelems, const int yelems, const int zelems,
                   const float xLen, const float yLen, const float zLen)
{
    xElemCt = xelems;
    yElemCt = yelems;
    zElemCt = zelems;

    // Calculate the number of nodes
    xNodeCt = xElemCt + 1;
    yNodeCt = yElemCt + 1;
    zNodeCt = zElemCt + 1;

    // Allocate space
    xNodes.resize(xNodeCt);
    yNodes.resize(yNodeCt);
    zNodes.resize(zNodeCt);

    dx.resize(xElemCt);
    dy.resize(yElemCt);
    dz.resize(zElemCt);

    vol.resize(xElemCt * yElemCt * zElemCt);

    m_maxX = xLen;
    m_maxY = yLen;
    m_maxZ = zLen;

    // The coordinates between mesh elements iterate the xMesh+1 to get the last bin
    for(unsigned int i = 0; i < xElemCt+1; i++)
        xNodes[i] = i * (xLen / xElemCt);

    for(unsigned int i = 0; i < yElemCt+1; i++)
        yNodes[i] = i * (yLen / yElemCt);

    for(unsigned int i = 0; i < zElemCt+1; i++)
        zNodes[i] = i * (zLen / zElemCt);

    // Calculate the segment sizes
    for(unsigned int i = 0; i < xElemCt; i++)
        dx[i] = xNodes[i+1] - xNodes[i];

    for(unsigned int i = 0; i < yElemCt; i++)
        dy[i] = yNodes[i+1] - yNodes[i];

    for(unsigned int i = 0; i < zElemCt; i++)
        dz[i] = zNodes[i+1] - zNodes[i];

    for(unsigned int xIndx = 0; xIndx < xElemCt; xIndx++)
        for(unsigned int yIndx = 0; yIndx < yElemCt; yIndx++)
            for(unsigned int zIndx = 0; zIndx < zElemCt; zIndx++)
                vol[xIndx * yElemCt*zElemCt + yIndx * zElemCt + zIndx] =
                                         dx[xIndx] * dy[yIndx] * dz[zIndx];

    zoneId.resize(xElemCt * yElemCt * zElemCt, 0);  // Initialize to all zeros

}

void Mesh::calcAreas(const Quadrature *quad, const int eGroups)
{
    Axy.resize(eGroups * quad->angleCount() * xElemCt * yElemCt);
    Ayz.resize(eGroups * quad->angleCount() * yElemCt * zElemCt);
    Axz.resize(eGroups * quad->angleCount() * xElemCt * zElemCt);

    // Calculate the cell face area for each angle as well as volume
    for(int eIndx = 0; eIndx < eGroups; eIndx++)
    {
        for(unsigned int qIndx = 0; qIndx < quad->angleCount(); qIndx++)
        {
            float vMu = fabs(quad->mu[qIndx]);
            float vEta = fabs(quad->eta[qIndx]);
            float vZi = fabs(quad->zi[qIndx]);

            for(unsigned int yIndx = 0; yIndx < yElemCt; yIndx++)
                for(unsigned int zIndx = 0; zIndx < zElemCt; zIndx++)
                {
                    Ayz[eIndx*quad->angleCount()*yElemCt*zElemCt +
                         qIndx*yElemCt*zElemCt + yIndx*zElemCt + zIndx] =
                                                  2 * vMu * dy[yIndx] * dz[zIndx];
                }

            for(unsigned int xIndx = 0; xIndx < xElemCt; xIndx++)
                for(unsigned int zIndx = 0; zIndx < zElemCt; zIndx++)
                {
                    Axz[eIndx*quad->angleCount()*xElemCt*zElemCt +
                         qIndx*xElemCt*zElemCt + xIndx*zElemCt + zIndx] =
                                                 2 * vZi * dx[xIndx] * dz[zIndx];
                }

            for(unsigned int xIndx = 0; xIndx < xElemCt; xIndx++)
                for(unsigned int yIndx = 0; yIndx < yElemCt; yIndx++)
                {
                    Axy[eIndx*quad->angleCount()*xElemCt*yElemCt +
                         qIndx*xElemCt*yElemCt + xIndx*yElemCt + yIndx] =
                                                 2 * vEta * dx[xIndx] * dy[yIndx];
                }
        }
    }
}

void Mesh::initCtVariables()
{
    assert(vol.size() >= 1 && "Cannot initialize CT variables before volume data!");

    ct.resize(vol.size());
    density.resize(vol.size());
    atomDensity.resize(vol.size());
}

int Mesh::getFlatIndex(unsigned int xindx, unsigned int yindx, unsigned int zindx) const
{
    assert(xindx < xElemCt && "x-index was too large!");
    assert(yindx < yElemCt && "y-index was too large!");
    assert(zindx < zElemCt && "z-index was too large!");

    return xindx * yElemCt * zElemCt + yindx * zElemCt + zindx;
}

int Mesh::getZoneIdAt(unsigned int xindx, unsigned int yindx, unsigned int zindx) const 
{
    assert(xindx < xElemCt && "x-index was too large!");
    assert(yindx < yElemCt && "y-index was too large!");
    assert(zindx < zElemCt && "z-index was too large!");

    return zoneId[getFlatIndex(xindx, yindx, zindx)];
}

float Mesh::getPhysicalDensityAt(unsigned int xindx, unsigned int yindx, unsigned int zindx) const 
{
    assert(xindx < xElemCt && "x-index was too large!");
    assert(yindx < yElemCt && "y-index was too large!");
    assert(zindx < zElemCt && "z-index was too large!");

    return density[getFlatIndex(xindx, yindx, zindx)];
}

float Mesh::getAtomDensityAt(unsigned int xindx, unsigned int yindx, unsigned int zindx) const 
{
    assert(xindx < xElemCt && "x-index was too large!");
    assert(yindx < yElemCt && "y-index was too large!");
    assert(zindx < zElemCt && "z-index was too large!");

    return atomDensity[getFlatIndex(xindx, yindx, zindx)];
}
