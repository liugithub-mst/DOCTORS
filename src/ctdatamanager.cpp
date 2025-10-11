#define _USE_MATH_DEFINES
#include <cmath>
#include <fstream>
#include <cstdint>
#include <iostream>

#include "ctdatamanager.h"
#include "materialutils.h"
#include "histogram.h"

CtDataManager::CtDataManager() : m_valid(false), m_mesh(nullptr)
{

}

CtDataManager::~CtDataManager()
{

}

Mesh* CtDataManager::parse16(int xbins, int ybins, int zbins, double xlen,
                             double ylen, double zlen, std::string filename)
{
    m_mesh = new Mesh;
    unsigned long tbins = xbins * ybins * zbins;  // total bins

    std::cout << "Reading CT data ..." << std::endl;
    std::cout << "xbins= " << xbins << std::endl;
    std::cout << "ybins= " << ybins << std::endl;
    std::cout << "zbins= " << zbins << std::endl;
    std::cout << "xlen = " << xlen << std::endl;
    std::cout << "ylen = " << ylen << std::endl;
    std::cout << "zlen = " << zlen << std::endl;

    m_mesh->uniform(xbins, ybins, zbins, xlen, ylen, zlen);  // uniform mesh grid
    m_mesh->initCtVariables();

    std::vector<unsigned short> zoneIds;
    zoneIds.resize(tbins);

    std::ifstream szChkFin(filename.c_str(), std::ios::binary|std::ios::in|std::ios::ate);
    unsigned long szFin = szChkFin.tellg();
    if(szFin != 2*tbins)
    {
        std::cout << "Mesh Data Size Mismatch!" << std::endl;
        std::cout << "The specified mesh requires " << tbins*2 << " bytes!" << std::endl;
        std::cout << "The binary file only reported " << szFin << " bytes!" << std::endl;

    }
    szChkFin.close();

    std::ifstream fin(filename.c_str(), std::ios::binary);

    if(fin.good())
    {
        for(unsigned int i = 0; i < tbins; i++) 
        {
            if(i % 10000 == 0)
            {
                std::cout << "Reading voxels ... " << i << " done" << std::endl;
            }
            fin.read((char*)&zoneIds[i], 2);
        }
    }
    else
    {
        std::cout << "CT data file reading error!" << std::endl;
        return nullptr;
    }

    int gindx = 0;
    int YX = xbins * ybins;
    int X = xbins;
    for(int i = 0; i < xbins; i++)
        for(int j = 0; j < ybins; j++)
            for(int k = 0; k < zbins; k++)
            {
                m_mesh->ct[gindx] = zoneIds[k*YX + j*X + i];
                gindx++;
            }

    return m_mesh;
}

void CtDataManager::rotateCTvol(Mesh* meshIn, Mesh* meshOut, float rotang)
{
    // Rotate input CT volume by rotang degree
    int NYZ = meshIn->xjmp();
    int NZ = meshIn->yjmp();
    float xold, yold, xnew, ynew;
    int xnearest, ynearest;

    // Initialize output mesh
    meshOut->uniform(meshIn->xElemCt, meshIn->yElemCt, meshIn->zElemCt,
                     meshIn->getMaxX(), meshIn->getMaxY(), meshIn->getMaxZ());
    meshOut->initCtVariables();

    int centerX = meshOut->xElemCt / 2;
    int centerY = meshOut->yElemCt / 2;

    std::cout << "centerX=" << centerX << std::endl;
    std::cout << "centerY=" << centerY << std::endl;

    meshOut->ct.resize(meshOut->voxelCount(),0);

    if(rotang == 0.0 or rotang == 360.0)
    {
       // No change
       meshOut->ct = meshIn->ct;
    }
    else
    {
        //Convert deg to radian
        rotang = rotang * M_PI / 180.0;
        for(int k = 0; k < meshIn->zElemCt; k++)
            for(int j = 0; j < meshIn->yElemCt; j++)
                for(int i = 0; i < meshIn->xElemCt; i++)
                {
                    // Move rotation to the origin
                    xnew = i - centerX;
                    ynew = j - centerY;

                    // Map the point of rotated volume to the old volume
                    xold = xnew * cos(rotang) + ynew * sin(rotang);
                    yold = ynew * cos(rotang) - xnew * sin(rotang);

                    // Round to the nearest point
                    xnearest = round(xold);
                    ynearest = round(yold);

                    // Move back to the original center
                    xnearest += centerX;
                    ynearest += centerY;

                    // Check over boundary
                    if(xnearest > meshIn->xElemCt -1)
                        xnearest = meshIn->xElemCt - 1;
                    else if(xnearest < 0)
                        xnearest = 0;

                    if(ynearest > meshIn->yElemCt -1)
                        ynearest = meshIn->yElemCt - 1;
                    else if(ynearest < 0)
                        ynearest = 0;

                    meshOut->ct[i*NYZ + j*NZ + k] = meshIn->ct[xnearest*NYZ + ynearest*NZ + k];
                }

    }

}

Mesh *CtDataManager::ctNumberToWater(Mesh *mesh)
{
    m_mesh->material = MaterialUtils::WATER;

    std::vector<float> atomPerG;
    std::vector<std::vector<float> > waterFractions;

    // Check the size of the material vector
    if(MaterialUtils::water.size() != MaterialUtils::waterWeights.size())
    {
        std::cout << "Water vectors mismatched: water: " << MaterialUtils::water.size() <<
                     ",  weight: " << MaterialUtils::waterWeights.size() << std::endl;
        return NULL;
    }
    for(unsigned int i = 0; i < MaterialUtils::waterWeights.size(); i++)
    {
        if(MaterialUtils::waterWeights[i].size() != MaterialUtils::waterElements.size())
        {
            std::cout << "Water weight vector mismatched: weight[i]: " <<
                          MaterialUtils::waterWeights[i].size()
                      << ",  elements: " << MaterialUtils::waterElements.size()
                      << " @ index " << i << std::endl;
            return NULL;
        }
    }

    // Convert the weight fractions to atom fractions for each material
    for(unsigned int i = 0; i < MaterialUtils::water.size(); i++)
        waterFractions.push_back(MaterialUtils::weightFracToAtomFrac(
                                 MaterialUtils::waterElements, MaterialUtils::waterWeights[i]));

    // Convert atom fraction to atom density for each material
    for(unsigned int i = 0; i < MaterialUtils::water.size(); i++)
        atomPerG.push_back(MaterialUtils::atomsPerGram(
                           MaterialUtils::waterElements, waterFractions[i]));

    const int offset = 1024;  //FixMe! hardcoded offset
    for(unsigned int i = 0; i < mesh->voxelCount(); i++)
    {
        if(mesh->ct[i] > 65500) // Artifact
            mesh->ct[i] = 0;

        // Raw data minus offset
        int ctv = mesh->ct[i] - offset;  // Underlying data is 0-2500 instead of -1000-1500

        // Determine the density (atom density is after the zoneId calculation)
        if(ctv <= -66)
        {
            mesh->density[i] = 0.001225f;  // Air
            
        }
        else if(ctv <= 60)
        {
            mesh->density[i] = 1.0f;  // Pure water
    
        }
        else
        {   //other mateirals, such as a detector
            mesh->density[i] = 0.451f; 
        }

        // Determine the material
        mesh->zoneId[i] = static_cast<unsigned short>(MaterialUtils::water.size() - 1);  // Last bin is "illegal"
        for(unsigned int j = 0; j < MaterialUtils::water.size()-1; j++)
        {
            if(ctv < MaterialUtils::water[j])
            {
                mesh->zoneId[i] = j;
                break;
            }
        }

        if(mesh->zoneId[i] >= (signed) atomPerG.size())
        {
            std::cout << "CT number to water: EXPLODE!" << std::endl;
        }

        mesh->atomDensity[i] = mesh->density[i] * atomPerG[mesh->zoneId[i]] * 1.0E-24f;  // g/cc * @/g * cm^2/b = @/cm-b

    }

    return mesh;
}

Mesh *CtDataManager::ctNumberToHumanPhantom(Mesh *mesh)
{
    m_mesh->material = MaterialUtils::HOUNSFIELD19;

    std::vector<float> atomPerG;
    std::vector<std::vector<float> > hounsfieldRangePhantom19Fractions;

    // Convert the weight fractions to atom fractions for each material
    for(unsigned int i = 0; i < MaterialUtils::hounsfieldRangePhantom19.size(); i++)
        hounsfieldRangePhantom19Fractions.push_back(MaterialUtils::weightFracToAtomFrac(
                              MaterialUtils::hounsfieldRangePhantom19Elements,
                              MaterialUtils::hounsfieldRangePhantom19Weights[i]));

    // Convert atom fraction to atom density for each material
    for(unsigned int i = 0; i < MaterialUtils::hounsfieldRangePhantom19.size(); i++)
        atomPerG.push_back(MaterialUtils::atomsPerGram(
                            MaterialUtils::hounsfieldRangePhantom19Elements,
                            hounsfieldRangePhantom19Fractions[i]));

    int offset = 1024;  //Fixme! hardcoded
   
    for(unsigned int i = 0; i < mesh->voxelCount(); i++)
    {
        // Raw data minus offset
        int ctv = mesh->ct[i] - offset;  

        // Determine the density (atom density is after the zoneId calculation)
        if(ctv <= 55)
        {
            if(ctv <= -1000)
            {
                mesh->density[i] = 0.001225f; //Air
            }
            else
            {
                mesh->density[i] = 0.0010186f * ctv + 1.013812f;
            }
        }
        else
        {
            mesh->density[i] = 0.000578402f * ctv + 1.103187f;
        }

        // Determine the material
        mesh->zoneId[i] = static_cast<unsigned short>(
                          MaterialUtils::hounsfieldRangePhantom19.size() - 1);  // Last bin is "illegal"
        for(unsigned int j = 0; j < MaterialUtils::hounsfieldRangePhantom19.size()-1; j++)
        {
            if(ctv < MaterialUtils::hounsfieldRangePhantom19[j])
            {
                mesh->zoneId[i] = j;
                break;
            }
        }

        if(mesh->zoneId[i] >= (signed) atomPerG.size())
        {
            std::cout << "CT number to human phantom: EXPLODE!" << std::endl;
        }

        if(mesh->zoneId[i] == 0) // Air
        {
            mesh->density[i] = 0.001225f;  // force air density because it is very sensitive to miscalibrations
        }

        if(mesh->zoneId[i] == 19) // FixMe! assume this material is metal such as iron or Titanium
        {
            mesh->density[i] = 5.24f;  //FixMe! hardcoded iron density
        }

        mesh->atomDensity[i] = mesh->density[i] * atomPerG[mesh->zoneId[i]] * 1.0E-24f;  // g/cc * @/g * cm^2/b = @/cm-b

        if(mesh->atomDensity[i] <= 1.0E-10f)
        {
            std::cout << " CT number to human phantom: Got zero density" << std::endl;
        }
    }

    Histogram h(0, 70000, 10000);
    for(unsigned int i = 0; i < mesh->voxelCount(); i++)
        h.insert(mesh->ct[i]);
    std::ofstream fout;
    fout.open("./ctunits.dat");
    fout << h;
    fout.close();
    std::cout << "CT number to human phantom: Finished histogram" << std::endl;

    return mesh;
}

Mesh *CtDataManager::ctNumberToQuickCheck(Mesh *mesh)
{
    float maxi = -1;
    float mini = 1E10;
    for(unsigned int i = 0; i < mesh->voxelCount(); i++)
    {
        int ctv = mesh->ct[i];

        // Determine the density
        if(ctv <= 55)
            mesh->density[i] = 0.001 * (1.02 * ctv - 7.65);
        else
            mesh->density[i] = 0.001 * (0.58 * ctv + 467.79);

        if(mesh->density[i] > maxi)
            maxi = mesh->density[i];
        if(mesh->density[i] < mini)
            mini = mesh->density[i];

        // Determine the material
        mesh->zoneId[i] = mesh->density[i]*1.68;
        if(mesh->zoneId[i] > 63)
            mesh->zoneId[i] = 63;
    }

    std::cout << "Max: " << maxi << "   Min: " << mini << std::endl;

    return mesh;
}
