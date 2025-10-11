#include <iostream>

#include "dtfrparser.h"
#include "sourceparams.h"

SourceParams::SourceParams(DtfrParser* p_parser) : sourceType(-1)
{
    // Initialize
    spectraIntensity.resize(p_parser->getGammaEnergyGroups(), 0);
    // Set group energy boundaires
    spectraEnergyLimits = p_parser->getGammaEnergy();

}

SourceParams::SourceParams(std::vector<float> energy_bounds): sourceType(-1)
{
    int num_groups =  energy_bounds.size()-1;
    // Initialize
    spectraIntensity.resize(num_groups, 0);
    // Set group energy boundaires
    spectraEnergyLimits = energy_bounds;
}

SourceParams::~SourceParams()
{

}

bool SourceParams::normalize()
{
    float t = 0; // total intensity

    for(unsigned int i = 0; i < spectraIntensity.size(); i++)
    {
        if(spectraIntensity[i] < 0)
            spectraIntensity[i] = 0;
    }

    for(unsigned int i = 0; i < spectraIntensity.size(); i++)
    {
        t += spectraIntensity[i];
    }

    if(t <= 1e-6)
    {
        std::cout << "The total energy spectra is zero!" << std::endl;
        return false;
    }

    for(unsigned int i = 0; i < spectraIntensity.size(); i++)
    {
        spectraIntensity[i] /= t;
    }
    return true;
}

bool SourceParams::update(std::vector<float> e, double x, double y, double z,
                          double phi, double theta, int n,
                          double d, double w, double h, int type)
{
    if(spectraIntensity.size() != e.size())
    {
        std::cout << "Source spectra size mismatch!" << std::endl;
        return false;
    }
    spectraIntensity = e;
    sourceX = x;
    sourceY = y;
    sourceZ = z;
    sourcePhi = phi;
    sourceTheta = theta;
    sourceN = n;
    isoX = d;
    isoY = w;
    isoZ = h;
    sourceType = type;
    return normalize();
}
