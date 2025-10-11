#ifndef SOURCEPARAMS_H_
#define SOURCEPARAMS_H_

#include <vector>

class DtfrParser;

class SourceParams
{
public:
    SourceParams(DtfrParser* p_parser);
    // VS: added additional constructor independent of parser
    SourceParams(std::vector<float> energy_bounds);
    ~SourceParams();

    int sourceType;

    float sourceX, sourceY, sourceZ;  // Position of the source
    float sourcePhi, sourceTheta;     // Cone beam
    int sourceN;                      // Number of source point around sample

    float isoX, isoY, isoZ;           // Rotation center

    std::vector<float> spectraEnergyLimits;
    std::vector<float> spectraIntensity;

    bool normalize();
    bool update(std::vector<float> e, double x, double y, double z,
                double phi, double theta, int n,
                double d, double w, double h, int type);

};

#endif // SOURCEPARAMS_H_
