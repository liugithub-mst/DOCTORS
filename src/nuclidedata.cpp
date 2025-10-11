#include <string>
#include <iostream>

#include "nuclidedata.h"

NuclideData::NuclideData() :
    gScatter1d(),
    gScatter2d(),
    tolxsection(),
    cohxsection(),
    incxsection(),
    phexsection(),
    parxsection()
{

}

NuclideData::~NuclideData()
{

}

bool NuclideData::parse(std::ifstream &binfile, int gGroups)
{
    // Read gamma cross section data
    std::string blankline, title;
    std::getline(binfile, blankline);  // The first line is a blank line

    if(!std::getline(binfile, title))  // Then the title line
    {
        std::cout << "DTFR Parse Data: Could not read the nuclide title!" << std::endl;
        return false;
    }
    else
    {
        std::cout << title << std::endl;
    }

    float tempdata;

    for(int i= 0; i < gGroups; i++)
    {
        binfile >> tempdata;
        cohxsection.push_back(tempdata);
        gScatter1d.push_back(tempdata);
        binfile >> tempdata;
        incxsection.push_back(tempdata);
        gScatter1d[i] += incxsection[i];  //add coherent and incoherent scatter
        binfile >> tempdata;
        phexsection.push_back(tempdata);
        binfile >> tempdata;    // not used
        binfile >> tempdata;    // not used
        binfile >> tempdata;    // total cross section
        tolxsection.push_back(tempdata);

        for(int j = 0; j < gGroups; j++)
        {
            binfile >> tempdata;
            gScatter2d.push_back(tempdata);
        }

    }

    std::string empline;
    if(std::getline(binfile, empline))  // DTFR an empty line between each (material, legendre-coefficient) pair.
    {
        if(empline.empty())
        {
            std::cout << "DTFR parse nuclidedata: Empty line reached" << std::endl;
        }
        else
        {
            std::cout << "DTFR parse nuclidedata : read scatter matrix error!" << std::endl;
            return false;
        }
    }

    return true;
}
