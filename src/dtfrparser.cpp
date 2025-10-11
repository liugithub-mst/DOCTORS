#include <ctime>
#include <cstdio>

#include <string>
#include <algorithm>
#include <iostream>

#include "dtfrparser.h"

DtfrParser::DtfrParser()
{

}

DtfrParser::~DtfrParser()
{
    closeFile();

    for(unsigned int i = 0; i < data.size(); i++)
        if(data[i] != nullptr)
            delete data[i];

}

bool DtfrParser::openFile(std::string filename)
{
    m_filename = filename;
    binfile.open(filename.c_str(), std::ios::in);

    if(!binfile.good())
    {
        std::cout << "Cannot open DTFR file!" << std::endl;
        return false;
    }
    else
        return true;
}

void DtfrParser::closeFile()
{
    if(binfile.is_open())
        binfile.close();
}

bool DtfrParser::parseFile(std::string filename)
{
    if(openFile(filename)) 
    {
        if(parseHeader()) 
        {
            if(parseData()) 
            {
                std::cout << "DTFR parser: data successfully parsed!" << std::endl;
                return true;
            }
            else
            {
                std::cout << "DTFR parser: parseData failed!" << std::endl;
                return false;
            }
        }
        else
        {
            std::cout << "DTFR parser: parseHeader failed!" << std::endl;
                return false;
        }
    }
    else
    {
        std::cout << "DTFR parser: open file failed!" << std::endl;
        return false;
    }

}

bool DtfrParser::parseHeader()
{
    if(!binfile.good())
    {
        std::cout << "DTFR parseHeader: File was not found!" << std::endl;
        return false;
    }

    std::string title;
    if(!std::getline(binfile, title))  // First line of header is title
    {
        std::cout << "DTFR parseHeader: Could not read the title!" << std::endl;
        return false;
    }

    std::cout << "DTFR title: " + title << std::endl;

    if(!(binfile >> gGroups))
    {
        std::cout << "DTFR parseHeader: Could not read the number of groups!" <<std::endl;
        return false;
    }

    std::cout << "DTFR number of groups = " << gGroups << std::endl;

    float tempdata;
    for(int i = 0; i < gGroups+1; i++)
    {
        if(!(binfile >> tempdata))
        {
            std::cout << "DTFR parseHeader: Could not read the energy bounds!" <<std::endl;
            return false;
        }
        gBounds.push_back(tempdata);  // Store energy group boundaries
    }
    std::reverse(gBounds.begin(), gBounds.end()); // Make energy boundary high to low

    for(int i = 0; i <gGroups+1; i++)
    {
        std::cout << "DTFR energy bins = " <<
        std::to_string(gBounds[i]) << std::endl;
    }

    if(!(binfile >> nPn))
    {
        std::cout << "DTFR parseHeader: Could not read the Legendre order!" <<std::endl;
        return false;
    }
    std::cout << "DTFR Legendre order = " << nPn << std::endl;

    if(!(binfile >> nNuclides))
    {
        std::cout << "DTFR parseHeader: Could not read number of nuclides!" <<std::endl;
        return false;
    }
    std::cout << "DTFR number of nuclides = " << nNuclides << std::endl;

    int tempid;
    for(int i = 0; i < nNuclides; i++)
    {
        if(!(binfile >> tempid))
        {
            std::cout << "DTFR parseHeader: Could not read nuclide zAid!" <<std::endl;
            return false;
        }
        zAids.push_back(tempid);
    }

    for(int i = 0; i < nNuclides; i++)
    {
        std::cout << "DTFR nuclides = " <<
        std::to_string(zAids[i]) << std::endl;
    }

    std::string empline;
    if(std::getline(binfile, empline)) // DTFR an empty line between header and data
    {
        if(empline.empty())
        {
            std::cout << "DTFR parseHeader: Empty line reached!" <<std::endl;
        }
        else
        {
            std::cout << "DTFR parseHeader: read header error!" <<std::endl;
            return false;
        }
    }

    return true;
}

bool DtfrParser::parseData()
{
    if(!binfile.good())
    {
        std::cout << "DTFR parseData: File was not good!" << std::endl;
        return false;
    }

    std::clock_t start;
    start = std::clock();

    for(int i = 0; i < nNuclides*(nPn+1); i++)
    {
        NuclideData* p_nextData = new NuclideData; // each Pn expansion is a new mat.
        p_nextData->parse(binfile, gGroups);
        data.push_back(p_nextData);
    }

    std::cout << "Time: " << (std::clock() - start)/(double)(CLOCKS_PER_SEC/1000)
                          << " ms" << std::endl;
    std::cout << "DTFR: finished parsing" << std::endl;

    return true;

}

int DtfrParser::getIndexbyZaid(int zaid) const
{
    int indx = -1;
    // Search for atomic number Z in the DTFR file
    for(int i = 0; i < nNuclides; i++)
    {
        if(zAids[i] == zaid)
        {
            if(indx == -1)
            {
                indx = i;
            }
            else
            {
                std::cout << "DtfrParser::getIndexByZaid(): multiple indices!"
                          << std::endl;
            }
        }
    }
    if(indx == -1)
    {
        std::cout << "DtfrParser::getIndexByZaid(): no index!" << std::endl;
    }
    return indx;
}

NuclideData* DtfrParser::getData(unsigned int indx, unsigned int nl) const
{
    if(indx*(nPn+1)+nl < data.size())
        return data[indx*(nPn+1)+nl]; //each Legendre expansion is a new material
    else
    {
        std::cout << "DtfrParser::getData(): Could not get nuclide number " +
                    std::to_string(indx) << std::endl;
        return nullptr;
    }
}

void DtfrParser::debugGammaEnergies()
{
    std::string log("");
    for(unsigned int i = 0; i < getGammaEnergy().size(); i++)
    {
        log += std::to_string(getGammaEnergy()[i]) + "\t";
    }
    std::cout << "Gamma Energy (" + std::to_string(getGammaEnergy().size()-1) +
                "): " + log << std::endl;
}
