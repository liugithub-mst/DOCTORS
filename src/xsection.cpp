#define _USE_MATH_DEFINES
#include <cmath>

#include <fstream>
#include <iostream>

#include "xsection.h"
#include "dtfrparser.h"
#include "materialutils.h"
#include "quadrature.h"

XSection::XSection() : m_groups(0), m_matsLoaded(0), m_elements(nullptr), m_weights(nullptr)
{

}

XSection::~XSection()
{

}

float XSection::scatXs2d(const int matid, const int gSrcIndx, const int gSinkIndx, const int n) const
{
    if(matid > m_matsLoaded || gSrcIndx >= m_groups || gSinkIndx >= m_groups || n > m_pn)
    {
        std::cout << "scatXs2d: Illegal index value" << std::endl;
        std::cout << "Violated matid > m_matsLoaded || gSrcIndx >= m_groups || gSinkIndx >= m_groups || n > m_pn" << std::endl;
        std::cout <<  matid << "> " << m_matsLoaded << " || " <<
        gSrcIndx << " >= " << m_groups << " || " << gSinkIndx << " >= " << m_groups
                 << " || " << n << " > " << m_pn;
        return -1.0;
    }
    const int pnCount = m_pn+1;
    return m_scat2d[matid*m_groups*m_groups*pnCount + gSrcIndx*m_groups*pnCount + gSinkIndx*pnCount + n];
}

bool XSection::allocateMemory(const unsigned int materialCount, const unsigned int groupCount, const unsigned int pn)
{
    m_tot1d.clear();
    m_scat1d.clear();
    m_scat2d.clear();

    m_mats = materialCount;
    m_groups = groupCount;
    m_pn = pn;
    unsigned long floats1d = materialCount * groupCount;
    unsigned long floats2d = materialCount * groupCount * groupCount * (m_pn+1);

    try
    {
        m_tot1d.resize(floats1d);
    }
    catch(std::bad_alloc &bad)
    {
        std::cout << "bad_alloc caught during XS initialization of the 1d total xs data, requested "
                  << floats1d * sizeof(float) << " bytes. Reported error: " << bad.what();

        return false;
    }

    try
    {
        m_scat1d.resize(floats1d);
    }
    catch(std::bad_alloc &bad)
    {
        std::cout << "bad_alloc caught during XS initialization of the 1d scatter xs data, requested "
                  << floats1d * sizeof(float) <<  " bytes. Reported error: " << bad.what();
        return false;
    }

    try
    {
        m_scat2d.resize(floats2d);
    }
    catch(std::bad_alloc &bad)
    {
        std::cout << "bad_alloc caught during XS initialization of the 2d scatter xs data, requested "
                  << floats2d * sizeof(float) << " bytes. Reported error: " << bad.what();

        return false;
    }

    return true;
}

bool XSection::allocateMemory(const unsigned int groupCount, const unsigned int pn)
{
    if(m_elements == nullptr || m_weights == nullptr)
        return false;

    return allocateMemory(m_weights->size()+1, groupCount, pn);
}

bool XSection::addMaterial(const std::vector<int> &z, const std::vector<float> &w, const DtfrParser *p)
{
    if(m_matsLoaded >= m_mats)
    {
        std::cout << "There was an internal error, the xs table is already full but another material was added to it. "
                  << "The table has room for " << m_mats << " materials." << std::endl;

        return false;
    }

    std::vector<float> atom_frac = MaterialUtils::weightFracToAtomFrac(z, w);

    // For each z
    for(unsigned int i = 0; i < z.size(); i++)
    {
        float afrac = atom_frac[i];
       // unsigned int elemIndx = z[i] - 1;

        // If there is a natural isotope in the library, use it
        //int naturalZAID = z[i]*1000;
        int naturalZAID = z[i]*100;  //In DTFR, atomic number is multiplied by 100

        int naturalIndx = p->getIndexbyZaid(naturalZAID);

        if(naturalIndx >= 0)
        {
            // Handle the scatter 1d xs
            NuclideData *d = p->getData(naturalIndx, 0); //1D xs data are stored in "0"th order Legendre expansion

            // Loop over (source) energy groups 
            for(int ei = 0; ei < m_groups; ei++)
            {
                m_scat1d[m_matsLoaded*m_groups + ei] += (d->getGammaScatter1D()).at(ei)*afrac;   //Get the total 1d scatter xs
                
                m_tot1d[m_matsLoaded*m_groups + ei] += (d->getGammaXs()).at(ei)*afrac;           //Get the total xs
            }


            // Handle the scatter 2d xs
            const int pnCount = m_pn+1;
            for(int n = 0; n < pnCount; n++)
            {
               NuclideData *d = p->getData(naturalIndx, n); //2D xs data are stored in the order of Legendre expansion
               
                for(int iesrc = 0; iesrc < m_groups; iesrc++)
                {
                    for(int iesnk = 0; iesnk < m_groups; iesnk++)
                    {
                        if(iesrc == iesnk)  //in-group scatter
                            m_scat2d[m_matsLoaded*m_groups*m_groups*pnCount + iesrc*m_groups*pnCount + iesnk*pnCount + n] += (d->getGammaScatter2D()).at(iesrc*m_groups)*afrac;
                        else if(iesrc < iesnk) //down scatter
                        {
                            m_scat2d[m_matsLoaded*m_groups*m_groups*pnCount + iesrc*m_groups*pnCount + iesnk*pnCount + n] += (d->getGammaScatter2D()).at(iesnk*m_groups+iesnk-iesrc)*afrac;
                        }
                        else
                            continue;  //assume no upscatter
                    }
                }

            }

        }

    } // for each z[i]

    // This material is now loaded
    m_matsLoaded++;

    std::cout << "Mat ID = " << m_matsLoaded << std::endl;

    return true;
}

bool XSection::addAll(DtfrParser *parser)
{
    if(m_elements == nullptr || m_weights == nullptr)
        return false;

    bool allPassed = true;

    // Add the materials to the xs library
    for(unsigned int i = 0; i < m_weights->size(); i++)
        if(allPassed)
            allPassed &= addMaterial(*m_elements, (*m_weights)[i], parser);

    // The last material is empty and should never be used
    if(allPassed)
        allPassed &= addMaterial(std::vector<int>{}, std::vector<float>{}, parser);

    return allPassed;
}

bool XSection::setElements(const std::vector<int> &elem, const std::vector<std::vector<float> > &wt)
{
   for(unsigned int i = 0; i < wt.size(); i++)
    {
        std::cout << "elem.size = " << elem.size() << std::endl;
        std::cout << "wt[i].size = " << wt[i].size() << std::endl;
        if(elem.size() != wt[i].size())
        {
            // Number of elements should match the number of weight fractions in material [i]
            std::cout << "setElements failed due to a size mismatch" << std::endl;

            return false;
        }
    }

    m_elements = &elem;
    m_weights = &wt;

    return true;
}
