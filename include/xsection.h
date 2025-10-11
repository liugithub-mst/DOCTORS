#ifndef XSECTION_H_
#define XSECTION_H_

//This class uses data from the DTFR cross sections

#include <vector>

class DtfrParser;
class Quadrature;

class XSection
{
public:
    XSection();
    ~XSection();

    inline unsigned int groupCount() const { return m_groups; }

    inline float scatXs1d(const int matid, const int g) const
    {
        return m_scat1d[matid*m_groups + g];
    }
    inline float totXs1d(const int matid, const int g) const
    {
        return m_tot1d[matid*m_groups + g];
    }
    float scatXs2d(const int matid, const int gSource, const int gSink, const int n) const;

    bool allocateMemory(const unsigned int materialCount, const unsigned int groupCount, const unsigned int PnCount);
    bool allocateMemory(const unsigned int groupCount, const unsigned int PnCount);

    // Needed if material cross-sections are computed directly in C++ based on MaterialUtils::
    bool addMaterial(const std::vector<int> &z, const std::vector<float> &w, const DtfrParser *p);
    bool addAll(DtfrParser *parser);
    bool setElements(const std::vector<int> &elem, const std::vector<std::vector<float> > &wt);

    // Added utilities to directly set material-wise cross-section data, w/o using MaterialUtils:: and XSection::addMaterial
    // Needed for python interface to this C++ backend
    inline void set_metadata(int num_groups, int num_mats, int pn)
    {
        m_groups = num_groups;
        m_mats = num_mats;
        m_pn = pn;
        m_matsLoaded = num_mats;
    }
    // set total cross-section directly.
    // you can use index-wise assignment rather than .push_back() since ::allocateMemory() is called prior to this
    inline void set_totalXs(float *totalXs)
    {
        int len = m_mats * m_groups;
        for(int i=0; i<len; i++)
            m_tot1d[i] = totalXs[i];
    }
    // set differential cross-section directly.
    // you can use index-wise assignment rather than .push_back() since ::allocateMemory() is called prior to this
    inline void set_diffXs(float *diffXs)
    {
        int len = m_mats * m_groups * m_groups * (m_pn+1);
        for(int i=0; i<len; i++)
            m_scat2d[i] = diffXs[i];
    }
    // we don't need this but just for completeness. contribution of coherent + incoherent scatter to total cross-section.
    inline void set_sumDualXs(float *sumDualXs)
    {    
        int len = m_mats * m_groups;
        for(int i=0; i<len; i++)
            m_scat1d[i] = sumDualXs[i];   
    }

private:
    int m_groups;        //energy groups
    int m_mats;          //materials
    int m_pn;            //order of Legendre expansion from zero
    int m_matsLoaded;    //materials loaded with cross sections
    std::vector<float> m_scat1d;
    std::vector<float> m_scat2d;
    std::vector<float> m_tot1d;

    // Needed if material cross-sections are computed directly in C++ based on MaterialUtils::
    const std::vector<int> *m_elements;                // Elements in mixture
    const std::vector<std::vector<float> > *m_weights; // Weight fractions
    
    std::vector<float> gbounds;                        // FixMe! not used?

};

#endif // XSECTION_H_
