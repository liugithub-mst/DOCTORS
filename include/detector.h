#ifndef DETECTOR_H_
#define DETECTOR_H_

#include <vector>
#include <iostream>
#include <tuple>

#include "mesh.h"
#include "sourceparams.h"
#include "xsection.h"
#include "quadrature.h"

//VS: macros.h included in .cpp

class Detector
{
public:
    Detector();
    ~Detector();

    float geoMag;  // Geometry magnification ratio ~SDD/SOD
    float scaFac;  // Scaling factor for the collided flux

    std::vector<float> xDetNodes;  // Coordinates between detector elements in x
    std::vector<float> zDetNodes;  // Coordinates between detector elements in z

    std::vector<float> detResfun;  // Detector response function

    std::vector<double> primaryImg;    // Detector primary image
    std::vector<double> scatterImg;    // Detector scatter image

    inline float getDetcenX () const { return m_detcenX; }
    inline float getDetcenY () const { return m_detcenY; }
    inline float getDetcenZ () const { return m_detcenZ; }
    inline unsigned int getDetxelem() const { return m_xElemDet; }
    inline unsigned int getDetzelem() const { return m_zElemDet; }
    inline unsigned int getDetxnode() const { return m_xNodeDet; }
    inline unsigned int getDetznode() const { return m_zNodeDet; }
    inline float getDetlowerX() const { return m_lowerX; }
    inline float getDetlowerZ() const { return m_lowerZ; }
    inline float getDetupperX() const { return m_upperX; }
    inline float getDetupperZ() const { return m_upperZ; }
    inline float getDetpxarea() const { return m_detDx * m_detDz; }

    // VS: Added
    inline float getDetDx() const { return m_detDx; }
    inline float getDetDz() const { return m_detDz; }

    // Setup detector elements and nodes
    void setDetector(const Mesh* ctmesh, const SourceParams* source);

    // Set detector response function
    void setDetResponse(std::string filename);

    // Get primary image from uncollided flux
    void getPrimaryImg(const std::vector<double>* uflux,
                         const Mesh* ctmesh,
                         const SourceParams* source);

    // If ray hit the detector, return 1D detector cell index
    std::tuple<int, float, float> hitDet(const float x, const float y, const float z,
                const float mu, const float eta, const float ksi);

    // VS: added below routines for forward-projector that samples multiple rays within detector-pixel
    float get_weight(float dist_from_pixel_center);
    float get_srcDetDist(int idx_det, const SourceParams* source);

    // Get scatter image from uncollided and collided flux
    float **getScatterImgIso(const std::vector<double>* uflux,
                             const std::vector<double>* cflux,
                             const Mesh* ctmesh,
                             const XSection* xs);

   // Get scatter image from uncollided and collided flux - anisotropic scattering
   float **getScatterImgAniso(const std::vector<double>* uflux,
                              const std::vector<double>* cflux,
                              const Mesh* ctmesh,
                              const XSection* xs,
                              const SourceParams* source,
                              const int Pn);

private:
    float m_detcenX;  // Physical detector center x position
    float m_detcenY;  // Physical detector center y position
    float m_detcenZ;  // Physical detector center z position

    unsigned int m_xElemDet;       // Number of pixels in x dimension
    unsigned int m_zElemDet;       // Number of pixels in z dimension

    unsigned int m_xNodeDet;       // Number of nodes in x dimension ~xElemDet+1
    unsigned int m_zNodeDet;       // Number of nodes in z dimension ~zElemDet+1

    float m_detDx;                 // Detector x pixel size
    float m_detDz;                 // Detector z pixel size

    float m_lowerX;                // Physical detector x boundary
    float m_lowerZ;                // Physical detector z boundary
    float m_upperX;                // Physical detector x boundary
    float m_upperZ;                // Physical detector z boundary

};

// TODO: this is a placeholder, replace with a common util
template <class T> float *vector_to_float_array(std::vector<T> v)
{
    float *array = new float[v.size()];
    for(unsigned int i=0; i<v.size(); i++)
        array[i]=v[i];
    return array;
}


#endif // DETECTOR_H_
