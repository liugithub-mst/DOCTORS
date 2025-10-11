#ifndef SOLVER_H_
#define SOLVER_H_

#include <cmath>

#include <vector>
#include <iostream>
#include <omp.h>


class Quadrature;
class Mesh;
class XSection;
class SourceParams;

// VS: Removed all CPU/GPU branching from here.
//     The GPU implementations are kept separate since they utilize certain different data structures for the metadata.
// VS: Consequently, I've removed some (wrapper) functions of class Solver. They are:
//     1) basicRaytraceGPU
//     2) raytraceIso
//     3) raytraceLegendre
// VS: removed following variables of class Solver.
//     1) gpu_accel

// VS: Changed return data-types for all gsSolver<Iso/Legendre> routines to std::vector<double>*

//VS: All member-functions of below class are ONLY CPU-based implementations
class Solver
{
public:
    Solver();
    ~Solver();

    const double m_pi = 3.1415926535897932384626433832795;
    const double m_4pi = 4.0 * m_pi;
    const double m_4pi_inv = 1.0 / m_4pi;

    unsigned int pn;   // Pn order

    // main-routine for ray-trace that other member functions such as raytraceIsoCPU() and raytraceLegendreCPU() wrap around.
    std::vector<double>* basicRaytraceCPU(const Quadrature*, const Mesh* mesh, const XSection* xs, const SourceParams* srcPar);

    // merely invokes basicRaytraceCPU()
    std::vector<double>* raytraceIsoCPU(const Quadrature* quad, const Mesh* mesh,
                            const XSection* xs,const SourceParams* srcPar);

    std::vector<double>* gsSolverIsoCPU(const Quadrature* quad, const Mesh* mesh,
                            const XSection* xs,const SourceParams* srcPar,
                            const std::vector<double>* uflux);

    // merely invokes basicRaytraceCPU() and then assigns direction to uncollided-flux.
    std::vector<double>* raytraceLegendreCPU(const Quadrature* quad, const Mesh* mesh,
                             const XSection* xs, const SourceParams* srcPar);

    std::vector<double>* gsSolverHarmonicCPU(const Quadrature* quad, const Mesh* mesh,
                             const XSection* xs, const SourceParams* srcPar,
                             const std::vector<double>* uflux);   
};



#endif // SOLVER_H_
