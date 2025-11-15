#ifndef MACROS_H_
#define MACROS_H_

#define SOL_T double //float or double - dataype for GPU anisotropic LBTE solver
#define RAY_T double //float or double - datatype for GPU-based ray-trace

// Approx to 0 and \infinity
#define TINY_T  1.0E-35f
#define HUGE_T  1.0E+35f

// method of computing legendre polynomials
#define LEGENDRE_DIRECT // LEGENDRE_DIRECT defined --> means the recursive method is used to compute legendre polynomials (fast)
                        // otherwise, the spherical harmonic addition theorem based on associated legendre basis is used to do the same.

// Block sizes used by various GPU kernels for LBTE solver and detector-ray trace
#define BLOCK_SIZE_PRECOMP 256
#define BLOCK_SIZE_LBTE 128
#define BLOCK_SIZE 128

#define MIN(a,b) (a<b)?a:b
#define MAX(a,b) (a>b)?a:b

#endif
