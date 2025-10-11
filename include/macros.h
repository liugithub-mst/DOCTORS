#ifndef MACROS_H_
#define MACROS_H_

// VS: MACROS: pre-defined flags, constants, compact functions used across various .cpp/.cu files
// VS: Later on move multiple MACROS and predefined constants to this header file

#define SOL_T double //float or double - dataype for GPU anisotropic LBTE solver
#define RAY_T double //float or double - datatype for GPU-based ray-trace

// VS: Approx to 0 and \infinity
#define TINY_T  1.0E-35f
#define HUGE_T  1.0E+35f
// VS: constant factor within forward-projector that affects sampling of rays within detector-pixel
#define SIGMA_MULT 0.5 //multiplier of pixel-size (along diagonal) for detector ray-trace
                       //limit SIGMA_MULT-->0 means nearest-neighbor sampling with respect to detector-pixel center

// vs: method of computing legendre polynomials
#define LEGENDRE_DIRECT // LEGENDRE_DIRECT defined --> means the recursive method is used to compute legendre polynomials (fast)
                        // otherwise, the spherical harmonic addition theorem based on associated legendre basis is used to do the same.

// VS: Block sizes used by various GPU kernels for LBTE solver and detector-ray trace
#define BLOCK_SIZE_PRECOMP 256
#define BLOCK_SIZE_LBTE 128
#define BLOCK_SIZE 128

// VS: The below 2 definitions are for multi-threaded CPU implementations
// VS: if NTEAMS is defined, it limits the max number of threads to NTEAMS*NTHREADS_PER_TEAM
//     if NTEAMS is undefined, maximum threads available are used.
/*
#define NTEAMS 8
*/
// VS: if NTHREADS_PER_TEAM > 1 we have nested parallelism for the CPU-based LBTE solver
//     threads are divided into NTEAMS different teams/groups (outer parallelism) each consisting of "NTHREADS_PER_TEAM" sub-threads (inner parallelism)
#define NTHREADS_PER_TEAM 1



#define MIN(a,b) (a<b)?a:b
#define MAX(a,b) (a>b)?a:b

#endif
