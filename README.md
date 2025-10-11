# DOCTORS
Computed tomogrpahy (CT) and radiogrpahy using ionizing radiations have been applied extensively in nondestructive testing and medical diagnosis. A complete description of the particle fluence distribution and energy transfer is essential for estimating radiation doses and scatter contamination in these imaging systems. Discrete Ordinates Computed TOmography and Radiography Simulator (i.e. DOCTORS) is a code package for producing neutral particle fluence distribution, including both primary and scattered, in a 3D voxelized object by solving the linear Boltzmann transport equation (LBTE) with GPU parallel computing techniques. 

## Overview
DOCTORS is build upon C++/CUDA on Linux platforms. To improve its usability, a graphicial user interface (GUI) is developed to help users generating an input file. Once the input file was generated, users can use that input file as a template future DOCTORS runs. DOCTORS consists of a set of modules, each performing a well-defined processing task. The major modules of the current vesion are listed in the following:
1. Geometry - load object volume and specify mesh grid, physical length and iso center.
2. Cross Section - load multi-group cross section data and specify material type.
3. Quandrature - select the angular discretizaton order (i.e. Sn order)
4. Source - select source type and specify source energy spectrum, position, and profile.
5. Detector - specify source to object distance (SOD) and object to detector distance (ODD).
6. Solver - select solver type and check GPU if using NVidia GPU parallel computing.
7. Input File - generate and edit the input file which can be read by the solver.
8. Launch Solver - simply execute the calculation

## Cautions and Warnings
As we continue to improve the functionality of DOCTORS, there are items which will be implemented in the future versions. The GPU parallel computing was tested using CUDA 12.4 and Nvidia Quadro 3000 and newer. Ig may not work on other platforms.

## Installation



