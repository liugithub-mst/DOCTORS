# DOCTORS
Computed tomography (CT) and radiography using ionizing radiation have been applied extensively in nondestructive testing and medical diagnosis. A complete description of the particle fluence distribution and energy transfer is essential for estimating radiation doses and scatter contamination in these imaging systems. Discrete Ordinates Computed TOmography and Radiography Simulator (i.e., DOCTORS) is a code package for producing neutral particle fluence distributions—including both primary and scattered components—in a 3D voxelized object by solving the linear Boltzmann transport equation (LBTE) with GPU parallel computing techniques. Its accuracy is close to that of a Monte Carlo simulation, but with much less computation time on a regular desktop computer.

## Overview
DOCTORS is built upon C++/CUDA on both Linux and Windows platforms. To improve its usability, a graphical user interface (GUI) is developed to help users generate an input file. Once the input file is generated, users can use it as a template for future DOCTORS runs. DOCTORS consists of a set of modules, each performing a well-defined processing task. The major modules of the current version are listed below:
- Geometry - load object volume and specify mesh grid, physical length and iso center.
- Cross Section - load multi-group cross section data and specify material type.
- Quandrature - select the angular discretizaton order (i.e. Sn order)
- Source - select source type and specify source energy spectrum, position, and profile.
- Detector - specify source to object distance (SOD) and object to detector distance (ODD).
- Solver - select solver type and check GPU if using NVIDIA GPU parallel computing.
- Input File - generate and edit the input file which can be read by the solver.
- Launch Solver - simply execute the calculation

## Cautions and Warnings
As we continue to improve the functionality of DOCTORS, there are items that will be implemented in future versions. The GPU parallel computing capability was tested using CUDA 12.4 and NVIDIA Quadro 3000 and newer GPUs. It may not work on other platforms.

## Installation
Users can build the DOCTORS executable file using CMake and the provided CMakeLists file. Users can also download the executable file directly and place it in any folder under their home directory.

## Required Files from User
In order to run a simulation using DOCTORS, several data files are required. These data files can be placed in any folder; however, it is recommended to put all the data files in the same directory where DOCTORS is located for good file organization.
The required fiels are:
1. CT volume data file. This is a 16-bit unsigned binary file, which contains CT number in every voxel of a 3D object. This data file can be easily generated from a series of DICOM images using ImageJ. It can also be generated numerically by a data processing tool, such as Matlab or Python.
2. Multi-group cross section data file. This data file is usually generated from ENDF/B data library using NJOY. NJOY2016 is an open source software, and is recommedn to teh user for the cross section data generation.
3. Source energy spectrum. This file contains normalzied enegy spectrum of the source particle.

In addition to these required files, the user needs to provide the material type of the object, which describes the chemical element composition of each voxel. Two prebuild material types are available, namely, water and human abdomen body part. Users can edit these examples for their applications.

## How to Run DOCTORS
To run DOCTORS, simply go to the directory where DOCTORS executalbe and the Python GUI is located, and type the following command:
~$python3 doctors_main.py

It is also possilbe to run DOCTORS without using the GUI. Simply go to the directory where both the DOCTORS executable and the input file are located, and type the following command: 
~$./DOCTORS your_input_file_name.txt

## How to Uninstall DOCTORS
To uninstall DOCTORS, simply delete the DOCTORS executalbe or the whole folder conataining DOCTORS and the Python GUI files.

## Example Results
As a simple demonstration of the accuracy of DOCTORS, we show below the effective dose using MCNP6 and DOCTORS of an abdoman body part. The relative differnece between DOCTORS and MCNP6 is less than 5%, and DOCTORS is over 10 times faster than MCNP6. Detailed results can be found in this paper,
- Edward T. Norris and Xin Liu, " Phtoton flunce and dose estimation in computed tomography using a discrete ordinates Boltzmann solver", Scientific Reports, 10.1 (2020): 11609.

## Authors and Contact
DOCTORS are developed by Edward Norris and Xin Liu, please contact Xin Liu (liux.xin@gmail.com) for any questions.

## License
DOCTORS is distributed under the terms of the MIT license. All new contributions must be made under this license. See LICENSE in this directory for the terms of the license. 

Please cite our work by referencing this github page and citing our paper:
- Edward T. Norris and Xin Liu, " Phtoton flunce and dose estimation in computed tomography using a discrete ordinates Boltzmann solver", Scientific Reports, 10.1 (2020): 11609.



