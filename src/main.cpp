#include <iostream>
#include <string>
#include <cstdio>
#include <ctime>
#include <limits>

#include "macros.h"
#include "dtfrparser.h"
#include "ctdatamanager.h"
#include "mesh.h"
#include "xsection.h"
#include "materialutils.h"
#include "sourceparams.h"
#include "solver.h"
#include "outwriter.h"
#include "quadrature.h"

#include "mcnpwriter.h"    
#include "legendre.h"    

#include "solver_metadata.h"
#include "solver_metadata_extract.h"
#include "solver_gpu.cuh"
#include "outwriter.h"

#include "device_utils.cuh"

double timeElapsed(struct timeval start, struct timeval end)
{
  double time_taken;
  time_taken = (end.tv_sec - start.tv_sec) * 1e6;
  time_taken = (time_taken + (end.tv_usec - start.tv_usec)) * 1e-6;
  return time_taken;
}

int main(int argc, char* argv[])
{

    std::cout << "Start Discrete Ordinates Solver ...." << std::endl;

    // Define variables
    std::string xsection_filename, spectrum_filename, ctdata_filename;
    float s_X, s_Y, s_Z, s_Phi, s_Theta;
    int s_N, s_Type;
    float isoX, isoY, isoZ;

    int xbins, ybins, zbins;
    float xlen, ylen, zlen;

    int mat_type, quad, GPU_flag, ISO_flag, LEG_order, num_angle;

    float geo_magnification, ang_range, start_angle;

    int mcnp_input_gen = 0;

   // Read the input file
    std::string doctors_input = argv[1]; //"doctors.inp";
    std::ifstream inputData;

    inputData.open(doctors_input);
    if (!inputData)
    {
       std::cout << "Doctors input file open error!" << std::endl;
       return 0;
    }

    while (inputData)
    {
       std::cout << "start " << std::endl;
       std::string header;
       std::getline(inputData, header, '=');
       if (header == "xsection")
       {
           inputData >> xsection_filename;
           inputData.ignore(std::numeric_limits<std::streamsize>::max(), '\n'); // Clear buffer
           std::cout << "cross section filename is " << xsection_filename << std::endl;
           continue;
       }

       if (header == "spectrum")
       {
           inputData >> spectrum_filename;
           inputData.ignore(std::numeric_limits<std::streamsize>::max(), '\n'); // Clear buffer
           std::cout << "X-ray spectrum filename is " << spectrum_filename << std::endl;
           continue;
       }

       if (header == "source_x")
       {
           inputData >> s_X;
           inputData.ignore(std::numeric_limits<std::streamsize>::max(), '\n'); // Clear buffer
           std::cout << "source_x = " << s_X << std::endl;
           continue;
       }

       if (header == "source_y")
       {
           inputData >> s_Y;
           inputData.ignore(std::numeric_limits<std::streamsize>::max(), '\n'); // Clear buffer
           std::cout << "source_y = " << s_Y << std::endl;
           continue;
       }

       if (header == "source_z")
       {
           inputData >> s_Z;
           inputData.ignore(std::numeric_limits<std::streamsize>::max(), '\n'); // Clear buffer
           std::cout << "source_z = " << s_Z << std::endl;
           continue;
       }

       if (header == "source_phi")
       {
           inputData >> s_Phi;
           inputData.ignore(std::numeric_limits<std::streamsize>::max(), '\n'); // Clear buffer
           std::cout << "source_phi = " << s_Phi << std::endl;
           continue;
       }

       if (header == "source_theta")
       {
           inputData >> s_Theta;
           inputData.ignore(std::numeric_limits<std::streamsize>::max(), '\n'); // Clear buffer
           std::cout << "source_theta = " << s_Theta << std::endl;
           continue;
       }

       if (header == "source_number")
       {
           inputData >> s_N;
           inputData.ignore(std::numeric_limits<std::streamsize>::max(), '\n'); // Clear buffer
           std::cout << "source_number = " << s_N << std::endl;
           continue;
       }

       if (header == "source_type")
       {
           inputData >> s_Type;
           inputData.ignore(std::numeric_limits<std::streamsize>::max(), '\n'); // Clear buffer
           std::cout << "source_type = " << s_Type << std::endl;
           continue;
       }

       if (header == "iso_x")
       {
           inputData >> isoX;
           inputData.ignore(std::numeric_limits<std::streamsize>::max(), '\n'); // Clear buffer
           std::cout << "iso_x = " << isoX << std::endl;
           continue;
       }

       if (header == "iso_y")
       {
           inputData >> isoY;
           inputData.ignore(std::numeric_limits<std::streamsize>::max(), '\n'); // Clear buffer
           std::cout << "iso_y = " << isoY << std::endl;
           continue;
       }

       if (header == "iso_z")
       {
           inputData >> isoZ;
           inputData.ignore(std::numeric_limits<std::streamsize>::max(), '\n'); // Clear buffer
           std::cout << "iso_z = " << isoZ << std::endl;
           continue;
       }

       if (header == "ctdata")
       {
           inputData >> ctdata_filename;
           inputData.ignore(std::numeric_limits<std::streamsize>::max(), '\n'); // Clear buffer
           std::cout << "ct data filename is " << ctdata_filename << std::endl;
           continue;
       }

       if (header == "xbins")
       {
           inputData >> xbins;
           inputData.ignore(std::numeric_limits<std::streamsize>::max(), '\n'); // Clear buffer
           std::cout << "xbins = " << xbins << std::endl;
           continue;
       }

       if (header == "ybins")
       {
           inputData >> ybins;
           inputData.ignore(std::numeric_limits<std::streamsize>::max(), '\n'); // Clear buffer
           std::cout << "ybins = " << ybins << std::endl;
           continue;
       }

       if (header == "zbins")
       {
           inputData >> zbins;
           inputData.ignore(std::numeric_limits<std::streamsize>::max(), '\n'); // Clear buffer
           std::cout << "zbins = " << zbins << std::endl;
           continue;
       }

       if (header == "xlen")
       {
           inputData >> xlen;
           inputData.ignore(std::numeric_limits<std::streamsize>::max(), '\n'); // Clear buffer
           std::cout << "xlen = " << xlen << std::endl;
           continue;
       }

       if (header == "ylen")
       {
           inputData >> ylen;
           inputData.ignore(std::numeric_limits<std::streamsize>::max(), '\n'); // Clear buffer
           std::cout << "ylen = " << ylen << std::endl;
           continue;
       }

       if (header == "zlen")
       {
           inputData >> zlen;
           inputData.ignore(std::numeric_limits<std::streamsize>::max(), '\n'); // Clear buffer
           std::cout << "zlen = " << zlen << std::endl;
           continue;
       }

       if (header == "geomag")
       {
           inputData >> geo_magnification;
           inputData.ignore(std::numeric_limits<std::streamsize>::max(), '\n'); // Clear buffer
           std::cout << "geometry magnification = " << geo_magnification << std::endl;
           continue;
       }

       if (header == "material_type")
       {
           inputData >> mat_type;
           inputData.ignore(std::numeric_limits<std::streamsize>::max(), '\n'); // Clear buffer
           std::cout << "material_type = " << mat_type << std::endl;
           continue;
       }

       if (header == "quadrature")
       {
           inputData >> quad;
           inputData.ignore(std::numeric_limits<std::streamsize>::max(), '\n'); // Clear buffer
           std::cout << "quadrature = " << quad << std::endl;
           continue;
       }

       if (header == "Legendre_order")
       {
           inputData >> LEG_order;
           inputData.ignore(std::numeric_limits<std::streamsize>::max(), '\n'); // Clear buffer
           std::cout << "Legendre order = " << LEG_order << std::endl;
           continue;
       }
       
       if (header == "GPU_flag")
       {
           inputData >> GPU_flag;
           inputData.ignore(std::numeric_limits<std::streamsize>::max(), '\n'); // Clear buffer
           std::cout << "GPU_flag = " << GPU_flag << std::endl;
           continue;
       }

       if (header == "ISO_flag")
       {
           inputData >> ISO_flag;
           inputData.ignore(std::numeric_limits<std::streamsize>::max(), '\n'); // Clear buffer
           std::cout << "ISO_flag = " << ISO_flag << std::endl;
           continue;
       }

       if (header == "start_angle")
       {
           inputData >> start_angle;
           inputData.ignore(std::numeric_limits<std::streamsize>::max(), '\n'); // Clear buffer
           std::cout << "start angle = " << start_angle << std::endl;
           continue;
       }

       if (header == "angular_range")
       {
           inputData >> ang_range;
           inputData.ignore(std::numeric_limits<std::streamsize>::max(), '\n'); // Clear buffer
           std::cout << "angular range = " << ang_range << std::endl;
           continue;
       }

       if (header == "num_angle")
       {
           inputData >> num_angle;
           inputData.ignore(std::numeric_limits<std::streamsize>::max(), '\n'); // Clear buffer
           std::cout << "number of angle = " << num_angle << std::endl;
           continue;
       }

       if (header == "mcnp_input_gen")
       {
           inputData >> mcnp_input_gen;
           inputData.ignore(std::numeric_limits<std::streamsize>::max(), '\n'); // Clear buffer
           std::cout << "mcnp_input_gen = " << mcnp_input_gen << std::endl;
           continue;
       }
    }

   inputData.close();

    // Construct the data readers
    DtfrParser* p_parser = new DtfrParser;         // DTFR parser
    CtDataManager* p_ctman = new CtDataManager;    // CT data manager
    Mesh* p_meshIn = new Mesh;                     // Input Mesh grid the LBTE build upon
    Mesh* p_meshOut = new Mesh;                    // Rotated mesh
    XSection* p_xs = new XSection;                 // Cross section of selected materials
    //Detector* p_det= new Detector;

    // Read cross section data
    //p_parser->parseFile("H_N_O_Ar_100keV_90gp2.dtfr");
    p_parser->parseFile(xsection_filename);

    // Setup source parameter
    SourceParams* p_source = new SourceParams(p_parser);  // Source parameter
    //std::string spec_name = "ninetygroup.spec";           // X-ray spectra file name
    std::string spec_name = spectrum_filename;            // X-ray spectra file name
    std::ifstream txtfile;
    txtfile.open(spec_name.c_str(), std::ios::in);
    float f;
    std::vector<float> values;  // temp vector to store source intensity
    if(!txtfile.good())
    {
        std::cout << "Cannot open X-ray spectra file!" << std::endl;
        return false;
    }
    else
    {
        while(txtfile >> f)
        {
            values.push_back(f);
        }
    }

    // Check to make sure the correct number of energy groups were loaded
    if((int) values.size() != p_parser->getGammaEnergyGroups())
    {
        std::cout << "X-ray spectra and cross section energy group mismatch! "
                  << " X-ray spectra energy group = " << values.size() << std::endl;
    }

    /*
    int s_Type = 0;             // Source type  0: Point source
                                //              1: Fan beam
                                //              2: Multi fan beam
                                //              3: Cone beam
                                //              4: Multi cone beam
    */

    if(!p_source->update(values, s_X, s_Y, s_Z, s_Phi, s_Theta,
                         s_N, isoX, isoY, isoZ, s_Type))
    {
        std::cout << "Source setup error!" << std::endl;
    }
    else
    {
         std::cout << "Source setup successfully!" << std::endl;
         std::cout << "Source Type: " << p_source->sourceType << std::endl;
         std::cout << "Source position X: " << p_source->sourceX << std::endl;
         std::cout << "Source position Y: " << p_source->sourceY << std::endl;
         std::cout << "Source position Z: " << p_source->sourceZ << std::endl;
         std::cout << "Iso position X: " << p_source->isoX << std::endl;
         std::cout << "Iso position Y: " << p_source->isoY << std::endl;
         std::cout << "Iso position Z: " << p_source->isoZ << std::endl;
         
         // Note the spectra start at highest energy
         for(unsigned int i = 0; i < p_source->spectraIntensity.size(); i++) 
         {
            //p_source->spectraIntensity[i] = p_source->spectraIntensity[i] * 1e3;  // Scale to 1e3 source particles 
            std::cout << "Energy bin #" << i << " = " << p_source->spectraIntensity[i] << std::endl;

         }
            
    }

    p_meshIn = p_ctman->parse16(xbins, ybins, zbins, xlen, ylen, zlen, ctdata_filename);

    // Construct the quadrature
    Quadrature* p_quad = new Quadrature(quad);  // Sn order

    // Construct the solver
    Solver* p_solver = new Solver;

    p_solver->pn = p_parser->getPnOrder();    // DTFR xsection Legendre order

    if(p_solver->pn != LEG_order) 
    {
        std::cout << "Warning: Legendre order in xsection is " << p_solver->pn
            << " but solver Legendre order is " << LEG_order << std::endl;
        p_solver->pn = LEG_order;  // Set user input Legendre order for the solver
    }

    // Loop over the number of projections here to generate rotated volume
    float delta_ang = 0.0;
    if (num_angle == 1)
        delta_ang = 0.0;
    else
        delta_ang = ang_range / (num_angle - 1);
    
    for (int iang = 0; iang < num_angle; iang++)
    {
        float current_ang = start_angle + iang * delta_ang;
        std::cout << "Projection # " << iang << " angle = " << current_ang << std::endl;
        
        p_ctman->rotateCTvol(p_meshIn, p_meshOut, current_ang);  // Rotate CT volume by degree

        std::cout << std::endl << "Generating materials ..." << std::endl;
        
        p_meshOut->material = mat_type;
        switch (p_meshOut->material)
        {
        case MaterialUtils::HOUNSFIELD19:
            std::cout << "Human Phantom is selected " << std::endl;
            p_meshOut = p_ctman->ctNumberToHumanPhantom(p_meshOut);  // Convert CT number to materials
            p_xs->setElements(MaterialUtils::hounsfieldRangePhantom19Elements,
                MaterialUtils::hounsfieldRangePhantom19Weights);
            break;
        case MaterialUtils::WATER:
            std::cout << "Water Phantom is selected " << std::endl;
            p_meshOut = p_ctman->ctNumberToWater(p_meshOut);  // Convert CT number to materials
            p_xs->setElements(MaterialUtils::waterElements,
                MaterialUtils::waterWeights);
            break;

        default:
            std::cout << "Material map " << p_meshOut->material
                << " was not found in the MaterialUtil database!" << std::endl;
            return false;
        }

        // Generate MCNP6 input file
        if (mcnp_input_gen > 0)
        {
            McnpWriter mcnpwriter;
            std::string mcnp_filename = "mcnp6_proj_" + std::to_string(iang) + ".inp";
            mcnpwriter.writeMcnp(mcnp_filename, p_meshOut, p_source, true);
        }
        
        // Construct materials' cross section
        if (!p_xs->allocateMemory(p_parser->getGammaEnergyGroups(), p_solver->pn))
            return false;
        else
            p_xs->addAll(p_parser);

        // Calculate cell areas
        p_meshOut->calcAreas(p_quad, p_parser->getGammaEnergyGroups());   

        // Run ray-tracing to estimate uncollided flux
        std::clock_t startTime = std::clock();
        std::vector<double>* uflux = p_solver->raytraceIsoCPU(p_quad, p_meshOut, p_xs, p_source);  // uncollided flux
        
        // Rename the output file name based on angle#
        std::string newname = "uncollided_flux_iso_" + std::to_string(iang) + ".dat";
        if (std::rename("uncollided_flux_iso.dat", newname.c_str()) != 0)
        {
            std::perror("Error renaming uncollided_flux_iso file!");
        }
        std::cout << "Time to complete raytracing on CPU: "
            << (std::clock() - startTime) / (double)(CLOCKS_PER_SEC / 1000.0)
            << " ms" << std::endl;
        
        // Define the collided flux vector for storage
        std::vector<double> *cflux = new std::vector<double>;

        if (!GPU_flag)  // CPU version
        {
            if (ISO_flag)
            {
                std::cout << "Isotropic solver ..." << std::endl;
                p_solver->gsSolverIsoCPU(p_quad, p_meshOut, p_xs, p_source,
                    p_solver->raytraceIsoCPU(p_quad, p_meshOut, p_xs, p_source));

            }
            else
            {
                //Testing Sn-Pn 
                std::cout << "Anisotropic solver using Spherical Harmonics ..." << std::endl;

                // Get scalar collided flux 
                cflux = p_solver->gsSolverHarmonicCPU(p_quad, p_meshOut, p_xs, p_source, uflux);

                // Rename the output file name 
                //newname = "collided_flux_harmonic_" + std::to_string(iang) + ".dat";
                //if (std::rename("collided_flux_harmonic.dat", newname.c_str()) != 0)
                //{
                //    std::perror("Error renaming collided_flux_harmonic file!");
                //}
            }

        }
        else   // GPU version
        {
            reportGpuData(); 
          // Run LBTE solver on GPU 
            struct solver_metadata metadata;
            SOL_T* cFlux;  //XL :make it double
            SOL_T* cMoments;
            RAY_T* uFlux; 
            // extract metadata for LBTE solver on GPU
            solver_metadata_extract(&metadata, p_quad, p_meshOut, p_xs, p_source, p_solver->pn);
            std::cout << "==Metadata extracted for LBTE solver==" << std::endl;

            // Run ray-tracing to estimate uncollided flux
            std::clock_t startTime = std::clock();
            uFlux = raytraceIsoGPU(metadata);  // uncollided flux
            std::cout << "Time to complete raytracing on GPU: "
                << (std::clock() - startTime) / (double)(CLOCKS_PER_SEC / 1000.0)
                << " ms" << std::endl;
             
            // Transform a double vector to a float*
            //float* uflux_float = new float[uflux->size()];
           // std::transform(uflux->begin(), uflux->end(), uflux_float, [](double val) { return static_cast<float>(val); });


            const Mesh* ptr_mesh = p_meshOut;
            const XSection* ptr_xs = p_xs;

            if (ISO_flag)
            {
                // run LBTE solver on GPU - Isotropic version
                std::cout << "==LBTE ISO solver on GPU==" << std::endl;
                std::clock_t startTime = std::clock();
                //cFlux = gsSolverIsoGPU(metadata, uflux);
                std::cout << "Time to complete Isotropic solver on GPU: "
                    << (std::clock() - startTime) / (double)(CLOCKS_PER_SEC / 1000.0)
                    << " ms" << std::endl;
                //Write collided flux to file. For now this file-writer requires a double vector
                std::vector<double> cFlux_vector;
                for (int i = 0; i < (metadata.num_energies * metadata.num_voxels); i++)
                    cFlux_vector.push_back((double)cFlux[i]);

                const std::vector<double>* ptr_cFlux = &cFlux_vector;

                std::string output_filename = "collided_flux_iso_gpu_" + std::to_string(iang) + ".dat";
                OutWriter::writeScalarFlux(output_filename, *ptr_xs, *ptr_mesh, *ptr_cFlux);
                
            }
            else
            {
                // run LBTE solver on GPU - AnIsotropic version
                std::cout << "==LBTE Harmonic solver on GPU==" << std::endl;
                std::clock_t startTime = std::clock();
                cFlux = gsSolverHarmonicGPU(metadata, uflux);
                //cMoments = gsSolverHarmonicGPU(metadata, uflux);

                std::cout << "Time to complete Anisotropic solver on GPU: "
                    << (std::clock() - startTime) / (double)(CLOCKS_PER_SEC / 1000.0)
                    << " ms" << std::endl;
                //Write collided flux to file. For now this file-writer requires a double vector
                std::vector<double> cFlux_vector;
                for (int i = 0; i < (metadata.num_energies * metadata.num_voxels); i++)
                    cFlux_vector.push_back((double)cFlux[i]);

                const std::vector<double>* ptr_cFlux = &cFlux_vector;

                std::string output_filename = "collided_flux_aniso_gpu_" + std::to_string(iang) + ".dat";

                OutWriter::writeScalarFlux(output_filename, *ptr_xs, *ptr_mesh, *ptr_cFlux);
            
                //Optional: Write collided flux moment to file. For now this file-writer requires a double vector
                //int mCount = (metadata.pn + 1) * (metadata.pn + 1);
                //std::vector<double> cMoments_vector;
                //for (int i = 0; i < (metadata.num_energies * metadata.num_voxels * mCount); i++)
                //    cMoments_vector.push_back((double)cMoments[i]);

                //const std::vector<double>* ptr_cMoments = &cMoments_vector;

                //OutWriter::writeScalarFlux("collided_flux_moment.dat", *ptr_xs, *ptr_mesh, *ptr_cMoments);
                
            }
        } 

    }
    
    
    delete p_parser;
    delete p_source;
    delete p_ctman;
    delete p_meshIn;
    delete p_meshOut;
    delete p_xs;
    delete p_solver;
    delete p_quad;

    return 0;
}



