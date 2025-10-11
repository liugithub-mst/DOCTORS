#define _USE_MATH_DEFINES
#include <cmath>
#include <ctime>

#include <iostream>
#include <algorithm>

#include "solver.h"
#include "quadrature.h"
#include "mesh.h"
#include "xsection.h"
#include "legendre.h"
#include "sourceparams.h"
#include "outwriter.h"
#include "sweep.h"  
#include "macros.h" 

Solver::Solver() : pn(0)
{

}

Solver::~Solver()
{

}

/************************************************************************************************
 * ========================================= CPU Code ========================================= *
 ************************************************************************************************/

 std::vector<double>* Solver::basicRaytraceCPU(const Quadrature*, const Mesh* mesh,
                                               const XSection* xs,
                                               const SourceParams* srcPar)

{
    unsigned int groups = xs->groupCount();

    const unsigned short DIRECTION_X = 1;
    const unsigned short DIRECTION_Y = 2;
    const unsigned short DIRECTION_Z = 3;

    std::vector<double> sx;  // Actual X position for each point source
    std::vector<double> sy;  // Actual Y position for each point source
    std::vector<double> sz;  // Actual Z position for each point source

    int sourceCt = 1;

    double xt;
    double yt;
    double dt;

    double focusX = srcPar->isoX; // Iso center X
    double focusY = srcPar->isoY; // Iso center y
    double focusZ = srcPar->isoZ; // Iso center z

    // Build a list of point sources
    switch(srcPar->sourceType)
    {
    case 0:  // Isotropic point source
    case 1:  // Fan
    case 3:  // Cone
        sx.push_back(static_cast<double>(srcPar->sourceX));
        sy.push_back(static_cast<double>(srcPar->sourceY));
        sz.push_back(static_cast<double>(srcPar->sourceZ));
        break;

    case 2:  // Multifan
    case 4:  // Multicone
        sx.push_back(static_cast<double>(srcPar->sourceX));
        sy.push_back(static_cast<double>(srcPar->sourceY));
        sz.push_back(static_cast<double>(srcPar->sourceZ));

        dt = 2*m_pi / srcPar->sourceN;  // Delta theta - angle interval

        for(int i = 0; i < srcPar->sourceN-1; i++)
        {
            xt = sx[sx.size()-1];
            yt = sy[sy.size()-1];
            sx.push_back((xt-focusX)*cos(dt) - (yt-focusY)*sin(dt) + focusX);
            sy.push_back((yt-focusY)*cos(dt) + (xt-focusX)*sin(dt) + focusY);
            sz.push_back(srcPar->sourceZ);
        }

        sourceCt = srcPar->sourceN;
        break;

    default:
        std::cout << "Error: Source type " << srcPar->sourceType
                  << " not understood" << std::endl;
        return nullptr; 
    }

    double phi;
    double theta;

    // Default of phi and theta are degrees. Convert degrees to radians
    phi = srcPar->sourcePhi * m_pi/180.0;
    theta = srcPar->sourceTheta * m_pi/180.0;

    // Allocate solution memory - Uncolided flux in each voxel
    std::vector<double> *uflux = new std::vector<double>;
    uflux->resize(groups * mesh->voxelCount(), static_cast<double>(0.0));

    unsigned int ejmp = mesh->voxelCount();
    unsigned int xjmp = mesh->xjmp();
    unsigned int yjmp = mesh->yjmp();

    std::cout << "Running raytracer ... " << std::endl;

    double tiny = 1.0E-35f;
    double huge = 1.0E35f;
    std::vector<double> meanFreePaths;
    meanFreePaths.resize(xs->groupCount());

    // Set the source energy distribution
    std::vector<double> srcStrength(groups, 0.0);
    for(unsigned int i = 0; i < groups; i++)
        srcStrength[i] = srcPar->spectraIntensity[i];

    int oldfraction = 10; // For raytracer progress purpose

    for(int is = 0; is < sourceCt; is++)
    {
        // Find the indices of the source index
        unsigned int srcIndxX = 0;
        unsigned int srcIndxY = 0;
        unsigned int srcIndxZ = 0;

        // If the source is outside the mesh then it can't have meaningful indices
        if(sx[is] < mesh->xNodes[0] || sx[is] > mesh->xNodes[mesh->xNodes.size()-1] ||
           sy[is] < mesh->yNodes[0] || sy[is] > mesh->yNodes[mesh->yNodes.size()-1] ||
           sz[is] < mesh->zNodes[0] || sz[is] > mesh->zNodes[mesh->zNodes.size()-1])
        {
            srcIndxX = srcIndxY = srcIndxZ = -1; 
        }
        else
        {
            while(mesh->xNodes[srcIndxX+1] < sx[is])
                srcIndxX++;

            while(mesh->yNodes[srcIndxY+1] < sy[is])
                srcIndxY++;

            while(mesh->zNodes[srcIndxZ+1] < sz[is])
                srcIndxZ++;
        }

        //unsigned int totalMissedVoxels = 0;

        unsigned int highestEnergy = 0;
        while(srcStrength[highestEnergy] <= 0 && highestEnergy < groups)
            highestEnergy++;

        if(highestEnergy >= groups)
            return nullptr; 

        double beamVectorX = sx[is] - focusX;
        double beamVectorY = sy[is] - focusY;
        double beamVectorZ = sz[is] - focusZ;
        double beamCenterMag = sqrt(beamVectorX*beamVectorX +
                                    beamVectorY*beamVectorY +
                                    beamVectorZ*beamVectorZ);
        beamVectorX /= beamCenterMag;
        beamVectorY /= beamCenterMag;
        beamVectorZ /= beamCenterMag;

        for(unsigned int zIndxStart = 0; zIndxStart < mesh->zElemCt; zIndxStart++)
            for(unsigned int yIndxStart = 0; yIndxStart < mesh->yElemCt; yIndxStart++)
                for(unsigned int xIndxStart = 0; xIndxStart < mesh->xElemCt; xIndxStart++)  // For every voxel
                {
                    double acceptance = 1.0;

                    std::vector<double> travlength;
                    std::vector<double> travmu;
                    std::vector<double> travoptic;

                    double x = mesh->xNodes[xIndxStart] + mesh->dx[xIndxStart]/2;  // X position of current voxel
                    double y = mesh->yNodes[yIndxStart] + mesh->dy[yIndxStart]/2;
                    double z = mesh->zNodes[zIndxStart] + mesh->dz[zIndxStart]/2;

                    if(xIndxStart == srcIndxX && yIndxStart == srcIndxY && zIndxStart == srcIndxZ)  // End condition if it is the source cell
                    {
                        double srcToCellDist = sqrt((x-sx[is])*(x-sx[is]) + (y-sy[is])*(y-sy[is]) + (z-sz[is])*(z-sz[is]));
                        unsigned int zid = mesh->zoneId[xIndxStart*xjmp + yIndxStart*yjmp + zIndxStart];   // 1D index of zone map
                        double xsval;
                        for(unsigned int ie = highestEnergy; ie < xs->groupCount(); ie++)
                        {
                            xsval = xs->totXs1d(zid, ie) *
                                    mesh->atomDensity[xIndxStart*xjmp + yIndxStart*yjmp + zIndxStart];

                            (*uflux)[ie*ejmp + xIndxStart*xjmp + yIndxStart*yjmp + zIndxStart] +=
                                                srcStrength[ie] * exp(-xsval*srcToCellDist) /
                                            (4 * m_pi * srcToCellDist * srcToCellDist * sourceCt);
                        }
                        
                        continue;
                    }

                    // Start raytracing through the geometry
                    unsigned int xIndx = xIndxStart;
                    unsigned int yIndx = yIndxStart;
                    unsigned int zIndx = zIndxStart;

                    double srcToCellX = sx[is] - x;
                    double srcToCellY = sy[is] - y;
                    double srcToCellZ = sz[is] - z;

                    double srcToCellDist = sqrt(srcToCellX*srcToCellX + srcToCellY*srcToCellY + srcToCellZ*srcToCellZ);
                    double srcToPtDist;

                    double xcos = srcToCellX/srcToCellDist;  // Fraction of direction biased in x-direction, unitless
                    double ycos = srcToCellY/srcToCellDist;
                    double zcos = srcToCellZ/srcToCellDist;

                    // Determine whether or not the ray will be outright rejected
                    double zetaTheta;
                    double zetaPhi;

                    switch(srcPar->sourceType)
                    {
                    case 0: // Isotropic (do nothing)
                        break;

                    case 1: // Fan beam
                    case 2: // Multifan
                        zetaTheta = acos(beamVectorZ) - acos(zcos);
                        if(std::abs(zetaTheta) > theta/2.0)
                        {
                            continue;
                        }

                        // Test phi rejection
                        zetaPhi = acos((beamVectorX*xcos + beamVectorY*ycos)/
                                       (sqrt(beamVectorX*beamVectorX + beamVectorY*beamVectorY) *
                                        sqrt(xcos*xcos + ycos*ycos)));
                        if(std::abs(zetaPhi) > phi/2.0)
                            continue;

                        acceptance = 1.0;
                        break;

                    case 3: // Cone beam
                    case 4: // Multicone
                        zetaTheta = acos(beamVectorX*xcos + beamVectorY*ycos + beamVectorZ*zcos);
                        if(std::abs(zetaTheta) > theta)
                        {
                            continue;
                        }

                        acceptance = 1.0;
                        break;
                    default:
                        std::cout << "This should never happen. basicRaytraceCPU got an illegal source type: "
                                  << srcPar->sourceType << std::endl;
                    } // End of switch

                    // Index of the boundary the particle is headed toward.
                    unsigned int xBoundIndx = (xcos >= 0 ? xIndx+1 : xIndx);
                    unsigned int yBoundIndx = (ycos >= 0 ? yIndx+1 : yIndx);
                    unsigned int zBoundIndx = (zcos >= 0 ? zIndx+1 : zIndx);

                    // Clear the MFP array to zeros
                    for(unsigned int i = 0; i < xs->groupCount(); i++)
                        meanFreePaths[i] = 0.0f;

                    bool exhaustedRay = false;
                    while(!exhaustedRay)
                    {

                        // recomputing the direction cosines ensures that roundoff error doesn't cause the ray to miss the source cell
                        srcToCellX = sx[is] - x;
                        srcToCellY = sy[is] - y;
                        srcToCellZ = sz[is] - z;
                        srcToPtDist = sqrt(srcToCellX*srcToCellX + srcToCellY*srcToCellY + srcToCellZ*srcToCellZ);
                        xcos = srcToCellX/srcToPtDist;  // Fraction of direction biased in x-direction, unitless
                        ycos = srcToCellY/srcToPtDist;
                        zcos = srcToCellZ/srcToPtDist;

                        // Determine the distance to cell boundaries
                        double tx = (fabs(xcos) < tiny ? huge : (mesh->xNodes[xBoundIndx] - x)/xcos);  // Distance traveled [cm] when next cell is
                        double ty = (fabs(ycos) < tiny ? huge : (mesh->yNodes[yBoundIndx] - y)/ycos);  // entered traveling in x direction
                        double tz = (fabs(zcos) < tiny ? huge : (mesh->zNodes[zBoundIndx] - z)/zcos);

                        // Determine the shortest distance traveled [cm] before _any_ surface is crossed
                        double tmin;
                        unsigned short dirHitFirst;

                        if(tx < ty && tx < tz)
                        {
                            tmin = tx;
                            dirHitFirst = DIRECTION_X;
                        }
                        else if(ty < tz)
                        {
                            tmin = ty;
                            dirHitFirst = DIRECTION_Y;
                        }
                        else
                        {
                            tmin = tz;
                            dirHitFirst = DIRECTION_Z;
                        }

                        if(tmin < -3E-8)  // Include a little padding for hitting an edge
                            std::cout << "basicRaytraceCPU: Reversed space!" << std::endl;

                        // Update mfp array
                        unsigned long int zid = mesh->zoneId[xIndx*xjmp + yIndx*yjmp + zIndx];
                      
                        for(unsigned int ie = highestEnergy; ie < xs->groupCount(); ie++)
                        {
                            meanFreePaths[ie] += tmin * xs->totXs1d(zid, ie) * mesh->atomDensity[xIndx*xjmp + yIndx*yjmp + zIndx];
                        }

                        travlength.push_back(tmin);
                        travmu.push_back(xs->totXs1d(zid, 0) * mesh->atomDensity[xIndx*xjmp + yIndx*yjmp + zIndx]);
                        travoptic.push_back(tmin * xs->totXs1d(zid, 0) * mesh->atomDensity[xIndx*xjmp + yIndx*yjmp + zIndx]);

                        // Update cell indices and positions
                        if(dirHitFirst == DIRECTION_X) // x direction
                        {
                            x = mesh->xNodes[xBoundIndx];
                            y += tmin*ycos;
                            z += tmin*zcos;

                            if(xcos >= 0)
                            {
                                if(xBoundIndx == mesh->xNodeCt-1) // If the mesh boundary is reached, jump to the source
                                    exhaustedRay = true;
                                xIndx++;
                                xBoundIndx++;
                            }
                            else
                            {
                                if(xBoundIndx == 0)
                                    exhaustedRay = true;
                                xIndx--;
                                xBoundIndx--;
                            }
                        }
                        else if(dirHitFirst == DIRECTION_Y) // y direction
                        {
                            x += tmin*xcos;
                            y = mesh->yNodes[yBoundIndx];
                            z += tmin*zcos;
                            if(ycos >= 0)
                            {
                                if(yBoundIndx == mesh->yNodeCt-1)
                                    exhaustedRay = true;
                                yIndx++;
                                yBoundIndx++;
                            }
                            else
                            {
                                if(yBoundIndx == 0)
                                    exhaustedRay = true;
                                yIndx--;
                                yBoundIndx--;
                            }
                        }
                        else if(dirHitFirst == DIRECTION_Z) // z direction
                        {
                            x += tmin*xcos;
                            y += tmin*ycos;
                            z = mesh->zNodes[zBoundIndx];
                            if(zcos >= 0)
                            {
                                if(zBoundIndx == mesh->zNodeCt-1)
                                    exhaustedRay = true;
                                zIndx++;
                                zBoundIndx++;
                            }
                            else
                            {
                                if(zBoundIndx == 0)
                                    exhaustedRay = true;
                                zIndx--;
                                zBoundIndx--;
                            }
                        }

                        if(xIndx == srcIndxX && yIndx == srcIndxY && zIndx == srcIndxZ)
                        {
                            // If I'm still in the mesh, compute the last voxel contribution
                            double finalDist = sqrt((x-sx[is])*(x-sx[is]) + (y-sy[is])*(y-sy[is]) + (z-sz[is])*(z-sz[is]));

                            for(unsigned int ie = highestEnergy; ie < xs->groupCount(); ie++)
                            {
                                //       [#]       = [cm] * [b] * [1/cm-b]
                                meanFreePaths[ie] += finalDist * xs->totXs1d(zid, ie) * mesh->atomDensity[xIndx*xjmp + yIndx*yjmp + zIndx];
                            }
                            //}

                            exhaustedRay = true;
                        }

                    } // End of while !exhausted loop

                    for(unsigned int ie = highestEnergy; ie < xs->groupCount(); ie++)
                    {

                        double flx = acceptance * srcStrength[ie] * exp(-meanFreePaths[ie]) / (m_4pi * srcToCellDist * srcToCellDist * sourceCt);

                        if(flx < 0)
                            std::cout << "solver.cpp: (291): Negative flux?" << std::endl;

                        if(flx > 1E6)
                            std::cout << "solver.cpp: (294): flux Too big!" << std::endl;

                        (*uflux)[ie*ejmp + xIndxStart*xjmp + yIndxStart*yjmp + zIndxStart] += static_cast<double>(flx);  //srcStrength * exp(-meanFreePaths[ie]) / (4 * M_PI * srcToCellDist * srcToCellDist);
                    }
                    // Following code for updateing mainProgressbar
                    int finishedfrac = (zIndxStart+1) * (yIndxStart+1) * (xIndxStart+1) * (is+1) * 100 / (mesh->xElemCt * mesh->yElemCt * mesh->zElemCt * sourceCt);
                    if(finishedfrac > oldfraction)
                    {
                        std::cout << '\r'<< "basicRaytraceCPU: finishedfrac = " << finishedfrac << std::flush;
                        oldfraction += 10;  // Update mainProgressbar every 10%
                    }

                } // End of each voxel
    } // End of each point source
    std::cout<<std::endl;
    return uflux;
}

std::vector<double>* Solver::raytraceIsoCPU(const Quadrature* quad, const Mesh* mesh, const XSection* xs,
                            const SourceParams* srcPar)
{
    std::clock_t startMoment = std::clock();

    std::vector<double>* uflux = basicRaytraceCPU(quad, mesh, xs, srcPar);

    std::cout << "rayTraceISOCPU: Time to complete raytracer: "
              << (std::clock() - startMoment)/(double)(CLOCKS_PER_SEC/1000)
              << " ms" << std::endl;
    //
    //Output the uncollided flux to a file
    OutWriter::writeScalarFlux("uncollided_flux_iso.dat", *xs, *mesh, *uflux);

    return uflux;
}

std::vector<double>* Solver::raytraceLegendreCPU(const Quadrature* quad, const Mesh* mesh,
                                 const XSection* xs, const SourceParams* params)
{
    std::clock_t startMoment = std::clock();

    unsigned int groups = xs->groupCount();

    //RAY_T sx = 25.3906f;
    //RAY_T sy = 50.0f - 46.4844f;
    //RAY_T sz = 6.8906f;
    const double sx = static_cast<double>(params->sourceX);
    const double sy = static_cast<double>(params->sourceY);
    const double sz = static_cast<double>(params->sourceZ);

    unsigned int ejmp = mesh->voxelCount();
    unsigned int xjmp = mesh->xjmp();
    unsigned int yjmp = mesh->yjmp();

    std::vector<double> *uflux = basicRaytraceCPU(quad, mesh, xs, params);

    std::cout << "basicRaytraceCPU: Time to complete raytracer: "
              << (std::clock() - startMoment)/(double)(CLOCKS_PER_SEC/1000)
              << " ms" << std::endl;

    std::vector<double> *ufluxAng = new std::vector<double>;
    ufluxAng->resize(groups * quad->angleCount() * mesh->voxelCount(), 0.0f);

    int eajmp = quad->angleCount() * mesh->voxelCount();
    int aajmp = mesh->voxelCount();
    int xajmp = mesh->yElemCt * mesh->zElemCt;
    int yajmp = mesh->zElemCt;

    for(unsigned int ie = 0; ie < groups; ie++)
        for(unsigned int iz = 0; iz < mesh->zElemCt; iz++)
            for(unsigned int iy = 0; iy < mesh->yElemCt; iy++)
                for(unsigned int ix = 0; ix < mesh->xElemCt; ix++)  // For every voxel
                {
                    float x = mesh->xNodes[ix] + mesh->dx[ix]/2;
                    float y = mesh->yNodes[iy] + mesh->dy[iy]/2;
                    float z = mesh->zNodes[iz] + mesh->dz[iz]/2;

                    float deltaX = x - sx;
                    float deltaY = y - sy;
                    float deltaZ = z - sz;

                    // normalize to unit vector
                    float mag = sqrt(deltaX*deltaX + deltaY*deltaY + deltaZ*deltaZ);
                    deltaX /= mag;
                    deltaY /= mag;
                    deltaZ /= mag;

                    unsigned int bestAngIndx = 0;
                    float bestCosT = -1.0f;

                    for(unsigned int ia = 0; ia < quad->angleCount(); ia++)
                    {
                        // Cosine of the angle between Sn direction and X-ray beam direction
                        float cosT = quad->mu[ia] * deltaX + quad->eta[ia] * deltaY + quad->zi[ia] * deltaZ;
                        if(cosT > bestCosT)
                        {
                            bestCosT = cosT;
                            bestAngIndx = ia;
                        }
                    }

                    (*ufluxAng)[ie*eajmp + bestAngIndx*aajmp + ix*xajmp + iy*yajmp + iz] = (*uflux)[ie*ejmp + ix*xjmp + iy*yjmp + iz]/quad->wt[bestAngIndx];
                }

    std::cout << "raytraceLegendreCPU: Raytracer + anisotropic mapping completed in "
              << (std::clock() - startMoment)/(double)(CLOCKS_PER_SEC/1000)
              << " ms" << std::endl;

    // Write uncollided flux to a file
    // OutWriter::writeAngularFlux("uncollided_angular_leg.dat", *xs, *quad, *mesh, *ufluxAng);

    // Return ufluxAng
    return ufluxAng;

}

std::vector<double>* Solver::gsSolverIsoCPU(const Quadrature* quad, const Mesh* mesh,
                            const XSection* xs, const SourceParams *srcPar,
                            const std::vector<double>* uFlux)
{
    std::cout << "gsSolverIsoCPU Solving "
              << mesh->voxelCount() * quad->angleCount() * xs->groupCount()
              << " elements in phase space" << std::endl;

    std::cout << "voxelCount:" << mesh->voxelCount() << std::endl;
    std::cout << "angleCount:" << quad->angleCount() << std::endl;
    std::cout << "groupCount:" << xs->groupCount() << std::endl;

    std::clock_t startTime = std::clock();

    const int maxIterations = 25;
    const double epsilon = 0.001f;

    std::vector<double>* cFlux = new std::vector<double>(xs->groupCount() * mesh->voxelCount(), 0.0f);
    std::vector<double> nextFlux(mesh->voxelCount(), 0.0f);
    std::vector<double> totalSource(mesh->voxelCount(), -100.0f);
    std::vector<double> outboundFluxX(mesh->voxelCount(), -100.0f);
    std::vector<double> outboundFluxY(mesh->voxelCount(), -100.0f);
    std::vector<double> outboundFluxZ(mesh->voxelCount(), -100.0f);

    std::vector<double> errMaxList;
    std::vector<std::vector<double> > errList;
    std::vector<std::vector<double> > errIntList;
    std::vector<int> converganceIters;
    std::vector<double> converganceTracker;

    errMaxList.resize(xs->groupCount());
    errList.resize(xs->groupCount());
    errIntList.resize(xs->groupCount());
    converganceIters.resize(xs->groupCount());
    converganceTracker.resize(xs->groupCount());

    double influxX = 0.0f;
    double influxY = 0.0f;
    double influxZ = 0.0f;

    double outx;
    double outy;
    double outz;

    int xjmp = mesh->xjmp();
    int yjmp = mesh->yjmp();

    if(uFlux == nullptr && srcPar == nullptr) 
    {
        std::cout << "uFlux and source params cannot both be NULL" << std::endl;
    }

    bool noDownscatterYet = true;
    unsigned int highestEnergy = 0;

    while(noDownscatterYet)
    {
        double dmax = 0.0;
        unsigned int vc = mesh->voxelCount();
        for(unsigned int ira = 0; ira < vc; ira++)
        {
            dmax = (dmax > (*uFlux)[highestEnergy*vc + ira]) ? dmax : (*uFlux)[highestEnergy*vc + ira];
        }
        if(dmax <= 0.0)
        {
            std::cout << "No external source or downscatter, skipping energy group "
                      << highestEnergy << std::endl;
            highestEnergy++;
        }
        else
        {
            noDownscatterYet = false;
        }

        if(highestEnergy >= xs->groupCount())
        {
            std::cout << "Zero flux everywhere from the raytracer" << std::endl;
            return nullptr;
        }
    }

    // for every energy group
    for(unsigned int ie = highestEnergy; ie < xs->groupCount(); ie++)
    {
        int iterNum = 1;
        double maxDiff = 1.0;
        double totDiff = 1.0E30;

        for(unsigned int i = 0; i < totalSource.size(); i++)
            totalSource[i] = 0;

        for(unsigned int iie = highestEnergy; iie < ie; iie++)
            for(unsigned int ir = 0; ir < mesh->voxelCount(); ir++)
            {
                //         [#]  += [#/cm^2]                                        * [b]                                        * [1/b-cm]              * [cm^3]
                totalSource[ir] += (*cFlux)[iie*mesh->voxelCount() + ir] * xs->scatXs2d(mesh->zoneId[ir], iie, ie, 0) * mesh->atomDensity[ir] * mesh->vol[ir];
            }

        for(unsigned int iie = highestEnergy; iie <= ie; iie++) // Source energy
            for(unsigned int ir = 0; ir < mesh->voxelCount(); ir++)
            {
                //         [#]  += [#/cm^2]                                        * [b]                                        *  [1/b-cm]             * [cm^3]
                totalSource[ir] += (*uFlux)[iie*mesh->voxelCount() + ir] * xs->scatXs2d(mesh->zoneId[ir], iie, ie, 0) * mesh->atomDensity[ir] * mesh->vol[ir];
            }


        while(iterNum <= maxIterations && maxDiff > epsilon)  // while not converged. 
        {
            // Clear for a new sweep
            for(unsigned int i = 0; i < nextFlux.size(); i++)
                nextFlux[i] = 0;

            // At a given energy-level, loop over angular-direction
            for(unsigned int iang = 0; iang < quad->angleCount(); iang++)  // for every angle
            {
                // Find the correct direction and starting voxel for spatial-sweep
                int izStart = 0;                  // Sweep start index
                int diz = 1;                      // Sweep direction
                if(quad->eta[iang] < 0)           // Condition to sweep backward
                {
                    izStart = mesh->zElemCt - 1;  // Start at the far end
                    diz = -1;                     // Sweep toward zero
                }

                int iyStart = 0;
                int diy = 1;
                if(quad->zi[iang] < 0)
                {
                    iyStart = mesh->yElemCt - 1;
                    diy = -1;
                }

                int ixStart = 0;
                int dix = 1;
                if(quad->mu[iang] < 0)
                {
                    ixStart = mesh->xElemCt - 1;
                    dix = -1;
                }

                int iz = izStart;
                while(iz < (signed) mesh->zElemCt && iz >= 0)
                {
                    int iy = iyStart;
                    while(iy < (signed) mesh->yElemCt && iy >= 0)
                    {
                        int ix = ixStart;
                        while(ix < (signed) mesh->xElemCt && ix >= 0)  // for every mesh element in the proper order
                        {
                            unsigned int ir = ix*xjmp + iy*yjmp + iz;
                            int zid = mesh->zoneId[ir];  // Get the zone id of this element

                            // Handle the x influx
                            if(quad->mu[iang] >= 0)                                       // Approach x = 0 -> xMesh
                            {
                                if(ix == 0)                                               // If this is a boundary cell
                                    influxX = 0.0f;                                       // then the in-flux is zero
                                else                                                      // otherwise
                                    influxX = outboundFluxX[(ix-1)*xjmp + iy*yjmp + iz];  // the in-flux is the out-flux from the previous cell
                            }
                            else                                                          // Approach x = xMesh-1 -> 0
                            {
                                if(ix == (signed) mesh->xElemCt-1)
                                    influxX = 0.0f;
                                else
                                    influxX = outboundFluxX[(ix+1)*xjmp + iy*yjmp + iz];
                            }

                            // Handle the y influx
                            if(quad->zi[iang] >= 0)                                       // Approach y = 0 -> yMesh
                            {
                                if(iy == 0)
                                    influxY = 0.0f;
                                else
                                    influxY = outboundFluxY[ix*xjmp + (iy-1)*yjmp + iz];
                            }
                            else                                                          // Approach y = yMesh-1 -> 0
                            {
                                if(iy == (signed) mesh->yElemCt-1)
                                    influxY = 0.0f;
                                else
                                    influxY = outboundFluxY[ix*xjmp + (iy+1)*yjmp + iz];
                            }

                            // Handle the z influx
                            if(quad->eta[iang] >= 0)
                            {
                                if(iz == 0)
                                    influxZ = 0.0f;
                                else
                                    influxZ = outboundFluxZ[ix*xjmp + iy*yjmp + iz-1];
                            }
                            else
                            {
                                if(iz == (signed) mesh->zElemCt-1)
                                    influxZ = 0.0f;
                                else
                                    influxZ = outboundFluxZ[ix*xjmp + iy*yjmp + iz+1];
                            }

                            double inscatter = (*cFlux)[ie*mesh->voxelCount() + ir] * xs->scatXs2d(zid, ie, ie, 0) * mesh->atomDensity[ir] * mesh->vol[ir];

                            double numer = totalSource[ir] +  inscatter +                            // [#]
                                    mesh->Ayz[ie*quad->angleCount()*mesh->yElemCt*mesh->zElemCt + iang*mesh->yElemCt*mesh->zElemCt + iy*mesh->zElemCt + iz] * influxX +  // [cm^2 * #/cm^2]  The 2x is already factored in
                                    mesh->Axz[ie*quad->angleCount()*mesh->xElemCt*mesh->zElemCt + iang*mesh->xElemCt*mesh->zElemCt + ix*mesh->zElemCt + iz] * influxY +
                                    mesh->Axy[ie*quad->angleCount()*mesh->xElemCt*mesh->yElemCt + iang*mesh->xElemCt*mesh->yElemCt + ix*mesh->yElemCt + iy] * influxZ;

                            double denom = mesh->vol[ir]*xs->totXs1d(zid, ie)*mesh->atomDensity[ir] +                               // [cm^3] * [b] * [1/b-cm]
                                    mesh->Ayz[ie*quad->angleCount()*mesh->yElemCt*mesh->zElemCt + iang*mesh->yElemCt*mesh->zElemCt + iy*mesh->zElemCt + iz] +            // [cm^2]
                                    mesh->Axz[ie*quad->angleCount()*mesh->xElemCt*mesh->zElemCt + iang*mesh->xElemCt*mesh->zElemCt + ix*mesh->zElemCt + iz] +
                                    mesh->Axy[ie*quad->angleCount()*mesh->xElemCt*mesh->yElemCt + iang*mesh->xElemCt*mesh->yElemCt + ix*mesh->yElemCt + iy];

                            //   [#/cm^2] = [#]  / [cm^2]

                            double angFlux = numer/denom;

                            outx = 2*angFlux - influxX;
                            outy = 2*angFlux - influxY;
                            outz = 2*angFlux - influxZ;

                            bool cflag = false;
                            if(outx < 0)
                            {
                                cflag = true;
                                outx = 0;
                            }
                            if(outy < 0)
                            {
                                cflag = true;
                                outy = 0;
                            }
                            if(outz < 0)
                            {
                                cflag = true;
                                outz = 0;
                            }
                            if(cflag)
                            {
                                angFlux = (influxX + influxY + influxZ + outx + outy + outz)/6.0;
                            }

                            outboundFluxX[ix*xjmp + iy*yjmp + iz] = outx;
                            outboundFluxY[ix*xjmp + iy*yjmp + iz] = outy;
                            outboundFluxZ[ix*xjmp + iy*yjmp + iz] = outz;

                            // Sum all the angular fluxes
                            nextFlux[ir] += quad->wt[iang]*angFlux;

                            ix += dix;
                        } // end of for ix

                        iy += diy;
                    } // end of for iy

                    iz += diz;
                } // end of for iz

                unsigned int xTracked = mesh->xElemCt/2;
                unsigned int yTracked = mesh->yElemCt/2;
                unsigned int zTracked = mesh->zElemCt/2;
                converganceTracker.push_back(nextFlux[xTracked*xjmp + yTracked*yjmp + zTracked]);

            } // end of all angles

            maxDiff = -1.0E35f;
            totDiff = 0.0f;
            for(unsigned int i = 0; i < nextFlux.size(); i++)
            {
                maxDiff = std::max(maxDiff, fabs((nextFlux[i] - (*cFlux)[ie*mesh->voxelCount() + i])/nextFlux[i]));
                totDiff += fabs(nextFlux[i] - (*cFlux)[ie*mesh->voxelCount() + i]);
            }

            for(unsigned int i = 0; i < nextFlux.size(); i++)
                (*cFlux)[ie*mesh->voxelCount() + i] = nextFlux[i];

            errList[ie].push_back(maxDiff);
            errIntList[ie].push_back(totDiff);
            errMaxList[ie] = maxDiff;
            converganceIters[ie] = iterNum;

            iterNum++;
        } // end not converged 

        if(!(iterNum <= maxIterations))
        {
            std::cout << "Max iterations hit" << std::endl;
        }
        else if(!(maxDiff > epsilon))
        {
            std::cout << "Converged on relative error" << std::endl;
        }
        else
        {
            std::cout << "Converged on precision bound" << std::endl;
        }
        // For updating Progress
        std::cout << (ie+1) * 100 / xs->groupCount() << "% finished" << std::endl;

    }  // end each energy group

    std::cout << "Time to complete: "
              << (std::clock() - startTime)/(double)(CLOCKS_PER_SEC/1000)
              << " ms" << std::endl;

    for(unsigned int i = 0; i < errList.size(); i++)
    {
        std::cout << "%Group: " << i << "   maxDiff: " << errMaxList[i]
                  << "   Iterations: " << converganceIters[i] << std::endl;
        std::cout << "cpu" << i << " = [";
        for(unsigned int j = 0; j < errList[i].size(); j++)
            std::cout << errList[i][j] << ",\t";
        std::cout << "];\ncpu" << i << "i = [";
        for(unsigned int j = 0; j < errIntList[i].size(); j++)
            std::cout << errIntList[i][j] << ",\t";
        std::cout << "];" << std::endl;
    }

    //Write collided flux to file
    OutWriter::writeScalarFlux("collided_flux_iso.dat", *xs, *mesh, *cFlux);

    return cFlux;
}


// ////////////////////////////////////////////////////////////////////////////////////////////// //
//                           Anisotropic versions of the above solvers                            //
// ////////////////////////////////////////////////////////////////////////////////////////////// //

std::vector<double>* Solver::gsSolverHarmonicCPU(const Quadrature* quad, const Mesh* mesh,
    const XSection* xs, const SourceParams* srcPar,
    const std::vector<double>* uflux)
{

    // Do some input checks
    if (this->pn > 8)
    {
        std::cout << "Pn check failed, Pn = " << this->pn << std::endl;
        std::cout << "Maximum Pn = 8 " << std::endl;
        // return nullptr;
    }

    std::cout << "Pn = " << this->pn << std::endl;

    std::clock_t startTime = std::clock();

    const int maxIterations = 50;
    const float epsilon = 0.001f;


    std::vector<SOL_T>* scalarFlux = new std::vector<SOL_T>(xs->groupCount() * mesh->voxelCount(), 0.0f);
    //std::vector<SOL_T> *moments = new std::vector<SOL_T>(xs->groupCount() * mesh->voxelCount() * momentCount, 0.0f);
    std::vector<SOL_T> tempScalarFlux(mesh->voxelCount());
    //std::vector<SOL_T> preFlux(mesh->voxelCount(), -100.0f);
    //std::vector<SOL_T> totalSource(mesh->voxelCount(), -100.0f);
    std::vector<SOL_T> outboundFluxX(mesh->voxelCount(), -100.0f);
    std::vector<SOL_T> outboundFluxY(mesh->voxelCount(), -100.0f);
    std::vector<SOL_T> outboundFluxZ(mesh->voxelCount(), -100.0f);
    //std::vector<SOL_T> extSource(xs->groupCount() * mesh->voxelCount(), 0.0f);
    // std::vector<SOL_T> momentToDiscrete(quad->angleCount() * momentCount);

    std::vector<SOL_T> errMaxList;
    std::vector<std::vector<SOL_T> > errList;
    std::vector<int> converganceIters;
    std::vector<SOL_T> converganceTracker;

    errMaxList.resize(xs->groupCount());
    errList.resize(xs->groupCount());
    converganceIters.resize(xs->groupCount());
    converganceTracker.resize(xs->groupCount());

    SOL_T outx;
    SOL_T outy;
    SOL_T outz;

    if (uflux == NULL && srcPar == NULL)
    {
        std::cout << "uFlux and params cannot both be NULL" << std::endl;
    }

    SOL_T influxX = 0.0f;
    SOL_T influxY = 0.0f;
    SOL_T influxZ = 0.0f;

    int ejmp = mesh->voxelCount();
    int xjmp = mesh->xjmp();
    int yjmp = mesh->yjmp();

    // Try to skip over any empty energy groups
    bool noDownscatterYet = true;
    unsigned int highestEnergy = 0;

    while (noDownscatterYet)
    {
        double dmax = 0.0;
        unsigned int vc = mesh->voxelCount();
        for (unsigned int ira = 0; ira < vc; ira++)
        {
            dmax = (dmax > (*uflux)[highestEnergy * vc + ira]) ? dmax : (*uflux)[highestEnergy * vc + ira];
        }
        if (dmax <= 0.0)
        {
            std::cout << "No external source or downscatter, skipping energy group "
                << highestEnergy << std::endl;
            highestEnergy++;
        }
        else
        {
            noDownscatterYet = false;
        }

        if (highestEnergy >= xs->groupCount())
        {
            std::cout << "Zero flux everywhere from the raytracer" << std::endl;
        }
    }

    std::vector<SOL_T>* umoments = new std::vector<SOL_T>;     // Uncollided flux moment
    std::vector<SOL_T>* cmoments = new std::vector<SOL_T>;     // Collided flux moment
    std::vector<SOL_T>* inmoments = new std::vector<SOL_T>;    // Summation of moments over groups

    const unsigned int momentCount = (this->pn + 1) * (this->pn + 1);
    std::cout << "MonmentCount = " << momentCount << std::endl;
    umoments->resize(xs->groupCount() * mesh->voxelCount() * momentCount, 0.0f);  // Uncollided flux moment
    cmoments->resize(xs->groupCount() * mesh->voxelCount() * momentCount, 0.0f);  // Collided flux moment
    inmoments->resize(mesh->voxelCount() * momentCount, 0.0f);                    // Summation of moments over groups

    std::vector<SOL_T> tempMoments(mesh->voxelCount() * momentCount);  // temporary moments storage

    int epjmp = mesh->voxelCount() * momentCount;
    int xpjmp = mesh->yElemCt * mesh->zElemCt * momentCount;
    int ypjmp = mesh->zElemCt * momentCount;
    int zpjmp = momentCount;

    const RAY_T sx = static_cast<RAY_T>(srcPar->sourceX);
    const RAY_T sy = static_cast<RAY_T>(srcPar->sourceY);
    const RAY_T sz = static_cast<RAY_T>(srcPar->sourceZ);

    SphericalHarmonic harmonic;
    AssocLegendre aLegendre;

    std::cout << "Solver::gsSolverHarmonic solving "
        << mesh->voxelCount() * momentCount * xs->groupCount()
        << " elements in phase space"
        << std::endl;

    // for every energy group
    for (unsigned int ie = highestEnergy; ie < xs->groupCount(); ie++)
    {
        std::cout << "Energy group #" << ie << std::endl;

        int iterNum = 1;
        SOL_T maxDiff = 1.0;

        // Initialize the inmoments to zero
        for (unsigned int i = 0; i < inmoments->size(); i++)
            (*inmoments)[i] = 0.0f;
            
        // Update uncollided flux moment and generate source in-moments for current energy group
        for (unsigned int iz = 0; iz < mesh->zElemCt; iz++)
            for (unsigned int iy = 0; iy < mesh->yElemCt; iy++)
                for (unsigned int ix = 0; ix < mesh->xElemCt; ix++)  // For every voxel
                {
                    const unsigned int ir = ix * xjmp + iy * yjmp + iz;  // global voxel index
                    int zid = mesh->zoneId[ir];  // Get the zone id of this voxel

                    SOL_T x = mesh->xNodes[ix] + mesh->dx[ix] / 2;
                    SOL_T y = mesh->yNodes[iy] + mesh->dy[iy] / 2;
                    SOL_T z = mesh->zNodes[iz] + mesh->dz[iz] / 2;

                    SOL_T deltaX = x - sx;
                    SOL_T deltaY = y - sy;
                    SOL_T deltaZ = z - sz;

                    // normalize to unit vector
                    SOL_T mag = sqrt(deltaX * deltaX + deltaY * deltaY + deltaZ * deltaZ);
                    deltaX /= mag;  // sin(theta)sin(phi)
                    deltaY /= mag;  // cos(theta)
                    deltaZ /= mag;  // sin(theta)cos(phi)

                    SOL_T phi = atan2(deltaX, deltaZ);
                    SOL_T theta = acos(deltaY);  // choose Y axis as the polar axis as X-ray is shooting down along Y axis

                    for (unsigned int il = 0; il <= this->pn; il++)
                    {
                        // Regular uncollided flux moment , i.e. m = 0
                        (*umoments)[ie * epjmp + ix * xpjmp + iy * ypjmp + iz * zpjmp + il * il] = (*uflux)[ie * ejmp + ix * xjmp + iy * yjmp + iz] * harmonic.yl0(il, theta);

                        // Generate regular source in-moments 
                        for (unsigned int iep = highestEnergy; iep <= ie; iep++)
                        {
                            double total_zeroMoments = (*umoments)[iep * epjmp + ix * xpjmp + iy * ypjmp + iz * zpjmp + il * il]
                                + (*cmoments)[iep * epjmp + ix * xpjmp + iy * ypjmp + iz * zpjmp + il * il];  // Note: cmoments=0 when iep=ie
                                //[cm^2] =                          [b]   *  [1/b-cm]             * [cm^3]
                            double coeff = xs->scatXs2d(zid, iep, ie, il) * mesh->atomDensity[ir] * mesh->vol[ir];
                            
                            (*inmoments)[ir * momentCount + il * il] += coeff * total_zeroMoments;  // zero source in-moments
                        }

                        for (unsigned int im = 1; im <= il; im++)
                        {
                            (*umoments)[ie * epjmp + ix * xpjmp + iy * ypjmp + iz * zpjmp + il * il + 2 * im - 1] = (*uflux)[ie * ejmp + ix * xjmp + iy * yjmp + iz] * harmonic.ylm_o(il, im, theta, phi);
                            (*umoments)[ie * epjmp + ix * xpjmp + iy * ypjmp + iz * zpjmp + il * il + 2 * im] = (*uflux)[ie * ejmp + ix * xjmp + iy * yjmp + iz] * harmonic.ylm_e(il, im, theta, phi);

                            // Generate cosine and sine source in-moments
                            for (unsigned int iep = highestEnergy; iep <= ie; iep++)
                            {
                                double total_sinMoments = (*umoments)[iep * epjmp + ix * xpjmp + iy * ypjmp + iz * zpjmp + il * il + 2 * im - 1]
                                    + (*cmoments)[iep * epjmp + ix * xpjmp + iy * ypjmp + iz * zpjmp + il * il + 2 * im - 1];
                                double total_cosMoments = (*umoments)[iep * epjmp + ix * xpjmp + iy * ypjmp + iz * zpjmp + il * il + 2 * im]
                                    + (*cmoments)[iep * epjmp + ix * xpjmp + iy * ypjmp + iz * zpjmp + il * il + 2 * im];

                                double coeff = xs->scatXs2d(zid, iep, ie, il) * mesh->atomDensity[ir] * mesh->vol[ir];
                                
                                (*inmoments)[ir * momentCount + il * il + 2 * im - 1] += coeff * total_sinMoments;  // sine source in-moments
                                (*inmoments)[ir * momentCount + il * il + 2 * im] += coeff * total_cosMoments;      // cosine source in-moments
                            }
                        }
                    }
                }

        // while not converged
        while (iterNum <= maxIterations && maxDiff > epsilon)
        {
            // Clear for a new sweep
            for (unsigned int i = 0; i < tempScalarFlux.size(); i++)
                tempScalarFlux[i] = 0;

            for (unsigned int i = 0; i < tempMoments.size(); i++)
                tempMoments[i] = 0;

            // For every angle
            for (unsigned int ia = 0; ia < quad->angleCount(); ia++)
            {
                // Find the angle (theta,phi)
                SOL_T theta = acos(quad->eta[ia]); // Choose Y axis as polar axis
                SOL_T phi = atan2(quad->mu[ia], quad->zi[ia]);

                // Find the correct direction to sweep
                int izStart = 0;                  // Sweep start index
                int diz = 1;                      // Sweep direction
                if (quad->zi[ia] < 0)
                    // Condition to sweep backward
                {
                    izStart = mesh->zElemCt - 1;  // Start at the far end
                    diz = -1;                     // Sweep toward zero
                }

                int iyStart = 0;
                int diy = 1;
                if (quad->eta[ia] < 0)
                {
                    iyStart = mesh->yElemCt - 1;
                    diy = -1;
                }

                int ixStart = 0;
                int dix = 1;
                if (quad->mu[ia] < 0)
                {
                    ixStart = mesh->xElemCt - 1;
                    dix = -1;
                }

                int iz = izStart;
                while (iz < (signed)mesh->zElemCt && iz >= 0)
                {
                    int iy = iyStart;
                    while (iy < (signed)mesh->yElemCt && iy >= 0)
                    {
                        int ix = ixStart;
                        while (ix < (signed)mesh->xElemCt && ix >= 0)  // for every mesh element in the proper order
                        {
                            const unsigned int ri = ix * xjmp + iy * yjmp + iz;  // global voxel index
                            int zid = mesh->zoneId[ri];  // Get the zone id of this element

                            // Handle the x influx
                            if (quad->mu[ia] >= 0)                                       // Approach x = 0 -> xMesh
                            {
                                if (ix == 0)                                               // If this is a boundary cell
                                    influxX = 0.0f;                                       // then the in-flux is zero
                                else                                                      // otherwise
                                    influxX = outboundFluxX[(ix - 1) * xjmp + iy * yjmp + iz];  // the in-flux is the out-flux from the previous cell
                            }
                            else                                                          // Approach x = xMesh-1 -> 0
                            {
                                if (ix == (signed)mesh->xElemCt - 1)
                                    influxX = 0.0f;
                                else
                                    influxX = outboundFluxX[(ix + 1) * xjmp + iy * yjmp + iz];
                            }

                            // Handle the y influx
                            if (quad->eta[ia] >= 0)                                       // Approach y = 0 -> yMesh
                            {
                                if (iy == 0)
                                    influxY = 0.0f;
                                else
                                    influxY = outboundFluxY[ix * xjmp + (iy - 1) * yjmp + iz];
                            }
                            else                                                          // Approach y = yMesh-1 -> 0
                            {
                                if (iy == (signed)mesh->yElemCt - 1)
                                    influxY = 0.0f;
                                else
                                    influxY = outboundFluxY[ix * xjmp + (iy + 1) * yjmp + iz];
                            }

                            // Handle the z influx
                            if (quad->zi[ia] >= 0)
                            {
                                if (iz == 0)
                                    influxZ = 0.0f;
                                else
                                    influxZ = outboundFluxZ[ix * xjmp + iy * yjmp + iz - 1];
                            }
                            else
                            {
                                if (iz == (signed)mesh->zElemCt - 1)
                                    influxZ = 0.0f;
                                else
                                    influxZ = outboundFluxZ[ix * xjmp + iy * yjmp + iz + 1];
                            }

                            // Calculate down-scattering source
                            double totalSource = 0.0;

                            // downscatter from the ucollided and collided flux including group inscatter   
                            double coeff = 0.0;
                            for (unsigned int il = 0; il <= this->pn; il++)
                            {
                                //[cm^2] =                          [b]            *  [1/b-cm]             * [cm^3]
                                coeff = xs->scatXs2d(mesh->zoneId[ri], ie, ie, il) * mesh->atomDensity[ri] * mesh->vol[ri];
                                

                                //[#/cm^2] += [#/cm^2]
                                totalSource += (2.0*il+1) * aLegendre(il, 0, cos(theta)) * (coeff * (*cmoments)[ie * epjmp + ix * xpjmp + iy * ypjmp + iz * zpjmp + il * il]
                                    + (*inmoments)[ri * momentCount + il * il]);
                                
                                float n_const = 2.0 * (2.0 * il + 1); // Normalization constant
                                for (unsigned int im = 1; im <= il; im++)
                                {
                                    totalSource += n_const * factorial(il-im) / factorial(il+im) * aLegendre(il, im, cos(theta)) * sin(im*phi) * (coeff * (*cmoments)[ie * epjmp + ix * xpjmp + iy * ypjmp + iz * zpjmp + il * il + 2 * im - 1]
                                        + (*inmoments)[ri * momentCount + il * il + 2 * im - 1]);
                                    totalSource += n_const * factorial(il-im) / factorial(il+im) * aLegendre(il, im, cos(theta)) * cos(im*phi) * (coeff * (*cmoments)[ie * epjmp + ix * xpjmp + iy * ypjmp + iz * zpjmp + il * il + 2 * im]
                                        + (*inmoments)[ri * momentCount + il * il + 2 * im]);
                                }

                            }

                            double numer = totalSource +                                                                                                               // [#]
                                mesh->Ayz[ie * quad->angleCount() * mesh->yElemCt * mesh->zElemCt + ia * mesh->yElemCt * mesh->zElemCt + iy * mesh->zElemCt + iz] * influxX +  // [cm^2 * #/cm^2]  
                                mesh->Axz[ie * quad->angleCount() * mesh->xElemCt * mesh->zElemCt + ia * mesh->xElemCt * mesh->zElemCt + ix * mesh->zElemCt + iz] * influxY +
                                mesh->Axy[ie * quad->angleCount() * mesh->xElemCt * mesh->yElemCt + ia * mesh->xElemCt * mesh->yElemCt + ix * mesh->yElemCt + iy] * influxZ;
                            double denom = mesh->vol[ix * xjmp + iy * yjmp + iz] * xs->totXs1d(zid, ie) * mesh->atomDensity[ix * xjmp + iy * yjmp + iz] +                               // [cm^3] * [b] * [1/b-cm]
                                mesh->Ayz[ie * quad->angleCount() * mesh->yElemCt * mesh->zElemCt + ia * mesh->yElemCt * mesh->zElemCt + iy * mesh->zElemCt + iz] +            // [cm^2]
                                mesh->Axz[ie * quad->angleCount() * mesh->xElemCt * mesh->zElemCt + ia * mesh->xElemCt * mesh->zElemCt + ix * mesh->zElemCt + iz] +
                                mesh->Axy[ie * quad->angleCount() * mesh->xElemCt * mesh->yElemCt + ia * mesh->xElemCt * mesh->yElemCt + ix * mesh->yElemCt + iy];

                            //           [#/cm^2] = [#]  / [cm^2]
                            double cellAvgAngFlux = numer / denom;

                            if (std::isnan(cellAvgAngFlux))
                            {
                                std::cout << "Found a nan!" << std::endl;
                                std::cout << "Vol = " << mesh->vol[ix * xjmp + iy * yjmp + iz] << std::endl;
                                std::cout << "xs = " << xs->totXs1d(zid, ie) << std::endl;
                                std::cout << "Ayz = " << mesh->Ayz[ia * mesh->yElemCt * mesh->zElemCt + iy * mesh->zElemCt + iz] << std::endl;
                                std::cout << "Axz = " << mesh->Axz[ia * mesh->xElemCt * mesh->zElemCt + ix * mesh->zElemCt + iz] << std::endl;
                                std::cout << "Axy = " << mesh->Axy[ia * mesh->xElemCt * mesh->yElemCt + ix * mesh->yElemCt + iy] << std::endl;
                            }

                            outx = 2 * cellAvgAngFlux - influxX;
                            outy = 2 * cellAvgAngFlux - influxY;
                            outz = 2 * cellAvgAngFlux - influxZ;

                            bool cflag = false;
                            if (outx < 0)
                            {
                                cflag = true;
                                outx = 0;
                            }
                            if (outy < 0)
                            {
                                cflag = true;
                                outy = 0;
                            }
                            if (outz < 0)
                            {
                                cflag = true;
                                outz = 0;
                            }
                            if (cflag)
                            {
                                cellAvgAngFlux = (influxX + influxY + influxZ + outx + outy + outz) / 6.0;
                            }

                            outboundFluxX[ix * xjmp + iy * yjmp + iz] = outx;
                            outboundFluxY[ix * xjmp + iy * yjmp + iz] = outy;
                            outboundFluxZ[ix * xjmp + iy * yjmp + iz] = outz;

                            // Sum all the angular fluxes
                            tempScalarFlux[ix * xjmp + iy * yjmp + iz] += quad->wt[ia] * cellAvgAngFlux;

                            // Update the collided flux moment for energy group ie
                            for (unsigned int il = 0; il <= this->pn; il++)
                            {
                                tempMoments[ix * xpjmp + iy * ypjmp + iz * zpjmp + il * il] += quad->wt[ia] * cellAvgAngFlux * aLegendre(il, 0, cos(theta));
                                for (unsigned int im = 1; im <= il; im++)
                                {
                                    tempMoments[ix * xpjmp + iy * ypjmp + iz * zpjmp + il * il + 2 * im - 1] += quad->wt[ia] * cellAvgAngFlux * aLegendre(il, im, cos(theta)) * sin(im * phi);
                                    tempMoments[ix * xpjmp + iy * ypjmp + iz * zpjmp + il * il + 2 * im] += quad->wt[ia] * cellAvgAngFlux * aLegendre(il, im, cos(theta)) * cos(im * phi);
                                }
                            }

                            ix += dix;
                        } // end of for ix

                        iy += diy;
                    } // end of for iy

                    iz += diz;
                } // end of for iz

            } // end of all angles

            unsigned int xTracked = mesh->xElemCt / 2;
            unsigned int yTracked = mesh->yElemCt / 2;
            unsigned int zTracked = mesh->zElemCt / 2;
            converganceTracker.push_back(tempScalarFlux[xTracked * xjmp + yTracked * yjmp + zTracked]);

            maxDiff = -1.0E35f;
            for (unsigned int i = 0; i < tempScalarFlux.size(); i++)
            {
                maxDiff = std::max(maxDiff, fabs((tempScalarFlux[i] - (*scalarFlux)[ie * mesh->voxelCount() + i]) / tempScalarFlux[i]));

                if (std::isnan(maxDiff))
                    std::cout << "Found a diff nan!" << std::endl;
            }
            std::cout << "Max diff = " << maxDiff << std::endl;

            errList[ie].push_back(maxDiff);
            errMaxList[ie] = maxDiff;
            converganceIters[ie] = iterNum;

            for (unsigned int i = 0; i < tempScalarFlux.size(); i++)
            {
                (*scalarFlux)[ie * mesh->voxelCount() + i] = tempScalarFlux[i];
            }

            for (unsigned int i = 0; i < tempMoments.size(); i++)
            {
                (*cmoments)[ie * epjmp + i] = tempMoments[i];
            }

            iterNum++;

        } // end not converged

        // For updating mainProgressbar
        std::cout << (ie + 1) * 100 / xs->groupCount() << "% finished" << std::endl;

    }  // end each energy group


    std::cout << "Time to complete: "
        << (std::clock() - startTime) / (double)(CLOCKS_PER_SEC / 1000.0)
        << " ms" << std::endl;

    for (unsigned int i = 0; i < errList.size(); i++)
    {
        std::cout << "Group: " << i << "   maxDiff: " << errMaxList[i] << std::endl;
        std::cout << "Iterations: " << converganceIters[i] << std::endl;
        for (unsigned int j = 0; j < errList[i].size(); j++)
            std::cout << errList[i][j] << "\t";
        std::cout << "\n" << std::endl;
    }

    //Write angular collided flux moment to file
    OutWriter::writeAngularFlux("collided_angular_harmonic_moment.dat", *xs, *quad, *mesh, *cmoments);

    //Write scalar collided flux to file
    OutWriter::writeScalarFlux("collided_flux_harmonic.dat", *xs, *mesh, *scalarFlux);

    return scalarFlux;
}
