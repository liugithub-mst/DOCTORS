
#include <iostream>
#include <stdio.h>
#include <stdlib.h>
#include "cuda.h"

// GPU profile
void reportGpuData()
{
    std::cout << "Reporting GPU resources ... " << std::endl;

    // Check the number of GPU resources
    int nDevices = 0;
    cudaGetDeviceCount(&nDevices);
    std::cout << "Found " << nDevices << " CUDA devices" << std::endl;

    // Console log
    int driverVersion = 0;
    int runtimeVersion = 0;
    cudaDriverGetVersion(&driverVersion);
    cudaRuntimeGetVersion(&runtimeVersion);
    printf("  CUDA Driver Version / Runtime Version          %d.%d / %d.%d\n",
        driverVersion / 1000, (driverVersion % 100) / 10,
        runtimeVersion / 1000, (runtimeVersion % 100) / 10);

    // Check free and total memory
    size_t free, total;

    for (unsigned int i = 0; i < nDevices; i++)
    {
        // Find a GPU
        cudaDeviceProp props;
        cudaGetDeviceProperties(&props, i);

        std::cout << "Device " << i << " : " << props.name << " with compute "
            << props.major << "." << props.minor << " capability" << std::endl;

        std::cout << "Max threads per block: " << props.maxThreadsPerBlock << std::endl;
        std::cout << "Max grid size: " << props.maxGridSize[0] << " x "
            << props.maxGridSize[1] << " x " << props.maxGridSize[2] << std::endl;

        std::cout << "Memory Clock Rate (KHz): " << props.memoryClockRate << std::endl;
        std::cout << "Memory Bus Width (bits): " << props.memoryBusWidth << std::endl;
        std::cout << "Peak Memory Bandwidth (GB/s): " << (props.memoryClockRate * (props.memoryBusWidth / 8) / 1.0e6) << std::endl;


        int mp = props.multiProcessorCount;
        int cores = 0;

        switch (props.major) {
        case 7: // Turing
            cores = mp * 64;
            break;
        default:
            printf("Unknown device type\n");
            break;
        }

        std::cout << "SMs: " << mp << std::endl;
        std::cout << "CUDA Cores: " << cores << std::endl;

        cudaMemGetInfo(&free, &total);
        std::cout << "GPU " << i << " memory: free= " << free / 1024 /1024 / 1024
                  << "GB , total= " << total / 1024 /1024 / 1024 << "GB\n" << std::endl;
    }
}


// Function to allocate GPU float array and copy data from CPU to GPU
float* alloc_gpuFloat(const int gpuID, const int elements, const float* cpuData)
{
    cudaError_t cudaerr;
    if ((cudaerr = cudaSetDevice(gpuID)) != cudaSuccess)
        std::cout << "alloc_gpuFloat failed to set the device with error code: "
                  << cudaerr << ":" << cudaGetErrorString(cudaerr) << std::endl;

    float* gpuData;
    if ((cudaerr = cudaMalloc(&gpuData, elements * sizeof(float))) != cudaSuccess)
        std::cout << "alloc_gpuFloat error while allocating CUDA memory with error code: "
                  << cudaerr << ":" << cudaGetErrorString(cudaerr) << std::endl;
    if (cpuData != NULL)
    {
        if ((cudaerr = cudaMemcpyAsync(gpuData, cpuData, elements * sizeof(float), cudaMemcpyHostToDevice)) != cudaSuccess)
            std::cout << "alloc_gpuFloat failed while copying data with error code: "
                      << cudaerr << ":" << cudaGetErrorString(cudaerr) << std::endl;
    }

    return gpuData;
}
