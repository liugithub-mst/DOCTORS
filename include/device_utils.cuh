
#ifndef DEVICE_INFO_CUH
#define DEVICE_INFO_CUH

void reportGpuData();
float* alloc_gpuFloat(const int gpuID, const int elements, const float* cpuData);

#endif
