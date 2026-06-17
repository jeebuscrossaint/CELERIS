#include "celeris/cuda/eigensolve.hpp"

// This file is PURE HOST CODE: it only calls cuSOLVER's host API and the CUDA
// runtime API. There are no __global__ kernels, so it compiles with the normal
// C++ compiler (MSVC) and just links the CUDA import libs — no nvcc, and no
// dependency on the CCCL headers (which scoop's CUDA package omits). Use the
// host-only <cuda_runtime_api.h>, NOT <cuda_runtime.h> (the latter drags in
// device intrinsics like cuda_fp16.h -> <nv/target>).
#include <cuComplex.h>
#include <cuda_runtime_api.h>
#include <cusolverDn.h>

#include <cstdlib>

namespace celeris::cuda {

bool available() {
    int count = 0;
    return cudaGetDeviceCount(&count) == cudaSuccess && count > 0;
}

// General complex eigendecomposition via cuSOLVER Xgeev (right eigenvectors).
// std::complex<double> and cuDoubleComplex are layout-compatible (two doubles).
bool geev(const std::complex<double>* A, int n,
          std::complex<double>* w, std::complex<double>* vr) {
    const int64_t N = n;
    cusolverDnHandle_t handle = nullptr;
    cusolverDnParams_t params = nullptr;
    if (cusolverDnCreate(&handle) != CUSOLVER_STATUS_SUCCESS) return false;
    cusolverDnCreateParams(&params);

    cuDoubleComplex *dA = nullptr, *dW = nullptr, *dVR = nullptr;
    int* dInfo = nullptr;
    void *workDev = nullptr, *workHost = nullptr;
    size_t lworkDev = 0, lworkHost = 0;
    bool ok = false;

    auto cleanup = [&] {
        if (dA) cudaFree(dA);
        if (dW) cudaFree(dW);
        if (dVR) cudaFree(dVR);
        if (dInfo) cudaFree(dInfo);
        if (workDev) cudaFree(workDev);
        if (workHost) std::free(workHost);
        if (params) cusolverDnDestroyParams(params);
        if (handle) cusolverDnDestroy(handle);
    };

    // Compiled as plain C++ (not nvcc), so cudaMalloc only has the void** form;
    // cast explicitly (no implicit T** -> void** in C++).
    if (cudaMalloc(reinterpret_cast<void**>(&dA),
                   sizeof(cuDoubleComplex) * N * N) != cudaSuccess ||
        cudaMalloc(reinterpret_cast<void**>(&dW),
                   sizeof(cuDoubleComplex) * N) != cudaSuccess ||
        cudaMalloc(reinterpret_cast<void**>(&dVR),
                   sizeof(cuDoubleComplex) * N * N) != cudaSuccess ||
        cudaMalloc(reinterpret_cast<void**>(&dInfo), sizeof(int)) !=
            cudaSuccess) {
        cleanup();
        return false;
    }

    cudaMemcpy(dA, A, sizeof(cuDoubleComplex) * N * N, cudaMemcpyHostToDevice);

    const cusolverEigMode_t novl = CUSOLVER_EIG_MODE_NOVECTOR;
    const cusolverEigMode_t vr_mode = CUSOLVER_EIG_MODE_VECTOR;

    if (cusolverDnXgeev_bufferSize(
            handle, params, novl, vr_mode, N, CUDA_C_64F, dA, N, CUDA_C_64F, dW,
            CUDA_C_64F, nullptr, N, CUDA_C_64F, dVR, N, CUDA_C_64F, &lworkDev,
            &lworkHost) != CUSOLVER_STATUS_SUCCESS) {
        cleanup();
        return false;
    }

    if (lworkDev && cudaMalloc(&workDev, lworkDev) != cudaSuccess) {  // workDev is void*
        cleanup();
        return false;
    }
    if (lworkHost) workHost = std::malloc(lworkHost);

    cusolverStatus_t st = cusolverDnXgeev(
        handle, params, novl, vr_mode, N, CUDA_C_64F, dA, N, CUDA_C_64F, dW,
        CUDA_C_64F, nullptr, N, CUDA_C_64F, dVR, N, CUDA_C_64F, workDev,
        lworkDev, workHost, lworkHost, dInfo);

    int info = 0;
    cudaMemcpy(&info, dInfo, sizeof(int), cudaMemcpyDeviceToHost);

    if (st == CUSOLVER_STATUS_SUCCESS && info == 0) {
        cudaMemcpy(w, dW, sizeof(cuDoubleComplex) * N, cudaMemcpyDeviceToHost);
        cudaMemcpy(vr, dVR, sizeof(cuDoubleComplex) * N * N,
                   cudaMemcpyDeviceToHost);
        ok = true;
    }

    cleanup();
    return ok;
}

} // namespace celeris::cuda
