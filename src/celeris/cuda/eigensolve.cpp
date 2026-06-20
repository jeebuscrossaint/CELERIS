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

#include <algorithm>
#include <cstdlib>
#include <cstring>
#include <vector>

namespace celeris::cuda {

// available()/device_name() are provided by the kernel TU (propagate.cu) when it
// is built. Define them here only for the cuSOLVER-without-kernel configuration,
// so the two TUs never both define them.
#ifndef CELERIS_USE_CUDA_KERNELS
bool available() {
    int count = 0;
    return cudaGetDeviceCount(&count) == cudaSuccess && count > 0;
}

const char* device_name() {
    static char name[256] = {0};
    int count = 0;
    if (cudaGetDeviceCount(&count) != cudaSuccess || count == 0) return "";
    cudaDeviceProp prop;
    if (cudaGetDeviceProperties(&prop, 0) != cudaSuccess) return "";
    std::strncpy(name, prop.name, sizeof(name) - 1);
    return name;
}
#endif

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

bool geev_batched(const std::complex<double>* As, int n, int batch,
                  std::complex<double>* ws, std::complex<double>* vrs,
                  int streams) {
    if (n <= 0 || batch <= 0) return false;
    const int64_t N = n;
    const size_t nn = static_cast<size_t>(N) * N;
    const int S = std::max(1, std::min(streams, batch));

    // One persistent "slot" per concurrent stream: handle + params + its own
    // device buffers + pinned host staging, all created ONCE and reused across
    // the whole batch (this is what kills the per-call overhead).
    struct Slot {
        cudaStream_t stream = nullptr;
        cusolverDnHandle_t handle = nullptr;
        cusolverDnParams_t params = nullptr;
        cuDoubleComplex* dA = nullptr;
        cuDoubleComplex* dW = nullptr;
        cuDoubleComplex* dVR = nullptr;
        int* dInfo = nullptr;
        void* workDev = nullptr;
        void* workHost = nullptr;
        cuDoubleComplex* hA = nullptr;   // pinned staging (async copies)
        cuDoubleComplex* hW = nullptr;
        cuDoubleComplex* hVR = nullptr;
    };
    std::vector<Slot> slot(S);
    size_t lworkDev = 0, lworkHost = 0;
    bool ok = true;

    const cusolverEigMode_t novl = CUSOLVER_EIG_MODE_NOVECTOR;
    const cusolverEigMode_t vr_mode = CUSOLVER_EIG_MODE_VECTOR;

    auto cleanup = [&] {
        for (Slot& s : slot) {
            if (s.dA) cudaFree(s.dA);
            if (s.dW) cudaFree(s.dW);
            if (s.dVR) cudaFree(s.dVR);
            if (s.dInfo) cudaFree(s.dInfo);
            if (s.workDev) cudaFree(s.workDev);
            if (s.workHost) std::free(s.workHost);
            if (s.hA) cudaFreeHost(s.hA);
            if (s.hW) cudaFreeHost(s.hW);
            if (s.hVR) cudaFreeHost(s.hVR);
            if (s.params) cusolverDnDestroyParams(s.params);
            if (s.handle) cusolverDnDestroy(s.handle);
            if (s.stream) cudaStreamDestroy(s.stream);
        }
    };

    for (int i = 0; i < S; ++i) {
        Slot& s = slot[i];
        if (cudaStreamCreate(&s.stream) != cudaSuccess ||
            cusolverDnCreate(&s.handle) != CUSOLVER_STATUS_SUCCESS) {
            cleanup();
            return false;
        }
        cusolverDnSetStream(s.handle, s.stream);
        cusolverDnCreateParams(&s.params);
        if (cudaMalloc(reinterpret_cast<void**>(&s.dA), sizeof(cuDoubleComplex) * nn) != cudaSuccess ||
            cudaMalloc(reinterpret_cast<void**>(&s.dW), sizeof(cuDoubleComplex) * N) != cudaSuccess ||
            cudaMalloc(reinterpret_cast<void**>(&s.dVR), sizeof(cuDoubleComplex) * nn) != cudaSuccess ||
            cudaMalloc(reinterpret_cast<void**>(&s.dInfo), sizeof(int)) != cudaSuccess ||
            cudaMallocHost(reinterpret_cast<void**>(&s.hA), sizeof(cuDoubleComplex) * nn) != cudaSuccess ||
            cudaMallocHost(reinterpret_cast<void**>(&s.hW), sizeof(cuDoubleComplex) * N) != cudaSuccess ||
            cudaMallocHost(reinterpret_cast<void**>(&s.hVR), sizeof(cuDoubleComplex) * nn) != cudaSuccess) {
            cleanup();
            return false;
        }
    }

    // Workspace size is identical for every solve (same n) — query once.
    if (cusolverDnXgeev_bufferSize(
            slot[0].handle, slot[0].params, novl, vr_mode, N, CUDA_C_64F,
            slot[0].dA, N, CUDA_C_64F, slot[0].dW, CUDA_C_64F, nullptr, N,
            CUDA_C_64F, slot[0].dVR, N, CUDA_C_64F, &lworkDev, &lworkHost) !=
        CUSOLVER_STATUS_SUCCESS) {
        cleanup();
        return false;
    }
    for (int i = 0; i < S; ++i) {
        if (lworkDev && cudaMalloc(&slot[i].workDev, lworkDev) != cudaSuccess) { cleanup(); return false; }
        if (lworkHost) slot[i].workHost = std::malloc(lworkHost);
    }

    // Process the batch in waves of S concurrent solves. Each wave enqueues all
    // S streams' work (H2D, Xgeev, D2H) before synchronizing, so the device
    // overlaps independent solves.
    for (int base = 0; base < batch; base += S) {
        int cnt = std::min(S, batch - base);
        for (int i = 0; i < cnt; ++i) {
            Slot& s = slot[i];
            int b = base + i;
            std::memcpy(s.hA, As + static_cast<size_t>(b) * nn,
                        sizeof(cuDoubleComplex) * nn);
            cudaMemcpyAsync(s.dA, s.hA, sizeof(cuDoubleComplex) * nn,
                            cudaMemcpyHostToDevice, s.stream);
            cusolverDnXgeev(s.handle, s.params, novl, vr_mode, N, CUDA_C_64F,
                            s.dA, N, CUDA_C_64F, s.dW, CUDA_C_64F, nullptr, N,
                            CUDA_C_64F, s.dVR, N, CUDA_C_64F, s.workDev, lworkDev,
                            s.workHost, lworkHost, s.dInfo);
            cudaMemcpyAsync(s.hW, s.dW, sizeof(cuDoubleComplex) * N,
                            cudaMemcpyDeviceToHost, s.stream);
            cudaMemcpyAsync(s.hVR, s.dVR, sizeof(cuDoubleComplex) * nn,
                            cudaMemcpyDeviceToHost, s.stream);
        }
        for (int i = 0; i < cnt; ++i) {
            Slot& s = slot[i];
            int b = base + i;
            cudaStreamSynchronize(s.stream);
            std::memcpy(ws + static_cast<size_t>(b) * N, s.hW,
                        sizeof(cuDoubleComplex) * N);
            std::memcpy(vrs + static_cast<size_t>(b) * nn, s.hVR,
                        sizeof(cuDoubleComplex) * nn);
        }
    }

    cleanup();
    return ok;
}

} // namespace celeris::cuda
