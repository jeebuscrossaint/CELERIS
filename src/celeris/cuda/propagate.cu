// CUDA kernel for Rayleigh-Sommerfeld focal-plane propagation. Compiled by nvcc
// (device code). One thread per output pixel; pillars are streamed through
// shared memory in tiles (the classic N-body tiling) so each pillar coordinate
// is read from global memory once per block instead of once per thread.
#include "celeris/cuda/propagate.hpp"

#include <cuda_runtime.h>
#include <math.h>

#include <vector>

namespace celeris::cuda {

namespace {
constexpr int kBlock = 256;  // threads per block == shared-memory tile size

__global__ void psf_kernel(const float* __restrict__ px,
                            const float* __restrict__ py,
                            const float* __restrict__ tr,
                            const float* __restrict__ ti, int npil, float cx,
                            float cy, float z, float k, int n, float W,
                            float step, float* __restrict__ out) {
    __shared__ float sx[kBlock];
    __shared__ float sy[kBlock];
    __shared__ float str[kBlock];
    __shared__ float sti[kBlock];

    const int idx = blockIdx.x * blockDim.x + threadIdx.x;
    const bool active = idx < n * n;
    const int i = active ? idx % n : 0;
    const int j = active ? idx / n : 0;
    const float fx = cx - W + i * step;
    const float fy = cy - W + j * step;

    float er = 0.0f, ei = 0.0f;

    for (int base = 0; base < npil; base += kBlock) {
        const int t = threadIdx.x;
        if (base + t < npil) {
            sx[t] = px[base + t];
            sy[t] = py[base + t];
            str[t] = tr[base + t];
            sti[t] = ti[base + t];
        }
        __syncthreads();

        const int lim = min(kBlock, npil - base);
        if (active) {
            for (int p = 0; p < lim; ++p) {
                const float dx = fx - sx[p];
                const float dy = fy - sy[p];
                const float R = sqrtf(dx * dx + dy * dy + z * z);
                float s, c;
                sincosf(k * R, &s, &c);
                const float inv = 1.0f / R;
                // (tr + i ti)(c + i s) / R
                er += (str[p] * c - sti[p] * s) * inv;
                ei += (str[p] * s + sti[p] * c) * inv;
            }
        }
        __syncthreads();
    }

    if (active) out[idx] = er * er + ei * ei;
}

__global__ void zx_kernel(const float* __restrict__ px,
                          const float* __restrict__ py,
                          const float* __restrict__ tr,
                          const float* __restrict__ ti, int npil, float xmin,
                          float xstep, int nx, float zmin, float zstep, int nz,
                          float k, float* __restrict__ out) {
    __shared__ float sx[kBlock];
    __shared__ float sy[kBlock];
    __shared__ float str[kBlock];
    __shared__ float sti[kBlock];

    const int idx = blockIdx.x * blockDim.x + threadIdx.x;
    const bool active = idx < nx * nz;
    const int ix = active ? idx % nx : 0;
    const int iz = active ? idx / nx : 0;
    const float fx = xmin + ix * xstep;
    const float z = zmin + iz * zstep;

    float er = 0.0f, ei = 0.0f;
    for (int base = 0; base < npil; base += kBlock) {
        const int t = threadIdx.x;
        if (base + t < npil) {
            sx[t] = px[base + t];
            sy[t] = py[base + t];
            str[t] = tr[base + t];
            sti[t] = ti[base + t];
        }
        __syncthreads();
        const int lim = min(kBlock, npil - base);
        if (active) {
            for (int p = 0; p < lim; ++p) {
                const float dx = fx - sx[p];
                const float dy = sy[p];  // observation point at y = 0
                const float R = sqrtf(dx * dx + dy * dy + z * z);
                float s, c;
                sincosf(k * R, &s, &c);
                const float inv = 1.0f / R;
                er += (str[p] * c - sti[p] * s) * inv;
                ei += (str[p] * s + sti[p] * c) * inv;
            }
        }
        __syncthreads();
    }
    if (active) out[idx] = er * er + ei * ei;
}
}  // namespace

bool propagate_psf(const double* px, const double* py,
                   const std::complex<double>* pt, int npil, double cx,
                   double cy, double z, double k, int n, double half_window,
                   double* out) {
    int devs = 0;
    if (cudaGetDeviceCount(&devs) != cudaSuccess || devs == 0) return false;
    if (npil <= 0 || n <= 0) return false;

    const size_t P = static_cast<size_t>(npil);
    const size_t NN = static_cast<size_t>(n) * n;
    const float step = static_cast<float>(2.0 * half_window / (n - 1));

    // Host float staging (kernel runs in single precision).
    std::vector<float> hx(P), hy(P), htr(P), hti(P);
    for (size_t p = 0; p < P; ++p) {
        hx[p] = static_cast<float>(px[p]);
        hy[p] = static_cast<float>(py[p]);
        htr[p] = static_cast<float>(pt[p].real());
        hti[p] = static_cast<float>(pt[p].imag());
    }

    float *dx = nullptr, *dy = nullptr, *dtr = nullptr, *dti = nullptr, *dout = nullptr;
    bool ok = false;
    auto cleanup = [&] {
        if (dx) cudaFree(dx);
        if (dy) cudaFree(dy);
        if (dtr) cudaFree(dtr);
        if (dti) cudaFree(dti);
        if (dout) cudaFree(dout);
    };

    if (cudaMalloc(&dx, P * sizeof(float)) != cudaSuccess ||
        cudaMalloc(&dy, P * sizeof(float)) != cudaSuccess ||
        cudaMalloc(&dtr, P * sizeof(float)) != cudaSuccess ||
        cudaMalloc(&dti, P * sizeof(float)) != cudaSuccess ||
        cudaMalloc(&dout, NN * sizeof(float)) != cudaSuccess) {
        cleanup();
        return false;
    }

    cudaMemcpy(dx, hx.data(), P * sizeof(float), cudaMemcpyHostToDevice);
    cudaMemcpy(dy, hy.data(), P * sizeof(float), cudaMemcpyHostToDevice);
    cudaMemcpy(dtr, htr.data(), P * sizeof(float), cudaMemcpyHostToDevice);
    cudaMemcpy(dti, hti.data(), P * sizeof(float), cudaMemcpyHostToDevice);

    const int blocks = static_cast<int>((NN + kBlock - 1) / kBlock);
    psf_kernel<<<blocks, kBlock>>>(dx, dy, dtr, dti, npil,
                                   static_cast<float>(cx), static_cast<float>(cy),
                                   static_cast<float>(z), static_cast<float>(k), n,
                                   static_cast<float>(half_window), step, dout);

    if (cudaGetLastError() == cudaSuccess && cudaDeviceSynchronize() == cudaSuccess) {
        std::vector<float> hout(NN);
        cudaMemcpy(hout.data(), dout, NN * sizeof(float), cudaMemcpyDeviceToHost);
        for (size_t i = 0; i < NN; ++i) out[i] = static_cast<double>(hout[i]);
        ok = true;
    }

    cleanup();
    return ok;
}

bool propagate_zx(const double* px, const double* py,
                  const std::complex<double>* pt, int npil, double xmin,
                  double xmax, int nx, double zmin, double zmax, int nz, double k,
                  double* out) {
    int devs = 0;
    if (cudaGetDeviceCount(&devs) != cudaSuccess || devs == 0) return false;
    if (npil <= 0 || nx <= 1 || nz <= 1) return false;

    const size_t P = static_cast<size_t>(npil);
    const size_t NN = static_cast<size_t>(nx) * nz;
    const float xstep = static_cast<float>((xmax - xmin) / (nx - 1));
    const float zstep = static_cast<float>((zmax - zmin) / (nz - 1));

    std::vector<float> hx(P), hy(P), htr(P), hti(P);
    for (size_t p = 0; p < P; ++p) {
        hx[p] = static_cast<float>(px[p]);
        hy[p] = static_cast<float>(py[p]);
        htr[p] = static_cast<float>(pt[p].real());
        hti[p] = static_cast<float>(pt[p].imag());
    }

    float *dx = nullptr, *dy = nullptr, *dtr = nullptr, *dti = nullptr, *dout = nullptr;
    bool ok = false;
    auto cleanup = [&] {
        if (dx) cudaFree(dx);
        if (dy) cudaFree(dy);
        if (dtr) cudaFree(dtr);
        if (dti) cudaFree(dti);
        if (dout) cudaFree(dout);
    };

    if (cudaMalloc(&dx, P * sizeof(float)) != cudaSuccess ||
        cudaMalloc(&dy, P * sizeof(float)) != cudaSuccess ||
        cudaMalloc(&dtr, P * sizeof(float)) != cudaSuccess ||
        cudaMalloc(&dti, P * sizeof(float)) != cudaSuccess ||
        cudaMalloc(&dout, NN * sizeof(float)) != cudaSuccess) {
        cleanup();
        return false;
    }
    cudaMemcpy(dx, hx.data(), P * sizeof(float), cudaMemcpyHostToDevice);
    cudaMemcpy(dy, hy.data(), P * sizeof(float), cudaMemcpyHostToDevice);
    cudaMemcpy(dtr, htr.data(), P * sizeof(float), cudaMemcpyHostToDevice);
    cudaMemcpy(dti, hti.data(), P * sizeof(float), cudaMemcpyHostToDevice);

    const int blocks = static_cast<int>((NN + kBlock - 1) / kBlock);
    zx_kernel<<<blocks, kBlock>>>(dx, dy, dtr, dti, npil,
                                  static_cast<float>(xmin), xstep, nx,
                                  static_cast<float>(zmin), zstep, nz,
                                  static_cast<float>(k), dout);

    if (cudaGetLastError() == cudaSuccess && cudaDeviceSynchronize() == cudaSuccess) {
        std::vector<float> hout(NN);
        cudaMemcpy(hout.data(), dout, NN * sizeof(float), cudaMemcpyDeviceToHost);
        for (size_t i = 0; i < NN; ++i) out[i] = static_cast<double>(hout[i]);
        ok = true;
    }
    cleanup();
    return ok;
}

}  // namespace celeris::cuda
