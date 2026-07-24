#!/usr/bin/env python3
"""Head-to-head cross-check: CELERIS vs Stanford S4 on an identical grating.

This reproduces the external-solver validation reported in the paper (§ validation):
CELERIS's from-scratch 1D RCWA solver is run against S4 (Liu & Fan, Comput. Phys.
Commun. 183, 2233, 2012) -- the field-standard Fourier-modal-method reference -- on the
SAME freestanding binary grating, for both polarizations, at increasing Fourier order M.

Requirements (both optional heavy deps, not needed for CELERIS itself):
    pip install .            # the `celeris` Python module (from the repo root)
    # S4: build a fork (e.g. github.com/phoebe-p/S4) with BLAS/LAPACK/FFTW/CHOLMOD/Boost:
    #   make S4_pyext BOOST_PREFIX=/usr BLAS_LIB=-lopenblas LAPACK_LIB=-lopenblas \
    #                 FFTW3_LIB=-lfftw3 CHOLMOD_LIB=-lcholmod \
    #                 CHOLMOD_INC=-I/usr/include/suitesparse PTHREAD_LIB=-lpthread

Representative result (RTX/GCC Linux, this repo):
    TE  M=40   CELERIS 0.933334   S4 0.933334   |diff| 7e-08   (== grcwa 0.93333)
    TM  M=40   CELERIS 0.962949   S4 0.962606   |diff| 3.4e-04  (TM converges slower in BOTH)
"""
import celeris as c
import S4

# Freestanding n=1.5 binary grating, normal incidence. Same case CELERIS locks
# against grcwa in selftest [5] (0th-order TE transmittance -> 0.93333).
LAMBDA_PERIOD, FILL, THICKNESS, WAVELENGTH = 0.3, 0.5, 0.5, 0.5
GLASS = c.Material.constant(complex(1.5, 0.0))
AIR = c.materials.air()


def celeris_t0(M, pol):
    g = c.BinaryGrating1D(GLASS, AIR, LAMBDA_PERIOD, FILL, THICKNESS)
    r = c.solve_rcwa_1d(AIR, g, AIR, WAVELENGTH, 0.0, M, pol)
    return r.de_t[len(r.orders) // 2]          # 0th order sits at the middle


def s4_t0(M, s_amp, p_amp):
    S = S4.New(Lattice=LAMBDA_PERIOD, NumBasis=M)
    S.SetMaterial(Name='Vac', Epsilon=1.0)
    S.SetMaterial(Name='Glass', Epsilon=2.25)   # n = 1.5
    S.AddLayer(Name='top', Thickness=0, Material='Vac')
    S.AddLayer(Name='grating', Thickness=THICKNESS, Material='Vac')
    S.SetRegionRectangle(Layer='grating', Material='Glass', Center=(0, 0),
                         Angle=0, Halfwidths=(FILL * LAMBDA_PERIOD / 2, 0))
    S.AddLayer(Name='bot', Thickness=0, Material='Vac')
    S.SetFrequency(1.0 / WAVELENGTH)
    S.SetExcitationPlanewave(IncidenceAngles=(0, 0), sAmplitude=s_amp, pAmplitude=p_amp)
    ft, _ = S.GetPowerFlux(Layer='top', zOffset=0)
    fb, _ = S.GetPowerFlux(Layer='bot', zOffset=0)
    return (fb / ft).real


if __name__ == "__main__":
    print("Freestanding n=1.5 grating: Lambda=0.3, fill=0.5, thickness=0.5, "
          "lambda=0.5, normal incidence")
    print(f"{'pol':>4} {'M':>4} {'CELERIS T0':>12} {'S4 T0':>12} {'|diff|':>10}")
    for pol_name, (s_amp, p_amp) in [('TE', (1, 0)), ('TM', (0, 1))]:
        pol = c.Pol.TE if pol_name == 'TE' else c.Pol.TM
        for M in (10, 20, 40):
            a, b = celeris_t0(M, pol), s4_t0(M, s_amp, p_amp)
            print(f"{pol_name:>4} {M:>4} {a:>12.6f} {b:>12.6f} {abs(a - b):>10.2e}")
