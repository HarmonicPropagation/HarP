import pytest
import numpy as np
from typing import Callable
from components import Main
from harmonic import Harmonic
from diagnostic import DiagEnergy, Diagnostics
from medium import Medium
from pulse import Pulse

# Some minimal tests to start from
def test_main_grid_initialization() -> None:
    """Verify that specifying grid points correctly computes steps."""
    main: Main = Main(Lr=10.0, Lz=100.0, nr=512, nz=256, simulation_time=50.0, nt=128, verbose=False)
    assert main.dr == 10.0 / 511.0
    assert main.dz == 100.0 / 255.0
    assert len(main.r) == 512
    assert len(main.z) == 256

def test_harmonic_missing_parameters() -> None:
    """Ensure Harmonic raises an error if mandatory parameters are missing."""
    with pytest.raises(ValueError):
        Harmonic(order=[3, 5]) # Missing table_path


def test_diagnostic_initialization() -> None:
    """Verify that diagnostic components initialize with default types."""
    diag: DiagEnergy = DiagEnergy(every_z=2, order=[1, 3])
    assert diag.every_z == 2
    assert 1 in diag.order
    assert 3 in diag.order


def make_medium(main: Main) -> Medium:
    # gaz_profile is expected to already return the density normalised to
    # N_r (see doc/source/theory/units.rst); 1e-3 is representative of a
    # real HarP input (e.g. input.py's own n0 works out to ~2e-3).
    medium = Medium(
        gaz_profile=lambda r, z: 1e-3,
        polarizability_0=1.8247e-40,
        gamma_kerr=1.2e-62,
        table_path="tables/form_factor/f12_Ar_henkegov.txt",
    )
    medium.initialize(main)
    return medium


def gaussian_envelope(tau: float, t0: float = 100.0) -> Callable[[float], float]:
    return lambda t: np.exp(-((t - t0) / tau) ** 2)


#def test_svea_aborts_for_too_short_pulse() -> None:
#    """A pulse FWHM well under 2 optical cycles must abort (Physical Constraint #6)."""
#    main = Main(Lr=10.0, Lz=10.0, nr=16, nz=8, simulation_time=200.0, nt=2000, verbose=False,approx_check=True)
#    medium = make_medium(main)
#    pulse = Pulse(radial_profile=lambda r: np.ones_like(r), time_envelope=gaussian_envelope(tau=2.0))
#    with pytest.raises(ValueError):
#        pulse.initialize(main, medium)


def test_svea_warns_but_does_not_abort_for_marginal_pulse(capsys: pytest.CaptureFixture[str]) -> None:
    """A pulse FWHM between 2 and 5 optical cycles should warn, not abort."""
    main = Main(Lr=10.0, Lz=10.0, nr=16, nz=8, simulation_time=200.0, nt=2000, verbose=True,approx_check=True)
    medium = make_medium(main)
    pulse = Pulse(radial_profile=lambda r: np.ones_like(r), time_envelope=gaussian_envelope(tau=16.0))
    pulse.initialize(main, medium)  # must not raise
    assert "WARNING" in capsys.readouterr().out


def test_svea_silent_for_well_resolved_pulse(capsys: pytest.CaptureFixture[str]) -> None:
    """A pulse FWHM well above 5 optical cycles should be silent."""
    main = Main(Lr=10.0, Lz=10.0, nr=16, nz=8, simulation_time=200.0, nt=2000, verbose=True,approx_check=True)
    medium = make_medium(main)
    pulse = Pulse(radial_profile=lambda r: np.ones_like(r), time_envelope=gaussian_envelope(tau=80.0))
    pulse.initialize(main, medium)
    assert "WARNING" not in capsys.readouterr().out


def test_svea_silent_for_flat_envelope_within_window(capsys: pytest.CaptureFixture[str]) -> None:
    """A flat/CW-like envelope never drops below half its peak within the
    simulated time window - its true FWHM is unmeasurable (at least as long
    as the window), so it must not be flagged as an SVEA violation."""
    main = Main(Lr=10.0, Lz=10.0, nr=16, nz=8, simulation_time=10.0, nt=128, verbose=True,approx_check=True)
    medium = make_medium(main)
    pulse = Pulse(radial_profile=lambda r: np.ones_like(r), time_envelope=lambda t: 1.0)
    pulse.initialize(main, medium)
    assert "WARNING" not in capsys.readouterr().out


#def test_svea_rejects_multi_peak_envelope() -> None:
#    """A time_envelope with two separate excursions above half its peak
#    intensity (e.g. a free-form pulse shape from an optimiser) isn't a
#    single-peaked pulse - "FWHM" isn't meaningful for it, so it must raise
#    rather than silently measure a bogus value."""
#    main = Main(Lr=10.0, Lz=10.0, nr=16, nz=8, simulation_time=100.0, nt=1000, verbose=False,approx_check=True)
#    medium = make_medium(main)
#    two_peaks = lambda t: np.exp(-((t - 30.0) / 5.0) ** 2) + np.exp(-((t - 70.0) / 5.0) ** 2)
#    pulse = Pulse(radial_profile=lambda r: np.ones_like(r), time_envelope=two_peaks)
#    with pytest.raises(ValueError):
#        pulse.initialize(main, medium)


def test_low_density_warns_for_high_density(capsys: pytest.CaptureFixture[str]) -> None:
    """A gas density approaching N_r should warn (Physical Constraint #8)."""
    main = Main(Lr=10.0, Lz=10.0, nr=16, nz=8, simulation_time=10.0, nt=4, verbose=True, approx_check=True)
    medium = Medium(
        gaz_profile=lambda r, z: 0.5,
        polarizability_0=1.8247e-40,
        gamma_kerr=1.2e-62,
        table_path="tables/form_factor/f12_Ar_henkegov.txt",
    )
    medium.initialize(main)
    assert "WARNING" in capsys.readouterr().out


def test_low_density_silent_for_realistic_density(capsys: pytest.CaptureFixture[str]) -> None:
    """A gas density well below N_r (typical of real HarP inputs) should be silent."""
    main = Main(Lr=10.0, Lz=10.0, nr=16, nz=8, simulation_time=10.0, nt=4, verbose=True,approx_check=True)
    make_medium(main)
    assert "WARNING" not in capsys.readouterr().out
