import pytest
import numpy as np
from simulation import Simulation
from factory import SimulationBuilder
from components import Main
from pulse import Pulse
from medium import Medium
from harmonic import Harmonic
from diagnostic import DiagEnergy

# Minimal tests just to see if testing works

def test_cutoff_energy_scaling() -> None:
    """Verify the analytical cutoff scaling formula used in documentation. No HarP run."""
    # Parameters from application example
    I0:         float = 2.5e14 # W/cm^2
    wavelength: float = 810e-9 # m
    omega:      float = 2 * np.pi * 299792458 / wavelength
    
    # Ponderomotive potential U_p (in atomic units scaled to eV)
    # 3.17 * (e^2 E0^2 / 4 m omega^2)
    E0_au: float = np.sqrt(I0 / 3.509445e16) # rough conversion to a.u.
    Up_eV: float = 9.33 * (I0 / 1e14) * (wavelength * 1e6)**2
    
    Ip_Ar: float = 15.76 # eV
    expected_cutoff: float = Ip_Ar + 3.17 * Up_eV
    
    assert expected_cutoff > 60.0 and expected_cutoff < 70.0
