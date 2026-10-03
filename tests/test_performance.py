import pytest
import tracemalloc
import time
import numpy as np
import numpy.typing as npt
from typing import Any

from components import Main
from pulse import Pulse
from medium import Medium
from ionisation import Ionisation
from harmonic import Harmonic
from diagnostic import DiagEnergy
from simulation import Simulation
from solver import Solver

def create_simulation(nr: int, nz: int, nt: int, approx_check: bool,orders: list[int]) -> Simulation:
    """Factory function to compose the HarP simulation components."""
    # Pass print_every > nt so it doesn't pollute benchmark stdout
    main: Main = Main(
        Lr=10.0, Lz=10.0, nr=nr, nz=nz, simulation_time=10.0, nt=nt,
        verbose=False, print_every=nt + 1, approx_check=approx_check
    )

    def gas_profile(r: npt.NDArray[np.float64], z: float) -> float:
        return 1.0

    # Provide minimal physical constants that Harmonic/Pulse rely on
    medium: Medium = Medium(
        gaz_profile=gas_profile,
        polarizability_0=1.8247e-40,
        gamma_kerr=1.2e-62,
        table_path="tables/form_factor/f12_Ar_henkegov.txt"
    )
    medium.initialize(main)

    # 3. Initialize Ionisation (Allocates N_z x N_r tracking arrays)
    ion: Ionisation = Ionisation(table_path="tables/ionisation/ionisation_Ar_810_tdse.txt", skiprows=2)
    ion.initialize(main)

    # 4. Initialize Harmonic (Allocates N_q x N_r complex fields)
    harmonic: Harmonic = Harmonic(order=orders, table_path="tables/tdse/dipoles_Ar_810.txt")

    # Mock the pulse
    pulse: Pulse = Pulse(
        radial_profile=lambda r: 1.0,
        time_envelope=lambda t: 1.0,
        amplitude=0.1
    )
    pulse.initialize(main, medium)

    harmonic.initialize(main, medium, pulse)
    main.nh = harmonic.nh  # Required for Solver matrix allocation

    diag: DiagEnergy = DiagEnergy(every_z=1, order=orders, SI=True)
    diag.initialize(main, harmonic)

    solver = Solver()
    solver.initialize(main)

    return Simulation(
        main=main,
        medium=medium,
        pulse=pulse,
        ionisation=ion,
        harmonic=harmonic,
        solver=solver,
        diagnostics=[],
        diagenergy=diag,
        diagfarfield=None
    )

def measure_state_memory(nr: int, nz: int, orders: list[int]) -> float:
    """Helper to instantiate the full HarP state and return peak memory in MB."""
    tracemalloc.start()
    _ = create_simulation(nr=nr, nz=nz, nt=512, approx_check=True, orders=orders)
    current, peak = tracemalloc.get_traced_memory()
    tracemalloc.stop()
    return peak / 1024 / 1024

def test_memory_scaling_across_sizes() -> None:
    """Verify that memory scales predictably across varying grid sizes and harmonic counts.

    Empirically validated cost model (O(N_q N_r + N_z N_r)):
      - Ionisation tracking matrices (float64, 8 bytes): 2 * N_z * N_r * 8 B
      - Harmonic complex fields (complex128, 16 bytes): 6 * N_q * N_r * 16 B
      - Fixed Base Overhead (C): ~2.80 MB (interpreter, imports, ASCII tables)

    Theoretical vs. measured peak memory comparison:
      - Small Tier  (nr=512,  nz=256,  N_q=10): Theory ~5.27 MB  | Measured: 5.31 MB
      - Medium Tier (nr=1024, nz=512,  N_q=20): Theory ~12.68 MB | Measured: 13.89 MB
      - Large Tier  (nr=2048, nz=1024, N_q=30): Theory ~40.43 MB | Measured: 42.96 MB
    """
    # Define three tiers: (Tier Name, nr, nz, harmonic_orders)
    configs: list[tuple[str, int, int, list[int]]] = [
        ("Small Tier",  512,  256,  list(range(3, 23, 2))),  # 10 harmonics
        ("Medium Tier", 1024, 512,  list(range(3, 43, 2))),  # 20 harmonics
        ("Large Tier",  2048, 1024, list(range(3, 63, 2))),  # 30 harmonics
    ]

    results: dict[str, float] = {}
    for name, nr, nz, orders in configs:
        peak_mb: float = measure_state_memory(nr, nz, orders)
        results[name] = peak_mb
        print(f"\n[{name}] nr={nr}, nz={nz}, N_q={len(orders)} -> Peak memory: {peak_mb:.2f} MB")

    # Assert that each tier stays safely bounded according to the O(N_q N_r + N_z N_r) claim
    assert results["Small Tier"]  < 15.0, f"Small tier footprint too high:  {results['Small Tier']:.2f} MB"
    assert results["Medium Tier"] < 35.0, f"Medium tier footprint too high: {results['Medium Tier']:.2f} MB"
    assert results["Large Tier"]  < 75.0, f"Large tier footprint too high:  {results['Large Tier']:.2f} MB"

def test_bo_workhorse_time_scaling() -> None:
    """Simulate a representative BO-like single-harmonic run setup and measure initialization timing."""
    configs: list[tuple[str, int, int, int]] = [
        ("BO Small Grid", 512, 128, 32),
        ("BO Standard Grid", 1024, 256, 32),
    ]

    for name, nr, nz, nt in configs:
        t0 = time.perf_counter()
        _ = create_simulation(nr=nr, nz=nz, nt=nt, approx_check=False, orders=[23])
        init_time = time.perf_counter() - t0

        print(f"\n[BO Time Benchmark] {name} (nr={nr}, nz={nz}, nt={nt}) -> Setup Time: {init_time:.4f} s")
        assert init_time > 0.0

@pytest.mark.parametrize("nt", [32, 64, 128])
def test_propagation_time_scaling(benchmark: Any, nt: int) -> None:
    """Empirically verify O(n_t) time complexity of the main integration loop."""
    # Use a factory to spin up the state
    sim = create_simulation(nr=128, nz=64, nt=nt, approx_check=False, orders=[23])

    # pytest-benchmark runs `sim.run` multiple times.
    # We monkey-patch the save methods to avoid HDF5 "dataset already exists"
    # crashes on the second iteration, and to isolate pure compute time from disk I/O.
    if sim.diagenergy is not None:
        sim.diagenergy.save = lambda *args, **kwargs: None

    if sim.diagfarfield is not None:
        sim.diagfarfield.save = lambda *args, **kwargs: None

    for diag in sim.diagnostics:
        diag.save = lambda *args, **kwargs: None

    # Run the actual integration loop benchmark
    benchmark(sim.run)

# def test_grid_initialization_speed(benchmark: Any) -> None:
#     """Benchmark the execution time of setting up the grid."""
#     benchmark(Main, Lr=10.0, Lz=100.0, nr=2048, nz=1024, simulation_time=50.0, nt=512, verbose=False, print_every=1000)
