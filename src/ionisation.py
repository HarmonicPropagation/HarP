from typing import Any
from components import HarpComponent
import h5py
import numpy as np
from numpy.typing import NDArray
from scipy.interpolate import interp1d

type RateField = NDArray[np.float64]
type TableData = NDArray[np.float64]

class Ionisation(HarpComponent):  # type: ignore[misc]
    table_path: str = ""
    skiprows: int = 0
    intensity_unit_Wcm2: float = 1.
    rate_unit_SI: float = 1.

    def __init__(self, **kwargs: Any) -> None:
        super(Ionisation, self).__init__(**kwargs)

        self.E_Table: TableData | None = None
        self.Rate_Table: TableData | None = None

        self.compute_rate: interp1d | None = None
        self.local_rate: RateField | None = None
        self.integrated_rate: RateField | None = None

    def initialize(self, main: Any) -> None:
        self.local_rate = np.zeros((main.nz, main.nr))
        self.integrated_rate = np.zeros((main.nz, main.nr))
        self.load_table(main)

        self.compute_rate = interp1d(self.E_Table, self.Rate_Table, kind='linear', bounds_error=False, fill_value=0.0)

    def load_table(self, main: Any) -> None:
        try:
            data = np.loadtxt(self.table_path, skiprows=self.skiprows)
        except Exception:
            raise Exception("   ERROR in the ionisation table path: file not found    ")

        self.E_Table    = np.sqrt(data[:, 0] * self.intensity_unit_Wcm2 / main.I_r)
        self.Rate_Table = data[:, 1] * self.rate_unit_SI * main.T_r

    def update(self, sim: Any, iz: int) -> None:
        assert self.local_rate is not None
        assert self.compute_rate is not None
        self.local_rate[iz, :] = self.compute_rate(np.abs(sim.pulse.E))[:]

    def compute_integrated_rate(self, sim: Any) -> None:
        assert self.integrated_rate is not None
        assert self.local_rate is not None
        self.integrated_rate[:, :] += self.local_rate[:, :] * sim.main.dt
