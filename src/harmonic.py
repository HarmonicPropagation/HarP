from typing import Any, Callable, TYPE_CHECKING, cast
from components import HarpComponent
import numpy as np
from numpy.typing import NDArray
from scipy.interpolate import interp1d,RegularGridInterpolator
import time

# Semantic type synonyms
type ComplexField = NDArray[np.complex128]
type RealField = NDArray[np.float64]

if TYPE_CHECKING:
    from components import Main
    from medium import Medium
    from pulse import Pulse
    from simulation import Simulation

class Harmonic(HarpComponent):  # type: ignore[misc]
    order: Any = None
    table_path: str = ""
    intensity_unit_Wcm2: float = 1e14
    statistical_factor: float = 2.
    smoothing_function: Callable[[RealField], RealField] | None = None

    def __init__(self, **kwargs: Any) -> None:
        super(Harmonic, self).__init__(**kwargs)
        required = ["order","table_path"]
        for key in required:
            if key not in kwargs:
                raise ValueError(f"Missing required parameter '{key}' in Harmonic")
        self.E: ComplexField | None = None
        self.S: ComplexField | None = None
        assert self.order is not None
        self.nh: int = len(self.order)
        self.k: ComplexField | None = None
        self.dk: ComplexField | None = None
        self.kold: ComplexField | None = None
        self.dkold: ComplexField | None = None
        self.compute_k: Callable[[int], Any] | None = None
        self.compute_dk: Callable[[int], Any] | None = None
        self.Phase: ComplexField | None = None
        self.chi_at: ComplexField | None = None
        self.medium: "Medium | None" = None
        self.pulse: "Pulse | None" = None
        self.timer: float = 0.
        self.E_Tables: list[RealField] = []
        self.Dipoleabs_Tables: list[RealField] = []
        self.Dipolearg_Tables: list[RealField] = []
        self.func_dipoleabs: Any = None
        self.func_dipolearg: Any = None
        self.compute_dipoleabs: list[Callable[[Any], Any]] = []
        self.compute_dipolearg: list[Callable[[Any], Any]] = []

        if self.order is not None:
            self.order.sort()

    def initialize(self, main: "Main", medium: "Medium", pulse: "Pulse") -> None:
        self.medium = medium
        self.pulse  = pulse
        self.timer  = 0.
        if main.anne_eq:
            self.compute_k  = self.compute_k_Anne
            self.compute_dk = self.compute_dk_Anne
        else:
            self.compute_k  = self.compute_k_AnneBar
            self.compute_dk = self.compute_dk_AnneBar

        self.E     = np.zeros(self.nh*main.nr, dtype = complex)

        self.kold   = np.ones(self.nh,dtype = complex)
        self.dkold  = np.zeros((self.nh,main.nr),dtype = complex)

        self.k      = np.ones(self.nh,dtype = complex)
        self.dk     = np.zeros((self.nh,main.nr),dtype = complex)

        self.Pqold  = np.zeros((self.nh,main.nr), dtype=complex)
        self.Pq     = np.zeros((self.nh,main.nr), dtype=complex)
        self.S      = np.zeros(self.nh*main.nr, dtype=complex)

        self.chi_at = np.zeros(self.nh, dtype=complex)
        for q in range (self.nh):
            self.compute_chi_at(q)
            self.k[q]  = self.compute_k(q)
            self.dk[q] = self.compute_dk(q)

        self.E_Tables = []
        self.Dipoleabs_Tables = []
        self.Dipolearg_Tables = []
        self.compute_dipoleabs = []
        self.compute_dipolearg = []
        for i in range (self.nh):
            self.load_table(main,int(self.order[i]))
            self.compute_dipoleabs.append(
                cast(Callable[[Any], Any], lambda absE, x=self.E_Tables[i][:], y=self.Dipoleabs_Tables[i][:]: np.interp(absE, x, y, left=0.0, right=self.Dipoleabs_Tables[i][-1]))
            )
            self.compute_dipolearg.append(
                cast(Callable[[Any], Any], lambda absE, x=self.E_Tables[i][:], y=self.Dipolearg_Tables[i][:]: np.interp(absE, x, y, left=0.0, right=self.Dipolearg_Tables[i][-1]))
            )

    def load_table(self, main: "Main", q: int) -> None:
        with open(self.table_path, 'r') as f:
            found = False
            Data = []
            for line in f:
                if "harmonic"==line.split()[0] and q==int(line.split()[1]):
                    found = True
                    continue
                if found:
                    if "harmonic" in line:
                        break
                    row = [float(val) for val in line.strip().split()]
                    Data.append(row)
            if not found:
                raise ValueError(f"Harmonic order {q} not found in file. ")
            self.E_Tables.append(np.sqrt(np.array(Data)[:,0]*self.intensity_unit_Wcm2/main.I_r))
            if self.smoothing_function is None:
                self.Dipoleabs_Tables.append(np.sqrt(np.array(Data)[:,1])*main.D_r_au * self.statistical_factor)
            else:
                self.Dipoleabs_Tables.append(self.smoothing_function(np.sqrt(np.array(Data)[:,1])*main.D_r_au*self.statistical_factor))
            self.Dipolearg_Tables.append(np.unwrap(np.array(Data)[:,2]))

    def compute_chi_at(self, q: int) -> None:
        assert self.medium is not None
        assert self.chi_at is not None
        self.chi_at[q]=(-self.medium.compute_f1(self.order[q]) + 1j*self.medium.compute_f2(self.order[q]))/self.order[q]

    def compute_k_Anne(self, q: int) -> Any:
        assert self.medium is not None
        assert self.chi_at is not None
        return self.order[q] + 1./2.*self.medium.n_gaz * self.chi_at[q]

    def compute_dk_Anne(self, q: int) -> Any:
        assert self.medium is not None
        assert self.pulse is not None
        assert self.chi_at is not None
        return 1/2.*( -self.chi_at[q] - 1./self.order[q] ) * self.medium.n_electron\
               +1./2.*self.medium.n_gaz * (self.chi_at[q]-self.order[q]*self.pulse.chi_at)

    def compute_k_AnneBar(self, q: int) -> Any:
        return self.order[q]+0.j

    def compute_dk_AnneBar(self, q: int) -> Any:
        assert self.medium is not None
        assert self.chi_at is not None
        return 1/2.*( self.chi_at[q] *  self.medium.n_atomic[:] - self.medium.n_electron[:]/self.order[q] )

    def reinitialize(self) -> None:
        assert self.E is not None
        assert self.k is not None
        assert self.dk is not None
        assert self.pulse is not None
        assert self.pulse.E is not None
        assert self.compute_k is not None
        assert self.compute_dk is not None
        assert self.medium is not None

        self.E[:] = 0.+0.j
        absE1 = np.abs(self.pulse.E[:])
        phiE1 = np.angle(self.pulse.E[:])
        for q in range (self.nh):
            self.compute_chi_at(q)
            self.k[q]  = self.compute_k(q)
            self.dk[q,:] = self.compute_dk(q)
            self.Pq[q,:] = - 2.*self.order[q]**2 * self.medium.n_atomic*self.compute_dipoleabs[q](absE1)*np.exp(1j*(self.compute_dipolearg[q](absE1)+self.order[q]*phiE1))

    def update(self, sim: "Simulation") -> None:
        assert self.pulse is not None
        assert self.pulse.E is not None
        assert self.k is not None
        assert self.kold is not None
        assert self.dk is not None
        assert self.dkold is not None
        assert self.E is not None
        assert self.compute_k is not None
        assert self.compute_dk is not None
        assert self.medium is not None
        assert self.S is not None

        self.timer  =time.time()
        absE1 = np.abs(self.pulse.E[:])
        phiE1 = np.angle(self.pulse.E[:])
        for q in range (self.nh):
            self.kold[q]     = self.k[q]
            self.dkold[q,:]  = self.dk[q,:]
            self.k[q]        = self.compute_k(q)
            self.dk[q,:]     = self.compute_dk(q)

            self.Pqold[q,:]  = self.Pq[q,:]
            self.Pq[q,:]     = - 2.*self.order[q]**2 * self.medium.n_atomic*self.compute_dipoleabs[q](absE1)*np.exp(1j*(self.order[q]*phiE1+self.compute_dipolearg[q](absE1)))
            self.S[q*sim.main.nr:(q+1)*sim.main.nr]      = - (self.Pqold[q,:]/self.kold[q]+self.Pq[q,:]/self.k[q])*sim.solver.alphakdr2
        self.E[:]      = sim.solver.solveS(self)
        sim.main.harmonic_timer += time.time()-self.timer
