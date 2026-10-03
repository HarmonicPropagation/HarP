from typing import Callable, Any, TYPE_CHECKING, cast
from components import HarpComponent
import numpy as np
from numpy.typing import NDArray
from scipy.sparse import diags
from scipy.sparse.linalg import spsolve
import scipy.linalg as la

type ComplexField = NDArray[np.complex128]
type MatrixBanded = NDArray[np.complex128]

if TYPE_CHECKING:
    from components import Main
    from pulse import Pulse
    from harmonic import Harmonic

class Solver(HarpComponent):  # type: ignore[misc]
    def __init__(self) -> None:
        self.Mold: MatrixBanded | None = None
        self.Mnew: MatrixBanded | None = None
        self.matrix: Any = None
        self.set_matrix: Callable[[float, Any, Any], tuple[ComplexField, ComplexField, ComplexField]] | None = None
        self.solve: Callable[["Pulse"], ComplexField] | None = None
        self.solveS: Callable[["Harmonic"], ComplexField] | None = None
        self.alphak: complex | None = None
        self.alphakdr2: complex | None = None
        self.main: "Main | None" = None
        self.nr: int | None = None

        self.Mold_times_E: ComplexField | None = None
        self.Mnew_q: MatrixBanded | None = None
        self.Mold_times_E_q: ComplexField | None = None
        super(Solver, self).__init__()

    def initialize(self, main: "Main") -> None:
        self.main = main
        assert self.main.nr is not None
        assert self.main.dz is not None
        assert self.main.dr is not None

        self.nr = self.main.nr

        if self.main.solver == "CrankNicolson":
            if self.main.boundary_condition == "absorbing":
                self.set_matrix = self.matrixCN_Abs
            elif self.main.boundary_condition == "reflective":
                self.set_matrix = self.matrixCN_Ref
            else:
                raise Exception(f"ERROR: boundary condition {self.main.boundary_condition} is not defined. Available boundary_condition : absorbing, reflective")
            
            self.solve  = self.solveCN
            self.solveS = self.solveCNS
            
            self.alphak = 1j * self.main.dz / 4. / self.main.dr**2
            self.alphakdr2 = 1j * self.main.dz / 4.

        self.Mold = np.zeros((3, self.nr), dtype=complex)
        self.Mnew = np.zeros((3, self.nr), dtype=complex)
        self.Mold_times_E = np.zeros(self.nr, dtype=complex)

        self.Mnew_q = np.zeros((3, self.main.nh * self.nr), dtype=complex)
        self.Mold_times_E_q = np.zeros(self.main.nh * self.nr, dtype=complex)

    def matrixCN_Abs(self, sign: float, k: Any, dk: Any) -> tuple[ComplexField, ComplexField, ComplexField]:
        assert self.main is not None
        assert self.alphak is not None
        assert self.alphakdr2 is not None
        assert self.main.dr is not None

        alpha = sign * self.alphak / k
        main_diag      = (1. + 2*alpha - 2*alpha * k * dk * self.main.dr**2)
        main_diag[-1]  = (1.+self.alphakdr2 * sign * ( k - 1j/self.main.r[-1]-2. * dk[-1]) )
        main_diag[0]  += 2. * alpha

        upper_diag     = np.zeros(self.main.nr - 1, dtype=complex)
        upper_diag[1:] = -alpha * (1 + self.main.dr / 2. / self.main.r[1:-1])
        upper_diag[0]  = -4 * alpha

        lower_diag     = -alpha * (1 - self.main.dr / 2. / self.main.r[1:])
        lower_diag[-1] = 0.+0j
        
        return upper_diag, main_diag, lower_diag
    
    def matrixCN_Ref(self, sign: float, k: Any, dk: Any) -> tuple[ComplexField, ComplexField, ComplexField]:
        assert self.main is not None
        assert self.alphak is not None
        assert self.main.dr is not None

        alpha = sign * self.alphak / k
        main_diag      = (1. + 2*alpha - 2*alpha * k * dk * self.main.dr**2)
        main_diag[-1]  = 1.
        main_diag[0]  += 2. * alpha

        upper_diag     = np.zeros(self.main.nr - 1, dtype=complex)
        upper_diag[1:] = -alpha * (1 + self.main.dr / 2. / self.main.r[1:-1])
        upper_diag[0]  = -4 * alpha

        lower_diag     = -alpha * (1 - self.main.dr / 2. / self.main.r[1:])
        lower_diag[-1] = 0.+0j
        
        return upper_diag, main_diag, lower_diag

    def solveCN(self, field: "Pulse") -> ComplexField:
        assert self.set_matrix is not None
        assert self.Mold is not None
        assert self.Mnew is not None
        assert self.Mold_times_E is not None
        assert field.E is not None

        self.Mold[0, 1:], self.Mold[1, :], self.Mold[2, :-1] = self.set_matrix(-1., field.kold, field.dkold)

        self.Mold_times_E[:] = 0.

        self.Mold_times_E[:]   += self.Mold[1, :] * field.E[:]             # Main diagonal contribution
        self.Mold_times_E[:-1] += self.Mold[0, 1:] * field.E[1:]           # Upper diagonal contribution
        self.Mold_times_E[1:]  += self.Mold[2, :-1] * field.E[:-1]

        self.Mnew[0, 1:], self.Mnew[1, :], self.Mnew[2, :-1] = self.set_matrix(1., field.k, field.dk)

        result = la.solve_banded((1, 1), self.Mnew, self.Mold_times_E)
        return cast(ComplexField, result)

    def solveCNS(self, field: "Harmonic") -> ComplexField:
        assert self.set_matrix is not None
        assert self.Mold_times_E_q is not None
        assert self.Mnew_q is not None
        assert self.nr is not None
        assert field.E is not None
        assert field.S is not None
        assert field.kold is not None
        assert field.dkold is not None
        assert field.k is not None
        assert field.dk is not None

        self.Mold_times_E_q[:] = 0.
        for q in range(field.nh):
            upper, main, lower = self.set_matrix(-1., field.kold[q], field.dkold[q])

            idx_start = q * self.nr
            idx_end = (q + 1) * self.nr

            self.Mold_times_E_q[idx_start:idx_end]   += main * field.E[idx_start:idx_end]              # Main diagonal
            self.Mold_times_E_q[idx_start:idx_end-1] += upper * field.E[idx_start+1:idx_end]           # Upper diagonal
            self.Mold_times_E_q[1+idx_start:idx_end] += lower * field.E[idx_start:idx_end-1]           # Lower diagonal

            self.Mnew_q[0, 1+idx_start:idx_end], self.Mnew_q[1, idx_start:idx_end], self.Mnew_q[2, idx_start:idx_end-1] = self.set_matrix(1., field.k[q], field.dk[q])

        RightEq = self.Mold_times_E_q + field.S

        result = la.solve_banded((1, 1), self.Mnew_q, RightEq)
        return cast(ComplexField, result)
