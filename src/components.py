from typing import Any, Callable, Type
import numpy as np
from numpy.typing import NDArray

type SpatialGrid = NDArray[np.float64]

def _add_metaclass(metaclass: Type[type]) -> Callable[[type], type]:
    def wrapper(cls: type) -> type:
        orig_vars = cls.__dict__.copy()
        slots = orig_vars.get('__slots__')
        if slots is not None:
            if isinstance(slots, str):
                slots = [slots]
            for slots_var in slots:
                orig_vars.pop(slots_var)
        orig_vars.pop('__dict__', None)
        orig_vars.pop('__weakref__', None)
        return metaclass(cls.__name__, cls.__bases__, orig_vars)
    return wrapper

class HarpComponentType(type):
    _list: list[Any]

    def __init__(cls, name: str, bases: tuple[type, ...], attrs: dict[str, Any]) -> None:
        cls._list = []
        super(HarpComponentType, cls).__init__(name, bases, attrs)

@_add_metaclass(HarpComponentType)
class HarpComponent(object):
    _list: list[Any]  # Provided by metaclass

    def __init__(self, **kwargs: Any) -> None:
        for k, v in kwargs.items():
            setattr(self, k, v)
        type(self)._list.append(self)

    @classmethod
    def get_instance(cls, mandatory: bool = True, singleton: bool = True) -> Any:
        if not cls._list:
            if mandatory:
                raise Exception(f"ERROR: {cls.__name__}() block is missing in input file")
            return None
        if singleton:
            if len(cls._list) > 1:
                raise Exception(f"ERROR: Multiple {cls.__name__}() blocks defined")
            return cls._list[0]
        else:
            return [clss for clss in cls._list]


class Main(HarpComponent):
    wavelength: float = 0.81 * 1e-6
    solver: str = "CrankNicolson"
    boundary_condition: str = "absorbing"
    dr: float | None = None
    dz: float | None = None
    nr: int | None = None
    nz: int | None = None
    Lr: float | None = None
    Lz: float | None = None
    dt: float | None = None
    nt: int | None = None
    simulation_time: float | None = None
    print_every: int | None = None
    verbose: bool = True
    anne_eq: bool = True
    predictor_corrector: bool = True
    approx_check: bool = False

    r: SpatialGrid
    z: SpatialGrid
    t: SpatialGrid
    timer: float
    T_r: float
    L_r: float
    E_r: float
    I_r: float
    N_r: float
    Polau_r: float
    Pol_r: float
    Kerau_r: float
    Ker_r: float
    D_r: float
    ev_r: float
    eps0c: float
    medium_timer: float
    pulse_timer: float
    diag_timer: float
    harmonic_timer: float
    nh: int

    def __init__(self, **kwargs: Any) -> None:

        super(Main, self).__init__(**kwargs)

        if self.solver != "CrankNicolson":
            raise Exception(f"ERROR: solver {self.solver} is not defined. Available solver : CrankNicolson")

        # TODO Some clarification about the cases may be needed.
        if self.nr is not None:
            assert self.Lr is not None
            self.dr = self.Lr / (self.nr - 1)
        if self.dr is not None and self.nr is None:
            assert self.Lr is not None
            self.nr = int(self.Lr / self.dr + 1)

        if self.nz is not None:
            assert self.Lz is not None
            self.dz = self.Lz / (self.nz - 1)
        if self.dz is not None and self.nz is None:
            assert self.Lz is not None
            self.nz = int(self.Lz / self.dz + 1)

        if self.nt is not None:
            assert self.simulation_time is not None
            self.dt = self.simulation_time / (self.nt - 1)
        if self.dt is not None and self.nt is None:
            assert self.simulation_time is not None
            self.nt = int(self.simulation_time / self.dt + 1)
        
        assert self.Lr is not None and self.nr is not None
        assert self.Lz is not None and self.nz is not None
        assert self.simulation_time is not None and self.nt is not None

        self.r = np.linspace(0, self.Lr, self.nr)
        self.z = np.linspace(0, self.Lz, self.nz)
        self.t = np.linspace(0., self.simulation_time, self.nt)

        if self.print_every is None:
            self.print_every = self.nt-1

        self.timer = 0.
        self.T_r = self.wavelength / (2 * np.pi * 299792458.)
        self.L_r = self.wavelength / (2 * np.pi)
        self.E_r = 3210701.09 / self.wavelength
        self.I_r = 1368167.2465547603 / self.wavelength**2
        self.N_r = 1114854216171690.4 / self.wavelength**2
        self.Polau_r = 2.076019128799422e-15 / self.wavelength**2
        self.Pol_r = 1.2591264605433468e+26 / self.wavelength**2
        self.Kerau_r = 8.09341907276034e-26 / self.wavelength**4
        self.Ker_r = 1.2979832996610822e39 / self.wavelength**4
        self.D_r_au = 3.3249184764493793e-10 / self.wavelength
        self.ev_r = 1.2398419843320026e-06 / self.wavelength
        self.eps0c = 0.0026544187279929563
        self.medium_timer = 0.0
        self.pulse_timer = 0.0
        self.diag_timer = 0.0
        self.harmonic_timer = 0.0
        self.nh = 0

    def printer(self, sentence: str) -> None:
        if self.verbose:
            print(sentence)
