from typing import Any, Callable, TYPE_CHECKING
from components import HarpComponent
import h5py
import numpy as np
from numpy.typing import NDArray
import time
from scipy.special import j0,jn_zeros

type RealField = NDArray[np.float64]

if TYPE_CHECKING:
    from components import Main
    from harmonic import Harmonic

class Diagnostics(HarpComponent):  # type: ignore[misc]
    diag_path: str = ""
    deposited_quantity: str = ""
    add_suffix: str  = ""
    every_z: int = 1
    every_r: int = 1
    every_t: int = 1
    save_every: int | None = None
    order: list[int]   = [1]
    SI: bool      = True

    def __init__(self, **kwargs: Any) -> None:
        super(Diagnostics, self).__init__(**kwargs)
        self.filename = self.diag_path+"Diag"+self.deposited_quantity+self.add_suffix+".h5"
        self.nzdiag: int | None = None
        self.nrdiag: int | None = None
        self.ntdiag: int | None = None
        self.dummy_r: RealField | None = None
        self.order  = list(self.order[:])
        if 1 not in self.order:
            self.order.append(1)
        self.order.sort()
        self.quantity: dict[int, list[Any]] = {o: [] for o in self.order}
        self.saved_t_count = 0
        self.type: str | None = ""
        self.unit: Any = 1.

    def initialize(self, main: "Main", harmonic: "Harmonic | None") -> None:
        self.nrdiag = len(main.r[::self.every_r])
        self.nzdiag = len(main.z[::self.every_z])
        self.ntdiag = len(main.t[::self.every_t])
        self.dummy_r = np.ones(main.nr)
        types = {
                "Field": "complex128",
                "RealField": "complex128",
                "Phase": "complex128",
                "S": "complex128",
                "LocalRate": "float64",
                "IntegratedRate": "float64",
                "ElectronDensity": "float64",
                "AtomicDensity": "float64",
                "GazDensity": "float64",
                "k": "complex128",
                "dk": "complex128",
        }
        self.type = types.get(self.deposited_quantity)
        units = {
                "Field": main.E_r,
                "RealField": main.E_r,
                "Phase": 1.,
                "S": main.E_r/main.L_r**2,
                "LocalRate": 1./main.T_r,
                "IntegratedRate": 1.,
                "ElectronDensity": main.N_r,
                "AtomicDensity": main.N_r,
                "GazDensity": main.N_r,
                "k": 1./main.L_r,
                "dk": 1./main.L_r,
        }
        if self.SI:
            self.unit = units.get(self.deposited_quantity)

        if self.deposited_quantity not in list(types.keys()):
            raise Exception(f"""   ERROR: {self.deposited_quantity} is not a valid Diagnostic quantity. Available quantities are:
                {list(types.keys())}  """)
        for o in self.order:
            if o != 1:
                assert harmonic is not None
                if o not in harmonic.order:
                    raise Exception(
                        f"""ERROR: {o} is not a valid harmonic order. Available orders are:
                        {1} for the fundamental and those in harmonic of the input file"""
                    )

        with h5py.File(self.filename, 'w') as f:
            if self.SI:
                f.create_dataset("r", data=main.r[::self.every_r]*main.L_r)
                f.create_dataset("z", data=main.z[::self.every_z]*main.L_r)
                f.create_dataset("t", data=main.t[::self.every_t]*main.T_r)
            else:
                f.create_dataset("r", data=main.r[::self.every_r])
                f.create_dataset("z", data=main.z[::self.every_z])
                f.create_dataset("t", data=main.t[::self.every_t])

    def apply(self, it: int, iz: int, sim: Any) -> None:
        if not (it % self.every_t == 0 and iz % self.every_z == 0):
            return

        for io, o in enumerate(self.order):

            data_map: dict[str, Callable[[], Any]] = {
                "Field":           lambda: sim.pulse.E[::self.every_r] if io == 0 else sim.harmonic.E[(io-1)*sim.main.nr:(io)*sim.main.nr:self.every_r],
                "RealField":       lambda: sim.pulse.E[::self.every_r]*np.exp(1j*o*sim.pulse.Phase[iz]) if io == 0 else sim.harmonic.E[(io-1)*sim.main.nr:(io)*sim.main.nr:self.every_r]*np.exp(1j*o*sim.pulse.Phase[iz]),
                "Phase":           lambda: (sim.pulse.Phase[iz]*self.dummy_r)[::self.every_r],
                "S":               lambda: np.zeros(self.nrdiag,dtype=complex) if io==0 else sim.harmonic.Pq[io-1, ::self.every_r],
                "LocalRate":       lambda: sim.ionisation.local_rate[iz, ::self.every_r],
                "IntegratedRate":  lambda: sim.ionisation.integrated_rate[iz, ::self.every_r],
                "ElectronDensity": lambda: sim.medium.n_electron[::self.every_r],
                "AtomicDensity":   lambda: sim.medium.n_atomic[::self.every_r],
                "GazDensity":      lambda: (sim.medium.n_gaz*self.dummy_r)[::self.every_r],
                "k":               lambda: (sim.pulse.k*self.dummy_r)[::self.every_r] if io == 0 else (sim.harmonic.k[io-1]*self.dummy_r)[::self.every_r],
                "dk":              lambda: sim.pulse.dk[::self.every_r] if io == 0 else sim.harmonic.dk[io-1, ::self.every_r],
            }

            func = data_map.get(self.deposited_quantity)
            if func:
                self.quantity[o].append(np.copy(func()*self.unit))

        if self.save_every is not None:
            assert self.nzdiag is not None
            steps_collected = len(self.quantity[self.order[0]]) // self.nzdiag
            if steps_collected >= self.save_every:
                if iz == (len(sim.main.z) - 1):
                    self.save()

    def save(self) -> None:
        assert self.nzdiag is not None
        assert self.nrdiag is not None
        with h5py.File(self.filename, 'a') as f:
            for o in self.order:
                grp: Any
                if self.deposited_quantity in ["LocalRate","IntegratedRate","ElectronDensity","AtomicDensity","GazDensity","Phase"]:
                    grp = f
                else:
                    grp_name = f"order_{o}"
                    if grp_name in f:
                        grp = f[grp_name]
                    else:
                        grp = f.create_group(grp_name)

                data_list = self.quantity[o]
                if not data_list:
                    continue

                data_array = np.array(data_list)
                n_t_in_buffer = len(data_list) // self.nzdiag
                data_array = data_array.reshape(n_t_in_buffer, self.nzdiag, self.nrdiag)

                for itt in range(n_t_in_buffer):
                    abs_t_idx = self.saved_t_count + itt
                    dataset_name = f"t_{abs_t_idx}"
                    if dataset_name not in grp:
                        grp.create_dataset(dataset_name, data=data_array[itt, :, :], dtype=self.type)

        self.saved_t_count += (len(self.quantity[self.order[0]]) // self.nzdiag)
        self.quantity = {o: [] for o in self.order}


class DiagEnergy(HarpComponent):  # type: ignore[misc]
    diag_path: str = ""
    add_suffix: str  = ""
    every_z: int = 1
    order: list[int] = [1]
    SI: bool = True

    def __init__(self, **kwargs: Any) -> None:
        super(DiagEnergy, self).__init__(**kwargs)
        self.filename = self.diag_path+"DiagEnergy"+self.add_suffix+".h5"
        self.nzdiag: int | None = None
        self.order = list(self.order[:])
        if 1 not in self.order:
            self.order.append(1)
        self.order.sort()
        self.energy: dict[int, RealField] | None = None
        self.unit: Any = 1.


    def initialize(self, main: "Main", harmonic: "Harmonic | None") -> None:
        self.nzdiag = len(main.z[::self.every_z])
        self.energy = {o: np.zeros(self.nzdiag) for o in self.order}

        for o in self.order:
            if o != 1:
                assert harmonic is not None
                if o not in harmonic.order:
                    raise Exception(
                        f"""ERROR: {o} is not a valid harmonic order. Available orders are:
                        {1} for the fundamental and those in harmonic of the input file"""
                    )

        if self.SI:
            self.unit = main.E_r**2*main.eps0c/2. * main.L_r**2 * main.T_r

        with h5py.File(self.filename, 'w') as f:
            if self.SI:
                f.create_dataset("z", data=main.z[::self.every_z]*main.L_r)
            else:
                f.create_dataset("z", data=main.z[::self.every_z])

    def apply(self, iz: int, sim: Any) -> None:
        assert self.energy is not None
        if not (iz % self.every_z == 0):
            return

        for io, o in enumerate(self.order):
            field: Callable[[], Any] = lambda: sim.pulse.E[:] if io == 0 else sim.harmonic.E[(io-1)*sim.main.nr:(io)*sim.main.nr]
            try:
                self.energy[o][iz // self.every_z] += np.trapz(sim.main.r*np.abs(field())**2,x=sim.main.r)*sim.main.dt  # type: ignore[attr-defined]
            except AttributeError:
                self.energy[o][iz // self.every_z] += np.trapezoid(sim.main.r*np.abs(field())**2,x=sim.main.r)*sim.main.dt

    def save(self) -> None:
        assert self.energy is not None
        with h5py.File(self.filename, 'a') as f:
            for o in self.order:
                data = self.energy[o]*self.unit* np.pi*2.
                data_name = f"order_{o}"
                if data_name not in f:
                    f.create_dataset(data_name,data=data[:], dtype="float64")


class DiagFarField(HarpComponent):  # type: ignore[misc]
    diag_path: str  = ""
    add_suffix: str  = ""
    positions: Any  = []
    order: list[int]= []
    paraxial_max_factor: float = 0.5
    every_r: int    = 1
    every_t: int    = 1
    SI: bool        = True

    def __init__(self, **kwargs: Any) -> None:
        super(DiagFarField, self).__init__(**kwargs)
        self.filename = self.diag_path+"DiagFarField"+self.add_suffix+".h5"
        self.nrdiag: int | None = None
        self.nzdiag: int | None = None
        self.ntdiag: int | None = None
        self.order = list(self.order[:])
        if 1 not in self.order:
            self.order.append(1)
        self.order.sort()
        self.farfield: dict[int, list[Any]] | None = None
        self.unit: Any   = 1.
        self.no = len(self.order)


    def initialize(self, main: "Main", harmonic: "Harmonic | None") -> None:
        self.nrdiag = len(main.r[::self.every_r])
        self.nzdiag = len(self.positions)
        self.ntdiag = len(main.t[::self.every_t])
        self.positions  = np.array(self.positions)

        self.farfield = {o: [] for o in self.order}

        for o in self.order:
            if o != 1:
                assert harmonic is not None
                if o not in harmonic.order:
                    raise Exception(
                        f"""ERROR: {o} is not a valid harmonic order. Available orders are:
                        {1} for the fundamental and those in harmonic of the input file"""
                    )

        if self.SI:
            self.unit = main.E_r

        with h5py.File(self.filename, 'w') as f:
            if self.SI:
                f.create_dataset("r", data=main.r[::self.every_r]*main.L_r)
                f.create_dataset("z", data=(self.positions[:]+main.Lz)*main.L_r)
                f.create_dataset("t", data=main.t[::self.every_t]*main.T_r)
            else:
                f.create_dataset("r", data=main.r[::self.every_r])
                f.create_dataset("z", data=self.positions[:]+main.Lz)
                f.create_dataset("t", data=main.t[::self.every_t])

    def apply(self, it: int, iz: int, sim: Any) -> None:
        assert self.farfield is not None
        if (iz != sim.main.nz-1) or (it % self.every_t != 0):
            return
        for io, o in enumerate(self.order):
            field: Callable[[], Any] = lambda: sim.pulse.E[:] if io == 0 else sim.harmonic.E[(io-1)*sim.main.nr:(io)*sim.main.nr]

            kmax_hankel   = np.minimum(np.pi/sim.main.dr,self.paraxial_max_factor * o) # min of Nyquist limit and Paraxial approx limit
            k_hankle      = np.linspace(0.,kmax_hankel,sim.main.nr)
            dk_hankle     = k_hankle[1]-k_hankle[0]

            bessel_jv0    = j0(np.outer(k_hankle,sim.main.r))
            j_hankel      = bessel_jv0*(sim.main.r*sim.main.dr)
            j_hankel_inv  = bessel_jv0.T*(k_hankle*dk_hankle)

            field_hankel  = j_hankel @ field()
            for zp in self.positions:
                farfield_hankel = field_hankel * np.exp(-1j*k_hankle**2/(2.*o)*(zp))
                self.farfield[o].append( (j_hankel_inv @ farfield_hankel)[::self.every_r] * self.unit )

    def save(self) -> None:
        assert self.farfield is not None
        assert self.nzdiag is not None
        assert self.nrdiag is not None
        with h5py.File(self.filename, 'a') as f:
            for o in self.order:
                grp_name = f"order_{o}"
                grp = f.create_group(grp_name)

                data_list = self.farfield[o]

                data_array = np.array(data_list)
                n_t_in_buffer = len(data_list) // self.nzdiag
                data_array = data_array.reshape(n_t_in_buffer, self.nzdiag, self.nrdiag)

                for itt in range(n_t_in_buffer):
                    dataset_name = f"t_{itt}"
                    if dataset_name not in grp:
                        grp.create_dataset(dataset_name, data=data_array[itt, :, :], dtype=complex)
