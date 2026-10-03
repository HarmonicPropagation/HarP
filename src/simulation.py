from typing import TYPE_CHECKING
import numpy as np
import time

if TYPE_CHECKING:
    from components import Main
    from medium import Medium
    from pulse import Pulse
    from ionisation import Ionisation
    from harmonic import Harmonic
    from solver import Solver
    from diagnostic import Diagnostics, DiagEnergy, DiagFarField

class Simulation:
    def __init__(self, main: "Main", medium: "Medium", pulse: "Pulse", ionisation: "Ionisation",
                 harmonic: "Harmonic | None", solver: "Solver", diagnostics: list["Diagnostics"],
                 diagenergy: "DiagEnergy | None", diagfarfield: "DiagFarField | None") -> None:
        self.main         = main
        self.pulse        = pulse
        self.medium       = medium
        self.ionisation   = ionisation
        self.harmonic     = harmonic
        self.solver       = solver
        self.diagnostics  = diagnostics
        self.diagenergy   = diagenergy
        self.diagfarfield = diagfarfield

        self.ndiags       = len(diagnostics)
        assert self.main.nt is not None
        self.w            = len(str(self.main.nt))

    def printer(self, it: int) -> None:
        assert self.main.print_every is not None
        if it % self.main.print_every == 0:
            self.main.printer(f"        Iteration {it:{self.w}d}/{(self.main.nt-1)}   |   Time: {time.time()-self.main.timer:>7.2f}s")


    # --- Diagnostics ---
    def applyAllDiag(self, it: int, iz: int) -> None:
        t0 = time.time()
        for id in range(self.ndiags):
            self.diagnostics[id].apply(it, iz, self)
        if self.diagenergy is not None:
            self.diagenergy.apply(iz, self)
        if self.diagfarfield is not None:
            self.diagfarfield.apply(it, iz, self)
        self.main.diag_timer += time.time() - t0

    def saveAllDiag(self) -> None:
        for id in range(self.ndiags):
            self.diagnostics[id].save()
        if self.diagenergy is not None:
            self.diagenergy.save()
        if self.diagfarfield is not None:
            self.diagfarfield.save()

    def run(self) -> None:
        assert self.main.nt is not None
        assert self.main.nz is not None

        self.main.printer(" ------------------------ HarP LOOP ------------------------ ")
        self.main.printer("                                                     ")
        self.main.printer("                                                     ")
        if self.harmonic is None:
            self.main.timer = time.time()

            for it in range (self.main.nt):
                self.medium.reinitialize(self)
                self.pulse.reinitialize(self, it)
                self.applyAllDiag(it, 0)
                for iz in range(self.main.nz - 1):
                    self.ionisation.update(self, iz)
                    self.medium.update(self, iz)
                    self.pulse.update(self, iz)
                    self.applyAllDiag(it, iz + 1)
                self.ionisation.compute_integrated_rate(self)
                self.printer(it)

        else:
            self.main.timer = time.time()
            for it in range (self.main.nt):
                self.medium.reinitialize(self)
                self.pulse.reinitialize(self, it)
                self.harmonic.reinitialize()
                self.applyAllDiag(it, 0)
                for iz in range(self.main.nz - 1):
                    self.ionisation.update(self, iz)
                    self.medium.update(self, iz)
                    self.pulse.update(self, iz)
                    self.harmonic.update(self)
                    self.applyAllDiag(it, iz + 1)
                self.ionisation.compute_integrated_rate(self)
                self.printer(it)

        self.main.printer("                                                     ")
        self.main.printer("                                                     ")
        self.main.printer(" ------------------------ SAVE DATA ------------------------ ")
        self.main.printer("                                                     ")
        self.saveAllDiag()
        self.main.printer(f"       Time in Diagnostics              {self.main.diag_timer:>7.2f}s")
        self.main.printer(f"       Time in Pulse                    {self.main.pulse_timer:>7.2f}s")
        self.main.printer(f"       Time in Harmonic                 {self.main.harmonic_timer:>7.2f}s")
        self.main.printer(f"       Time in Medium                   {self.main.medium_timer:>7.2f}s")
        self.main.printer("                                                     ")
        self.main.printer(" --------------------------- END --------------------------- ")
