from components import Main
from pulse import Pulse
from medium import Medium
from ionisation import Ionisation
from harmonic import Harmonic
from diagnostic import Diagnostics, DiagEnergy, DiagFarField
from simulation import Simulation
from solver import Solver


class SimulationBuilder:

    @staticmethod
    def build() -> Simulation:
        main = Main.get_instance()
        main.printer("""
              ██╗  ██╗  █████╗  ██████╗  ██████╗
              ██║  ██║ ██╔══██╗ ██╔══██╗ ██╔══██╗
              ███████║ ███████║ ██████╔╝ ██████╔╝
              ██╔══██║ ██╔══██║ ██╔══██╗ ██╔═══╝
              ██║  ██║ ██║  ██║ ██║  ██║ ██║
              ╚═╝  ╚═╝ ╚═╝  ╚═╝ ╚═╝  ╚═╝ ╚═╝
        """)
        main.printer("                                                     ")

        main.printer(" -------------------------- BUILD -------------------------- ")
        main.printer("                                                     ")


        medium = Medium.get_instance()
        medium.initialize(main)

        ionisation = Ionisation.get_instance()
        ionisation.initialize(main)

        pulse = Pulse.get_instance()
        pulse.initialize(main, medium)

        harmonic = Harmonic.get_instance(mandatory=False)
        if harmonic is not None:
            harmonic.initialize(main, medium, pulse)
            main.nh = harmonic.nh

        diagnostics = [
            d for d in (Diagnostics._list)
        ]
        for diag in diagnostics:
            diag.initialize(main, harmonic)

        diagenergy = DiagEnergy.get_instance(mandatory=False)
        if diagenergy is not None:
            diagenergy.initialize(main, harmonic)

        diagfarfield = DiagFarField.get_instance(mandatory=False)
        if diagfarfield is not None:
            diagfarfield.initialize(main, harmonic)

        solver = Solver()
        solver.initialize(main)

        sim = Simulation(
            main=main,
            medium=medium,
            pulse=pulse,
            ionisation=ionisation,
            harmonic=harmonic,
            solver=solver,
            diagnostics=diagnostics,
            diagenergy=diagenergy,
            diagfarfield=diagfarfield
        )

        main.printer(f"       Radial cells                     {main.nr}")
        main.printer(f"       Transversal cells                {main.nz}")
        main.printer(f"       Timesteps                        {main.nt}")
        main.printer("                                                             ")


        return sim
