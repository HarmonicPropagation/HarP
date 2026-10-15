from components import Main
from pulse import Pulse
from medium import Medium
from ionisation import Ionisation
from diagnostic import Diagnostics, DiagEnergy, DiagFarField
from harmonic import Harmonic
import numpy as np
import os
from scipy import constants



######################### Physical constants #########################
Temperature         = 293                       # gaz temperature, K
lambda0             = 0.81e-6                   # laser wavelength, m
c                   = constants.c               # lightspeed, m/s
omega0              = 2*np.pi*c/lambda0         # laser angular frequency, rad/s
eps0                = constants.epsilon_0       # Vacuum permittivity, F/m
e                   = constants.e               # Elementary charge, C
me                  = constants.m_e             # Electron mass, kg
kb                  = constants.Boltzmann       # Boltmann constant, J/K

######################### Variables for unit conversions #########################
T_r                 = 1./omega0                 # reference time, s
L_r                 = lambda0/2./np.pi          # reference length, m
n_r                 = eps0*omega0**2*me/e**2    # reference number density, m-3
E_r                 = me*omega0*c/e             # reference electric field, V/m
I_r                 = eps0*c*E_r**2/2. *1E-4    # reference intensity, W/cm2
Pol_r               = me*omega0**2/e**2         # reference polarisabiltiy, 
mm                  = 1.e-3/L_r                 # 1 mm                      in normalized units
fs                  = 1.e-15/T_r                # 1 femtosecond             in normalized units


######################### Laser Parameters #########################
z_R    = 20.*mm                                # rayleigh lenght
w0     = np.sqrt(z_R*(lambda0/L_r)/np.pi)       # laser waist
tfwhm  = 22*fs                                  # laser FWHM in intensity
I0     = 2.5 * 1e14/I_r                         # laser maximal peak intensity if propagating in vaccum (corresponding to 2.5 W/cm2)
E0     = np.sqrt(I0)                            # laser maximal electric field if propagating in vaccum

######################### Gaz Parameters  #########################
p_in_mbar = 20.                # set pressure in mbar assuming p*zr=xxx mbar*cm 
p         = p_in_mbar*1E2                       # pressure, Pa
n0        = p/kb/Temperature/n_r                # gaz density in code unit
Lgaz      = 10*mm                       # gaz lenght

def gazprofile(r,z):
    """if depend on r need to set anne_eq=False"""
    return n0*np.heaviside(z-0.,1.)*np.heaviside(Lgaz-z,1.)


######################### Simulations paramaters #########################
nr  = 512                                       # number of mesh points in the transverse direction
nz  = 256                                       # number of mesh points in the longitudinal direction 
nt  = 256                                       # number of mesh points in time

Lr    = 4.*w0                                   # radial size of the box
Lz    = Lgaz                                    # longitudinal size of the box
Tsim  = 4.*tfwhm                                # total simulated time

dr = Lr/nr                                      # radial resolution 
dz = Lz/nz                                      # longitudinal resolution 
dt = Tsim/nt                                    # time resolution 

######################### Laser profile #########################
tfocus  = Tsim/2.                               # focal time-position of the laser
zfocus  = Lz/2.                                 # focal longitudinal-position of the laser

def tenv(t):
    return np.exp(-2*np.log(2)*(t-tfocus)**2/tfwhm**2)

def eprofile(r):                                # radial profile of the pulse at z=0
    wz       = w0 * np.sqrt(1+ (zfocus/z_R)**2 )
    R        = zfocus*(1.+(z_R/zfocus)**2)
    Phi_gouy = np.arctan2(zfocus,z_R)
    return w0/wz * np.exp(-(r/wz)**2)*np.exp(-1j*(r**2/2./R-Phi_gouy))


##### Paths
tables_path = 'tables/'                         # path for the TDSE tables (relative to where HarP.py is run from)
diag_path   = os.path.dirname(os.path.abspath(__file__)) + os.sep  # diagnostics are written next to this namelist by default


######################### Main #########################
Main(
    wavelength                  = lambda0,
    dr                          = dr,
    dz                          = dz,
    Lr                          = Lr,
    Lz                          = Lz, 
    dt                          = dt,
    simulation_time             = Tsim,
    solver                      = "CrankNicolson",
    boundary_condition          = "absorbing",
    anne_eq                     = True,
    predictor_corrector         = False,
    print_every                 = 10,
    approx_check                = True
)

#########################  Laser #########################
Pulse(
    amplitude                   = E0,
    radial_profile              = eprofile,
    time_envelope               = tenv
)

######################### Medium #########################
Medium(
    gaz_profile                 = gazprofile,
    polarizability_0            = 1.8247e-40,      # in SI from P. Soldan et al. 2001 PCCP 
    gamma_kerr                  = 7.2767118e-62/6. ,  # in SI from P. Soldan et al. 2001 PCCP 
    table_path                  = tables_path+"form_factor/f12_Ar_henkegov_lhuillier90.txt",
    unit_ev                     = 1.,         # photon energy in ev
    unit_f1                     = 1.,        # polarisability as form factor f1
    unit_f2                     = 1.,         # absorption cross section as form factor f2
    skiprows                    = 1
)

Ionisation(
    table_path          = tables_path+"ionisation/ionisation_Ar_810_tdse.txt",
    skiprows            = 2,
    intensity_unit_Wcm2 = 1.,
    rate_unit_SI        = 2.         #add statistical factor of 2 on the rate due to double occupancy state SAE approx
)


order = np.arange(3,60,2,dtype=int)#[23]                                  
Harmonic(
    order = order,                             #photon energy corresponding to order*hbar*omega
    table_path = tables_path+'tdse/dipoles_Ar_810.txt',
    statistical_factor = 2.,        #add statistical factor of 2 on the dipole due to double occupancy state SAE approx
    intensity_unit_Wcm2 = 1e14
)

every_r=1
every_z=1
every_t=1
save_every=50

DiagEnergy(
    diag_path = diag_path,
    every_z = 1,
    order   = order,
    SI      = True
)

Diagnostics(
    diag_path          = diag_path,
    deposited_quantity = "RealField",
    every_r            = every_r,
    every_z            = 64,
    every_t            = every_t,
    save_every         = save_every,
    order              = order,
    SI                 = True
)

#Diagnostics(
#    diag_path          = diag_path,
#    deposited_quantity = "Field",
#    every_r            = every_r,
#    every_z            = every_z,
#    every_t            = every_t,
#    save_every         = save_every,
#    order              = order,
#    SI                 = True
#)

#DiagFarField(
#    diag_path = diag_path,
#    name      = "HankelField",
#    positions = [15.*mm],
#    order     = order,
#    paraxial_max_factor = 0.5,
#    every_r   = 1,
#    every_t   = 1,
#    SI        = True
#)

Diagnostics(
    diag_path          = diag_path,
    deposited_quantity = "AtomicDensity",
    every_r            = every_r,
    every_z            = every_z,
    every_t            = 64,
    save_every         = save_every,
    SI                 = True
)

#Diagnostics(
#    diag_path          = diag_path,
#    deposited_quantity = "GazDensity",
#    every_r            = every_r,
#    every_z            = every_z,
#    every_t            = 64,
#    save_every         = save_every,
#    SI                 = True
#)


Diagnostics(
    deposited_quantity = "ElectronDensity",
    every_r            = every_r,
    every_z            = every_z,
    every_t            = 64,
    save_every         = save_every,
    SI                 = True
)

#Diagnostics(
#    deposited_quantity = "LocalRate",
#    every_r            = every_r,
#    every_z            = every_z,
#    every_t            = every_t
#)

#Diagnostics(
#    deposited_quantity = "IntegratedRate",
#    every_r            = every_r,
#    every_z            = every_z,
#    every_t            = every_t
#)

