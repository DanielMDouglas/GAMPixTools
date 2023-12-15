import numpy as np
"""
Set physics constants
"""

## Physical params
#: Recombination :math:`\alpha` constant for the Box model
BOX_ALPHA = 0.93
#: Recombination :math:`\beta` value for the Box model in :math:`(kV/cm)(g/cm^2)/MeV`
BOX_BETA = 0.207 #0.3 (MeV/cm)^-1 * 1.383 (g/cm^3)* 0.5 (kV/cm), R. Acciarri et al JINST 8 (2013) P08005
#: Recombination :math:`A_b` value for the Birks Model
BIRKS_Ab = 0.800
#: Recombination :math:`k_b` value for the Birks Model in :math:`(kV/cm)(g/cm^2)/MeV`
BIRKS_kb = 0.0486 # g/cm2/MeV Amoruso, et al NIM A 523 (2004) 275
#: Electron charge in Coulomb
E_CHARGE = 1.602e-19
#: Average energy expended per ion pair in LAr in :math:`MeV` from Phys. Rev. A 10, 1452
W_ION = 23.6e-6

## Quenching parameters
BOX = 1
BIRKS = 2




"""
Set detector constants
"""
# only the following are currently used #
#: Liquid argon density in :math:`g/cm^3`
LAR_DENSITY = 1.38 # g/cm^3
#: Electric field magnitude in :math:`kV/cm`
E_FIELD = 0.50 # kV/cm

# unused #
#: Detector temperature in K
TEMPERATURE = 87.17
#: Drift velocity in :math:`cm/\mu s`
V_DRIFT = 0.1648 # cm / us,
#: Electron lifetime in :math:`\mu s`
ELECTRON_LIFETIME = 2.2e3 # us,
#: Time sampling in :math:`\mu s`
TIME_SAMPLING = 0.1 # us
#: Drift time window in :math:`\mu s`
TIME_INTERVAL = (0, 200.) # us
#: Signal time window padding in :math:`\mu s`
TIME_PADDING = 10
#: Number of sampled points for each segment slice
SAMPLED_POINTS = 40
#: Longitudinal diffusion coefficient in :math:`cm^2/\mu s`
LONG_DIFF = 4.0e-6 # cm * cm / us
#: Transverse diffusion coefficient in :math:`cm^2/\mu s`
TRAN_DIFF = 8.8e-6 # cm * cm / us
#: Numpy array containing all the time ticks in the drift time window
TIME_TICKS = np.linspace(TIME_INTERVAL[0],
                         TIME_INTERVAL[1],
                         int(round(TIME_INTERVAL[1]-TIME_INTERVAL[0])/TIME_SAMPLING)+1)
#: Time window of current response in :math:`\mu s`
TIME_WINDOW = 8.9 # us
#: TPC drift length in :math:`cm`
DRIFT_LENGTH = 0
#: Time sampling in the pixel response file in :math:`\mu s`
RESPONSE_SAMPLING = 0.1
#: Spatial sampling in the pixel reponse file in :math:`cm`
RESPONSE_BIN_SIZE = 0.04434
#: Borders of each TPC volume in :math:`cm`
TPC_BORDERS = np.zeros((0, 3, 2))
#: TPC offsets wrt the origin in :math:`cm`
TPC_OFFSETS = np.zeros((0, 3, 2))
#: Pixel tile borders in :math:`cm`
TILE_BORDERS = np.zeros((2,2))
#: Default value for pixel_plane, to indicate out-of-bounds edep
DEFAULT_PLANE_INDEX = 0x0000BEEF
#: Total number of pixels
N_PIXELS = 0, 0
#: Number of pixels in each tile
N_PIXELS_PER_TILE = 0, 0
#: Dictionary between pixel ID and its position in the pixel array
PIXEL_CONNECTION_DICT = {}
#: Pixel pitch in :math:`cm`
PIXEL_PITCH = 0.4434
#: Tile position wrt the center of the anode in :math:`cm`
TILE_POSITIONS = {}
#: Tile orientations in each anode
TILE_ORIENTATIONS = {}
#: Map of tiles in each anode
TILE_MAP = ()
#: Association between chips and io channels
TILE_CHIP_TO_IO = {}
#: Association between modules and io groups
MODULE_TO_IO_GROUPS = {}
#: Association between modules and tpcs
MODULE_TO_TPCS = {}
TPC_TO_MODULE = {}

ELECTRON_MOBILITY_PARAMS = 551.6, 7158.3, 4440.43, 4.29, 43.63, 0.2053





"""
Sets ligth-related constants
"""

#: Number of true segments to track for each time tick (`MAX_MC_TRUTH_IDS=0` to disable complete truth tracking)
MAX_MC_TRUTH_IDS = 0 #256
#: Threshold for propogating truth information on a given SiPM
MC_TRUTH_THRESHOLD = 0.1 # pe/us
ENABLE_LUT_SMEARING = False

N_OP_CHANNEL = 0
LIGHT_SIMULATED = True
OP_CHANNEL_EFFICIENCY = np.zeros(0)
OP_CHANNEL_TO_TPC = np.zeros(0)
TPC_TO_OP_CHANNEL = np.zeros((0,0))

#: Prescale factor analogous to ScintPreScale in LArSoft FIXME
SCINT_PRESCALE = 1
#: Ion + excitation work function in `MeV`
W_PH = 19.5e-6 # MeV

#: Step size for light simulation [microseconds]
LIGHT_TICK_SIZE = 0.001 # us
#: Pre- and post-window for light simulation [microseconds]
LIGHT_WINDOW = (1, 10) # us

#: Fraction of total light emitted from singlet state
SINGLET_FRACTION = 0.3
#: Singlet decay time [microseconds]
TAU_S = 0.001 # us
#: Triplet decay time [microseconds]
TAU_T = 1.530

#: Conversion from PE/microsecond to ADC
DEFAULT_LIGHT_GAIN = -2.30 # ADC * us/PE
LIGHT_GAIN = np.zeros((0,))
#: Set response model type (0=RLC response, 1=arbitrary input)
SIPM_RESPONSE_MODEL = 0
#: Response RC time [microseconds]
LIGHT_RESPONSE_TIME = 0.055
#: Reponse oscillation period [microseconds]
LIGHT_OSCILLATION_PERIOD = 0.095
#: Sample rate for input noise spectrum [microseconds]
LIGHT_DET_NOISE_SAMPLE_SPACING = 0.01 # us
#: Arbitrary input model (normalized to sum of 1)
IMPULSE_MODEL = np.array([1,0])
#: Arbitrary input model tick size [microseconds]
IMPULSE_TICK_SIZE = 0.001

#: Number of SiPMs per detector (used by trigger)
OP_CHANNEL_PER_TRIG = 6
#: Total detector light threshold [ADC] (one value for every OP_CHANNEL_PER_TRIG detector sum)
LIGHT_TRIG_THRESHOLD = np.zeros((0,))
#: Light digitization window [microseconds]
LIGHT_TRIG_WINDOW = (0.9, 1.66) # us
#: Light waveform sample rate [microseconds]
LIGHT_DIGIT_SAMPLE_SPACING = 0.01 # us
#: Light digitizer bits
LIGHT_NBIT = 10
