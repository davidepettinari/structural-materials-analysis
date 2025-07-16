import openmc
import openmc.deplete
import numpy as np
from pathlib import Path
import structural_material_analysis as sma
import tokamak_radiation_environment as tre

openmc.config['chain_file'] = r'/home/davide/openmc_models/CROSS_SECTIONS/decay/endfb-viii.0-chain/chain_endfb80_sfr.xml'