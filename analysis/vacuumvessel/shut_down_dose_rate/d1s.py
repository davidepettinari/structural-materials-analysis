import openmc
from openmc.deplete import d1s
import numpy as np
import structural_material_analysis as sma
import tokamak_radiation_environment as tre

openmc.config['chain_file'] = r'/home/davide/openmc_models/CROSS_SECTIONS/decay/endfb-viii.0-chain/chain_endfb80_sfr.xml'