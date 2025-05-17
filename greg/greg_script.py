import openmc 
from openmc.deplete import d1s
from openmc.data import effective_dose
from openmc_fusion_benchmarks.neutron_sources import fng_source
import subprocess

model = openmc.Model()

# Define the materials
material = openmc.Materials.from_xml('materials.xml')
model.materials = material

# Define the geometry
geom = openmc.Geometry.from_xml('geometry.xml')
model.geometry = geom

# Define Settings
# source definition
fng_center = (0, 0, 0)
fng_uvw = (0., 1., 0)

source = fng_source(center=fng_center,
                    reference_uvw=fng_uvw, beam_energy=230)

# settings
settings = openmc.Settings(run_mode='fixed source')
settings.batches = 500
settings.particles = 10000
settings.source = source
settings.photon_transport = True
settings.verbosity = 5
settings.output = {'tallies': False}
settings.use_decay_photons = True

model.settings = settings

############################################################################
# Specify Tallies
model.tallies = openmc.Tallies()

# particle filters
neutron_filter = openmc.ParticleFilter(['neutron'])
photon_filter = openmc.ParticleFilter(['photon'])

# cell filters
cell_filter_9 = openmc.CellFilter([621, 622, 624, 625])
cell_filter_10 = openmc.CellFilter([178, 179, 180, 181, 182, 230, 231, 232, 233, 234, 235, 236, 237, 238, 239, 636, 637, 638])
cell_filter_photon = openmc.CellFilter([227])

# Energy filters
energies, coeffs = effective_dose.dose.dose_coefficients('photon', data_source='icrp74')
dose_filter = openmc.EnergyFunctionFilter(energies, coeffs)

tally = openmc.Tally(name="sdr_dose")
tally.filters = [cell_filter_photon, photon_filter, dose_filter]
tally.scores = ['flux']
model.tallies.append(tally)

tally = openmc.Tally(name="foil_activation")
tally.filters = [cell_filter_10, neutron_filter]
tally.scores = ['(n,gamma)', '(n,p)', '(n,2n)', '(n,3n)']
model.tallies.append(tally)

tally = openmc.Tally(name="neutron_flux")
tally.filters = [cell_filter_photon, neutron_filter]
tally.scores = ['flux']
model.tallies.append(tally)

nuclides = d1s.prepare_tallies(model, chain_file='/home/davide/openmc_models/CROSS_SECTIONS/decay/endfb-viii.0-chain/chain_endfb80_sfr.xml')
model.export_to_xml()

model.run()

output_path = 'statepoint.500.h5'

# campaign_1
time_steps = [19440, 61680, 32940, 54840, 15720, 6360, 8940]
source_rates = [2.32e10, 0.0, 2.87e10, 0.0, 1.90e10, 0.0, 1.36e10]

factors = d1s.time_correction_factors(nuclides, time_steps, source_rates)
with openmc.StatePoint(output_path) as sp:
    tally_dose = sp.get_tally(name='sdr_dose')
    tally_rr = sp.get_tally(name='foil_activation')
    tally_neutron = sp.get_tally(name='neutron_flux')

print(tally_neutron.get_values())
tally_dose = d1s.apply_time_correction(tally_dose, factors)
print(tally_dose.get_values())
print(tally_rr.get_values())