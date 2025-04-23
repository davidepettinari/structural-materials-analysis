import openmc
from openmc.deplete import d1s
import numpy as np
import structural_material_analysis as sma
import tokamak_radiation_environment as tre

openmc.config['chain_file'] = r'/home/davide/openmc_models/CROSS_SECTIONS/decay/endfb-viii.0-chain/chain_endfb80_sfr.xml'

def reactor_model(first_wall_material):
    """Returns an openmc model of the ARC-class model.

    Returns:
        the openmc model
    """
    dt_plasma = sma.materials.dt_plasma
    inconel718 = sma.materials.inconel718
    flibe = sma.materials.flibe
    beryllium = sma.materials.beryllium
    boroncarbide = sma.materials.boroncarbide

    materials = openmc.Materials([
        dt_plasma, 
        flibe, 
        boroncarbide, 
        beryllium, 
        first_wall_material, 
        inconel718])

    (root, geometry) = sma.geometry.reactor_geometry(dt_plasma, first_wall_material, inconel718, flibe, beryllium, boroncarbide)
    
    ############################################################################
    # Define Settings
    source = openmc.IndependentSource()
    source.particle = "neutron"
    radius = openmc.stats.Discrete([285], [1])
    z_values = openmc.stats.Discrete([0], [1])
    angle = openmc.stats.Uniform(a=np.radians(-10), b=np.radians(10))
    source.space = openmc.stats.CylindricalIndependent(
        r=radius, phi=angle, z=z_values, origin=(0.0, 0.0, 0.0)
    )
    source.angle = openmc.stats.Isotropic()
    source.energy = openmc.stats.muir(e0=14080000.0, m_rat=5.0, kt=20000.0)

    settings = openmc.Settings()
    settings.batches = 100
    settings.particles = 10000
    settings.inactive = 0
    settings.run_mode = "fixed source"
    settings.source = source
    settings.photon_transport = True  # This line is required to switch on photons tracking
    settings.output = {"tallies": False}

    # FILTERS
    # adds a tally to record neutron spectra on STR1
    energy_bins = openmc.mgxs.GROUP_STRUCTURES["CCFE-709"]
    energy_filter = openmc.EnergyFilter(energy_bins)

    plasma_filter = openmc.CellFilter(root[0])
    sol_filter = openmc.CellFilter(root[1])
    fw_filter = openmc.CellFilter(root[2])
    str1_filter = openmc.CellFilter(root[3])
    cc_filter = openmc.CellFilter(root[4])
    nm_filter = openmc.CellFilter(root[5])
    str2_filter = openmc.CellFilter(root[6])
    blanket_filter = openmc.CellFilter(root[7])
    shield_filter = openmc.CellFilter(root[8])

    neutron_particle_filter = openmc.ParticleFilter(["neutron"])
    photon_particle_filter = openmc.ParticleFilter(["photon"])
    electron_particle_filter = openmc.ParticleFilter(["electron"])
    positron_particle_filter = openmc.ParticleFilter(["positron"])

    # creates a mesh that covers the geometry
    mesh = openmc.RegularMesh()
    mesh.dimension = [300, 100, 200]
    mesh.lower_left = [
        100,
        -100,
        -250,
    ]  # x,y,z coordinates start at 0 as this is a sector model
    mesh.upper_right = [500, 0, 250]
    mesh_filter = openmc.MeshFilter(mesh)  # creating a mesh

    tallies = openmc.Tallies()

    energies, pSv_cm2 = openmc.data.dose_coefficients(particle="photon", geometry="ISO")
    dose_filter = openmc.EnergyFunctionFilter(
        energies, pSv_cm2, interpolation="cubic"  # interpolation method recommended by ICRP
    )
    flux_tally = openmc.Tally()
    flux_tally.filters = [fw_filter, dose_filter, photon_particle_filter]
    flux_tally.scores = ["flux"]
    flux_tally.name = "photon_dose_on_mesh"
    tallies.append(flux_tally)

    model = openmc.Model(
        materials=materials, geometry=geometry, settings=settings, tallies=tallies
    )
    model.settings.photon_transport = True
    model.settings.use_decay_photons = True
    d1s.prepare_tallies(model)
    
    return model, root

############################################################################

available_materials = {
    "tungsten": sma.materials.tungsten,
    "zrnbc": sma.materials.zrnbc,
    "tinbc": sma.materials.tinbc,
    "vmoc": sma.materials.vmoc,
    "vtac": sma.materials.vtac,
    "zrhfc": sma.materials.zrhfc,
    "zrtac": sma.materials.zrtac,
    "nbhfc": sma.materials.nbhfc,
    "hftac": sma.materials.hftac,
}

if __name__ == "__main__":
    import argparse
    import os

    parser = argparse.ArgumentParser(description="Run OpenMC model with selected first wall material")
    parser.add_argument(
        "-m",
        type=str,
        default="tungsten",
        help="Material to assign to first wall (e.g. tungsten, zrnbc, tinbc, etc.)"
    )
    args = parser.parse_args()
    

    if args.m not in available_materials:
        raise ValueError(f"Material '{args.m}' not recognized. Available: {list(available_materials)}")

    first_wall_material = available_materials[args.m]
    print(first_wall_material)

    model, _ = reactor_model(first_wall_material)    
    model.run(geometry_debug=False)

    material_tag = args.m
    original_sp = "statepoint.100.h5"
    renamed_sp = f"{material_tag}.h5"
    if os.path.exists(original_sp):
        os.rename(original_sp, renamed_sp)
        print(f"✔️ Renamed: {original_sp} → {renamed_sp}")

    summary_file = "summary.h5"
    if os.path.exists(summary_file):
        os.remove(summary_file)
        print(f"🗑️ Deleted: {summary_file}")