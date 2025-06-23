import openmc
import numpy as np
import structural_material_analysis as sma
import tokamak_radiation_environment as tre


def reactor_model(inner_vessel_material):
    """Returns an openmc model of the ARC-class model.

    Returns:
        the openmc model
    """
    dt_plasma = sma.materials.dt_plasma
    tungsten = sma.materials.tungsten
    flibe = sma.materials.flibe
    beryllium = sma.materials.beryllium
    boroncarbide = sma.materials.boroncarbide
    inconel718 = sma.materials.inconel718

    materials = openmc.Materials([
        dt_plasma,
        tungsten,
        inner_vessel_material,
        flibe, 
        boroncarbide, 
        beryllium, 
        ])

    (root, geometry) = sma.geometry.reactor_geometry(dt_plasma, tungsten, inner_vessel_material, flibe, beryllium, inconel718, boroncarbide)
    
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
    # TBR
    # makes a tally to distinguish tritium production from Li6 and Li7
    tbr_cc_tally = openmc.Tally(name="TBR channel")
    tbr_cc_tally.filters = [cc_filter]
    tbr_cc_tally.scores = ["(n,Xt)"]
    tbr_cc_tally.nuclides = [openmc.Nuclide("Li6"), openmc.Nuclide("Li7")]
    tallies.append(tbr_cc_tally)
    tbr_blanket_tally = openmc.Tally(name="TBR blanket")
    tbr_blanket_tally.filters = [blanket_filter]
    tbr_blanket_tally.scores = ["(n,Xt)"]
    tbr_blanket_tally.nuclides = [openmc.Nuclide("Li6"), openmc.Nuclide("Li7")]
    tallies.append(tbr_blanket_tally)
    tbr_tally = openmc.Tally(name="TBR")
    tbr_tally.scores = ["(n,Xt)"]
    tbr_tally.nuclides = [openmc.Nuclide("Li6"), openmc.Nuclide("Li7")]
    tallies.append(tbr_tally)
    tbr_energy_tally = openmc.Tally(name="TBR(energy)")
    tbr_energy_tally.filters = [energy_filter]
    tbr_energy_tally.scores = ["(n,Xt)"]
    tallies.append(tbr_energy_tally)

    # NEUTRON SPECTRA
    # FW
    spectrafw_cell_tally = openmc.Tally(name="neutron_spectra_fw")
    spectrafw_cell_tally.filters = [fw_filter, neutron_particle_filter, energy_filter]
    spectrafw_cell_tally.scores = ["flux"]
    tallies.append(spectrafw_cell_tally)
    # STR1
    spectra_str1_cell_tally = openmc.Tally(name="neutron_spectra_str1")
    spectra_str1_cell_tally.filters = [str1_filter, neutron_particle_filter, energy_filter]
    spectra_str1_cell_tally.scores = ["flux"]
    tallies.append(spectra_str1_cell_tally)
    # STR2
    spectra_str2_cell_tally = openmc.Tally(name="neutron_spectra_str2")
    spectra_str2_cell_tally.filters = [str2_filter, neutron_particle_filter, energy_filter]
    spectra_str2_cell_tally.scores = ["flux"]
    tallies.append(spectra_str2_cell_tally)

    # HEATING
    # FW
    heating_fw_tally = openmc.Tally(name="heating_fw")
    heating_fw_tally.scores = ["heating"]
    heating_fw_tally.filters = [fw_filter]
    tallies.append(heating_fw_tally)
    heating_fw_n_tally = openmc.Tally(name="heating_fw_n")
    heating_fw_n_tally.scores = ["heating"]
    heating_fw_n_tally.filters = [fw_filter, neutron_particle_filter]
    tallies.append(heating_fw_n_tally)
    heating_fw_e_tally = openmc.Tally(name="heating_fw_electron")
    heating_fw_e_tally.scores = ["heating"]
    heating_fw_e_tally.filters = [fw_filter, electron_particle_filter]
    tallies.append(heating_fw_e_tally)
    heating_fw_g_tally = openmc.Tally(name="heating_fw_photon")
    heating_fw_g_tally.scores = ["heating"]
    heating_fw_g_tally.filters = [fw_filter, photon_particle_filter]
    tallies.append(heating_fw_g_tally)
    heating_fw_p_tally = openmc.Tally(name="heating_fw_positron")
    heating_fw_p_tally.scores = ["heating"]
    heating_fw_p_tally.filters = [fw_filter, positron_particle_filter]
    tallies.append(heating_fw_p_tally)
    # STR1
    heating_str1_tally = openmc.Tally(name="heating_str1")
    heating_str1_tally.scores = ["heating"]
    heating_str1_tally.filters = [str1_filter]
    tallies.append(heating_str1_tally)
    heating_str1_n_tally = openmc.Tally(name="heating_str1_n")
    heating_str1_n_tally.scores = ["heating"]
    heating_str1_n_tally.filters = [str1_filter, neutron_particle_filter]
    tallies.append(heating_str1_n_tally)
    heating_str1_e_tally = openmc.Tally(name="heating_str1_electron")
    heating_str1_e_tally.scores = ["heating"]
    heating_str1_e_tally.filters = [str1_filter, electron_particle_filter]
    tallies.append(heating_str1_e_tally)
    heating_str1_g_tally = openmc.Tally(name="heating_str1_photon")
    heating_str1_g_tally.scores = ["heating"]
    heating_str1_g_tally.filters = [str1_filter, photon_particle_filter]
    tallies.append(heating_str1_g_tally)
    heating_str1_p_tally = openmc.Tally(name="heating_str1_positron")
    heating_str1_p_tally.scores = ["heating"]
    heating_str1_p_tally.filters = [str1_filter, positron_particle_filter]
    tallies.append(heating_str1_p_tally)
    # STR2
    heating_str2_tally = openmc.Tally(name="heating_str2")
    heating_str2_tally.scores = ["heating"]
    heating_str2_tally.filters = [str2_filter]
    tallies.append(heating_str2_tally)
    heating_str2_n_tally = openmc.Tally(name="heating_str2_n")
    heating_str2_n_tally.scores = ["heating"]
    heating_str2_n_tally.filters = [str2_filter, neutron_particle_filter]
    tallies.append(heating_str2_n_tally)
    heating_str2_e_tally = openmc.Tally(name="heating_str2_electron")
    heating_str2_e_tally.scores = ["heating"]
    heating_str2_e_tally.filters = [str2_filter, electron_particle_filter]
    tallies.append(heating_str2_e_tally)
    heating_str2_g_tally = openmc.Tally(name="heating_str2_photon")
    heating_str2_g_tally.scores = ["heating"]
    heating_str2_g_tally.filters = [str2_filter, photon_particle_filter]
    tallies.append(heating_str2_g_tally)
    heating_str2_p_tally = openmc.Tally(name="heating_str2_positron")
    heating_str2_p_tally.scores = ["heating"]
    heating_str2_p_tally.filters = [str2_filter, positron_particle_filter]
    tallies.append(heating_str2_p_tally)

    # DPA
    # FW
    damage_fw_tally = openmc.Tally(name="dpa_fw")
    damage_fw_tally.scores = ["damage-energy"]
    damage_fw_tally.filters = [fw_filter]
    tallies.append(damage_fw_tally)
    # STR1
    damage_str1_tally = openmc.Tally(name="dpa_str1")
    damage_str1_tally.scores = ["damage-energy"]
    damage_str1_tally.filters = [str1_filter]
    tallies.append(damage_str1_tally)
    # STR2
    damage_str2_tally = openmc.Tally(name="dpa_str2")
    damage_str2_tally.scores = ["damage-energy"]
    damage_str2_tally.filters = [str2_filter]
    tallies.append(damage_str2_tally)

    model = openmc.Model(
        materials=materials, geometry=geometry, settings=settings, tallies=tallies
    )
    
    return model, root

############################################################################

available_materials = {
    "inconel718": sma.materials.inconel718,
    "ss316l": sma.materials.ss316l,
}

if __name__ == "__main__":
    import argparse
    import os

    parser = argparse.ArgumentParser(description="Run OpenMC model with selected structural material")
    parser.add_argument(
        "-m",
        type=str,
        default="inconel718",
        help="Material to assign to first wall (e.g. inconel718, ss316l, eurofer, etc.)"
    )
    args = parser.parse_args()

    if args.m not in available_materials:
        raise ValueError(f"Material '{args.m}' not recognized. Available: {list(available_materials)}")

    vacuum_vessel_material = available_materials[args.m]

    model, _ = reactor_model(vacuum_vessel_material)
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