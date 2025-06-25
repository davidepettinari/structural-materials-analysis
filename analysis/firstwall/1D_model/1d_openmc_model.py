import openmc
import numpy as np
import structural_material_analysis as sma

def create_model(fw_material, plasma_thickness=1.0, fw_thickness=1.5, str1_thickness=1.5,
                 cc_thickness=2.0, nmul_thickness=1.0, str2_thickness=3.0, blanket_thickness=50.0):

    #############################
    #          GEOMETRY         #
    #############################

    delta_yz = 1
    source_plane_position = 0.5

    x0 = openmc.XPlane(x0=0.0, boundary_type='vacuum')
    x1 = openmc.XPlane(x0=plasma_thickness)
    x2 = openmc.XPlane(x0=plasma_thickness + fw_thickness)
    x3 = openmc.XPlane(x0=plasma_thickness + fw_thickness + str1_thickness)
    x4 = openmc.XPlane(x0=plasma_thickness + fw_thickness + str1_thickness + cc_thickness)
    x5 = openmc.XPlane(x0=plasma_thickness + fw_thickness + str1_thickness + cc_thickness + nmul_thickness)
    x6 = openmc.XPlane(x0=plasma_thickness + fw_thickness + str1_thickness + cc_thickness + nmul_thickness + str2_thickness)
    x7 = openmc.XPlane(x0=plasma_thickness + fw_thickness + str1_thickness + cc_thickness + nmul_thickness + str2_thickness + blanket_thickness)
    x7.boundary_type = 'vacuum'

    y0 = openmc.YPlane(y0=-delta_yz/2, boundary_type='reflective')
    y1 = openmc.YPlane(y0=+delta_yz/2, boundary_type='reflective')
    z0 = openmc.ZPlane(z0=-delta_yz/2, boundary_type='reflective')
    z1 = openmc.ZPlane(z0=+delta_yz/2, boundary_type='reflective')

    plas_cell = openmc.Cell(fill=sma.materials.dt_plasma, region=+x0 & -x1 & +y0 & -y1 & +z0 & -z1)
    fw_cell   = openmc.Cell(fill=fw_material, region=+x1 & -x2 & +y0 & -y1 & +z0 & -z1)
    str1_cell = openmc.Cell(fill=sma.materials.inconel718, region=+x2 & -x3 & +y0 & -y1 & +z0 & -z1)
    cc_cell   = openmc.Cell(fill=sma.materials.flibe, region=+x3 & -x4 & +y0 & -y1 & +z0 & -z1)
    nmul_cell = openmc.Cell(fill=sma.materials.beryllium, region=+x4 & -x5 & +y0 & -y1 & +z0 & -z1)
    str2_cell = openmc.Cell(fill=sma.materials.inconel718, region=+x5 & -x6 & +y0 & -y1 & +z0 & -z1)
    blk_cell  = openmc.Cell(fill=sma.materials.flibe, region=+x6 & -x7 & +y0 & -y1 & +z0 & -z1)

    root_univ = openmc.Universe(cells=[plas_cell, fw_cell, str1_cell, cc_cell, nmul_cell, str2_cell, blk_cell])
    geom = openmc.Geometry(root_univ)

    materials = openmc.Materials([
        sma.materials.dt_plasma,
        fw_material,
        sma.materials.inconel718,
        sma.materials.flibe,
        sma.materials.beryllium
    ])

    #############################
    #           SOURCE          #
    #############################

    source = openmc.IndependentSource()
    source.particle = "neutron"
    x_values = openmc.stats.Discrete([source_plane_position], [1.0])
    y_values = openmc.stats.Uniform(a=-delta_yz/2, b=delta_yz/2)
    z_values = openmc.stats.Uniform(a=-delta_yz/2, b=delta_yz/2)
    source.space = openmc.stats.CartesianIndependent(x=x_values, y=y_values, z=z_values)
    source.angle = openmc.stats.Monodirectional(reference_uvw=[1, 0, 0])
    source.energy = openmc.stats.muir(e0=14080000.0, m_rat=5.0, kt=20000.0)

    #############################
    #           TALLIES         #
    #############################

    energy_bins = openmc.mgxs.GROUP_STRUCTURES["CCFE-709"]
    energy_filter = openmc.EnergyFilter(energy_bins)

    fw_filter = openmc.CellFilter(fw_cell)
    cc_filter = openmc.CellFilter(cc_cell)
    blk_filter = openmc.CellFilter(blk_cell)

    neutron_filter = openmc.ParticleFilter(["neutron"])
    photon_filter = openmc.ParticleFilter(["photon"])
    electron_filter = openmc.ParticleFilter(["electron"])
    positron_filter = openmc.ParticleFilter(["positron"])

    tallies = openmc.Tallies()

    def add_tally(name, scores, filters, nuclides=None):
        tally = openmc.Tally(name=name)
        tally.scores = scores
        tally.filters = filters
        if nuclides: tally.nuclides = nuclides
        tallies.append(tally)

    # TBR tallies
    add_tally("TBR channel", ["(n,Xt)"], [cc_filter], [openmc.Nuclide("Li6"), openmc.Nuclide("Li7")])
    add_tally("TBR blanket", ["(n,Xt)"], [blk_filter], [openmc.Nuclide("Li6"), openmc.Nuclide("Li7")])
    add_tally("TBR", ["(n,Xt)"], [], [openmc.Nuclide("Li6"), openmc.Nuclide("Li7")])

    # Spectra
    add_tally("neutron_spectra_fw", ["flux"], [fw_filter, neutron_filter, energy_filter])

    # Heating tallies
    add_tally("heating_fw", ["heating"], [fw_filter])
    add_tally("heating_fw_n", ["heating"], [fw_filter, neutron_filter])
    add_tally("heating_fw_electron", ["heating"], [fw_filter, electron_filter])
    add_tally("heating_fw_photon", ["heating"], [fw_filter, photon_filter])
    add_tally("heating_fw_positron", ["heating"], [fw_filter, positron_filter])

    # DPA
    add_tally("dpa_fw", ["damage-energy"], [fw_filter])

    # Absorption
    add_tally("absorption_fw", ["absorption"], [fw_filter])

    #############################
    #          SETTINGS         #
    #############################

    settings = openmc.Settings()
    settings.batches = 100
    settings.particles = int(1e3)
    settings.inactive = 0
    settings.run_mode = "fixed source"
    settings.source = source
    settings.photon_transport = True
    settings.output = {"tallies": False}

    return openmc.Model(materials=materials, geometry=geom, settings=settings, tallies=tallies)


if __name__ == "__main__":
    model = create_model(sma.materials.tungsten)
    model.export_to_xml()
    model.run(geometry_debug=True)
