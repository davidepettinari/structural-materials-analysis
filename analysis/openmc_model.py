import openmc
import tokamak_radiation_environment as tre
import matplotlib.pyplot as plt
import numpy as np

############################################################################
# Define Materials
# Source: PNNL Materials Compendium April 2021
# PNNL-15870, Rev. 2
inconel718 = openmc.Material(name="Inconel 718")
inconel718.add_element("B", 0.000055, "wo")
inconel718.add_element("C", 0.000728, "wo")
inconel718.add_element("Al", 0.005000, "wo")
inconel718.add_element("Si", 0.003184, "wo")
inconel718.add_element("P", 0.000136, "wo")
inconel718.add_element("S", 0.000136, "wo")
inconel718.add_element("Ti", 0.009000, "wo")
inconel718.add_element("Cr", 0.190000, "wo")
inconel718.add_element("Mn", 0.003184, "wo")
inconel718.add_element("Fe", 0.170000, "wo")
inconel718.add_element("Ni", 0.525000, "wo")
inconel718.add_element("Co", 0.009098, "wo")
inconel718.add_element("Cu", 0.002729, "wo")
inconel718.add_element("Nb", 0.051250, "wo")
inconel718.add_element("Mo", 0.030500, "wo")
inconel718.set_density("g/cm3", 8.19)

tungsten = openmc.Material(name="Tungsten")
tungsten.add_element("W", 1.0, "wo")
tungsten.set_density("g/cm3", 19.3)

beryllium = openmc.Material(name="Beryllium")
beryllium.add_element("Be", 1.0, "wo")
beryllium.set_density("g/cm3", 1.848)

boroncarbide = openmc.Material(name="Boron Carbide")
boroncarbide.add_element("B", 0.782671, "wo")
boroncarbide.add_element("C", 0.217329, "wo")
boroncarbide.set_density("g/cm3", 2.52)

flibe = openmc.Material(name="FLiBe")
flibe.add_element("F", 4.0)
flibe.add_element("Be", 1.0)
flibe.add_nuclide("Li6", 1.8)
flibe.add_nuclide("Li7", 0.2)
flibe.set_density("g/cm3", 1.94)
flibe.temperature = 900.0

dt_plasma = openmc.Material(name="dt_plasma")
dt_plasma.add_nuclide("H2", 1.0)
dt_plasma.add_nuclide("H3", 1.0)
dt_plasma.set_density("g/cm3", 1e-5)

materials = openmc.Materials(
    [dt_plasma, flibe, boroncarbide, beryllium, tungsten, inconel718]
)

# Create a list of all nodes
plasma_nodes = (
    [443, 0],
    [435, 35],
    [417, 73],
    [392, 113],
    [356, 156],
    [312, 193],
    [295, 199],
    [281, 196],
    [272, 185],
    [260, 156],
    [240, 92],
    [223, 36],
    [217, 0],
    [223, -36],
    [240, -92],
    [260, -156],
    [272, -185],
    [281, -196],
    [295, -199],
    [312, -193],
    [356, -156],
    [392, -113],
    [417, -73],
    [435, -35]
)

fw_nodes = (
    [454, 0],
    [453, 22],
    [435, 89],
    [334, 197],
    [322, 207],
    [287, 227],
    [248, 225],
    [247, 204],
    [243, 183],
    [223, 119],
    [210, 36],
    [209, 0],
    [210, -36],
    [223, -119],
    [243, -183],
    [247, -204],
    [248, -225],
    [287, -227],
    [322, -207],
    [334, -197],
    [435, -89],
    [453, -22]  # Closing vertex to complete the hexagon
)

blanket_nodes = (
    [577, 0],
    [563, 81],
    [536, 156],
    [507, 183],
    [470, 182],
    [470, 352],
    [416, 384],
    [253, 384],
    [197, 330],
    [196, 277],
    [155, 276],
    [155, 0],
    [155, -276],
    [196, -277],
    [197, -330],
    [253, -384],
    [416, -384],
    [470, -352],
    [470, -182],
    [507, -183],
    [536, -156],
    [563, -81],
)

# geometry

# boundary surfaces
angle = (-10, 10)

# Plasma
plasma = tre.components.Plasma(
    outer_nodes=plasma_nodes, material=None, surf_offset=0.0, angle=angle
)

# Vessel

first_wall = tre.components.FirstWall(
    inner_nodes=fw_nodes, thickness=1.0, material=tungsten, angle=angle
)

vessel_inner_structure = tre.components.VesselInnerStructure(
    first_wall=first_wall, thickness=1.5, material=inconel718, angle=angle
)

vessel_cooling_channel = tre.components.VesselCoolingChannel(
    vessel_inner_structure=vessel_inner_structure,
    thickness=2.0,
    material=flibe,
    angle=angle,
)

vessel_neutron_multiplier = tre.components.VesselNeutronMultiplier(
    vessel_cooling_channel=vessel_cooling_channel,
    thickness=1.0,
    material=beryllium,
    angle=angle,
)

vessel_outer_structure = tre.components.VesselOuterStructure(
    vessel_neutron_multiplier=vessel_neutron_multiplier,
    thickness=3.5,
    material=inconel718,
    angle=angle,
)

# SOL vacuum

sol = tre.components.SOLVacuum(
    plasma=plasma, first_wall=first_wall, material=None, angle=angle
)

# Blanket

blanket = tre.components.Blanket(
    vacuum_vessel=vessel_outer_structure,
    thickness=40,
    material=flibe,
    nodes=blanket_nodes,
    angle=angle,
)


# Shield

shield = tre.components.Shield(
    blanket=blanket,
    thickness=10,
    material=boroncarbide,
    nodes=None,
    angle=angle,
)

reactor_components = [
    plasma,
    first_wall,
    vessel_inner_structure,
    vessel_cooling_channel,
    vessel_neutron_multiplier,
    vessel_outer_structure,
    sol,
    blanket,
    shield,
]

# building enclosure
enclosure_surf = openmc.Sphere(r=5000, boundary_type="vacuum")
enclosure_region = -enclosure_surf
for rc in reactor_components:
    enclosure_region = enclosure_region & ~(rc.region)
enclosure_cell = openmc.Cell(region=enclosure_region, fill=None)

root = [
    plasma.cell,
    sol.cell,
    first_wall.cell,
    vessel_inner_structure.cell,
    vessel_cooling_channel.cell,
    vessel_neutron_multiplier.cell,
    vessel_outer_structure.cell,
    blanket.cell,
    shield.cell,
    enclosure_cell,
]

geometry = openmc.Geometry(root=root)
geometry.merge_surfaces = True

# universe = openmc.Universe(cells=root)

# colors={root[0]:'magenta', root[1]:'white', root[2]:'red', root[3]:'grey', root[4]:'aqua', 
#         root[5]:'yellow', root[6]:'grey', root[7]:'aqua', root[8]:'darkred'}

# fig, (ax1, ax2) = plt.subplots(ncols=2, figsize=(14,7))
# universe.plot(origin=(0, 0, 0), basis='xz', width=[1300, 800], pixels=[600, 600], axes=ax1, 
#               color_by='cell', colors=colors, outline=True)
# universe.plot(origin=(455, 0, 0), basis='xz', width=[30, 15], pixels=[400, 400], axes=ax2, 
#               color_by='cell', colors=colors, outline=True)
# plt.show()

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
settings.particles = 100000
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

# NEUTRON SPECTRA
# FW
spectrafw_cell_tally = openmc.Tally(name="neutron_spectra_fw")
spectrafw_cell_tally.filters = [fw_filter, neutron_particle_filter, energy_filter]
spectrafw_cell_tally.scores = ["flux"]
tallies.append(spectrafw_cell_tally)

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

# DPA
# FW
damage_fw_tally = openmc.Tally(name="dpa_fw")
damage_fw_tally.scores = ["damage-energy"]
damage_fw_tally.filters = [fw_filter]
tallies.append(damage_fw_tally)

model = openmc.Model(
    materials=materials, geometry=geometry, settings=settings, tallies=tallies
)

model.run(geometry_debug=False)