import openmc
import numpy as np
import tokamak_radiation_environment as tre
import os 

os.environ["OMP_NUM_THREADS"] = "12"

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

materials = openmc.Materials([inconel718, tungsten, beryllium, boroncarbide, flibe, dt_plasma]) 

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

# geometry

# boundary surfaces
angle = (-10, 10)

# Plasma
plasma = tre.components.Plasma(
    outer_nodes=plasma_nodes, material=dt_plasma, surf_offset=0.0, angle=angle
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
    thickness=3.0,
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
    thickness=50,
    material=flibe,
    nodes=None,
    # nodes=blanket_nodes,
    angle=angle,
)

# Shield

shield = tre.components.Shield(
    blanket=blanket,
    thickness=30,
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
geometry.export_to_xml()

############################################################################
# Define Settings
lower_left = (100, -100, -350)
upper_right = (600, 100, 350)

vol_calc = openmc.VolumeCalculation([root[2]], 1000, lower_left, upper_right)

vol_calc.set_trigger(1e-01, 'std_dev')    

settings = openmc.Settings()
settings.particles = 100
settings.batches = 10
settings.volume_calculations = [vol_calc]
settings.run_mode = "volume"

model = openmc.Model(
    materials=materials, geometry=geometry, settings=settings
)