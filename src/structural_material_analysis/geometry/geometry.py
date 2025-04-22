import openmc
import tokamak_radiation_environment as tre
import structural_material_analysis as sma

def reactor_geometry(plasma_mat, firstwall_mat, structure_mat, breeder_mat, multiplier_mat, shield_mat):
    """Create the reactor geometry."""
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
        outer_nodes=plasma_nodes, material=plasma_mat, surf_offset=0.0, angle=angle
    )

    # Vessel

    first_wall = tre.components.FirstWall(
        inner_nodes=fw_nodes, thickness=1.0, material=firstwall_mat, angle=angle
    )

    vessel_inner_structure = tre.components.VesselInnerStructure(
        first_wall=first_wall, thickness=1.5, material=structure_mat, angle=angle
    )

    vessel_cooling_channel = tre.components.VesselCoolingChannel(
        vessel_inner_structure=vessel_inner_structure,
        thickness=2.0,
        material=breeder_mat,
        angle=angle,
    )

    vessel_neutron_multiplier = tre.components.VesselNeutronMultiplier(
        vessel_cooling_channel=vessel_cooling_channel,
        thickness=1.0,
        material=multiplier_mat,
        angle=angle,
    )

    vessel_outer_structure = tre.components.VesselOuterStructure(
        vessel_neutron_multiplier=vessel_neutron_multiplier,
        thickness=3.0,
        material=structure_mat,
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
        material=breeder_mat,
        nodes=None,
        # nodes=blanket_nodes,
        angle=angle,
    )

    # Shield

    shield = tre.components.Shield(
        blanket=blanket,
        thickness=30,
        material=shield_mat,
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
    
    return (
        root,
        geometry
    )