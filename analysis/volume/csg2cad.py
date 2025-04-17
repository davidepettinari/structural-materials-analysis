import geouned

geo = geouned.CsgToCad()
print("Read model")
geo.read_csg_file('geometry.xml', 'openmc_xml')
print("Start conversion to CAD")
geo.build_universe()
print("Write stp file")
geo.export_cad()