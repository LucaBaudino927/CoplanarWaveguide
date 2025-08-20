import gmsh

gmsh.initialize()
gmsh.model.add("coplanar_waveguide")

lc = 1e+2

def create_rectangle(x0, y0, width, height):
    p = [
        gmsh.model.occ.addPoint(x0, y0, 0, lc),
        gmsh.model.occ.addPoint(x0 + width, y0, 0, lc),
        gmsh.model.occ.addPoint(x0 + width, y0 + height, 0, lc),
        gmsh.model.occ.addPoint(x0, y0 + height, 0, lc)
    ]
    l = [gmsh.model.occ.addLine(p[i], p[(i+1)%4]) for i in range(4)]
    cl = gmsh.model.occ.addCurveLoop(l)
    return gmsh.model.occ.addPlaneSurface([cl])

# Surfaces
surf_kapton = create_rectangle(0, -1000, 20000, 2000)
#surf_air    = create_rectangle(0, -1000, 20000, 2000)
surf_signal = create_rectangle(0, 0, 20000, 100)
surf_gnd1   = create_rectangle(0, -200, 20000, 100)
surf_gnd2   = create_rectangle(0, 200, 20000, 100)

gmsh.model.occ.synchronize()

# Extrusions
out_kapton = gmsh.model.occ.extrude([(2, surf_kapton)], 0, 0, -100)
out_air    = gmsh.model.occ.extrude([(2, surf_kapton)], 0, 0,  100)
out_signal = gmsh.model.occ.extrude([(2, surf_signal)], 0, 0,  20)
out_gnd1   = gmsh.model.occ.extrude([(2, surf_gnd1)],   0, 0,  20)
out_gnd2   = gmsh.model.occ.extrude([(2, surf_gnd2)],   0, 0,  20)

gmsh.model.occ.synchronize()

# Physical Volumes
gmsh.model.addPhysicalGroup(3, [out_kapton[1][1]], tag=3001)
gmsh.model.setPhysicalName(3, 3001, "kapton")

gmsh.model.addPhysicalGroup(3, [out_air[1][1]], tag=3002)
gmsh.model.setPhysicalName(3, 3002, "air")

# Helper for boundary assignment
def assign_signal_or_gnd_surfaces(out, vol_tag, lateral_tag, port_in_tag, port_out_tag,
                                   vol_name, lat_name, port_in_name, port_out_name):
    vol = out[1][1]
    lat1, port1 = out[2][1], out[3][1]
    lat2, port2 = out[4][1], out[5][1]

    gmsh.model.addPhysicalGroup(3, [vol], tag=vol_tag)
    gmsh.model.setPhysicalName(3, vol_tag, vol_name)

    gmsh.model.addPhysicalGroup(2, [lat1, lat2], tag=lateral_tag)
    gmsh.model.setPhysicalName(2, lateral_tag, lat_name)

    gmsh.model.addPhysicalGroup(2, [port1], tag=port_in_tag)
    gmsh.model.setPhysicalName(2, port_in_tag, port_in_name)

    gmsh.model.addPhysicalGroup(2, [port2], tag=port_out_tag)
    gmsh.model.setPhysicalName(2, port_out_tag, port_out_name)

# Signal and GND volumes with detailed surfaces
assign_signal_or_gnd_surfaces(out_signal, 3003, 2001, 2004, 2005, "signal", "boundary_signal", "signal_port_in", "signal_port_out")
assign_signal_or_gnd_surfaces(out_gnd1,   3004, 2002, 2006, 2007, "gnd1",   "boundary_gnd1",   "gnd1_port_in",   "gnd1_port_out")
assign_signal_or_gnd_surfaces(out_gnd2,   3005, 2003, 2008, 2009, "gnd2",   "boundary_gnd2",   "gnd2_port_in",   "gnd2_port_out")

gmsh.model.mesh.generate(3)
gmsh.option.setNumber("Mesh.MshFileVersion", 2.2)
gmsh.option.setNumber("Mesh.Binary", 0)
gmsh.write("coplanar_waveguide.msh2")
gmsh.finalize()

