import gmsh
import math
import os
import meshio
import pyvista as pv
import matplotlib.pyplot as plt
import numpy as np
from scipy.spatial import cKDTree
from mpmath import ellipk as K

"""
Funzione che importa un file geo o un msh2 generato con gmsh,
crea una mesh, usa pyvista per elaborare una sezione della linea
e calcola l'impedenza caratteristica.

- filename: nome del file geo o msh2. In caso venga fornito un file geo, il file msh2 
  verrà generato nella cartella dello script. Entrambi i file devono contenere i physical group
  "trace" e "trace2" per identificare le superfici dei conduttori, "kapton" per il substrato.
- conductor_thickness_um: spessore del conduttore in micrometri. Usato per definire
  la minima dimensione della linea per fare la dimensione della mesh minore di questo valore.
- substrate_dielectric_constant: la costante dielettrica del substrato.
- scan_coordinate_um: coordinata X in micrometri del piano della scansione.
"""
def calculate_impedance_of(
    filename = "",
    conductor_thickness_um = 20.0,
    substrate_dielectric_constant=3.3,
    scan_coordinate_um =10.0
):
    assert filename != '', "Filename must be provided"
    gmsh.initialize()
    gmsh.option.setNumber("General.Verbosity", 5)
    # Add model
    if "MeshfromImpedanceCalculator" in gmsh.model.list():
        gmsh.model.setCurrent("MeshfromImpedanceCalculator")
        gmsh.model.remove()
    gmsh.model.add("MeshfromImpedanceCalculator")

    desired_physical_groups = []
    if filename.endswith(".geo") or filename.endswith(".geo_unrolled"):
        
        gmsh.open(filename)
        gmsh.model.occ.synchronize()

        trace_total = []
        for pg in gmsh.model.getPhysicalGroups():
            dim, tag = pg
            name = gmsh.model.getPhysicalName(dim, tag)
            print("Physical Group: dim =", dim, ", tag =", tag, ", name =", name)
            if name == "trace":
                desired_physical_groups.append(("trace", tag))
                trace_tags = gmsh.model.getEntitiesForPhysicalGroup(2, tag)
                trace_total.extend(trace_tags)
            if name == "trace2":
                desired_physical_groups.append(("trace2", tag))
                extruded_trace_tags = gmsh.model.getEntitiesForPhysicalGroup(2, tag)
                trace_total.extend(extruded_trace_tags)
            if name == "kapton":
                desired_physical_groups.append(("kapton", tag))

        # Generate mesh
        gmsh.option.setNumber("Mesh.MeshSizeMin", conductor_thickness_um / 2.0)
        gmsh.option.setNumber("Mesh.MeshSizeMax", conductor_thickness_um * 10.0)
        gmsh.option.setNumber("Mesh.MeshSizeFromPoints", 0)
        gmsh.option.setNumber("Mesh.MeshSizeFromCurvature", 0)
        gmsh.option.setNumber("Mesh.MeshSizeExtendFromBoundary", 0)

        for field in gmsh.model.mesh.field.list():
            gmsh.model.mesh.field.remove(field)
        gmsh.model.mesh.field.add("Distance", 1)
        gmsh.model.mesh.field.setNumbers(1, "SurfacesList", trace_total)
        #gmsh.model.mesh.field.setNumber(1, "Sampling", math.ceil(length_μm / l_trace))
        gmsh.model.mesh.field.setNumber(1, "Sampling", 100)

        gmsh.model.mesh.field.add("Threshold", 2)
        gmsh.model.mesh.field.setNumber(2, "InField", 1)
        gmsh.model.mesh.field.setNumber(2, "SizeMin", conductor_thickness_um / 1.1)  # ##prima l_trace
        gmsh.model.mesh.field.setNumber(2, "SizeMax", conductor_thickness_um * 10.0)  # ##prima l_farfield
        gmsh.model.mesh.field.setNumber(2, "DistMin", 1.0 * conductor_thickness_um)  # ##prima 1.0 * trace_width_μm
        gmsh.model.mesh.field.setNumber(2, "DistMax", 2.0 * conductor_thickness_um)  # ##prima 0.7 * sep_dy

        gmsh.model.mesh.field.add("MinAniso", 101)
        gmsh.model.mesh.field.setNumbers(101, "FieldsList", [2])
        gmsh.model.mesh.field.setAsBackgroundMesh(101)
        gmsh.option.setNumber("Mesh.Algorithm", 6)
        gmsh.option.setNumber("Mesh.Algorithm3D", 10)

        gmsh.model.mesh.generate(3)
        gmsh.model.mesh.setOrder(1)

        # Save mesh
        gmsh.option.setNumber("Mesh.MshFileVersion", 2.2)
        gmsh.option.setNumber("Mesh.Binary", 0)
        if filename.endswith(".geo"):
            meshFilename = filename.replace(".geo", ".msh2")
        elif filename.endswith(".geo_unrolled"):
            meshFilename = filename.replace(".geo_unrolled", ".msh2")
        directory = os.path.dirname(__file__) if '__file__' in globals() else os.getcwd()
        gmsh.write(os.path.join(directory, meshFilename))
        gmsh.finalize()
        
    elif filename.endswith(".msh2"):
        gmsh.merge(filename)
        gmsh.model.occ.synchronize()
        meshFilename = filename
        directory = os.path.dirname(__file__) if '__file__' in globals() else os.getcwd()
        for pg in gmsh.model.getPhysicalGroups():
            dim, tag = pg
            name = gmsh.model.getPhysicalName(dim, tag)
            print("Physical Group: dim =", dim, ", tag =", tag, ", name =", name)
            if name == "trace":
                desired_physical_groups.append(("trace", tag))
                trace_tags = gmsh.model.getEntitiesForPhysicalGroup(2, tag)
            if name == "trace2":
                desired_physical_groups.append(("trace2", tag))
                extruded_trace_tags = gmsh.model.getEntitiesForPhysicalGroup(2, tag)
            if name == "kapton":
                desired_physical_groups.append(("kapton", tag))

    mesh = meshio.read(os.path.join(directory, meshFilename), file_format="gmsh")
    pv_mesh = pv.wrap(mesh)   # Converts meshio object to PyVista mesh
    #pv_mesh.plot()
    single_slice = pv_mesh.slice(normal=[1, 0, 0], origin=[scan_coordinate_um/3, 0, 0])

    substrate_thickness, conductor_thickness, gap_width, trace_width = process_slice(single_slice, desired_physical_groups)
    Z0 = calculate_Z0_from_conductor_losses(
        trace_width=trace_width*1e-6,
        trace_thickness=conductor_thickness*1e-6,
        gap_width=gap_width*1e-6,
        substrate_thickness=substrate_thickness*1e-6,
        dielectric_constant=substrate_dielectric_constant
    )
    print(f'Calculated Z0: {Z0} Ohm')

    # Visualize
    p = pv.Plotter()
    p.add_mesh(pv_mesh.outline(), color='k')
    cmap = plt.get_cmap('viridis', 4)
    p.add_mesh(single_slice, cmap=cmap)
    p.show()

    return




def process_slice(single_slice, desired_physical_groups):
    
    for label, tag in desired_physical_groups:
        if label == "kapton":
            kapton_pg = tag
            # 1) Extract only the slice points belonging to the desired physical groups
            substrate_pg_pts = extract_substrate_thickness(single_slice, [kapton_pg])
            # 2) Split those PG points by min-Z, max-Z, and remaining
            substrate_pts_zmin, substrate_pts_zmax, _ = split_by_z(substrate_pg_pts)
            substrate_thickness = np.max(substrate_pts_zmax[:,2]) - np.min(substrate_pts_zmin[:,2])
            print("Extracted substrate thickness: ", substrate_thickness)
    desired_physical_groups.remove(("kapton", kapton_pg))

    # 1) Extract only the slice points belonging to the desired physical groups
    conductor_pgs = []
    for label, tag in desired_physical_groups:
        conductor_pgs.append(tag)
    pg_pts = extract_physical_points(single_slice, conductor_pgs)
    # 2) Split those PG points by min-Z, max-Z, and remaining
    pts_zmin, pts_zmax, remaining = split_by_z(pg_pts)
    conductor_thickness = np.max(pts_zmax[:,2]) - np.min(pts_zmin[:,2])
    print("Extracted conductor thickness: ", conductor_thickness)
    # 3) Group remaining points by identical Y coordinate
    groups_y = group_by_y(remaining)
    y_values = sorted(groups_y.keys())
    trace_width = y_values[1] - y_values[0]
    print("Extracted trace width: ", trace_width)
    gap_width = y_values[2] - y_values[1]
    print("Extracted gap width: ", gap_width)

    return substrate_thickness, conductor_thickness, gap_width, trace_width

def extract_substrate_thickness(mesh, physical_groups, tol=1e-9):
    """
    Returns all points belonging to the specified physical groups.
    """
    # Physical group ID *per cell*
    pg_cell = mesh["gmsh:physical"]  # size = n_cells

    # Get unique physical-group cells
    selected_cell_ids = np.where(np.isin(pg_cell, physical_groups))[0]

    # Extract connectivity for each selected cell
    cell_point_ids = []
    for cid in selected_cell_ids:
        cell = mesh.get_cell(cid)  # PyVista Cell object
        cell_point_ids.extend(cell.point_ids)

    point_ids = np.unique(cell_point_ids)
    return mesh.points[point_ids]

def extract_physical_points(mesh, physical_groups, tol=1e-9):
    """
    Returns all points belonging to the specified physical groups.
    """
    # Physical group ID *per cell*
    pg_cell = mesh["gmsh:physical"]  # size = n_cells

    # Get unique physical-group cells
    selected_cell_ids = np.where(np.isin(pg_cell, physical_groups))[0]

    # Extract connectivity for each selected cell
    cell_point_ids = []
    for cid in selected_cell_ids:
        cell = mesh.get_cell(cid)  # PyVista Cell object
        cell_point_ids.extend(cell.point_ids)

    point_ids = np.unique(cell_point_ids)
    return mesh.points[point_ids]

# ----------------------------------------------------------
# Helper: group by max/min Z and remove them
# ----------------------------------------------------------
def split_by_z(pts, tol=1e-9):
    #print("pts.shape: ", pts.shape, "pts: ", pts)
    z = pts[:, 2]
    zmin = z.min()
    zmax = z.max()

    pts_zmin = pts[np.abs(z - zmin) < tol]
    pts_zmax = pts[np.abs(z - zmax) < tol]

    # Remaining points
    mask = (np.abs(z - zmin) >= tol) & (np.abs(z - zmax) >= tol)
    remaining = pts[mask]
    return pts_zmin, pts_zmax, remaining

# ----------------------------------------------------------
# Helper: group remaining points by identical Y
# ----------------------------------------------------------
def group_by_y(pts, tol=1e-9):
    #print("pts.shape: ", pts.shape, "pts: ", pts)
    yvals = pts[:, 1]
    
    # Unique y values within tolerance
    # We cluster Y values using rounding:
    y_key = np.round(yvals / tol).astype(int)

    groups = {}
    for key, p in zip(y_key, pts):
        groups.setdefault(key, []).append(p)

    # Convert keys back to approx Y value
    result = {np.mean([p[1] for p in arr]): np.array(arr) for arr in groups.values()}
    return result

def calculate_Z0_from_conductor_losses(trace_width=90e-6, trace_thickness=20e-6, gap_width=100e-6, substrate_thickness=25e-6, dielectric_constant=3.3):
    """
    Calculate Z0 from the predicted conductor loss using following paper IPC formula and comparing it with the usual conductor loss equation.
    "A new analytical, cad-oriented model for the ohmic and radiation loesses of asymmetric coplanar waveguides in
    hybrid and monolithic mic's", G. Ghione, C.U. Naldi
    """
    a = trace_width / 2
    b = a + gap_width
    ks = a / b
    ksprime = np.sqrt(1 - ks*ks)
    k2 = np.sinh(np.pi * a / (2 * substrate_thickness)) / np.sinh(np.pi * b / (2 * substrate_thickness))
    k2prime = np.sqrt(1 - k2*k2)
    e_eff = 1 + ((dielectric_constant - 1) / 2) * (K(ksprime)/K(ks)) * (K(k2)/K(k2prime))
    #Ac = 8.686 * Rs.real * np.sqrt(dielectric_constant) / (480 * np.pi * K(ks) * K(ksprime) * (1 - ks*ks)) * (1/a * (np.pi + np.log(8 * np.pi * a * (1 - ks) / (trace_thickness * (1 + ks)))) + 
    #                                                                                            1/b * (np.pi + np.log(8 * np.pi * b * (1 - ks) / (trace_thickness * (1 + ks)))))
    #print(f'CPW impedance v1: {float((30 * np.pi) / np.sqrt(dielectric_constant) * (K(ksprime)/K(ks))):.2f} Ohm')
    Z0 = (1/trace_width) * (240 * np.pi * K(ks) * K(ksprime) * (1 - ks*ks)) / (np.sqrt(dielectric_constant) * (1/a * (np.pi + np.log(8 * np.pi * a * (1 - ks) / (trace_thickness * (1 + ks)))) + 1/b * (np.pi + np.log(8 * np.pi * b * (1 - ks) / (trace_thickness * (1 + ks))))))

    #versione approssimata, infatti non funziona
    #print(f'CPW impedance v1: {float((30 * np.pi) / np.sqrt(3.3) * (K(ksprime)/K(ks))):.2f} Ohm')
    return float(Z0)




# Example usage:
calculate_impedance_of(
    filename="cpw_lumped_mesh22.5um_port75.0um_length1000.0um.msh2",
    conductor_thickness_um=20.0,
    substrate_dielectric_constant=3.3,
    scan_coordinate_um=10.0
)
