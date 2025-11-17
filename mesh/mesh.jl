# Palace: Coplanar Waveguides Mesh Generation
import Elliptic
using Gmsh: gmsh

"""
    function generate_cpw_wave_mesh(;
        filename::AbstractString,
        refinement::Integer       = 0,
        order::Integer            = 1,
        trace_width_μm::Real      = 30.0,
        gap_width_μm::Real        = 18.0,
        separation_width_μm::Real = 200.0,
        ground_width_μm::Real     = 800.0,
        substrate_height_μm::Real = 500.0,
        metal_height_μm::Real     = 0.0,
        remove_metal_vol::Bool    = true,
        length_μm::Real           = 4000.0,
        coax_ports::Bool          = false,
        verbose::Integer          = 5,
        gui::Bool                 = false
    )

Generate a mesh for the coplanar waveguide with wave ports using Gmsh

# Arguments

  - filename - the filename to use for the generated mesh
  - refinement - measure of how many elements to include, 0 is least
  - order - the polynomial order of the approximation, minimum 1
  - trace_width_μm - width of the coplanar waveguide trace, in μm
  - gap_width_μm - width of the coplanar waveguide gap, in μm
  - separation_width_μm - separation distance between the two waveguides, in μm
  - ground_width_μm - width of the ground plane, in μm
  - substrate_height_μm - height of the substrate, in μm
  - metal_height_μm - metal thickness, in μm
  - remove_metal_vol - for positive metal thickness, whether to remove the metal domain
  - length_μm - length of the waveguides, in μm
  - coax_ports - flag to use coaxial lumped ports instead of wave ports
  - verbose - flag to dictate the level of print to REPL, passed to Gmsh
  - gui - whether to launch the Gmsh GUI on mesh generation
"""
function generate_cpw_wave_mesh(;
    filename::AbstractString,
    refinement::Integer        = 0,
    order::Integer             = 1,
    trace_width_μm::Real      = 30.0,
    gap_width_μm::Real        = 18.0,
    boundary_distance_um::Real = 200.0,
    ground_width_μm::Real     = 800.0,
    substrate_height_μm::Real = 500.0,
    metal_height_μm::Real     = 0.0,
    remove_metal_vol::Bool     = false,
    length_μm::Real           = 4000.0,
    air_distance_um::Real      = 300.0,
    coax_ports::Bool           = false,
    verbose::Integer           = 5,
    gui::Bool                  = false
)
    @assert refinement >= 0
    @assert order > 0
    @assert trace_width_μm > 0
    @assert gap_width_μm > 0
    @assert ground_width_μm > 0
    @assert substrate_height_μm > 0
    @assert metal_height_μm >= 0
    @assert length_μm > 0

    kernel = gmsh.model.occ

    gmsh.initialize()
    gmsh.option.setNumber("General.Verbosity", verbose)

    # Add model
    if "cpw" in gmsh.model.list()
        gmsh.model.setCurrent("cpw")
        gmsh.model.remove()
    end
    gmsh.model.add("cpw")

    sep_dz = air_distance_um
    sep_dy = air_distance_um/3.0

    # Mesh parameters
    l_trace = 0.25 * trace_width_μm * (2.0^-refinement)
    l_farfield = 1.0 * sep_dz * (2.0^-refinement)

    # Chip pattern
    dy = 0.0
    n0 = kernel.addRectangle(0.0, dy, 0.0, length_μm, boundary_distance_um)
    dy += boundary_distance_um
    g1 = kernel.addRectangle(0.0, dy, 0.0, length_μm, ground_width_μm)
    dy += ground_width_μm
    n1 = kernel.addRectangle(0.0, dy, 0.0, length_μm, gap_width_μm)
    dy += gap_width_μm
    t1 = kernel.addRectangle(0.0, dy, 0.0, length_μm, trace_width_μm)
    dy += trace_width_μm
    n2 = kernel.addRectangle(0.0, dy, 0.0, length_μm, gap_width_μm)
    dy += gap_width_μm
    g2 = kernel.addRectangle(0.0, dy, 0.0, length_μm, ground_width_μm)
    dy += ground_width_μm
    n3 = kernel.addRectangle(0.0, dy, 0.0, length_μm, boundary_distance_um)
    dy += boundary_distance_um

    # Metal thickness
    metal_boundary = [g1, g2, t1]
    metal_boundary_top = typeof(metal_boundary)(undef, 0)
    metal = typeof(metal_boundary)(undef, 0)
    if metal_height_μm > 0
        metal_dimtags =
            kernel.extrude([(2, x) for x in metal_boundary], 0.0, 0.0, metal_height_μm)
        metal = [x[2] for x in filter(x -> x[1] == 3, metal_dimtags)]
        for domain in metal
            _, boundary = kernel.getSurfaceLoops(domain)
            @assert length(boundary) == 1
            append!(metal_boundary_top, first(boundary))
        end
        filter!(x -> !(x in metal_boundary), metal_boundary_top)
    end

    # Substrate
    substrate =
        kernel.addBox(0.0, 0.0, -substrate_height_μm, length_μm, dy, substrate_height_μm)
    ##
    substrate_boundary = typeof(metal_boundary)(undef, 0)
    for domain in substrate
        _, boundary = kernel.getSurfaceLoops(domain)
        @assert length(boundary) == 1
        append!(substrate_boundary, first(boundary))
    end
    filter!(x -> !(x in metal_boundary) || !(x in metal_boundary_top), substrate_boundary)
    ##

    # Exterior box
    domain =
        kernel.addBox(0.0, -sep_dy, -sep_dz, length_μm, dy + 2.0 * sep_dy, 2.0 * sep_dz)
    _, domain_boundary = kernel.getSurfaceLoops(domain)
    @assert length(domain_boundary) == 1
    domain_boundary = first(domain_boundary)

    # Ports
    cy1 = boundary_distance_um + ground_width_μm + gap_width_μm + 0.5 * trace_width_μm
    if coax_ports
        ra = 0.5 * trace_width_μm
        rb = 0.5 * trace_width_μm + gap_width_μm
        let da, db
            da = kernel.addDisk(0.0, cy1, 0.0, ra, ra, -1, [1, 0, 0], [])
            db = kernel.addDisk(0.0, cy1, 0.0, rb, rb, -1, [1, 0, 0], [])
            global p1 = last(first(first(kernel.cut((2, db), (2, da)))))
        end
        let da, db
            da = kernel.addDisk(length_μm, cy1, 0.0, ra, ra, -1, [1, 0, 0], [])
            db = kernel.addDisk(length_μm, cy1, 0.0, rb, rb, -1, [1, 0, 0], [])
            global p2 = last(first(first(kernel.cut((2, db), (2, da)))))
        end
    else
        dyp1 = 4.0 * trace_width_μm
        dyp2 = dyp1
        dzp1 = 2 * trace_width_μm #0.5 * (sep_dz + substrate_height_μm)
        dzp2 = 3 * trace_width_μm
        let pa, pb, l
            pa = kernel.addPoint(0.0, cy1 - dyp2, -dzp1)
            pb = kernel.addPoint(0.0, cy1 + dyp1, -dzp1)
            l = kernel.addLine(pa, pb)
            global p1 = first(
                filter(x -> x[1] == 2, kernel.extrude([1, l], 0.0, 0.0, dzp1 + dzp2))
            )[2]
        end
        let pa, pb, l
            pa = kernel.addPoint(length_μm, cy1 - dyp2, -dzp1)
            pb = kernel.addPoint(length_μm, cy1 + dyp1, -dzp1)
            l = kernel.addLine(pa, pb)
            global p2 = first(
                filter(x -> x[1] == 2, kernel.extrude([1, l], 0.0, 0.0, dzp1 + dzp2))
            )[2]
        end
    end
    let pa, pb, l
        pa = kernel.addPoint(0.0, -sep_dy, -sep_dz)
        pb = kernel.addPoint(0.0, dy + sep_dy, -sep_dz)
        l = kernel.addLine(pa, pb)
        global p5 =
            first(filter(x -> x[1] == 2, kernel.extrude([1, l], 0.0, 0.0, 2.0 * sep_dz)))[2]
    end
    let pa, pb, l
        pa = kernel.addPoint(length_μm, -sep_dy, -sep_dz)
        pb = kernel.addPoint(length_μm, dy + sep_dy, -sep_dz)
        l = kernel.addLine(pa, pb)
        global p6 =
            first(filter(x -> x[1] == 2, kernel.extrude([1, l], 0.0, 0.0, 2.0 * sep_dz)))[2]
    end

    # Embedding
    geom_dimtags = filter(x -> x[1] == 2 || x[1] == 3, kernel.getEntities())
    _, geom_map = kernel.fragment(geom_dimtags, [])
    kernel.synchronize()

    # Add physical groups
    metal_domains = last.(
        collect(
            Iterators.flatten(
                geom_map[findall(x -> x[1] == 3 && x[2] in metal, geom_dimtags)]
            )
        )
    )

    si_domain = last.(geom_map[findfirst(x -> x == (3, substrate), geom_dimtags)])
    @assert length(si_domain) == 1
    si_domain = first(si_domain)

    air_domain = last.(geom_map[findfirst(x -> x == (3, domain), geom_dimtags)])
    filter!(x -> !(x == si_domain || x in metal_domains), air_domain)
    @assert length(air_domain) == 1
    air_domain = first(air_domain)

    if length(metal_domains) > 0 && remove_metal_vol
        remove_dimtags = [(3, x) for x in metal_domains]
        for tag in last.(
            filter(
                x -> x[1] == 2,
                gmsh.model.getBoundary(
                    [(3, z) for z in metal_domains],
                    false,
                    false,
                    false
                )
            )
        )
            normal = gmsh.model.getNormal(tag, [0, 0])
            if abs(normal[1]) == 1.0
                push!(remove_dimtags, (2, tag))
            end
        end
        kernel.remove(remove_dimtags)
        kernel.synchronize()
        filter!.(x -> !(x in remove_dimtags), geom_map)
        empty!(metal_domains)
    end

    air_domain_group = gmsh.model.addPhysicalGroup(3, [air_domain], -1, "air")
    si_domain_group = gmsh.model.addPhysicalGroup(3, [si_domain], -1, "si")
    metal_domain_group = gmsh.model.addPhysicalGroup(3, metal_domains, -1, "metal")

    port1 = last.(geom_map[findfirst(x -> x == (2, p1), geom_dimtags)])
    port2 = last.(geom_map[findfirst(x -> x == (2, p2), geom_dimtags)])

    port1_group = gmsh.model.addPhysicalGroup(2, port1, -1, "port1")
    port2_group = gmsh.model.addPhysicalGroup(2, port2, -1, "port2")

    end1 = last.(geom_map[findfirst(x -> x == (2, p5), geom_dimtags)])
    end2 = last.(geom_map[findfirst(x -> x == (2, p6), geom_dimtags)])
    filter!(x -> !(x in port1), end1)
    filter!(x -> !(x in port2), end2)

    end1_group = gmsh.model.addPhysicalGroup(2, end1, -1, "end1")
    end2_group = gmsh.model.addPhysicalGroup(2, end2, -1, "end2")

    farfield = last.(
        collect(
            Iterators.flatten(
                geom_map[findall(x -> x[1] == 2 && x[2] in domain_boundary, geom_dimtags)]
            )
        )
    )
    filter!(
        x -> !(
            x in port1 || x in port2 || x in end1 || x in end2
        ),
        farfield
    )

    farfield_group = gmsh.model.addPhysicalGroup(2, farfield, -1, "farfield")

    trace = last.(
        collect(
            Iterators.flatten(
                geom_map[findall(x -> x[1] == 2 && x[2] in metal_boundary, geom_dimtags)]
            )
        )
    )
    gap = last.(
        collect(
            Iterators.flatten(
                geom_map[findall(x -> x[1] == 2 && x[2] in [n0, n1, n2, n3], geom_dimtags)]
            )
        )
    )

    trace_group = gmsh.model.addPhysicalGroup(2, trace, -1, "trace")
    gap_group = gmsh.model.addPhysicalGroup(2, gap, -1, "gap")

    trace_top = last.(
        collect(
            Iterators.flatten(
                geom_map[findall(
                    x -> x[1] == 2 && x[2] in metal_boundary_top,
                    geom_dimtags
                )]
            )
        )
    )
    filter!(
        x -> !(
            x in port1 || x in port2 || x in end1 || x in end2
        ),
        trace_top
    )

    trace_top_group = gmsh.model.addPhysicalGroup(2, trace_top, -1, "trace2")

    ##
    filter!(
        x -> !(
            x in port1 || x in port2 || x in end1 || x in end2
        ),
        substrate_boundary
    )
    substrate_group = gmsh.model.addPhysicalGroup(2, substrate_boundary, -1, "substrate")
    ##

    # Generate mesh
    gmsh.option.setNumber("Mesh.MeshSizeMin", l_trace)
    gmsh.option.setNumber("Mesh.MeshSizeMax", l_farfield)
    gmsh.option.setNumber("Mesh.MeshSizeFromPoints", 0)
    gmsh.option.setNumber("Mesh.MeshSizeFromCurvature", 0)
    gmsh.option.setNumber("Mesh.MeshSizeExtendFromBoundary", 0)

    gap_points = last.(
        filter(
            x -> x[1] == 0,
            gmsh.model.getBoundary([(2, z) for z in gap], false, true, true)
        )
    )
    gap_curves = last.(
        filter(
            x -> x[1] == 1,
            gmsh.model.getBoundary([(2, z) for z in gap], false, false, false)
        )
    )

    gmsh.model.mesh.field.add("Distance", 1)
    #gmsh.model.mesh.field.setNumbers(1, "PointsList", gap_points)
    #gmsh.model.mesh.field.setNumbers(1, "CurvesList", gap_curves)
    gmsh.model.mesh.field.setNumbers(1, "SurfacesList", trace)
    gmsh.model.mesh.field.setNumber(1, "Sampling", ceil(length_μm / l_trace))

    gmsh.model.mesh.field.add("Threshold", 2)
    gmsh.model.mesh.field.setNumber(2, "InField", 1)
    gmsh.model.mesh.field.setNumber(2, "SizeMin", l_trace)
    gmsh.model.mesh.field.setNumber(2, "SizeMax", l_farfield)
    gmsh.model.mesh.field.setNumber(2, "DistMin", trace_width_μm)
    gmsh.model.mesh.field.setNumber(2, "DistMax", 0.7 * sep_dz)

    gmsh.model.mesh.field.add("Min", 101)
    gmsh.model.mesh.field.setNumbers(101, "FieldsList", [2])
    gmsh.model.mesh.field.setAsBackgroundMesh(101)

    gmsh.option.setNumber("Mesh.Algorithm", 6)
    gmsh.option.setNumber("Mesh.Algorithm3D", 10)
    for tag in Iterators.flatten((gap, trace, trace_top))
        gmsh.model.mesh.setAlgorithm(2, tag, 8)
    end

    gmsh.model.mesh.generate(3)
    gmsh.model.mesh.setOrder(order)

    # Save mesh
    gmsh.option.setNumber("Mesh.MshFileVersion", 2.2)
    gmsh.option.setNumber("Mesh.Binary", 0)
    meshFilename = replace(filename, ".msh2" => string("_wave_", l_trace, "um", ".msh2"))
    gmsh.write(joinpath(@__DIR__, meshFilename))

    # Save geo file (geometry script)
    geo_file = replace(meshFilename, ".msh2" => ".geo_unrolled")
    gmsh.write(joinpath(@__DIR__, geo_file))

    # Print some information
    if verbose > 0
        println("\nFinished generating mesh. Physical group tags:")
        println("Air domain: ", air_domain_group)
        println("Si domain: ", si_domain_group)
        if length(metal_domains) > 0
            println("Metal domain: ", metal_domain_group)
        end
        println("Near-end boundaries: ", end1_group)
        println("Far-end boundaries: ", end2_group)
        println("Farfield boundaries: ", farfield_group)
        println("Metal boundaries: ", trace_group)
        if length(trace_top) > 0
            println("Extruded metal boundaries: ", trace_top_group)
        end
        println("Negative trace boundaries: ", gap_group)

        println("\nPorts:")
        println("Port 1: ", port1_group)
        println("Port 2: ", port2_group)
        println()
    end

    # Optionally launch GUI
    if gui
        gmsh.fltk.run()
    end

    return gmsh.finalize()

end




"""
    function generate_cpw_lumped_mesh(;
        filename::AbstractString,
        refinement::Integer       = 0,
        order::Integer            = 1,
        trace_width_μm::Real      = 30.0,
        gap_width_μm::Real        = 18.0,
        separation_width_μm::Real = 200.0,
        ground_width_μm::Real     = 800.0,
        substrate_height_μm::Real = 500.0,
        metal_height_μm::Real     = 0.0,
        remove_metal_vol::Bool    = true,
        length_μm::Real           = 4000.0,
        verbose::Integer          = 5,
        gui::Bool                 = false
    )

Generate a mesh for the coplanar waveguide with lumped ports using Gmsh

# Arguments

  - filename - the filename to use for the generated mesh
  - refinement - measure of how many elements to include, 0 is least
  - order - the polynomial order of the approximation, minimum 1
  - trace_width_μm - width of the coplanar waveguide trace, in μm
  - gap_width_μm - width of the coplanar waveguide gap, in μm
  - separation_width_μm - separation distance between the two waveguides, in μm
  - ground_width_μm - width of the ground plane, in μm
  - substrate_height_μm - height of the substrate, in μm
  - metal_height_μm - metal thickness, in μm
  - remove_metal_vol - for positive metal thickness, whether to remove the metal domain
  - length_μm - length of the waveguides, in μm
  - verbose - flag to dictate the level of print to REPL, passed to Gmsh
  - gui - whether to launch the Gmsh GUI on mesh generation
"""
function generate_cpw_lumped_mesh(;
    filename::AbstractString,
    refinement::Integer        = 0,
    order::Integer             = 1,
    trace_width_μm::Real       = 90.0,
    gap_width_μm::Real         = 100.0,
    boundary_distance_um::Real = 300.0,
    substrate_height_μm::Real  = 25.0,
    metal_height_μm::Real      = 20.0,
    remove_metal_vol::Bool     = false,
    length_μm::Real            = 8000.0,
    air_distance_um::Real      = 200.0,
    port_factor::Real          = 0.075,
    verbose::Integer           = 5,
    gui::Bool                  = false
)
    @assert refinement >= 0
    @assert order > 0
    @assert trace_width_μm > 0
    @assert gap_width_μm > 0
    @assert boundary_distance_um > 0
    @assert substrate_height_μm > 0
    @assert metal_height_μm >= 0
    @assert length_μm > 0

    kernel = gmsh.model.occ

    gmsh.initialize()
    gmsh.option.setNumber("General.Verbosity", verbose)

    # Add model
    if "coplanar_waveguide_driven_lumped" in gmsh.model.list()
        gmsh.model.setCurrent("coplanar_waveguide_driven_lumped")
        gmsh.model.remove()
    end
    gmsh.model.add("coplanar_waveguide_driven_lumped")

    sep_dz = air_distance_um
    sep_dy = air_distance_um/3.0
    #@assert sep_dy >= mesh_gap_factor * gap_width_μm string("sep_dy (", sep_dy, ") is lower than ", mesh_gap_factor, " gap_width_μm (", mesh_gap_factor *gap_width_μm, ")! Mesh will not connect the two adjacent conductors!!!")
    #@assert sep_dy > trace_width_μm string("sep_dy (", sep_dy, ") is lower than trace_width_μm (", trace_width_μm, ")! This is not good for mesh construction!!!")

    # Mesh parameters
    l_trace = (1.0/4.0) * trace_width_μm * (2.0^-refinement)
    l_farfield = 1.0 * sep_dz * (2.0^-refinement)

    # Chip pattern
    dy = 0.0
    n1 = kernel.addRectangle(0.0, dy, 0.0, length_μm, boundary_distance_um)
    dy += boundary_distance_um
    gnd1 = kernel.addRectangle(0.0, dy, 0.0, length_μm, trace_width_μm)
    dy += trace_width_μm
    n2 = kernel.addRectangle(0.0, dy, 0.0, length_μm, gap_width_μm)
    dy += gap_width_μm
    signal = kernel.addRectangle(0.0, dy, 0.0, length_μm, trace_width_μm)
    dy += trace_width_μm
    n3 = kernel.addRectangle(0.0, dy, 0.0, length_μm, gap_width_μm)
    dy += gap_width_μm
    gnd2 = kernel.addRectangle(0.0, dy, 0.0, length_μm, trace_width_μm)
    dy += trace_width_μm
    n4 = kernel.addRectangle(0.0, dy, 0.0, length_μm, boundary_distance_um)
    dy += boundary_distance_um

    # Metal thickness
    metal_boundary = [gnd1, gnd2, signal]
    metal = typeof(metal_boundary)(undef, 0)
    if metal_height_μm > 0
        metal_dimtags =
            kernel.extrude([(2, x) for x in metal_boundary], 0.0, 0.0, metal_height_μm)
        metal = [x[2] for x in filter(x -> x[1] == 3, metal_dimtags)]
        ##
        extruded_metal_boundary = [x[2] for x in filter(x -> x[1] == 2, metal_dimtags)]
        metal_volumes = [x[2] for x in filter(x -> x[1] == 3, metal_dimtags)]
        filter!(x -> !(x == metal_boundary), extruded_metal_boundary)
        println("\nmetal_boundary:", metal_boundary)
        println("\nmetal_dimtags:", metal_dimtags)
        println("\nextruded_metal_boundary:", extruded_metal_boundary)
        println("\nmetal_volumes:", metal_volumes)
        ##
        for domain in metal
            _, boundary = kernel.getSurfaceLoops(domain)
            @assert length(boundary) == 1
        end
    end

    # Substrate
    substrate = kernel.addBox(0.0, 0.0, -substrate_height_μm, length_μm, dy, substrate_height_μm)
        
    # Exterior box
    domain = kernel.addBox(-0.5 * sep_dy, -0.5 * sep_dy, -sep_dz, length_μm + 1.0 * sep_dy, dy + 1.0 * sep_dy, 2.0 * sep_dz)
    _, domain_boundary = kernel.getSurfaceLoops(domain)
    @assert length(domain_boundary) == 1
    domain_boundary = first(domain_boundary)

    # Ports
    dy = boundary_distance_um + trace_width_μm
    p1a = kernel.addRectangle(0.0, dy, 0.0, length_μm*port_factor, gap_width_μm)
    p2a = kernel.addRectangle(length_μm - length_μm*port_factor, dy, 0.0, length_μm*port_factor, gap_width_μm)
    dy += gap_width_μm + trace_width_μm
    p1b = kernel.addRectangle(0.0, dy, 0.0, length_μm*port_factor, gap_width_μm)
    p2b = kernel.addRectangle(length_μm - length_μm*port_factor, dy, 0.0, length_μm*port_factor, gap_width_μm)


    """
    ## test section of the volume at fixed x, evaluation of Z0 before meshing
    x_plane = length_μm*0.3
    board_width = trace_width_μm + gap_width_μm + trace_width_μm + gap_width_μm + trace_width_μm
    plane = kernel.addPlaneSurface([kernel.addRectangle(x_plane, boundary_distance_um, -sep_dz * 0.9, 2.0 * sep_dz * 0.9, board_width)])
    kernel.rotate([(2, plane)], x_plane, -0.5 * sep_dy * 0.9, -sep_dz * 0.9, 0.0, 1.0, 0.0, 3*pi/2)
    section = kernel.intersect([(3, metal_volumes[1]), (3, metal_volumes[2]), (3, metal_volumes[3]), (3, substrate)], 
                               [(2, plane)], 
                               false, false)
    

    println("ok, section: ", section)
    """


    # Embedding
    geom_dimtags = filter(x -> x[1] == 2 || x[1] == 3, kernel.getEntities())
    _, geom_map = kernel.fragment(geom_dimtags, [])
    kernel.synchronize()

    # Add physical groups
    metal_domains = last.(
        collect(
            Iterators.flatten(
                geom_map[findall(x -> x[1] == 3 && x[2] in metal, geom_dimtags)]
            )
        )
    )

    si_domain = last.(geom_map[findfirst(x -> x == (3, substrate), geom_dimtags)])
    @assert length(si_domain) == 1
    si_domain = first(si_domain)

    air_domain = last.(geom_map[findfirst(x -> x == (3, domain), geom_dimtags)])
    filter!(x -> !(x == si_domain || x in metal_domains), air_domain)
    @assert length(air_domain) == 1
    air_domain = first(air_domain)

    if length(metal_domains) > 0 && remove_metal_vol
        remove_dimtags = [(3, x) for x in metal_domains]
        println("boundaries of metal domains: ", gmsh.model.getBoundary([(3, z) for z in metal_domains],false,false,false))
        for tag in last.(
            filter(
                x -> x[1] == 2,
                gmsh.model.getBoundary(
                    [(3, z) for z in metal_domains],
                    false,
                    false,
                    false
                )
            )
        )
            normal = gmsh.model.getNormal(tag, [0, 0])
            if abs(normal[1]) == 1.0
                push!(remove_dimtags, (2, tag))
            end
            println("normal for tag ", tag, ": ", normal)
        end
        kernel.remove(remove_dimtags)
        println("removed dimtags: ", remove_dimtags)
        kernel.synchronize()
        filter!.(x -> !(x in remove_dimtags), geom_map)
        empty!(metal_domains)
    end

    air_domain_group = gmsh.model.addPhysicalGroup(3, [air_domain], -1, "air")
    si_domain_group = gmsh.model.addPhysicalGroup(3, [si_domain], -1, "kapton")
    metal_domain_group = gmsh.model.addPhysicalGroup(3, metal_domains, -1, "metal")

    farfield = last.(
        collect(
            Iterators.flatten(
                geom_map[findall(x -> x[1] == 2 && x[2] in domain_boundary, geom_dimtags)]
            )
        )
    )

    farfield_group = gmsh.model.addPhysicalGroup(2, farfield, -1, "farfield")

    port1a = last.(geom_map[findfirst(x -> x == (2, p1a), geom_dimtags)])
    port2a = last.(geom_map[findfirst(x -> x == (2, p2a), geom_dimtags)])
    port1b = last.(geom_map[findfirst(x -> x == (2, p1b), geom_dimtags)])
    port2b = last.(geom_map[findfirst(x -> x == (2, p2b), geom_dimtags)])

    port1a_group = gmsh.model.addPhysicalGroup(2, port1a, -1, "port1a")
    port2a_group = gmsh.model.addPhysicalGroup(2, port2a, -1, "port2a")
    port1b_group = gmsh.model.addPhysicalGroup(2, port1b, -1, "port1b")
    port2b_group = gmsh.model.addPhysicalGroup(2, port2b, -1, "port2b")

    trace = last.(
        collect(
            Iterators.flatten(
                geom_map[findall(x -> x[1] == 2 && x[2] in metal_boundary, geom_dimtags)]
            )
        )
    )
    ##
    extruded_trace = last.(
        collect(
            Iterators.flatten(
                geom_map[findall(x -> x[1] == 2 && x[2] in extruded_metal_boundary, geom_dimtags)]
            )
        )
    )
    trace_total = unique(vcat(trace, extruded_trace))
    println("trace:", trace)
    println("extruded_trace:", extruded_trace)
    println("trace_total:", trace_total)
    ##
    gap = last.(
        collect(
            Iterators.flatten(
                geom_map[findall(x -> x[1] == 2 && x[2] in [n2, n3], geom_dimtags)]
            )
        )
    )
    filter!(
        x -> !(
            x in port1a ||
            x in port2a ||
            x in port1b ||
            x in port2b
        ),
        gap
    )
    traceAndGap = unique(vcat(trace_total, gap))

    trace_group = gmsh.model.addPhysicalGroup(2, trace, -1, "trace")
    gap_group = gmsh.model.addPhysicalGroup(2, gap, -1, "gap")
    trace_top_group = gmsh.model.addPhysicalGroup(2, extruded_trace, -1, "trace2")

    trace_points = last.(
        filter(
            x -> x[1] == 0,
            gmsh.model.getBoundary([(2, z) for z in trace], false, true, true)
        )
    )
    trace_curves = last.(
        filter(
            x -> x[1] == 1,
            gmsh.model.getBoundary([(2, z) for z in trace], false, false, false)
        )
    )

    finer_mesh_refinement = false

    if finer_mesh_refinement

        #new attemps to further refine mesh

        #per fare dei test modifico i valori qui
        #l_trace = 2.0
        #l_farfield = 200.0
        #sep_dy = 100.0

        # Generate mesh
        gmsh.option.setNumber("Mesh.MeshSizeMin", l_trace)
        gmsh.option.setNumber("Mesh.MeshSizeMax", l_farfield)
        gmsh.option.setNumber("Mesh.MeshSizeFromPoints", 0)
        gmsh.option.setNumber("Mesh.MeshSizeFromCurvature", 0)
        gmsh.option.setNumber("Mesh.MeshSizeExtendFromBoundary", 0)

        gmsh.model.mesh.field.add("Distance", 1)
        #gmsh.model.mesh.field.setNumbers(1, "PointsList", trace_points)
        #gmsh.model.mesh.field.setNumbers(1, "CurvesList", trace_curves)
        gmsh.model.mesh.field.setNumbers(1, "SurfacesList", trace_total) ##traceAndGap
        #gmsh.model.mesh.field.setNumber(1, "Sampling", ceil(length_μm / l_trace))
        gmsh.model.mesh.field.setNumber(1, "Sampling", 100)

        gmsh.model.mesh.field.add("Constant", 2)
        gmsh.model.mesh.field.setNumber(2, "VIn", 1)
        gmsh.model.mesh.field.setNumber(2, "VOut", 0)
        gmsh.model.mesh.field.setNumbers(2, "VolumesList", metal_volumes)

        gmsh.model.mesh.field.add("MathEval", 3)
        gmsh.model.mesh.field.setString(3, "F",
            # Gmsh MathEval syntax: use F1 (distance) and F2 (volume)
            # F2*(l_trace+x*F1) refers to regions inside the volume,
            # (1-F2)*(l_trace+y*F1) to regions outside the volume
            # Se x = y non c'è differenza nello scaling della mesh tra interno e esterno.
            # Maggiore è il coefficiente moltiplicativo di F1, maggiori saranno le dimensioni delle mesh prodotte (e viceversa)
            string("F2*(", 0, "+1.4*F1) + (1-F2)*(", 0, "+1.6*F1)")
        )

        #Thr interna
        gmsh.model.mesh.field.add("Threshold", 4)
        gmsh.model.mesh.field.setNumber(4, "InField", 3)
        gmsh.model.mesh.field.setNumber(4, "SizeMin", l_trace) ##versione migliore: l_trace
        gmsh.model.mesh.field.setNumber(4, "SizeMax", 1.0 * metal_height_μm) ##versione migliore: metal_height_μm
        gmsh.model.mesh.field.setNumber(4, "DistMin", 1.0 * l_trace) ##versione migliore: 1.0 * l_trace
        gmsh.model.mesh.field.setNumber(4, "DistMax", 2.0 * metal_height_μm) ##versione migliore: 1.0 * metal_height_μm

        #Thr esterna
        gmsh.model.mesh.field.add("Threshold", 5)
        gmsh.model.mesh.field.setNumber(5, "InField", 3)
        gmsh.model.mesh.field.setNumber(5, "SizeMin", l_trace) ##versione migliore: l_trace
        gmsh.model.mesh.field.setNumber(5, "SizeMax", 572) ##versione migliore: l_farfield
        gmsh.model.mesh.field.setNumber(5, "DistMin", 80) ##versione migliore: 1.0 * trace_width_μm
        gmsh.model.mesh.field.setNumber(5, "DistMax", 200) ##versione migliore: 1.0 * sep_dy

        gmsh.model.mesh.field.add("Max", 101) ##versione migliore: Max
        gmsh.model.mesh.field.setNumbers(101, "FieldsList", [4, 5])
        gmsh.model.mesh.field.setAsBackgroundMesh(101)

        gmsh.option.setNumber("Mesh.Algorithm", 6)
        gmsh.option.setNumber("Mesh.Algorithm3D", 10)

    else

        #raw refinement

        # Generate mesh
        gmsh.option.setNumber("Mesh.MeshSizeMin", l_trace)
        gmsh.option.setNumber("Mesh.MeshSizeMax", l_farfield)
        gmsh.option.setNumber("Mesh.MeshSizeFromPoints", 0)
        gmsh.option.setNumber("Mesh.MeshSizeFromCurvature", 0)
        gmsh.option.setNumber("Mesh.MeshSizeExtendFromBoundary", 0)

        gmsh.model.mesh.field.add("Distance", 1)
        #gmsh.model.mesh.field.setNumbers(1, "PointsList", trace_points)
        #gmsh.model.mesh.field.setNumbers(1, "CurvesList", trace_curves)
        gmsh.model.mesh.field.setNumbers(1, "SurfacesList", trace_total)
        gmsh.model.mesh.field.setNumber(1, "Sampling", ceil(length_μm / l_trace))
        #gmsh.model.mesh.field.setNumber(1, "Sampling", 1000)

        #distMin = 1.0 * trace_width_μm
        #distMax = 3.0 * trace_width_μm
        #if distMax <= distMin
        #    distMin = 0.7 * distMin
        #    distMax = 1.3 * distMin
        #end
        #println("\ndistMin: ", distMin, " distMax: ", distMax)
        gmsh.model.mesh.field.add("Threshold", 2)
        gmsh.model.mesh.field.setNumber(2, "InField", 1)
        gmsh.model.mesh.field.setNumber(2, "SizeMin", l_trace) ##prima l_trace
        gmsh.model.mesh.field.setNumber(2, "SizeMax", l_farfield) ##prima l_farfield
        gmsh.model.mesh.field.setNumber(2, "DistMin", 1.0 * trace_width_μm) ##prima 1.0 * trace_width_μm
        gmsh.model.mesh.field.setNumber(2, "DistMax", 1.3 * trace_width_μm) ##prima 0.7 * sep_dy
        
        #=
        ##
        gmsh.model.mesh.field.add("Threshold", 2)
        gmsh.model.mesh.field.setNumber(2, "InField", 1)
        gmsh.model.mesh.field.setNumber(2, "SizeMin", l_trace) #0.5um
        gmsh.model.mesh.field.setNumber(2, "SizeMax", l_farfield) #200um
        gmsh.model.mesh.field.setNumber(2, "DistMin", trace_width_μm / 90) #1um
        gmsh.model.mesh.field.setNumber(2, "DistMax", trace_width_μm / 9) # 10um

        gmsh.model.mesh.field.add("Threshold", 3)
        gmsh.model.mesh.field.setNumber(3, "InField", 1)
        gmsh.model.mesh.field.setNumber(3, "SizeMin", 50 * l_trace) #25um
        gmsh.model.mesh.field.setNumber(3, "SizeMax", l_farfield)
        gmsh.model.mesh.field.setNumber(3, "DistMin", trace_width_μm / 4.5) #20um
        gmsh.model.mesh.field.setNumber(3, "DistMax", trace_width_μm) #90um
        ##
        =#

        gmsh.model.mesh.field.add("MinAniso", 101)
        gmsh.model.mesh.field.setNumbers(101, "FieldsList", [2])
        gmsh.model.mesh.field.setAsBackgroundMesh(101)

        gmsh.option.setNumber("Mesh.Algorithm", 6)
        gmsh.option.setNumber("Mesh.Algorithm3D", 10)
    end


    gmsh.model.mesh.generate(3)
    gmsh.model.mesh.setOrder(order)

    # Save mesh
    gmsh.option.setNumber("Mesh.MshFileVersion", 2.2)
    gmsh.option.setNumber("Mesh.Binary", 0)
    meshFilename = replace(filename, ".msh2" => string("_lumped_mesh", l_trace, "um_port", length_μm*port_factor, "um_length", length_μm,"um.msh2"))
    gmsh.write(joinpath(@__DIR__, meshFilename))

    # Save geo file (geometry script)
    geo_file = replace(meshFilename, ".msh2" => ".geo_unrolled")
    gmsh.write(joinpath(@__DIR__, geo_file))

    # Print some information
    if verbose > 0
        println("\nFinished generating mesh. Physical group tags:")
        println("Air domain: ", air_domain_group)
        println("Si domain: ", si_domain_group)
        if length(metal_domains) > 0
            println("Metal domain: ", metal_domain_group)
        end
        println("Farfield boundaries: ", farfield_group)
        println("Metal boundaries: ", trace_group)
        if length(trace_top_group) > 0
            println("Extruded metal boundaries: ", trace_top_group)
        end
        println("Negative trace boundaries: ", gap_group)

        println("\nMultielement lumped ports:")
        println("Port 1: ", port1a_group, ", ", port1b_group)
        println("Port 2: ", port2a_group, ", ", port2b_group)
        println()
    end

    # Optionally launch GUI
    if gui
        gmsh.fltk.run()
    end

    return gmsh.finalize()
end





"""
    function generate_cpw_double_traces_lumped_mesh(;
        filename::AbstractString,
        refinement::Integer       = 0,
        order::Integer            = 1,
        trace_width_μm::Real      = 30.0,
        gap_width_μm::Real        = 18.0,
        separation_width_μm::Real = 200.0,
        ground_width_μm::Real     = 800.0,
        substrate_height_μm::Real = 500.0,
        metal_height_μm::Real     = 0.0,
        remove_metal_vol::Bool    = true,
        length_μm::Real           = 4000.0,
        verbose::Integer          = 5,
        gui::Bool                 = false
    )

Same as generate_cpw_lumped_mesh but with one more trace and one more gnd
# Arguments

  - filename - the filename to use for the generated mesh
  - refinement - measure of how many elements to include, 0 is least
  - order - the polynomial order of the approximation, minimum 1
  - trace_width_μm - width of the coplanar waveguide trace, in μm
  - gap_width_μm - width of the coplanar waveguide gap, in μm
  - separation_width_μm - separation distance between the two waveguides, in μm
  - ground_width_μm - width of the ground plane, in μm
  - substrate_height_μm - height of the substrate, in μm
  - metal_height_μm - metal thickness, in μm
  - remove_metal_vol - for positive metal thickness, whether to remove the metal domain
  - length_μm - length of the waveguides, in μm
  - verbose - flag to dictate the level of print to REPL, passed to Gmsh
  - gui - whether to launch the Gmsh GUI on mesh generation
"""
function generate_cpw_double_traces_lumped_mesh(;
    filename::AbstractString,
    refinement::Integer        = 0,
    order::Integer             = 1,
    trace_width_μm::Real       = 90.0,
    gap_width_μm::Real         = 100.0,
    boundary_distance_um::Real = 765.0,
    substrate_height_μm::Real  = 25.0,
    metal_height_μm::Real      = 20.0,
    remove_metal_vol::Bool     = false,
    length_μm::Real            = 8000.0,
    verbose::Integer           = 5,
    gui::Bool                  = false
)
    @assert refinement >= 0
    @assert order > 0
    @assert trace_width_μm > 0
    @assert gap_width_μm > 0
    @assert boundary_distance_um > 0
    @assert substrate_height_μm > 0
    @assert metal_height_μm >= 0
    @assert length_μm > 0

    kernel = gmsh.model.occ

    gmsh.initialize()
    gmsh.option.setNumber("General.Verbosity", verbose)

    # Add model
    if "coplanar_waveguide_double_traces_driven_lumped" in gmsh.model.list()
        gmsh.model.setCurrent("coplanar_waveguide_double_traces_driven_lumped")
        gmsh.model.remove()
    end
    gmsh.model.add("coplanar_waveguide_double_traces_driven_lumped")

    sep_dz = 400.0
    sep_dy = 0.5 * sep_dz

    # Mesh parameters
    l_trace = 0.5 * trace_width_μm * (2.0^-refinement)
    l_farfield = 2.0 * air_distance_um * (2.0^-refinement)

    # Chip pattern
    dy = 0.0
    n1 = kernel.addRectangle(0.0, dy, 0.0, length_μm, boundary_distance_um)
    dy += boundary_distance_um
    gnd1 = kernel.addRectangle(0.0, dy, 0.0, length_μm, trace_width_μm)
    dy += trace_width_μm
    n2 = kernel.addRectangle(0.0, dy, 0.0, length_μm, gap_width_μm)
    dy += gap_width_μm
    signal1 = kernel.addRectangle(0.0, dy, 0.0, length_μm, trace_width_μm)
    dy += trace_width_μm
    n3 = kernel.addRectangle(0.0, dy, 0.0, length_μm, gap_width_μm)
    dy += gap_width_μm
    gnd2 = kernel.addRectangle(0.0, dy, 0.0, length_μm, trace_width_μm)
    dy += trace_width_μm
    n4 = kernel.addRectangle(0.0, dy, 0.0, length_μm, gap_width_μm)
    dy += gap_width_μm
    signal2 = kernel.addRectangle(0.0, dy, 0.0, length_μm, trace_width_μm)
    dy += trace_width_μm
    n5 = kernel.addRectangle(0.0, dy, 0.0, length_μm, gap_width_μm)
    dy += gap_width_μm
    gnd3 = kernel.addRectangle(0.0, dy, 0.0, length_μm, trace_width_μm)
    dy += trace_width_μm
    n6 = kernel.addRectangle(0.0, dy, 0.0, length_μm, boundary_distance_um)
    dy += boundary_distance_um

    # Metal thickness
    metal_boundary = [gnd1, gnd2, signal1, gnd3, signal2]
    metal = typeof(metal_boundary)(undef, 0)
    if metal_height_μm > 0
        metal_dimtags =
            kernel.extrude([(2, x) for x in metal_boundary], 0.0, 0.0, metal_height_μm)
        metal = [x[2] for x in filter(x -> x[1] == 3, metal_dimtags)]
        for domain in metal
            _, boundary = kernel.getSurfaceLoops(domain)
            @assert length(boundary) == 1
        end
        
    end

    # Substrate
    substrate = kernel.addBox(0.0, 0.0, -substrate_height_μm, length_μm, dy, substrate_height_μm)

    # Exterior box
    domain = kernel.addBox(
        -0.5 * sep_dy,
        -sep_dy,
        -sep_dz,
        length_μm + sep_dy,
        dy + 2.0 * sep_dy,
        2.0 * sep_dz
    )
    _, domain_boundary = kernel.getSurfaceLoops(domain)
    @assert length(domain_boundary) == 1
    domain_boundary = first(domain_boundary)

    # Ports
    dy = boundary_distance_um + trace_width_μm
    p1a = kernel.addRectangle(0.0, dy, 0.0, gap_width_μm, gap_width_μm)
    p2a = kernel.addRectangle(length_μm - gap_width_μm, dy, 0.0, gap_width_μm, gap_width_μm)
    dy += gap_width_μm + trace_width_μm
    p1b = kernel.addRectangle(0.0, dy, 0.0, gap_width_μm, gap_width_μm)
    p2b = kernel.addRectangle(length_μm - gap_width_μm, dy, 0.0, gap_width_μm, gap_width_μm)
    dy += gap_width_μm + trace_width_μm
    p1c = kernel.addRectangle(0.0, dy, 0.0, gap_width_μm, gap_width_μm)
    p2c = kernel.addRectangle(length_μm - gap_width_μm, dy, 0.0, gap_width_μm, gap_width_μm)
    dy += gap_width_μm + trace_width_μm
    p1d = kernel.addRectangle(0.0, dy, 0.0, gap_width_μm, gap_width_μm)
    p2d = kernel.addRectangle(length_μm - gap_width_μm, dy, 0.0, gap_width_μm, gap_width_μm)

    # Embedding
    geom_dimtags = filter(x -> x[1] == 2 || x[1] == 3, kernel.getEntities())
    _, geom_map = kernel.fragment(geom_dimtags, [])
    kernel.synchronize()

    # Add physical groups
    metal_domains = last.(
        collect(
            Iterators.flatten(
                geom_map[findall(x -> x[1] == 3 && x[2] in metal, geom_dimtags)]
            )
        )
    )

    si_domain = last.(geom_map[findfirst(x -> x == (3, substrate), geom_dimtags)])
    @assert length(si_domain) == 1
    si_domain = first(si_domain)

    air_domain = last.(geom_map[findfirst(x -> x == (3, domain), geom_dimtags)])
    filter!(x -> !(x == si_domain || x in metal_domains), air_domain)
    @assert length(air_domain) == 1
    air_domain = first(air_domain)

    if length(metal_domains) > 0 && remove_metal_vol
        remove_dimtags = [(3, x) for x in metal_domains]
        for tag in last.(
            filter(
                x -> x[1] == 2,
                gmsh.model.getBoundary(
                    [(3, z) for z in metal_domains],
                    false,
                    false,
                    false
                )
            )
        )
            normal = gmsh.model.getNormal(tag, [0, 0])
            if abs(normal[1]) == 1.0
                push!(remove_dimtags, (2, tag))
            end
        end
        kernel.remove(remove_dimtags)
        kernel.synchronize()
        filter!.(x -> !(x in remove_dimtags), geom_map)
        empty!(metal_domains)
    end

    air_domain_group = gmsh.model.addPhysicalGroup(3, [air_domain], -1, "air")
    si_domain_group = gmsh.model.addPhysicalGroup(3, [si_domain], -1, "kapton")
    metal_domain_group = gmsh.model.addPhysicalGroup(3, metal_domains, -1, "metal")

    farfield = last.(
        collect(
            Iterators.flatten(
                geom_map[findall(x -> x[1] == 2 && x[2] in domain_boundary, geom_dimtags)]
            )
        )
    )

    farfield_group = gmsh.model.addPhysicalGroup(2, farfield, -1, "farfield")

    port1a = last.(geom_map[findfirst(x -> x == (2, p1a), geom_dimtags)])
    port2a = last.(geom_map[findfirst(x -> x == (2, p2a), geom_dimtags)])
    port1b = last.(geom_map[findfirst(x -> x == (2, p1b), geom_dimtags)])
    port2b = last.(geom_map[findfirst(x -> x == (2, p2b), geom_dimtags)])
    port1c = last.(geom_map[findfirst(x -> x == (2, p1c), geom_dimtags)])
    port2c = last.(geom_map[findfirst(x -> x == (2, p2c), geom_dimtags)])
    port1d = last.(geom_map[findfirst(x -> x == (2, p1d), geom_dimtags)])
    port2d = last.(geom_map[findfirst(x -> x == (2, p2d), geom_dimtags)])

    port1a_group = gmsh.model.addPhysicalGroup(2, port1a, -1, "port1a")
    port2a_group = gmsh.model.addPhysicalGroup(2, port2a, -1, "port2a")
    port1b_group = gmsh.model.addPhysicalGroup(2, port1b, -1, "port1b")
    port2b_group = gmsh.model.addPhysicalGroup(2, port2b, -1, "port2b")
    port1c_group = gmsh.model.addPhysicalGroup(2, port1c, -1, "port1c")
    port2c_group = gmsh.model.addPhysicalGroup(2, port2c, -1, "port2c")
    port1d_group = gmsh.model.addPhysicalGroup(2, port1d, -1, "port1d")
    port2d_group = gmsh.model.addPhysicalGroup(2, port2d, -1, "port2d")

    trace = last.(
        collect(
            Iterators.flatten(
                geom_map[findall(x -> x[1] == 2 && x[2] in metal_boundary, geom_dimtags)]
            )
        )
    )
    gap = last.(
        collect(
            Iterators.flatten(
                geom_map[findall(x -> x[1] == 2 && x[2] in [n1, n2, n3, n4, n5, n6], geom_dimtags)]
            )
        )
    )
    filter!(
        x -> !(
            x in port1a ||
            x in port2a ||
            x in port1b ||
            x in port2b ||
            x in port1c ||
            x in port2c ||
            x in port1d ||
            x in port2d
        ),
        gap
    )

    trace_group = gmsh.model.addPhysicalGroup(2, trace, -1, "trace")
    gap_group = gmsh.model.addPhysicalGroup(2, gap, -1, "gap")


    # Generate mesh
    gmsh.option.setNumber("Mesh.MeshSizeMin", l_trace)
    gmsh.option.setNumber("Mesh.MeshSizeMax", l_farfield)
    gmsh.option.setNumber("Mesh.MeshSizeFromPoints", 0)
    gmsh.option.setNumber("Mesh.MeshSizeFromCurvature", 0)
    gmsh.option.setNumber("Mesh.MeshSizeExtendFromBoundary", 0)

    trace_points = last.(
        filter(
            x -> x[1] == 0,
            gmsh.model.getBoundary([(2, z) for z in trace], false, true, true)
        )
    )
    trace_curves = last.(
        filter(
            x -> x[1] == 1,
            gmsh.model.getBoundary([(2, z) for z in trace], false, false, false)
        )
    )

    gmsh.model.mesh.field.add("Distance", 1)
    gmsh.model.mesh.field.setNumbers(1, "PointsList", trace_points)
    gmsh.model.mesh.field.setNumbers(1, "CurvesList", trace_curves)
    gmsh.model.mesh.field.setNumbers(1, "SurfacesList", trace)
    #gmsh.model.mesh.field.setNumber(1, "Sampling", ceil(length_μm / l_trace))
    gmsh.model.mesh.field.setNumber(1, "Sampling", 100)

    gmsh.model.mesh.field.add("Threshold", 2)
    gmsh.model.mesh.field.setNumber(2, "InField", 1)
    gmsh.model.mesh.field.setNumber(2, "SizeMin", l_trace)
    gmsh.model.mesh.field.setNumber(2, "SizeMax", l_farfield)
    gmsh.model.mesh.field.setNumber(2, "DistMin", trace_width_μm)
    gmsh.model.mesh.field.setNumber(2, "DistMax", 0.9 * sep_dz)

    gmsh.model.mesh.field.add("Min", 101)
    gmsh.model.mesh.field.setNumbers(101, "FieldsList", [2])
    gmsh.model.mesh.field.setAsBackgroundMesh(101)

    gmsh.option.setNumber("Mesh.Algorithm", 6)
    gmsh.option.setNumber("Mesh.Algorithm3D", 10)

    #gmsh.model.mesh.generate(3)
    #gmsh.model.mesh.setOrder(order)

    # Save mesh
    gmsh.option.setNumber("Mesh.MshFileVersion", 2.2)
    gmsh.option.setNumber("Mesh.Binary", 0)
    gmsh.write(joinpath(@__DIR__, filename))

    # Save geo file (geometry script)
    geo_file = replace(filename, ".msh" => ".geo_unrolled")
    gmsh.write(joinpath(@__DIR__, geo_file))

    # Print some information
    if verbose > 0
        println("\nFinished generating mesh. Physical group tags:")
        println("Air domain: ", air_domain_group)
        println("Si domain: ", si_domain_group)
        if length(metal_domains) > 0
            println("Metal domain: ", metal_domain_group)
        end
        println("Farfield boundaries: ", farfield_group)
        println("Metal boundaries: ", trace_group)
        println("Negative trace boundaries: ", gap_group)

        println("\nMultielement lumped ports:")
        println("Port 1: ", port1a_group, ", ", port1b_group)
        println("Port 2: ", port2a_group, ", ", port2b_group)
        println()
    end

    # Optionally launch GUI
    if gui
        gmsh.fltk.run()
    end

    return gmsh.finalize()
end









"""
Generate a mesh for the open/short coplanar waveguide with lumped ports using Gmsh

# Arguments

  - filename - the filename to use for the generated mesh
  - refinement - measure of how many elements to include, 0 is least
  - order - the polynomial order of the approximation, minimum 1
  - trace_width_μm - width of the coplanar waveguide trace, in μm
  - gap_width_μm - width of the coplanar waveguide gap, in μm
  - separation_width_μm - separation distance between the two waveguides, in μm
  - ground_width_μm - width of the ground plane, in μm
  - substrate_height_μm - height of the substrate, in μm
  - metal_height_μm - metal thickness, in μm
  - remove_metal_vol - for positive metal thickness, whether to remove the metal domain
  - length_μm - length of the waveguides, in μm
  - open - whether to make the second port open (true) or shorted (false)
  - verbose - flag to dictate the level of print to REPL, passed to Gmsh
  - gui - whether to launch the Gmsh GUI on mesh generation
"""
function generate_open_short_cpw_lumped_mesh(;
    filename::AbstractString,
    refinement::Integer        = 0,
    order::Integer             = 1,
    trace_width_μm::Real       = 90.0,
    ground_width_μm::Real      = 300.0,
    gap_width_μm::Real         = 100.0,
    boundary_distance_um::Real = 300.0,
    substrate_height_μm::Real  = 25.0,
    metal_height_μm::Real      = 20.0,
    remove_metal_vol::Bool     = false,
    length_μm::Real            = 8000.0,
    air_distance_um::Real      = 200.0,
    port_factor::Real          = 0.075,
    open::Bool                 = true,
    verbose::Integer           = 5,
    gui::Bool                  = false
)
    @assert refinement >= 0
    @assert order > 0
    @assert trace_width_μm > 0
    @assert gap_width_μm > 0
    @assert boundary_distance_um > 0
    @assert substrate_height_μm > 0
    @assert metal_height_μm >= 0
    @assert length_μm > 0

    kernel = gmsh.model.occ

    gmsh.initialize()
    gmsh.option.setNumber("General.Verbosity", verbose)

    # Add model
    if open
        if "open_coplanar_waveguide_driven_lumped" in gmsh.model.list()
            gmsh.model.setCurrent("open_coplanar_waveguide_driven_lumped")
            gmsh.model.remove()
        end
        gmsh.model.add("open_coplanar_waveguide_driven_lumped")
    else
        if "short_coplanar_waveguide_driven_lumped" in gmsh.model.list()
            gmsh.model.setCurrent("short_coplanar_waveguide_driven_lumped")
            gmsh.model.remove()
        end
        gmsh.model.add("short_coplanar_waveguide_driven_lumped")
    end
    
    sep_dz = air_distance_um
    sep_dy = air_distance_um/2.0
    #@assert sep_dy >= mesh_gap_factor * gap_width_μm string("sep_dy (", sep_dy, ") is lower than ", mesh_gap_factor, " gap_width_μm (", mesh_gap_factor *gap_width_μm, ")! Mesh will not connect the two adjacent conductors!!!")
    #@assert sep_dy > trace_width_μm string("sep_dy (", sep_dy, ") is lower than trace_width_μm (", trace_width_μm, ")! This is not good for mesh construction!!!")

    # Mesh parameters
    l_trace = (1.0/4.0) * trace_width_μm * (2.0^-refinement)
    l_farfield = 1.0 * sep_dz * (2.0^-refinement)

    # Chip pattern
    dy = 0.0
    n1 = kernel.addRectangle(0.0, dy, 0.0, length_μm + ground_width_μm, boundary_distance_um)
    dy += boundary_distance_um
    gnd1 = kernel.addRectangle(0.0, dy, 0.0, length_μm + ground_width_μm, ground_width_μm)
    dy += ground_width_μm
    n2 = kernel.addRectangle(0.0, dy, 0.0, length_μm, gap_width_μm)
    dy += gap_width_μm
    signal = kernel.addRectangle(0.0, dy, 0.0, length_μm, trace_width_μm)
    dy += trace_width_μm
    n3 = kernel.addRectangle(0.0, dy, 0.0, length_μm, gap_width_μm)
    dy += gap_width_μm
    gnd2 = kernel.addRectangle(0.0, dy, 0.0, length_μm + ground_width_μm, ground_width_μm)
    dy += ground_width_μm
    n4 = kernel.addRectangle(0.0, dy, 0.0, length_μm + ground_width_μm, boundary_distance_um)
    dy += boundary_distance_um

    if open
        open_gap = kernel.addRectangle(length_μm, 
                                dy - boundary_distance_um - ground_width_μm - gap_width_μm - trace_width_μm - gap_width_μm, 
                                0.0, 
                                ground_width_μm,
                                gap_width_μm + trace_width_μm + gap_width_μm)
        gap_fuse = kernel.fuse([(2, n2), (2, n3)], [(2, open_gap)])
    else
        short_wire = kernel.addRectangle(length_μm, 
                                   dy - boundary_distance_um - ground_width_μm - gap_width_μm - trace_width_μm - gap_width_μm, 
                                   0.0,
                                   ground_width_μm,
                                   gap_width_μm + trace_width_μm + gap_width_μm)
        conductors_fuse = kernel.fuse([(2, gnd1), (2, gnd2), (2, signal)], [(2, short_wire)])
    end

    # Metal thickness
    metal_boundary = []
    if open
        push!(metal_boundary, gnd1, gnd2, signal)
    else
        push!(metal_boundary, conductors_fuse[1][1][2])
    end
    metal = typeof(metal_boundary)(undef, 0)
    if metal_height_μm > 0
        metal_dimtags = kernel.extrude([(2, x) for x in metal_boundary], 0.0, 0.0, metal_height_μm)
        ##
        extruded_metal_boundary = [x[2] for x in filter(x -> x[1] == 2, metal_dimtags)]
        metal_volumes = [x[2] for x in filter(x -> x[1] == 3, metal_dimtags)]
        filter!(x -> !(x == metal_boundary), extruded_metal_boundary)
        println("metal_boundary:", metal_boundary)
        println("metal_dimtags:", metal_dimtags)
        println("extruded_metal_boundary:", extruded_metal_boundary)
        println("metal_volume:", metal_volumes)
        println("metal:", metal)
        ##
        for domain in metal
            _, boundary = kernel.getSurfaceLoops(domain)
            @assert length(boundary) == 1
        end
    end

    # Substrate
    substrate = kernel.addBox(0.0, 0.0, -substrate_height_μm, length_μm + ground_width_μm, dy, substrate_height_μm)
        
    # Exterior box
    domain = 
        kernel.addBox(-0.5 * sep_dy, -0.5 * sep_dy, -substrate_height_μm-sep_dz, length_μm + ground_width_μm + 1.0 * sep_dy, dy + 1.0 * sep_dy, substrate_height_μm + metal_height_μm + 2.0 * sep_dz)
    _, domain_boundary = kernel.getSurfaceLoops(domain)
    @assert length(domain_boundary) == 1
    domain_boundary = first(domain_boundary)

    # Ports
    dy = boundary_distance_um + ground_width_μm
    p1a = kernel.addRectangle(0.0, dy, 0.0, length_μm*port_factor, gap_width_μm)
    dy += gap_width_μm + trace_width_μm
    p1b = kernel.addRectangle(0.0, dy, 0.0, length_μm*port_factor, gap_width_μm)

    # Embedding
    geom_dimtags = filter(x -> x[1] == 2 || x[1] == 3, kernel.getEntities())
    _, geom_map = kernel.fragment(geom_dimtags, [])
    kernel.synchronize()

    # Add physical groups
    metal_domains = last.(
        collect(
            Iterators.flatten(
                geom_map[findall(x -> x[1] == 3 && x[2] in metal_volumes, geom_dimtags)]
            )
        )
    )
    println("metal_domains:", metal_domains)

    si_domain = last.(geom_map[findfirst(x -> x == (3, substrate), geom_dimtags)])
    @assert length(si_domain) == 1
    si_domain = first(si_domain)

    air_domain = last.(geom_map[findfirst(x -> x == (3, domain), geom_dimtags)])
    filter!(x -> !(x == si_domain || x in metal_domains), air_domain)
    @assert length(air_domain) == 1
    air_domain = first(air_domain)

    if length(metal_domains) > 0 && remove_metal_vol
        remove_dimtags = [(3, x) for x in metal_domains]
        println("boundaries of metal domains: ", gmsh.model.getBoundary([(3, z) for z in metal_domains],false,false,false))
        for tag in last.(
            filter(
                x -> x[1] == 2,
                gmsh.model.getBoundary(
                    [(3, z) for z in metal_domains],
                    false,
                    false,
                    false
                )
            )
        )
            normal = gmsh.model.getNormal(tag, [0, 0])
            if abs(normal[1]) == 1.0
                push!(remove_dimtags, (2, tag))
            end
            println("normal for tag ", tag, ": ", normal)
        end
        kernel.remove(remove_dimtags)
        println("removed dimtags: ", remove_dimtags)
        kernel.synchronize()
        filter!.(x -> !(x in remove_dimtags), geom_map)
        empty!(metal_domains)
    end
    

    air_domain_group = gmsh.model.addPhysicalGroup(3, [air_domain], -1, "air")
    si_domain_group = gmsh.model.addPhysicalGroup(3, [si_domain], -1, "kapton")
    metal_domain_group = gmsh.model.addPhysicalGroup(3, metal_domains, -1, "metal")
    println("metal_domain_group:", metal_domain_group)

    farfield = last.(
        collect(
            Iterators.flatten(
                geom_map[findall(x -> x[1] == 2 && x[2] in domain_boundary, geom_dimtags)]
            )
        )
    )

    farfield_group = gmsh.model.addPhysicalGroup(2, farfield, -1, "farfield")

    port1a = last.(geom_map[findfirst(x -> x == (2, p1a), geom_dimtags)])
    port1b = last.(geom_map[findfirst(x -> x == (2, p1b), geom_dimtags)])

    port1a_group = gmsh.model.addPhysicalGroup(2, port1a, -1, "port1a")
    port1b_group = gmsh.model.addPhysicalGroup(2, port1b, -1, "port1b")

    trace = last.(
        collect(
            Iterators.flatten(
                geom_map[findall(x -> x[1] == 2 && x[2] in metal_boundary, geom_dimtags)]
            )
        )
    )
    ##
    extruded_trace = last.(
        collect(
            Iterators.flatten(
                geom_map[findall(x -> x[1] == 2 && x[2] in extruded_metal_boundary, geom_dimtags)]
            )
        )
    )
    trace_total = unique(vcat(trace, extruded_trace))
    println("trace:", trace)
    println("extruded_trace:", extruded_trace)
    println("trace_total:", trace_total)
    ##
    if open
        gap = last.(
            collect(
                Iterators.flatten(
                    geom_map[findall(x -> x[1] == 2 && x[2] in [gap_fuse[1][1][2]], geom_dimtags)]
                )
            )
        )
    else
        gap = last.(
            collect(
                Iterators.flatten(
                    geom_map[findall(x -> x[1] == 2 && x[2] in [n2, n3], geom_dimtags)]
                )
            )
        )
    end
    filter!(
        x -> !(
            x in port1a ||
            x in port1b
        ),
        gap
    )
    traceAndGap = unique(vcat(trace_total, gap))

    trace_group = gmsh.model.addPhysicalGroup(2, trace, -1, "trace")
    gap_group = gmsh.model.addPhysicalGroup(2, gap, -1, "gap")
    trace_top_group = gmsh.model.addPhysicalGroup(2, extruded_trace, -1, "trace2")

    trace_points = last.(
        filter(
            x -> x[1] == 0,
            gmsh.model.getBoundary([(2, z) for z in trace], false, true, true)
        )
    )
    trace_curves = last.(
        filter(
            x -> x[1] == 1,
            gmsh.model.getBoundary([(2, z) for z in trace], false, false, false)
        )
    )

    finer_mesh_refinement = false

    if finer_mesh_refinement

        #new attemps to further refine mesh

        #per fare dei test modifico i valori qui
        #l_trace = 2.0
        #l_farfield = 200.0
        #sep_dy = 100.0

        # Generate mesh
        gmsh.option.setNumber("Mesh.MeshSizeMin", l_trace)
        gmsh.option.setNumber("Mesh.MeshSizeMax", l_farfield)
        gmsh.option.setNumber("Mesh.MeshSizeFromPoints", 0)
        gmsh.option.setNumber("Mesh.MeshSizeFromCurvature", 0)
        gmsh.option.setNumber("Mesh.MeshSizeExtendFromBoundary", 0)

        gmsh.model.mesh.field.add("Distance", 1)
        #gmsh.model.mesh.field.setNumbers(1, "PointsList", trace_points)
        #gmsh.model.mesh.field.setNumbers(1, "CurvesList", trace_curves)
        gmsh.model.mesh.field.setNumbers(1, "SurfacesList", trace_total) ##traceAndGap
        #gmsh.model.mesh.field.setNumber(1, "Sampling", ceil(length_μm / l_trace))
        gmsh.model.mesh.field.setNumber(1, "Sampling", 100)

        gmsh.model.mesh.field.add("Constant", 2)
        gmsh.model.mesh.field.setNumber(2, "VIn", 1)
        gmsh.model.mesh.field.setNumber(2, "VOut", 0)
        gmsh.model.mesh.field.setNumbers(2, "VolumesList", metal_volumes)

        gmsh.model.mesh.field.add("MathEval", 3)
        gmsh.model.mesh.field.setString(3, "F",
            # Gmsh MathEval syntax: use F1 (distance) and F2 (volume)
            # F2*(l_trace+x*F1) refers to regions inside the volume,
            # (1-F2)*(l_trace+y*F1) to regions outside the volume
            # Se x = y non c'è differenza nello scaling della mesh tra interno e esterno.
            # Maggiore è il coefficiente moltiplicativo di F1, maggiori saranno le dimensioni delle mesh prodotte (e viceversa)
            string("F2*(", 0, "+1.4*F1) + (1-F2)*(", 0, "+1.6*F1)")
        )

        #Thr interna
        gmsh.model.mesh.field.add("Threshold", 4)
        gmsh.model.mesh.field.setNumber(4, "InField", 3)
        gmsh.model.mesh.field.setNumber(4, "SizeMin", l_trace) ##versione migliore: l_trace
        gmsh.model.mesh.field.setNumber(4, "SizeMax", 1.0 * metal_height_μm) ##versione migliore: metal_height_μm
        gmsh.model.mesh.field.setNumber(4, "DistMin", 1.0 * l_trace) ##versione migliore: 1.0 * l_trace
        gmsh.model.mesh.field.setNumber(4, "DistMax", 2.0 * metal_height_μm) ##versione migliore: 1.0 * metal_height_μm

        #Thr esterna
        gmsh.model.mesh.field.add("Threshold", 5)
        gmsh.model.mesh.field.setNumber(5, "InField", 3)
        gmsh.model.mesh.field.setNumber(5, "SizeMin", l_trace) ##versione migliore: l_trace
        gmsh.model.mesh.field.setNumber(5, "SizeMax", 572) ##versione migliore: l_farfield
        gmsh.model.mesh.field.setNumber(5, "DistMin", 80) ##versione migliore: 1.0 * trace_width_μm
        gmsh.model.mesh.field.setNumber(5, "DistMax", 200) ##versione migliore: 1.0 * sep_dy

        gmsh.model.mesh.field.add("Max", 101) ##versione migliore: Max
        gmsh.model.mesh.field.setNumbers(101, "FieldsList", [4, 5])
        gmsh.model.mesh.field.setAsBackgroundMesh(101)

        gmsh.option.setNumber("Mesh.Algorithm", 6)
        gmsh.option.setNumber("Mesh.Algorithm3D", 10)

    else

        #raw refinement

        # Generate mesh
        gmsh.option.setNumber("Mesh.MeshSizeMin", l_trace)
        gmsh.option.setNumber("Mesh.MeshSizeMax", l_farfield)
        gmsh.option.setNumber("Mesh.MeshSizeFromPoints", 0)
        gmsh.option.setNumber("Mesh.MeshSizeFromCurvature", 0)
        gmsh.option.setNumber("Mesh.MeshSizeExtendFromBoundary", 0)

        gmsh.model.mesh.field.add("Distance", 1)
        #gmsh.model.mesh.field.setNumbers(1, "PointsList", trace_points)
        #gmsh.model.mesh.field.setNumbers(1, "CurvesList", trace_curves)
        gmsh.model.mesh.field.setNumbers(1, "SurfacesList", trace_total)
        gmsh.model.mesh.field.setNumber(1, "Sampling", ceil(length_μm / l_trace))
        #gmsh.model.mesh.field.setNumber(1, "Sampling", 1000)

        gmsh.model.mesh.field.add("Threshold", 2)
        gmsh.model.mesh.field.setNumber(2, "InField", 1)
        gmsh.model.mesh.field.setNumber(2, "SizeMin", l_trace) ##prima l_trace
        gmsh.model.mesh.field.setNumber(2, "SizeMax", l_farfield) ##prima l_farfield
        gmsh.model.mesh.field.setNumber(2, "DistMin", 1.0 * trace_width_μm) ##prima 1.0 * trace_width_μm
        gmsh.model.mesh.field.setNumber(2, "DistMax", 1.3 * trace_width_μm) ##prima 0.7 * sep_dy

        gmsh.model.mesh.field.add("MinAniso", 101)
        gmsh.model.mesh.field.setNumbers(101, "FieldsList", [2])
        gmsh.model.mesh.field.setAsBackgroundMesh(101)

        gmsh.option.setNumber("Mesh.Algorithm", 6)
        gmsh.option.setNumber("Mesh.Algorithm3D", 10)
    end


    gmsh.model.mesh.generate(3)
    gmsh.model.mesh.setOrder(order)

    # Save mesh
    gmsh.option.setNumber("Mesh.MshFileVersion", 2.2)
    gmsh.option.setNumber("Mesh.Binary", 0)
    meshFilename = replace(filename, ".msh2" => string("_lumped_mesh", l_trace, "um_port", length_μm*port_factor, "um_length", length_μm,"um.msh2"))
    gmsh.write(joinpath(@__DIR__, meshFilename))

    # Save geo file (geometry script)
    geo_file = replace(meshFilename, ".msh2" => ".geo_unrolled")
    gmsh.write(joinpath(@__DIR__, geo_file))

    # Print some information
    if verbose > 0
        println("\nFinished generating mesh. Physical group tags:")
        println("Air domain: ", air_domain_group)
        println("Si domain: ", si_domain_group)
        if length(metal_domains) > 0
            println("Metal domain: ", metal_domain_group)
        end
        println("Farfield boundaries: ", farfield_group)
        println("Metal boundaries: ", trace_group)
        if length(trace_top_group) > 0
            println("Extruded metal boundaries: ", trace_top_group)
        end
        println("Negative trace boundaries: ", gap_group)

        println("\nMultielement lumped ports:")
        println("Port 1: ", port1a_group, ", ", port1b_group)
        println()
    end

    # Optionally launch GUI
    if gui
        gmsh.fltk.run()
    end

    return gmsh.finalize()

end






"""
Generate a mesh for the gap capacitor coplanar waveguide with lumped ports using Gmsh

# Arguments

  - filename - the filename to use for the generated mesh
  - refinement - measure of how many elements to include, 0 is least
  - order - the polynomial order of the approximation, minimum 1
  - trace_width_μm - width of the coplanar waveguide trace, in μm
  - gap_width_μm - width of the coplanar waveguide gap, in μm
  - separation_width_μm - separation distance between the two waveguides, in μm
  - ground_width_μm - width of the ground plane, in μm
  - substrate_height_μm - height of the substrate, in μm
  - metal_height_μm - metal thickness, in μm
  - remove_metal_vol - for positive metal thickness, whether to remove the metal domain
  - length_μm - length of the waveguides, in μm
  - gap_capacitor_length_um - length of the gap capacitor section, in μm
  - verbose - flag to dictate the level of print to REPL, passed to Gmsh
  - gui - whether to launch the Gmsh GUI on mesh generation
"""
function generate_gap_capacitor_cpw_lumped_mesh(;
    filename::AbstractString,
    refinement::Integer        = 0,
    order::Integer             = 1,
    trace_width_μm::Real       = 90.0,
    gap_width_μm::Real         = 100.0,
    boundary_distance_um::Real = 300.0,
    substrate_height_μm::Real  = 25.0,
    metal_height_μm::Real      = 20.0,
    remove_metal_vol::Bool     = false,
    length_μm::Real            = 8000.0,
    air_distance_um::Real      = 200.0,
    port_factor::Real          = 0.075,
    gap_capacitor_length_um::Real = 20.0,
    verbose::Integer           = 5,
    gui::Bool                  = false
)
    @assert refinement >= 0
    @assert order > 0
    @assert trace_width_μm > 0
    @assert gap_width_μm > 0
    @assert boundary_distance_um > 0
    @assert substrate_height_μm > 0
    @assert metal_height_μm >= 0
    @assert length_μm > 0

    kernel = gmsh.model.occ

    gmsh.initialize()
    gmsh.option.setNumber("General.Verbosity", verbose)

    # Add model
    if "gap_capacitor_coplanar_waveguide_driven_lumped" in gmsh.model.list()
        gmsh.model.setCurrent("gap_capacitor_coplanar_waveguide_driven_lumped")
        gmsh.model.remove()
    end
    gmsh.model.add("gap_capacitor_coplanar_waveguide_driven_lumped")

    sep_dz = air_distance_um
    sep_dy = air_distance_um/2.0
    #@assert sep_dy >= mesh_gap_factor * gap_width_μm string("sep_dy (", sep_dy, ") is lower than ", mesh_gap_factor, " gap_width_μm (", mesh_gap_factor *gap_width_μm, ")! Mesh will not connect the two adjacent conductors!!!")
    #@assert sep_dy > trace_width_μm string("sep_dy (", sep_dy, ") is lower than trace_width_μm (", trace_width_μm, ")! This is not good for mesh construction!!!")

    # Mesh parameters
    l_trace = (1.0/4.0) * trace_width_μm * (2.0^-refinement)
    l_farfield = 1.0 * sep_dz * (2.0^-refinement)

    # Chip pattern
    dy = 0.0
    n1 = kernel.addRectangle(0.0, dy, 0.0, length_μm, boundary_distance_um)
    dy += boundary_distance_um
    gnd1 = kernel.addRectangle(0.0, dy, 0.0, length_μm, trace_width_μm)
    dy += trace_width_μm
    n2 = kernel.addRectangle(0.0, dy, 0.0, length_μm, gap_width_μm)
    dy += gap_width_μm
    signal1 = kernel.addRectangle(0.0, dy, 0.0, 0.5 * length_μm - 0.5 * gap_capacitor_length_um, trace_width_μm)
    n_gap_cap = kernel.addRectangle(0.5 * length_μm - 0.5 * gap_capacitor_length_um, dy, 0.0, gap_capacitor_length_um, trace_width_μm)
    signal2 = kernel.addRectangle(0.5 * length_μm + 0.5 * gap_capacitor_length_um, dy, 0.0, 0.5 * length_μm - 0.5 * gap_capacitor_length_um, trace_width_μm)
    dy += trace_width_μm
    n3 = kernel.addRectangle(0.0, dy, 0.0, length_μm, gap_width_μm)
    dy += gap_width_μm
    gnd2 = kernel.addRectangle(0.0, dy, 0.0, length_μm, trace_width_μm)
    dy += trace_width_μm
    n4 = kernel.addRectangle(0.0, dy, 0.0, length_μm, boundary_distance_um)
    dy += boundary_distance_um

    # Metal thickness
    metal_boundary = [gnd1, gnd2, signal1, signal2]
    metal = typeof(metal_boundary)(undef, 0)
    if metal_height_μm > 0
        metal_dimtags =
            kernel.extrude([(2, x) for x in metal_boundary], 0.0, 0.0, metal_height_μm)
        metal = [x[2] for x in filter(x -> x[1] == 3, metal_dimtags)]
        ##
        extruded_metal_boundary = [x[2] for x in filter(x -> x[1] == 2, metal_dimtags)]
        metal_volumes = [x[2] for x in filter(x -> x[1] == 3, metal_dimtags)]
        filter!(x -> !(x == metal_boundary), extruded_metal_boundary)
        println("\nmetal_boundary:", metal_boundary)
        println("\nmetal_dimtags:", metal_dimtags)
        println("\nextruded_metal_boundary:", extruded_metal_boundary)
        println("\nmetal_volume:", metal_volumes)
        ##
        for domain in metal
            _, boundary = kernel.getSurfaceLoops(domain)
            @assert length(boundary) == 1
        end
    end

    # Substrate
    substrate = kernel.addBox(0.0, 0.0, -substrate_height_μm, length_μm, dy, substrate_height_μm)
        
    # Exterior box
    domain = 
        kernel.addBox(-0.5 * sep_dy, -0.5 * sep_dy, -sep_dz, length_μm + 1.0 * sep_dy, dy + 1.0 * sep_dy, 2.0 * sep_dz)
    _, domain_boundary = kernel.getSurfaceLoops(domain)
    @assert length(domain_boundary) == 1
    domain_boundary = first(domain_boundary)

    # Ports
    dy = boundary_distance_um + trace_width_μm
    p1a = kernel.addRectangle(0.0, dy, 0.0, length_μm*port_factor, gap_width_μm)
    p2a = kernel.addRectangle(length_μm - length_μm*port_factor, dy, 0.0, length_μm*port_factor, gap_width_μm)
    dy += gap_width_μm + trace_width_μm
    p1b = kernel.addRectangle(0.0, dy, 0.0, length_μm*port_factor, gap_width_μm)
    p2b = kernel.addRectangle(length_μm - length_μm*port_factor, dy, 0.0, length_μm*port_factor, gap_width_μm)

    # Embedding
    geom_dimtags = filter(x -> x[1] == 2 || x[1] == 3, kernel.getEntities())
    _, geom_map = kernel.fragment(geom_dimtags, [])
    kernel.synchronize()

    # Add physical groups
    metal_domains = last.(
        collect(
            Iterators.flatten(
                geom_map[findall(x -> x[1] == 3 && x[2] in metal, geom_dimtags)]
            )
        )
    )

    si_domain = last.(geom_map[findfirst(x -> x == (3, substrate), geom_dimtags)])
    @assert length(si_domain) == 1
    si_domain = first(si_domain)

    air_domain = last.(geom_map[findfirst(x -> x == (3, domain), geom_dimtags)])
    filter!(x -> !(x == si_domain || x in metal_domains), air_domain)
    @assert length(air_domain) == 1
    air_domain = first(air_domain)

    if length(metal_domains) > 0 && remove_metal_vol
        remove_dimtags = [(3, x) for x in metal_domains]
        for tag in last.(
            filter(
                x -> x[1] == 2,
                gmsh.model.getBoundary(
                    [(3, z) for z in metal_domains],
                    false,
                    false,
                    false
                )
            )
        )
            normal = gmsh.model.getNormal(tag, [0, 0])
            if abs(normal[1]) == 1.0
                push!(remove_dimtags, (2, tag))
            end
        end
        kernel.remove(remove_dimtags)
        kernel.synchronize()
        filter!.(x -> !(x in remove_dimtags), geom_map)
        empty!(metal_domains)
    end

    air_domain_group = gmsh.model.addPhysicalGroup(3, [air_domain], -1, "air")
    si_domain_group = gmsh.model.addPhysicalGroup(3, [si_domain], -1, "kapton")
    metal_domain_group = gmsh.model.addPhysicalGroup(3, metal_domains, -1, "metal")

    farfield = last.(
        collect(
            Iterators.flatten(
                geom_map[findall(x -> x[1] == 2 && x[2] in domain_boundary, geom_dimtags)]
            )
        )
    )

    farfield_group = gmsh.model.addPhysicalGroup(2, farfield, -1, "farfield")

    port1a = last.(geom_map[findfirst(x -> x == (2, p1a), geom_dimtags)])
    port2a = last.(geom_map[findfirst(x -> x == (2, p2a), geom_dimtags)])
    port1b = last.(geom_map[findfirst(x -> x == (2, p1b), geom_dimtags)])
    port2b = last.(geom_map[findfirst(x -> x == (2, p2b), geom_dimtags)])

    port1a_group = gmsh.model.addPhysicalGroup(2, port1a, -1, "port1a")
    port2a_group = gmsh.model.addPhysicalGroup(2, port2a, -1, "port2a")
    port1b_group = gmsh.model.addPhysicalGroup(2, port1b, -1, "port1b")
    port2b_group = gmsh.model.addPhysicalGroup(2, port2b, -1, "port2b")

    trace = last.(
        collect(
            Iterators.flatten(
                geom_map[findall(x -> x[1] == 2 && x[2] in metal_boundary, geom_dimtags)]
            )
        )
    )
    ##
    extruded_trace = last.(
        collect(
            Iterators.flatten(
                geom_map[findall(x -> x[1] == 2 && x[2] in extruded_metal_boundary, geom_dimtags)]
            )
        )
    )
    trace_total = unique(vcat(trace, extruded_trace))
    println("trace:", trace)
    println("extruded_trace:", extruded_trace)
    println("trace_total:", trace_total)
    ##
    gap = last.(
        collect(
            Iterators.flatten(
                geom_map[findall(x -> x[1] == 2 && x[2] in [n2, n3], geom_dimtags)]
            )
        )
    )
    filter!(
        x -> !(
            x in port1a ||
            x in port2a ||
            x in port1b ||
            x in port2b
        ),
        gap
    )
    traceAndGap = unique(vcat(trace_total, gap))

    trace_group = gmsh.model.addPhysicalGroup(2, trace, -1, "trace")
    gap_group = gmsh.model.addPhysicalGroup(2, gap, -1, "gap")
    trace_top_group = gmsh.model.addPhysicalGroup(2, extruded_trace, -1, "trace2")

    trace_points = last.(
        filter(
            x -> x[1] == 0,
            gmsh.model.getBoundary([(2, z) for z in trace], false, true, true)
        )
    )
    trace_curves = last.(
        filter(
            x -> x[1] == 1,
            gmsh.model.getBoundary([(2, z) for z in trace], false, false, false)
        )
    )

    finer_mesh_refinement = false

    if finer_mesh_refinement

        #new attemps to further refine mesh

        #per fare dei test modifico i valori qui
        #l_trace = 2.0
        #l_farfield = 200.0
        #sep_dy = 100.0

        # Generate mesh
        gmsh.option.setNumber("Mesh.MeshSizeMin", l_trace)
        gmsh.option.setNumber("Mesh.MeshSizeMax", l_farfield)
        gmsh.option.setNumber("Mesh.MeshSizeFromPoints", 0)
        gmsh.option.setNumber("Mesh.MeshSizeFromCurvature", 0)
        gmsh.option.setNumber("Mesh.MeshSizeExtendFromBoundary", 0)

        gmsh.model.mesh.field.add("Distance", 1)
        #gmsh.model.mesh.field.setNumbers(1, "PointsList", trace_points)
        #gmsh.model.mesh.field.setNumbers(1, "CurvesList", trace_curves)
        gmsh.model.mesh.field.setNumbers(1, "SurfacesList", trace_total) ##traceAndGap
        #gmsh.model.mesh.field.setNumber(1, "Sampling", ceil(length_μm / l_trace))
        gmsh.model.mesh.field.setNumber(1, "Sampling", 100)

        gmsh.model.mesh.field.add("Constant", 2)
        gmsh.model.mesh.field.setNumber(2, "VIn", 1)
        gmsh.model.mesh.field.setNumber(2, "VOut", 0)
        gmsh.model.mesh.field.setNumbers(2, "VolumesList", metal_volumes)

        gmsh.model.mesh.field.add("MathEval", 3)
        gmsh.model.mesh.field.setString(3, "F",
            # Gmsh MathEval syntax: use F1 (distance) and F2 (volume)
            # F2*(l_trace+x*F1) refers to regions inside the volume,
            # (1-F2)*(l_trace+y*F1) to regions outside the volume
            # Se x = y non c'è differenza nello scaling della mesh tra interno e esterno.
            # Maggiore è il coefficiente moltiplicativo di F1, maggiori saranno le dimensioni delle mesh prodotte (e viceversa)
            string("F2*(", 0, "+1.4*F1) + (1-F2)*(", 0, "+1.6*F1)")
        )

        #Thr interna
        gmsh.model.mesh.field.add("Threshold", 4)
        gmsh.model.mesh.field.setNumber(4, "InField", 3)
        gmsh.model.mesh.field.setNumber(4, "SizeMin", l_trace) ##versione migliore: l_trace
        gmsh.model.mesh.field.setNumber(4, "SizeMax", 1.0 * metal_height_μm) ##versione migliore: metal_height_μm
        gmsh.model.mesh.field.setNumber(4, "DistMin", 1.0 * l_trace) ##versione migliore: 1.0 * l_trace
        gmsh.model.mesh.field.setNumber(4, "DistMax", 2.0 * metal_height_μm) ##versione migliore: 1.0 * metal_height_μm

        #Thr esterna
        gmsh.model.mesh.field.add("Threshold", 5)
        gmsh.model.mesh.field.setNumber(5, "InField", 3)
        gmsh.model.mesh.field.setNumber(5, "SizeMin", l_trace) ##versione migliore: l_trace
        gmsh.model.mesh.field.setNumber(5, "SizeMax", 572) ##versione migliore: l_farfield
        gmsh.model.mesh.field.setNumber(5, "DistMin", 80) ##versione migliore: 1.0 * trace_width_μm
        gmsh.model.mesh.field.setNumber(5, "DistMax", 200) ##versione migliore: 1.0 * sep_dy

        gmsh.model.mesh.field.add("Max", 101) ##versione migliore: Max
        gmsh.model.mesh.field.setNumbers(101, "FieldsList", [4, 5])
        gmsh.model.mesh.field.setAsBackgroundMesh(101)

        gmsh.option.setNumber("Mesh.Algorithm", 6)
        gmsh.option.setNumber("Mesh.Algorithm3D", 10)

    else

        #raw refinement

        # Generate mesh
        gmsh.option.setNumber("Mesh.MeshSizeMin", l_trace)
        gmsh.option.setNumber("Mesh.MeshSizeMax", l_farfield)
        gmsh.option.setNumber("Mesh.MeshSizeFromPoints", 0)
        gmsh.option.setNumber("Mesh.MeshSizeFromCurvature", 0)
        gmsh.option.setNumber("Mesh.MeshSizeExtendFromBoundary", 0)

        gmsh.model.mesh.field.add("Distance", 1)
        #gmsh.model.mesh.field.setNumbers(1, "PointsList", trace_points)
        #gmsh.model.mesh.field.setNumbers(1, "CurvesList", trace_curves)
        gmsh.model.mesh.field.setNumbers(1, "SurfacesList", trace_total)
        gmsh.model.mesh.field.setNumber(1, "Sampling", ceil(length_μm / l_trace))
        #gmsh.model.mesh.field.setNumber(1, "Sampling", 1000)

        #distMin = 1.0 * trace_width_μm
        #distMax = 3.0 * trace_width_μm
        #if distMax <= distMin
        #    distMin = 0.7 * distMin
        #    distMax = 1.3 * distMin
        #end
        #println("\ndistMin: ", distMin, " distMax: ", distMax)
        gmsh.model.mesh.field.add("Threshold", 2)
        gmsh.model.mesh.field.setNumber(2, "InField", 1)
        gmsh.model.mesh.field.setNumber(2, "SizeMin", l_trace) ##prima l_trace
        gmsh.model.mesh.field.setNumber(2, "SizeMax", l_farfield) ##prima l_farfield
        gmsh.model.mesh.field.setNumber(2, "DistMin", 1.0 * trace_width_μm) ##prima 1.0 * trace_width_μm
        gmsh.model.mesh.field.setNumber(2, "DistMax", 1.3 * trace_width_μm) ##prima 0.7 * sep_dy
        
        #=
        ##
        gmsh.model.mesh.field.add("Threshold", 2)
        gmsh.model.mesh.field.setNumber(2, "InField", 1)
        gmsh.model.mesh.field.setNumber(2, "SizeMin", l_trace) #0.5um
        gmsh.model.mesh.field.setNumber(2, "SizeMax", l_farfield) #200um
        gmsh.model.mesh.field.setNumber(2, "DistMin", trace_width_μm / 90) #1um
        gmsh.model.mesh.field.setNumber(2, "DistMax", trace_width_μm / 9) # 10um

        gmsh.model.mesh.field.add("Threshold", 3)
        gmsh.model.mesh.field.setNumber(3, "InField", 1)
        gmsh.model.mesh.field.setNumber(3, "SizeMin", 50 * l_trace) #25um
        gmsh.model.mesh.field.setNumber(3, "SizeMax", l_farfield)
        gmsh.model.mesh.field.setNumber(3, "DistMin", trace_width_μm / 4.5) #20um
        gmsh.model.mesh.field.setNumber(3, "DistMax", trace_width_μm) #90um
        ##
        =#

        gmsh.model.mesh.field.add("MinAniso", 101)
        gmsh.model.mesh.field.setNumbers(101, "FieldsList", [2])
        gmsh.model.mesh.field.setAsBackgroundMesh(101)

        gmsh.option.setNumber("Mesh.Algorithm", 6)
        gmsh.option.setNumber("Mesh.Algorithm3D", 10)
    end


    gmsh.model.mesh.generate(3)
    gmsh.model.mesh.setOrder(order)

    # Save mesh
    gmsh.option.setNumber("Mesh.MshFileVersion", 2.2)
    gmsh.option.setNumber("Mesh.Binary", 0)
    meshFilename = replace(filename, ".msh2" => string("_lumped_mesh", l_trace, "um_port", length_μm*port_factor, "um_length", length_μm,"um.msh2"))
    gmsh.write(joinpath(@__DIR__, meshFilename))

    # Save geo file (geometry script)
    geo_file = replace(meshFilename, ".msh2" => ".geo_unrolled")
    gmsh.write(joinpath(@__DIR__, geo_file))

    # Print some information
    if verbose > 0
        println("\nFinished generating mesh. Physical group tags:")
        println("Air domain: ", air_domain_group)
        println("Si domain: ", si_domain_group)
        if length(metal_domains) > 0
            println("Metal domain: ", metal_domain_group)
        end
        println("Farfield boundaries: ", farfield_group)
        println("Metal boundaries: ", trace_group)
        if length(trace_top_group) > 0
            println("Extruded metal boundaries: ", trace_top_group)
        end
        println("Negative trace boundaries: ", gap_group)

        println("\nMultielement lumped ports:")
        println("Port 1: ", port1a_group, ", ", port1b_group)
        println("Port 2: ", port2a_group, ", ", port2b_group)
        println()
    end

    # Optionally launch GUI
    if gui
        gmsh.fltk.run()
    end

    return gmsh.finalize()
end





"""
Generate a mesh for the differential cpw with ground layer with lumped ports using Gmsh

# Arguments

  - filename - the filename to use for the generated mesh
  - refinement - measure of how many elements to include, 0 is least
  - order - the polynomial order of the approximation, minimum 1
  - trace_width_μm - width of the coplanar waveguide trace, in μm
  - gap_width_μm - width of the coplanar waveguide gap, in μm
  - separation_width_μm - separation distance between the two waveguides, in μm
  - ground_width_μm - width of the ground plane, in μm
  - substrate_height_μm - height of the substrate, in μm
  - metal_height_μm - metal thickness, in μm
  - remove_metal_vol - for positive metal thickness, whether to remove the metal domain
  - length_μm - length of the waveguides, in μm
  - verbose - flag to dictate the level of print to REPL, passed to Gmsh
  - gui - whether to launch the Gmsh GUI on mesh generation
"""
function generate_diff_cpw_lumped_mesh(;
    filename::AbstractString,
    refinement::Integer        = 0,
    order::Integer             = 1,
    trace_width_μm::Real       = 90.0,
    gap_width_μm::Real         = 100.0,
    boundary_distance_um::Real = 300.0,
    substrate_height_μm::Real  = 25.0,
    metal_height_μm::Real      = 20.0,
    remove_metal_vol::Bool     = false,
    length_μm::Real            = 8000.0,
    air_distance_um::Real      = 200.0,
    port_factor::Real          = 0.075,
    verbose::Integer           = 5,
    gui::Bool                  = false
)
    @assert refinement >= 0
    @assert order > 0
    @assert trace_width_μm > 0
    @assert gap_width_μm > 0
    @assert boundary_distance_um > 0
    @assert substrate_height_μm > 0
    @assert metal_height_μm >= 0
    @assert length_μm > 0

    kernel = gmsh.model.occ

    gmsh.initialize()
    gmsh.option.setNumber("General.Verbosity", verbose)

    # Add model
    if "gap_capacitor_coplanar_waveguide_driven_lumped" in gmsh.model.list()
        gmsh.model.setCurrent("gap_capacitor_coplanar_waveguide_driven_lumped")
        gmsh.model.remove()
    end
    gmsh.model.add("gap_capacitor_coplanar_waveguide_driven_lumped")

    sep_dz = air_distance_um
    sep_dy = air_distance_um/2.0
    #@assert sep_dy >= mesh_gap_factor * gap_width_μm string("sep_dy (", sep_dy, ") is lower than ", mesh_gap_factor, " gap_width_μm (", mesh_gap_factor *gap_width_μm, ")! Mesh will not connect the two adjacent conductors!!!")
    #@assert sep_dy > trace_width_μm string("sep_dy (", sep_dy, ") is lower than trace_width_μm (", trace_width_μm, ")! This is not good for mesh construction!!!")

    # Mesh parameters
    l_trace = (1.0/4.0) * trace_width_μm * (2.0^-refinement)
    l_farfield = 1.0 * sep_dz * (2.0^-refinement)

    # Chip pattern
    dy = 0.0
    n1 = kernel.addRectangle(0.0, dy, 0.0, length_μm, boundary_distance_um)
    dy += boundary_distance_um
    signal1 = kernel.addRectangle(0.0, dy, 0.0, length_μm, trace_width_μm)
    dy += trace_width_μm
    n2 = kernel.addRectangle(0.0, dy, 0.0, length_μm, gap_width_μm)
    dy += gap_width_μm
    signal2 = kernel.addRectangle(0.0, dy, 0.0, length_μm, trace_width_μm)
    dy += trace_width_μm
    n3 = kernel.addRectangle(0.0, dy, 0.0, length_μm, boundary_distance_um)
    dy += boundary_distance_um

    gnd = kernel.addRectangle(0.0, 0.0, -substrate_height_μm - metal_height_μm, length_μm, dy)

    # Metal thickness
    metal_boundary = [gnd, signal1, signal2]
    metal = typeof(metal_boundary)(undef, 0)
    if metal_height_μm > 0
        metal_dimtags =
            kernel.extrude([(2, x) for x in metal_boundary], 0.0, 0.0, metal_height_μm)
        metal = [x[2] for x in filter(x -> x[1] == 3, metal_dimtags)]
        ##
        extruded_metal_boundary = [x[2] for x in filter(x -> x[1] == 2, metal_dimtags)]
        metal_volumes = [x[2] for x in filter(x -> x[1] == 3, metal_dimtags)]
        filter!(x -> !(x == metal_boundary), extruded_metal_boundary)
        println("\nmetal_boundary:", metal_boundary)
        println("\nmetal_dimtags:", metal_dimtags)
        println("\nextruded_metal_boundary:", extruded_metal_boundary)
        println("\nmetal_volume:", metal_volumes)
        ##
        for domain in metal
            _, boundary = kernel.getSurfaceLoops(domain)
            @assert length(boundary) == 1
        end
    end

    # Substrate
    substrate = kernel.addBox(0.0, 0.0, -substrate_height_μm, length_μm, dy, substrate_height_μm)
        
    # Exterior box
    domain = 
        kernel.addBox(-0.5 * sep_dy, -0.5 * sep_dy, -sep_dz, length_μm + 1.0 * sep_dy, dy + 1.0 * sep_dy, 2.0 * sep_dz)
    _, domain_boundary = kernel.getSurfaceLoops(domain)
    @assert length(domain_boundary) == 1
    domain_boundary = first(domain_boundary)

    # Ports
    dy = boundary_distance_um
    p1 = kernel.addRectangle(0.0, dy, 0.0, substrate_height_μm, trace_width_μm)
    kernel.rotate([(2, p1)], 0, 0, 0, 0, 1, 0, pi/2)
    p2 = kernel.addRectangle(length_μm, dy, 0.0, substrate_height_μm, trace_width_μm)
    kernel.rotate([(2, p2)], length_μm, 0, 0, 0, 1, 0, pi/2)
    dy += trace_width_μm + gap_width_μm 
    p3 = kernel.addRectangle(0.0, dy, 0.0, substrate_height_μm, trace_width_μm)
    kernel.rotate([(2, p3)], 0, 0, 0, 0, 1, 0, pi/2)
    p4 = kernel.addRectangle(length_μm, dy, 0.0, substrate_height_μm, trace_width_μm)
    kernel.rotate([(2, p4)], length_μm, 0, 0, 0, 1, 0, pi/2)

    # Embedding
    geom_dimtags = filter(x -> x[1] == 2 || x[1] == 3, kernel.getEntities())
    _, geom_map = kernel.fragment(geom_dimtags, [])
    kernel.synchronize()

    # Add physical groups
    metal_domains = last.(
        collect(
            Iterators.flatten(
                geom_map[findall(x -> x[1] == 3 && x[2] in metal, geom_dimtags)]
            )
        )
    )

    si_domain = last.(geom_map[findfirst(x -> x == (3, substrate), geom_dimtags)])
    @assert length(si_domain) == 1
    si_domain = first(si_domain)

    air_domain = last.(geom_map[findfirst(x -> x == (3, domain), geom_dimtags)])
    filter!(x -> !(x == si_domain || x in metal_domains), air_domain)
    @assert length(air_domain) == 1
    air_domain = first(air_domain)

    if length(metal_domains) > 0 && remove_metal_vol
        remove_dimtags = [(3, x) for x in metal_domains]
        for tag in last.(
            filter(
                x -> x[1] == 2,
                gmsh.model.getBoundary(
                    [(3, z) for z in metal_domains],
                    false,
                    false,
                    false
                )
            )
        )
            normal = gmsh.model.getNormal(tag, [0, 0])
            if abs(normal[1]) == 1.0
                push!(remove_dimtags, (2, tag))
            end
        end
        kernel.remove(remove_dimtags)
        kernel.synchronize()
        filter!.(x -> !(x in remove_dimtags), geom_map)
        empty!(metal_domains)
    end

    air_domain_group = gmsh.model.addPhysicalGroup(3, [air_domain], -1, "air")
    si_domain_group = gmsh.model.addPhysicalGroup(3, [si_domain], -1, "kapton")
    metal_domain_group = gmsh.model.addPhysicalGroup(3, metal_domains, -1, "metal")

    farfield = last.(
        collect(
            Iterators.flatten(
                geom_map[findall(x -> x[1] == 2 && x[2] in domain_boundary, geom_dimtags)]
            )
        )
    )

    farfield_group = gmsh.model.addPhysicalGroup(2, farfield, -1, "farfield")

    port1 = last.(geom_map[findfirst(x -> x == (2, p1), geom_dimtags)])
    port2 = last.(geom_map[findfirst(x -> x == (2, p2), geom_dimtags)])
    port3 = last.(geom_map[findfirst(x -> x == (2, p3), geom_dimtags)])
    port4 = last.(geom_map[findfirst(x -> x == (2, p4), geom_dimtags)])

    port1_group = gmsh.model.addPhysicalGroup(2, port1, -1, "port1")
    port2_group = gmsh.model.addPhysicalGroup(2, port2, -1, "port2")
    port3_group = gmsh.model.addPhysicalGroup(2, port3, -1, "port3")
    port4_group = gmsh.model.addPhysicalGroup(2, port4, -1, "port4")

    trace = last.(
        collect(
            Iterators.flatten(
                geom_map[findall(x -> x[1] == 2 && x[2] in metal_boundary, geom_dimtags)]
            )
        )
    )
    ##
    extruded_trace = last.(
        collect(
            Iterators.flatten(
                geom_map[findall(x -> x[1] == 2 && x[2] in extruded_metal_boundary, geom_dimtags)]
            )
        )
    )
    trace_total = unique(vcat(trace, extruded_trace))
    println("trace:", trace)
    println("extruded_trace:", extruded_trace)
    println("trace_total:", trace_total)
    ##
    gap = last.(
        collect(
            Iterators.flatten(
                geom_map[findall(x -> x[1] == 2 && x[2] in [n2], geom_dimtags)]
            )
        )
    )
    filter!(
        x -> !(
            x in port1 ||
            x in port2 ||
            x in port3 ||
            x in port4
        ),
        gap
    )
    traceAndGap = unique(vcat(trace_total, gap))

    trace_group = gmsh.model.addPhysicalGroup(2, trace, -1, "trace")
    gap_group = gmsh.model.addPhysicalGroup(2, gap, -1, "gap")
    trace_top_group = gmsh.model.addPhysicalGroup(2, extruded_trace, -1, "trace2")

    trace_points = last.(
        filter(
            x -> x[1] == 0,
            gmsh.model.getBoundary([(2, z) for z in trace], false, true, true)
        )
    )
    trace_curves = last.(
        filter(
            x -> x[1] == 1,
            gmsh.model.getBoundary([(2, z) for z in trace], false, false, false)
        )
    )

    finer_mesh_refinement = false

    if finer_mesh_refinement

        #new attemps to further refine mesh

        #per fare dei test modifico i valori qui
        #l_trace = 2.0
        #l_farfield = 200.0
        #sep_dy = 100.0

        # Generate mesh
        gmsh.option.setNumber("Mesh.MeshSizeMin", l_trace)
        gmsh.option.setNumber("Mesh.MeshSizeMax", l_farfield)
        gmsh.option.setNumber("Mesh.MeshSizeFromPoints", 0)
        gmsh.option.setNumber("Mesh.MeshSizeFromCurvature", 0)
        gmsh.option.setNumber("Mesh.MeshSizeExtendFromBoundary", 0)

        gmsh.model.mesh.field.add("Distance", 1)
        #gmsh.model.mesh.field.setNumbers(1, "PointsList", trace_points)
        #gmsh.model.mesh.field.setNumbers(1, "CurvesList", trace_curves)
        gmsh.model.mesh.field.setNumbers(1, "SurfacesList", trace_total) ##traceAndGap
        #gmsh.model.mesh.field.setNumber(1, "Sampling", ceil(length_μm / l_trace))
        gmsh.model.mesh.field.setNumber(1, "Sampling", 100)

        gmsh.model.mesh.field.add("Constant", 2)
        gmsh.model.mesh.field.setNumber(2, "VIn", 1)
        gmsh.model.mesh.field.setNumber(2, "VOut", 0)
        gmsh.model.mesh.field.setNumbers(2, "VolumesList", metal_volumes)

        gmsh.model.mesh.field.add("MathEval", 3)
        gmsh.model.mesh.field.setString(3, "F",
            # Gmsh MathEval syntax: use F1 (distance) and F2 (volume)
            # F2*(l_trace+x*F1) refers to regions inside the volume,
            # (1-F2)*(l_trace+y*F1) to regions outside the volume
            # Se x = y non c'è differenza nello scaling della mesh tra interno e esterno.
            # Maggiore è il coefficiente moltiplicativo di F1, maggiori saranno le dimensioni delle mesh prodotte (e viceversa)
            string("F2*(", 0, "+1.4*F1) + (1-F2)*(", 0, "+1.6*F1)")
        )

        #Thr interna
        gmsh.model.mesh.field.add("Threshold", 4)
        gmsh.model.mesh.field.setNumber(4, "InField", 3)
        gmsh.model.mesh.field.setNumber(4, "SizeMin", l_trace) ##versione migliore: l_trace
        gmsh.model.mesh.field.setNumber(4, "SizeMax", 1.0 * metal_height_μm) ##versione migliore: metal_height_μm
        gmsh.model.mesh.field.setNumber(4, "DistMin", 1.0 * l_trace) ##versione migliore: 1.0 * l_trace
        gmsh.model.mesh.field.setNumber(4, "DistMax", 2.0 * metal_height_μm) ##versione migliore: 1.0 * metal_height_μm

        #Thr esterna
        gmsh.model.mesh.field.add("Threshold", 5)
        gmsh.model.mesh.field.setNumber(5, "InField", 3)
        gmsh.model.mesh.field.setNumber(5, "SizeMin", l_trace) ##versione migliore: l_trace
        gmsh.model.mesh.field.setNumber(5, "SizeMax", 572) ##versione migliore: l_farfield
        gmsh.model.mesh.field.setNumber(5, "DistMin", 80) ##versione migliore: 1.0 * trace_width_μm
        gmsh.model.mesh.field.setNumber(5, "DistMax", 200) ##versione migliore: 1.0 * sep_dy

        gmsh.model.mesh.field.add("Max", 101) ##versione migliore: Max
        gmsh.model.mesh.field.setNumbers(101, "FieldsList", [4, 5])
        gmsh.model.mesh.field.setAsBackgroundMesh(101)

        gmsh.option.setNumber("Mesh.Algorithm", 6)
        gmsh.option.setNumber("Mesh.Algorithm3D", 10)

    else

        #raw refinement

        # Generate mesh
        gmsh.option.setNumber("Mesh.MeshSizeMin", l_trace)
        gmsh.option.setNumber("Mesh.MeshSizeMax", l_farfield)
        gmsh.option.setNumber("Mesh.MeshSizeFromPoints", 0)
        gmsh.option.setNumber("Mesh.MeshSizeFromCurvature", 0)
        gmsh.option.setNumber("Mesh.MeshSizeExtendFromBoundary", 0)

        gmsh.model.mesh.field.add("Distance", 1)
        #gmsh.model.mesh.field.setNumbers(1, "PointsList", trace_points)
        #gmsh.model.mesh.field.setNumbers(1, "CurvesList", trace_curves)
        gmsh.model.mesh.field.setNumbers(1, "SurfacesList", trace_total)
        gmsh.model.mesh.field.setNumber(1, "Sampling", ceil(length_μm / l_trace))
        #gmsh.model.mesh.field.setNumber(1, "Sampling", 1000)

        #distMin = 1.0 * trace_width_μm
        #distMax = 3.0 * trace_width_μm
        #if distMax <= distMin
        #    distMin = 0.7 * distMin
        #    distMax = 1.3 * distMin
        #end
        #println("\ndistMin: ", distMin, " distMax: ", distMax)
        gmsh.model.mesh.field.add("Threshold", 2)
        gmsh.model.mesh.field.setNumber(2, "InField", 1)
        gmsh.model.mesh.field.setNumber(2, "SizeMin", l_trace) ##prima l_trace
        gmsh.model.mesh.field.setNumber(2, "SizeMax", l_farfield) ##prima l_farfield
        gmsh.model.mesh.field.setNumber(2, "DistMin", 1.0 * trace_width_μm) ##prima 1.0 * trace_width_μm
        gmsh.model.mesh.field.setNumber(2, "DistMax", 1.3 * trace_width_μm) ##prima 0.7 * sep_dy
        
        #=
        ##
        gmsh.model.mesh.field.add("Threshold", 2)
        gmsh.model.mesh.field.setNumber(2, "InField", 1)
        gmsh.model.mesh.field.setNumber(2, "SizeMin", l_trace) #0.5um
        gmsh.model.mesh.field.setNumber(2, "SizeMax", l_farfield) #200um
        gmsh.model.mesh.field.setNumber(2, "DistMin", trace_width_μm / 90) #1um
        gmsh.model.mesh.field.setNumber(2, "DistMax", trace_width_μm / 9) # 10um

        gmsh.model.mesh.field.add("Threshold", 3)
        gmsh.model.mesh.field.setNumber(3, "InField", 1)
        gmsh.model.mesh.field.setNumber(3, "SizeMin", 50 * l_trace) #25um
        gmsh.model.mesh.field.setNumber(3, "SizeMax", l_farfield)
        gmsh.model.mesh.field.setNumber(3, "DistMin", trace_width_μm / 4.5) #20um
        gmsh.model.mesh.field.setNumber(3, "DistMax", trace_width_μm) #90um
        ##
        =#

        gmsh.model.mesh.field.add("MinAniso", 101)
        gmsh.model.mesh.field.setNumbers(101, "FieldsList", [2])
        gmsh.model.mesh.field.setAsBackgroundMesh(101)

        gmsh.option.setNumber("Mesh.Algorithm", 6)
        gmsh.option.setNumber("Mesh.Algorithm3D", 10)
    end


    gmsh.model.mesh.generate(3)
    gmsh.model.mesh.setOrder(order)

    # Save mesh
    gmsh.option.setNumber("Mesh.MshFileVersion", 2.2)
    gmsh.option.setNumber("Mesh.Binary", 0)
    meshFilename = replace(filename, ".msh2" => string("_lumped_mesh", l_trace, "um_port", length_μm*port_factor, "um_length", length_μm,"um.msh2"))
    gmsh.write(joinpath(@__DIR__, meshFilename))

    # Save geo file (geometry script)
    geo_file = replace(meshFilename, ".msh2" => ".geo_unrolled")
    gmsh.write(joinpath(@__DIR__, geo_file))

    # Print some information
    if verbose > 0
        println("\nFinished generating mesh. Physical group tags:")
        println("Air domain: ", air_domain_group)
        println("Si domain: ", si_domain_group)
        if length(metal_domains) > 0
            println("Metal domain: ", metal_domain_group)
        end
        println("Farfield boundaries: ", farfield_group)
        println("Metal boundaries: ", trace_group)
        if length(trace_top_group) > 0
            println("Extruded metal boundaries: ", trace_top_group)
        end
        println("Negative trace boundaries: ", gap_group)

        println("\nMultielement lumped ports:")
        println("Port 1: ", port1_group)
        println("Port 2: ", port2_group)
        println("Port 3: ", port3_group)
        println("Port 4: ", port4_group)
        println()
    end

    # Optionally launch GUI
    if gui
        gmsh.fltk.run()
    end

    return gmsh.finalize()
end









function generate_stepOut_cpw_mesh(;
    filename::AbstractString,
    refinement::Integer        = 0,
    order::Integer             = 1,
    trace_width_μm::Real       = 90.0,
    gap_width_μm::Real         = 100.0,
    boundary_distance_um::Real = 300.0,
    substrate_height_μm::Real  = 25.0,
    metal_height_μm::Real      = 20.0,
    remove_metal_vol::Bool     = false,
    length_μm::Real            = 100.0,
    air_distance_um::Real      = 100.0,
    port_factor::Real          = 0.075,
    gnd_reduction::Real        = 0.5,
    verbose::Integer           = 5,
    gui::Bool                  = false
)

    @assert refinement >= 0
    @assert order > 0
    @assert trace_width_μm > 0
    @assert gap_width_μm > 0
    @assert boundary_distance_um > 0
    @assert substrate_height_μm > 0
    @assert metal_height_μm >= 0
    @assert length_μm > 0
    @assert gnd_reduction >= 0.0 && gnd_reduction <= 1.0

    kernel = gmsh.model.occ

    gmsh.initialize()
    gmsh.option.setNumber("General.Verbosity", verbose)

    # Add model
    if "stepOut_coplanar_waveguide_driven_lumped" in gmsh.model.list()
        gmsh.model.setCurrent("stepOut_coplanar_waveguide_driven_lumped")
        gmsh.model.remove()
    end
    gmsh.model.add("stepOut_coplanar_waveguide_driven_lumped")

    sep_dz = air_distance_um
    sep_dy = air_distance_um/3.0
    #@assert sep_dy >= mesh_gap_factor * gap_width_μm string("sep_dy (", sep_dy, ") is lower than ", mesh_gap_factor, " gap_width_μm (", mesh_gap_factor *gap_width_μm, ")! Mesh will not connect the two adjacent conductors!!!")
    #@assert sep_dy > trace_width_μm string("sep_dy (", sep_dy, ") is lower than trace_width_μm (", trace_width_μm, ")! This is not good for mesh construction!!!")

    # Mesh parameters
    l_trace = (1.0/4.0) * trace_width_μm * (2.0^-refinement)
    l_farfield = 1.0 * sep_dz * (2.0^-refinement)

    # Chip pattern
    dy = 0.0
    n1 = kernel.addRectangle(0.0, dy, 0.0, length_μm, boundary_distance_um)
    dy += boundary_distance_um
    gnd1 = kernel.addRectangle(0.0, dy, 0.0, length_μm/2, trace_width_μm)
    gnd1_stepOut = kernel.addRectangle(length_μm/2, dy, 0.0, length_μm/2, trace_width_μm * gnd_reduction)
    n_gnd1_stepOut = kernel.addRectangle(length_μm/2, dy + trace_width_μm * gnd_reduction, 0.0, length_μm/2 - length_μm*port_factor, trace_width_μm * (1 - gnd_reduction))
    gnd1_fuse = kernel.fuse([(2, gnd1)], [(2, gnd1_stepOut)])
    dy += trace_width_μm
    n2 = kernel.addRectangle(0.0, dy, 0.0, length_μm - length_μm*port_factor, gap_width_μm)
    n2_fuse = kernel.fuse([(2, n2)], [(2, n_gnd1_stepOut)])
    dy += gap_width_μm
    signal = kernel.addRectangle(0.0, dy, 0.0, length_μm, trace_width_μm)
    dy += trace_width_μm
    n3 = kernel.addRectangle(0.0, dy, 0.0, length_μm - length_μm*port_factor, gap_width_μm)
    dy += gap_width_μm
    gnd2 = kernel.addRectangle(0.0, dy, 0.0, length_μm/2, trace_width_μm)
    n_gnd2_stepOut = kernel.addRectangle(length_μm/2, dy, 0.0, length_μm/2 - length_μm*port_factor, trace_width_μm * (1 - gnd_reduction))
    n3_fuse = kernel.fuse([(2, n3)], [(2, n_gnd2_stepOut)])
    gnd2_stepOut = kernel.addRectangle(length_μm/2, dy + trace_width_μm * (1 - gnd_reduction), 0.0, length_μm/2, trace_width_μm * gnd_reduction)
    gnd2_fuse = kernel.fuse([(2, gnd2)], [(2, gnd2_stepOut)])
    dy += trace_width_μm
    n4 = kernel.addRectangle(0.0, dy, 0.0, length_μm, boundary_distance_um)
    dy += boundary_distance_um

    # Metal thickness
    metal_boundary = [gnd1_fuse[1][1][2], gnd2_fuse[1][1][2], signal]
    metal = typeof(metal_boundary)(undef, 0)
    if metal_height_μm > 0
        metal_dimtags = kernel.extrude([(2, x) for x in metal_boundary], 0.0, 0.0, metal_height_μm)
        metal = [x[2] for x in filter(x -> x[1] == 3, metal_dimtags)]
        ##
        extruded_metal_boundary = [x[2] for x in filter(x -> x[1] == 2, metal_dimtags)]
        metal_volumes = [x[2] for x in filter(x -> x[1] == 3, metal_dimtags)]
        filter!(x -> !(x == metal_boundary), extruded_metal_boundary)
        println("\nmetal_boundary:", metal_boundary)
        println("\nmetal_dimtags:", metal_dimtags)
        println("\nextruded_metal_boundary:", extruded_metal_boundary)
        println("\nmetal_volume:", metal_volumes)
        ##
        for domain in metal
            _, boundary = kernel.getSurfaceLoops(domain)
            @assert length(boundary) == 1
        end
    end

    # Substrate
    substrate = kernel.addBox(0.0, 0.0, -substrate_height_μm, length_μm, dy, substrate_height_μm)
        
    # Exterior box
    domain = 
        kernel.addBox(-0.5 * sep_dy, -0.5 * sep_dy, -sep_dz, length_μm + 1.0 * sep_dy, dy + 1.0 * sep_dy, 2.0 * sep_dz)
    _, domain_boundary = kernel.getSurfaceLoops(domain)
    @assert length(domain_boundary) == 1
    domain_boundary = first(domain_boundary)

    # Ports
    dy = boundary_distance_um + trace_width_μm
    p1a = kernel.addRectangle(0.0, dy, 0.0, length_μm*port_factor, gap_width_μm)
    dy = boundary_distance_um + (trace_width_μm * gnd_reduction)
    p2a = kernel.addRectangle(length_μm - length_μm*port_factor, dy, 0.0, length_μm*port_factor, gap_width_μm + trace_width_μm * (1 - gnd_reduction))
    dy = boundary_distance_um + trace_width_μm + gap_width_μm + trace_width_μm
    p1b = kernel.addRectangle(0.0, dy, 0.0, length_μm*port_factor, gap_width_μm)
    p2b = kernel.addRectangle(length_μm - length_μm*port_factor, dy, 0.0, length_μm*port_factor, gap_width_μm + trace_width_μm * (1 - gnd_reduction))

    # Embedding
    geom_dimtags = filter(x -> x[1] == 2 || x[1] == 3, kernel.getEntities())
    _, geom_map = kernel.fragment(geom_dimtags, [])
    kernel.synchronize()

    # Add physical groups
    metal_domains = last.(
        collect(
            Iterators.flatten(
                geom_map[findall(x -> x[1] == 3 && x[2] in metal, geom_dimtags)]
            )
        )
    )

    si_domain = last.(geom_map[findfirst(x -> x == (3, substrate), geom_dimtags)])
    @assert length(si_domain) == 1
    si_domain = first(si_domain)

    air_domain = last.(geom_map[findfirst(x -> x == (3, domain), geom_dimtags)])
    filter!(x -> !(x == si_domain || x in metal_domains), air_domain)
    @assert length(air_domain) == 1
    air_domain = first(air_domain)

    if length(metal_domains) > 0 && remove_metal_vol
        remove_dimtags = [(3, x) for x in metal_domains]
        for tag in last.(
            filter(
                x -> x[1] == 2,
                gmsh.model.getBoundary(
                    [(3, z) for z in metal_domains],
                    false,
                    false,
                    false
                )
            )
        )
            normal = gmsh.model.getNormal(tag, [0, 0])
            if abs(normal[1]) == 1.0
                push!(remove_dimtags, (2, tag))
            end
        end
        kernel.remove(remove_dimtags)
        kernel.synchronize()
        filter!.(x -> !(x in remove_dimtags), geom_map)
        empty!(metal_domains)
    end

    air_domain_group = gmsh.model.addPhysicalGroup(3, [air_domain], -1, "air")
    si_domain_group = gmsh.model.addPhysicalGroup(3, [si_domain], -1, "kapton")
    metal_domain_group = gmsh.model.addPhysicalGroup(3, metal_domains, -1, "metal")

    farfield = last.(
        collect(
            Iterators.flatten(
                geom_map[findall(x -> x[1] == 2 && x[2] in domain_boundary, geom_dimtags)]
            )
        )
    )

    farfield_group = gmsh.model.addPhysicalGroup(2, farfield, -1, "farfield")

    port1a = last.(geom_map[findfirst(x -> x == (2, p1a), geom_dimtags)])
    port2a = last.(geom_map[findfirst(x -> x == (2, p2a), geom_dimtags)])
    port1b = last.(geom_map[findfirst(x -> x == (2, p1b), geom_dimtags)])
    port2b = last.(geom_map[findfirst(x -> x == (2, p2b), geom_dimtags)])

    port1a_group = gmsh.model.addPhysicalGroup(2, port1a, -1, "port1a")
    port2a_group = gmsh.model.addPhysicalGroup(2, port2a, -1, "port2a")
    port1b_group = gmsh.model.addPhysicalGroup(2, port1b, -1, "port1b")
    port2b_group = gmsh.model.addPhysicalGroup(2, port2b, -1, "port2b")

    trace = last.(
        collect(
            Iterators.flatten(
                geom_map[findall(x -> x[1] == 2 && x[2] in metal_boundary, geom_dimtags)]
            )
        )
    )
    ##
    extruded_trace = last.(
        collect(
            Iterators.flatten(
                geom_map[findall(x -> x[1] == 2 && x[2] in extruded_metal_boundary, geom_dimtags)]
            )
        )
    )
    trace_total = unique(vcat(trace, extruded_trace))
    println("trace:", trace)
    println("extruded_trace:", extruded_trace)
    println("trace_total:", trace_total)
    ##
    gap = last.(
        collect(
            Iterators.flatten(
                geom_map[findall(x -> x[1] == 2 && x[2] in [n2, n3], geom_dimtags)]
            )
        )
    )
    filter!(
        x -> !(
            x in port1a ||
            x in port2a ||
            x in port1b ||
            x in port2b
        ),
        gap
    )
    traceAndGap = unique(vcat(trace_total, gap))

    trace_group = gmsh.model.addPhysicalGroup(2, trace, -1, "trace")
    gap_group = gmsh.model.addPhysicalGroup(2, gap, -1, "gap")
    trace_top_group = gmsh.model.addPhysicalGroup(2, extruded_trace, -1, "trace2")

    trace_points = last.(
        filter(
            x -> x[1] == 0,
            gmsh.model.getBoundary([(2, z) for z in trace], false, true, true)
        )
    )
    trace_curves = last.(
        filter(
            x -> x[1] == 1,
            gmsh.model.getBoundary([(2, z) for z in trace], false, false, false)
        )
    )

    finer_mesh_refinement = false

    if finer_mesh_refinement

        #new attemps to further refine mesh

        #per fare dei test modifico i valori qui
        #l_trace = 2.0
        #l_farfield = 200.0
        #sep_dy = 100.0

        # Generate mesh
        gmsh.option.setNumber("Mesh.MeshSizeMin", l_trace)
        gmsh.option.setNumber("Mesh.MeshSizeMax", l_farfield)
        gmsh.option.setNumber("Mesh.MeshSizeFromPoints", 0)
        gmsh.option.setNumber("Mesh.MeshSizeFromCurvature", 0)
        gmsh.option.setNumber("Mesh.MeshSizeExtendFromBoundary", 0)

        gmsh.model.mesh.field.add("Distance", 1)
        #gmsh.model.mesh.field.setNumbers(1, "PointsList", trace_points)
        #gmsh.model.mesh.field.setNumbers(1, "CurvesList", trace_curves)
        gmsh.model.mesh.field.setNumbers(1, "SurfacesList", trace_total) ##traceAndGap
        #gmsh.model.mesh.field.setNumber(1, "Sampling", ceil(length_μm / l_trace))
        gmsh.model.mesh.field.setNumber(1, "Sampling", 100)

        gmsh.model.mesh.field.add("Constant", 2)
        gmsh.model.mesh.field.setNumber(2, "VIn", 1)
        gmsh.model.mesh.field.setNumber(2, "VOut", 0)
        gmsh.model.mesh.field.setNumbers(2, "VolumesList", metal_volumes)

        gmsh.model.mesh.field.add("MathEval", 3)
        gmsh.model.mesh.field.setString(3, "F",
            # Gmsh MathEval syntax: use F1 (distance) and F2 (volume)
            # F2*(l_trace+x*F1) refers to regions inside the volume,
            # (1-F2)*(l_trace+y*F1) to regions outside the volume
            # Se x = y non c'è differenza nello scaling della mesh tra interno e esterno.
            # Maggiore è il coefficiente moltiplicativo di F1, maggiori saranno le dimensioni delle mesh prodotte (e viceversa)
            string("F2*(", 0, "+1.4*F1) + (1-F2)*(", 0, "+1.6*F1)")
        )

        #Thr interna
        gmsh.model.mesh.field.add("Threshold", 4)
        gmsh.model.mesh.field.setNumber(4, "InField", 3)
        gmsh.model.mesh.field.setNumber(4, "SizeMin", l_trace) ##versione migliore: l_trace
        gmsh.model.mesh.field.setNumber(4, "SizeMax", 1.0 * metal_height_μm) ##versione migliore: metal_height_μm
        gmsh.model.mesh.field.setNumber(4, "DistMin", 1.0 * l_trace) ##versione migliore: 1.0 * l_trace
        gmsh.model.mesh.field.setNumber(4, "DistMax", 2.0 * metal_height_μm) ##versione migliore: 1.0 * metal_height_μm

        #Thr esterna
        gmsh.model.mesh.field.add("Threshold", 5)
        gmsh.model.mesh.field.setNumber(5, "InField", 3)
        gmsh.model.mesh.field.setNumber(5, "SizeMin", l_trace) ##versione migliore: l_trace
        gmsh.model.mesh.field.setNumber(5, "SizeMax", 572) ##versione migliore: l_farfield
        gmsh.model.mesh.field.setNumber(5, "DistMin", 80) ##versione migliore: 1.0 * trace_width_μm
        gmsh.model.mesh.field.setNumber(5, "DistMax", 200) ##versione migliore: 1.0 * sep_dy

        gmsh.model.mesh.field.add("Max", 101) ##versione migliore: Max
        gmsh.model.mesh.field.setNumbers(101, "FieldsList", [4, 5])
        gmsh.model.mesh.field.setAsBackgroundMesh(101)

        gmsh.option.setNumber("Mesh.Algorithm", 6)
        gmsh.option.setNumber("Mesh.Algorithm3D", 10)

    else

        #raw refinement

        # Generate mesh
        gmsh.option.setNumber("Mesh.MeshSizeMin", l_trace)
        gmsh.option.setNumber("Mesh.MeshSizeMax", l_farfield)
        gmsh.option.setNumber("Mesh.MeshSizeFromPoints", 0)
        gmsh.option.setNumber("Mesh.MeshSizeFromCurvature", 0)
        gmsh.option.setNumber("Mesh.MeshSizeExtendFromBoundary", 0)

        gmsh.model.mesh.field.add("Distance", 1)
        #gmsh.model.mesh.field.setNumbers(1, "PointsList", trace_points)
        #gmsh.model.mesh.field.setNumbers(1, "CurvesList", trace_curves)
        gmsh.model.mesh.field.setNumbers(1, "SurfacesList", trace_total)
        gmsh.model.mesh.field.setNumber(1, "Sampling", ceil(length_μm / l_trace))
        #gmsh.model.mesh.field.setNumber(1, "Sampling", 1000)

        #distMin = 1.0 * trace_width_μm
        #distMax = 3.0 * trace_width_μm
        #if distMax <= distMin
        #    distMin = 0.7 * distMin
        #    distMax = 1.3 * distMin
        #end
        #println("\ndistMin: ", distMin, " distMax: ", distMax)
        gmsh.model.mesh.field.add("Threshold", 2)
        gmsh.model.mesh.field.setNumber(2, "InField", 1)
        gmsh.model.mesh.field.setNumber(2, "SizeMin", l_trace) ##prima l_trace
        gmsh.model.mesh.field.setNumber(2, "SizeMax", l_farfield) ##prima l_farfield
        gmsh.model.mesh.field.setNumber(2, "DistMin", 1.0 * trace_width_μm) ##prima 1.0 * trace_width_μm
        gmsh.model.mesh.field.setNumber(2, "DistMax", 1.3 * trace_width_μm) ##prima 0.7 * sep_dy
        
        #=
        ##
        gmsh.model.mesh.field.add("Threshold", 2)
        gmsh.model.mesh.field.setNumber(2, "InField", 1)
        gmsh.model.mesh.field.setNumber(2, "SizeMin", l_trace) #0.5um
        gmsh.model.mesh.field.setNumber(2, "SizeMax", l_farfield) #200um
        gmsh.model.mesh.field.setNumber(2, "DistMin", trace_width_μm / 90) #1um
        gmsh.model.mesh.field.setNumber(2, "DistMax", trace_width_μm / 9) # 10um

        gmsh.model.mesh.field.add("Threshold", 3)
        gmsh.model.mesh.field.setNumber(3, "InField", 1)
        gmsh.model.mesh.field.setNumber(3, "SizeMin", 50 * l_trace) #25um
        gmsh.model.mesh.field.setNumber(3, "SizeMax", l_farfield)
        gmsh.model.mesh.field.setNumber(3, "DistMin", trace_width_μm / 4.5) #20um
        gmsh.model.mesh.field.setNumber(3, "DistMax", trace_width_μm) #90um
        ##
        =#

        gmsh.model.mesh.field.add("MinAniso", 101)
        gmsh.model.mesh.field.setNumbers(101, "FieldsList", [2])
        gmsh.model.mesh.field.setAsBackgroundMesh(101)

        gmsh.option.setNumber("Mesh.Algorithm", 6)
        gmsh.option.setNumber("Mesh.Algorithm3D", 10)
    end


    gmsh.model.mesh.generate(3)
    gmsh.model.mesh.setOrder(order)

    # Save mesh
    gmsh.option.setNumber("Mesh.MshFileVersion", 2.2)
    gmsh.option.setNumber("Mesh.Binary", 0)
    meshFilename = replace(filename, ".msh2" => string("_lumped_mesh", l_trace, "um_port", length_μm*port_factor, "um_length", length_μm,"um.msh2"))
    gmsh.write(joinpath(@__DIR__, meshFilename))

    # Save geo file (geometry script)
    geo_file = replace(meshFilename, ".msh2" => ".geo_unrolled")
    gmsh.write(joinpath(@__DIR__, geo_file))

    # Print some information
    if verbose > 0
        println("\nFinished generating mesh. Physical group tags:")
        println("Air domain: ", air_domain_group)
        println("Si domain: ", si_domain_group)
        if length(metal_domains) > 0
            println("Metal domain: ", metal_domain_group)
        end
        println("Farfield boundaries: ", farfield_group)
        println("Metal boundaries: ", trace_group)
        if length(trace_top_group) > 0
            println("Extruded metal boundaries: ", trace_top_group)
        end
        println("Negative trace boundaries: ", gap_group)

        println("\nMultielement lumped ports:")
        println("Port 1: ", port1a_group, ", ", port1b_group)
        println("Port 2: ", port2a_group, ", ", port2b_group)
        println()
    end

    # Optionally launch GUI
    if gui
        gmsh.fltk.run()
    end

    return gmsh.finalize()
end





"""
Generate a mesh for the coplanar waveguide with lumped ports using Gmsh

# Arguments

  - filename - the filename to use for the generated mesh
  - refinement - measure of how many elements to include, 0 is least
  - order - the polynomial order of the approximation, minimum 1
  - trace_width_μm - width of the coplanar waveguide trace, in μm
  - gap_width_μm - width of the coplanar waveguide gap, in μm
  - separation_width_μm - separation distance between the two waveguides, in μm
  - ground_width_μm - width of the ground plane, in μm
  - substrate_height_μm - height of the substrate, in μm
  - metal_height_μm - metal thickness, in μm
  - remove_metal_vol - for positive metal thickness, whether to remove the metal domain
  - length_μm - length of the waveguides, in μm
  - verbose - flag to dictate the level of print to REPL, passed to Gmsh
  - gui - whether to launch the Gmsh GUI on mesh generation
"""
function generate_trapezoidal_cpw_lumped_mesh(;
    filename::AbstractString,
    refinement::Integer        = 0,
    order::Integer             = 1,
    trace_width_μm::Real       = 90.0,
    gap_width_μm::Real         = 100.0,
    boundary_distance_um::Real = 300.0,
    substrate_height_μm::Real  = 25.0,
    metal_height_μm::Real      = 20.0,
    remove_metal_vol::Bool     = false,
    length_μm::Real            = 8000.0,
    air_distance_um::Real      = 200.0,
    port_factor::Real          = 0.075,
    scale_factor::Real         = 0.75,
    verbose::Integer           = 5,
    gui::Bool                  = false
)
    @assert refinement >= 0
    @assert order > 0
    @assert trace_width_μm > 0
    @assert gap_width_μm > 0
    @assert boundary_distance_um > 0
    @assert substrate_height_μm > 0
    @assert metal_height_μm >= 0
    @assert length_μm > 0
    @assert scale_factor >= 0.0 && scale_factor <= 1.0

    calculate_Z0(trace_width=90e-6, trace_thickness=20e-6, gap_width=100e-6, dielectric_constant=3.3)

    kernel = gmsh.model.occ

    gmsh.initialize()
    gmsh.option.setNumber("General.Verbosity", verbose)

    # Add model
    if "trapezoidal_coplanar_waveguide_driven_lumped" in gmsh.model.list()
        gmsh.model.setCurrent("trapezoidal_coplanar_waveguide_driven_lumped")
        gmsh.model.remove()
    end
    gmsh.model.add("trapezoidal_coplanar_waveguide_driven_lumped")

    sep_dz = air_distance_um
    sep_dy = air_distance_um/3.0
    #@assert sep_dy >= mesh_gap_factor * gap_width_μm string("sep_dy (", sep_dy, ") is lower than ", mesh_gap_factor, " gap_width_μm (", mesh_gap_factor *gap_width_μm, ")! Mesh will not connect the two adjacent conductors!!!")
    #@assert sep_dy > trace_width_μm string("sep_dy (", sep_dy, ") is lower than trace_width_μm (", trace_width_μm, ")! This is not good for mesh construction!!!")

    # Mesh parameters
    l_trace = (1.0/4.0) * trace_width_μm * (2.0^-refinement)
    l_farfield = 1.0 * sep_dz * (2.0^-refinement)

    # Chip pattern
    dy = 0.0
    n1 = kernel.addRectangle(0.0, dy, 0.0, length_μm, boundary_distance_um)
    dy += boundary_distance_um
    gnd1 = kernel.addRectangle(0.0, dy, 0.0, length_μm, trace_width_μm)
    top_gnd1 = kernel.copy([(2, gnd1)])
    kernel.dilate(top_gnd1, 0.0 + length_μm/2, dy + trace_width_μm/2, 0.0, 1.0, scale_factor, 1.0)
    kernel.translate(top_gnd1, 0.0, 0.0, metal_height_μm)
    dy += trace_width_μm
    n2 = kernel.addRectangle(0.0, dy, 0.0, length_μm, gap_width_μm)
    dy += gap_width_μm
    signal = kernel.addRectangle(0.0, dy, 0.0, length_μm, trace_width_μm)
    top_signal = kernel.copy([(2, signal)])
    kernel.dilate(top_signal, 0.0 + length_μm/2, dy + trace_width_μm/2, 0.0, 1.0, scale_factor, 1.0)
    kernel.translate(top_signal, 0.0, 0.0, metal_height_μm)
    dy += trace_width_μm
    n3 = kernel.addRectangle(0.0, dy, 0.0, length_μm, gap_width_μm)
    dy += gap_width_μm
    gnd2 = kernel.addRectangle(0.0, dy, 0.0, length_μm, trace_width_μm)
    top_gnd2 = kernel.copy([(2, gnd2)])
    kernel.dilate(top_gnd2, 0.0 + length_μm/2, dy + trace_width_μm/2, 0.0, 1.0, scale_factor, 1.0)
    kernel.translate(top_gnd2, 0.0, 0.0, metal_height_μm)
    dy += trace_width_μm
    n4 = kernel.addRectangle(0.0, dy, 0.0, length_μm, boundary_distance_um)
    dy += boundary_distance_um

    # Metal thickness
    metal_boundary = [gnd1, signal, gnd2]
    metal = typeof(metal_boundary)(undef, 0)
    metal_dimtags = [] #kernel.extrude([(2, x) for x in metal_boundary], 0.0, 0.0, metal_height_μm)
    gnd1_vol = kernel.addThruSections([gnd1, top_gnd1[1][2]])[1]
    signal_vol = kernel.addThruSections([signal, top_signal[1][2]])[1]
    gnd2_vol = kernel.addThruSections([gnd2, top_gnd2[1][2]])[1]
    push!(metal_dimtags, (3, gnd1_vol[2]))
    push!(metal_dimtags, (3, signal_vol[2]))
    push!(metal_dimtags, (3, gnd2_vol[2]))
    metal = [x[2] for x in filter(x -> x[1] == 3, metal_dimtags)]
    for domain in metal
        _, boundary = kernel.getSurfaceLoops(domain)
        @assert length(boundary) == 1
    end
    
    # Substrate
    substrate = kernel.addBox(0.0, 0.0, -substrate_height_μm, length_μm, dy, substrate_height_μm)
        
    # Exterior box
    domain = 
        kernel.addBox(-0.5 * sep_dy, -0.5 * sep_dy, -sep_dz, length_μm + 1.0 * sep_dy, dy + 1.0 * sep_dy, 2.0 * sep_dz)
    _, domain_boundary = kernel.getSurfaceLoops(domain)
    @assert length(domain_boundary) == 1
    domain_boundary = first(domain_boundary)

    # Ports
    dy = boundary_distance_um + trace_width_μm
    p1a = kernel.addRectangle(0.0, dy, 0.0, length_μm*port_factor, gap_width_μm)
    p2a = kernel.addRectangle(length_μm - length_μm*port_factor, dy, 0.0, length_μm*port_factor, gap_width_μm)
    dy += gap_width_μm + trace_width_μm
    p1b = kernel.addRectangle(0.0, dy, 0.0, length_μm*port_factor, gap_width_μm)
    p2b = kernel.addRectangle(length_μm - length_μm*port_factor, dy, 0.0, length_μm*port_factor, gap_width_μm)

    # Embedding
    geom_dimtags = filter(x -> x[1] == 2 || x[1] == 3, kernel.getEntities())
    _, geom_map = kernel.fragment(geom_dimtags, [])
    kernel.synchronize()

    push!(metal_dimtags, gmsh.model.getBoundary([gnd1_vol, signal_vol, gnd2_vol], true, false, false)...)
    metal = [x[2] for x in filter(x -> x[1] == 3, metal_dimtags)]
    geom_dimtags = filter(x -> x[1] == 2 || x[1] == 3, kernel.getEntities())
    _, geom_map = kernel.fragment(geom_dimtags, [])
    kernel.synchronize()

    extruded_metal_boundary = [x[2] for x in filter(x -> x[1] == 2, metal_dimtags)]
    metal_volumes = [x[2] for x in filter(x -> x[1] == 3, metal_dimtags)]
    filter!(x -> !(x in metal_boundary), extruded_metal_boundary)
    println("\nmetal_boundary:", metal_boundary)
    println("\nmetal_dimtags:", metal_dimtags)
    println("\nextruded_metal_boundary:", extruded_metal_boundary)
    println("\nmetal_volume:", metal_volumes)


    # Add physical groups
    metal_domains = last.(
        collect(
            Iterators.flatten(
                geom_map[findall(x -> x[1] == 3 && x[2] in metal, geom_dimtags)]
            )
        )
    )

    si_domain = last.(geom_map[findfirst(x -> x == (3, substrate), geom_dimtags)])
    @assert length(si_domain) == 1
    si_domain = first(si_domain)

    air_domain = last.(geom_map[findfirst(x -> x == (3, domain), geom_dimtags)])
    filter!(x -> !(x == si_domain || x in metal_domains), air_domain)
    @assert length(air_domain) == 1
    air_domain = first(air_domain)

    if length(metal_domains) > 0 && remove_metal_vol
        remove_dimtags = [(3, x) for x in metal_domains]
        for tag in last.(
            filter(
                x -> x[1] == 2,
                gmsh.model.getBoundary(
                    [(3, z) for z in metal_domains],
                    false,
                    false,
                    false
                )
            )
        )
            normal = gmsh.model.getNormal(tag, [0, 0])
            if abs(normal[1]) == 1.0
                push!(remove_dimtags, (2, tag))
            end
        end
        kernel.remove(remove_dimtags)
        kernel.synchronize()
        filter!.(x -> !(x in remove_dimtags), geom_map)
        empty!(metal_domains)
    end

    air_domain_group = gmsh.model.addPhysicalGroup(3, [air_domain], -1, "air")
    si_domain_group = gmsh.model.addPhysicalGroup(3, [si_domain], -1, "kapton")
    metal_domain_group = gmsh.model.addPhysicalGroup(3, metal_domains, -1, "metal")

    farfield = last.(
        collect(
            Iterators.flatten(
                geom_map[findall(x -> x[1] == 2 && x[2] in domain_boundary, geom_dimtags)]
            )
        )
    )

    farfield_group = gmsh.model.addPhysicalGroup(2, farfield, -1, "farfield")

    port1a = last.(geom_map[findfirst(x -> x == (2, p1a), geom_dimtags)])
    port2a = last.(geom_map[findfirst(x -> x == (2, p2a), geom_dimtags)])
    port1b = last.(geom_map[findfirst(x -> x == (2, p1b), geom_dimtags)])
    port2b = last.(geom_map[findfirst(x -> x == (2, p2b), geom_dimtags)])

    port1a_group = gmsh.model.addPhysicalGroup(2, port1a, -1, "port1a")
    port2a_group = gmsh.model.addPhysicalGroup(2, port2a, -1, "port2a")
    port1b_group = gmsh.model.addPhysicalGroup(2, port1b, -1, "port1b")
    port2b_group = gmsh.model.addPhysicalGroup(2, port2b, -1, "port2b")

    trace = last.(
        collect(
            Iterators.flatten(
                geom_map[findall(x -> x[1] == 2 && x[2] in metal_boundary, geom_dimtags)]
            )
        )
    )
    ##
    extruded_trace = last.(
        collect(
            Iterators.flatten(
                geom_map[findall(x -> x[1] == 2 && x[2] in extruded_metal_boundary, geom_dimtags)]
            )
        )
    )
    trace_total = unique(vcat(trace, extruded_trace))
    println("trace:", trace)
    println("extruded_trace:", extruded_trace)
    println("trace_total:", trace_total)
    ##
    gap = last.(
        collect(
            Iterators.flatten(
                geom_map[findall(x -> x[1] == 2 && x[2] in [n2, n3], geom_dimtags)]
            )
        )
    )
    filter!(
        x -> !(
            x in port1a ||
            x in port2a ||
            x in port1b ||
            x in port2b
        ),
        gap
    )
    traceAndGap = unique(vcat(trace_total, gap))

    trace_group = gmsh.model.addPhysicalGroup(2, trace, -1, "trace")
    gap_group = gmsh.model.addPhysicalGroup(2, gap, -1, "gap")
    trace_top_group = gmsh.model.addPhysicalGroup(2, extruded_trace, -1, "trace2")

    trace_points = last.(
        filter(
            x -> x[1] == 0,
            gmsh.model.getBoundary([(2, z) for z in trace], false, true, true)
        )
    )
    trace_curves = last.(
        filter(
            x -> x[1] == 1,
            gmsh.model.getBoundary([(2, z) for z in trace], false, false, false)
        )
    )

    finer_mesh_refinement = false

    if finer_mesh_refinement

        #new attemps to further refine mesh

        #per fare dei test modifico i valori qui
        #l_trace = 2.0
        #l_farfield = 200.0
        #sep_dy = 100.0

        # Generate mesh
        gmsh.option.setNumber("Mesh.MeshSizeMin", l_trace)
        gmsh.option.setNumber("Mesh.MeshSizeMax", l_farfield)
        gmsh.option.setNumber("Mesh.MeshSizeFromPoints", 0)
        gmsh.option.setNumber("Mesh.MeshSizeFromCurvature", 0)
        gmsh.option.setNumber("Mesh.MeshSizeExtendFromBoundary", 0)

        gmsh.model.mesh.field.add("Distance", 1)
        #gmsh.model.mesh.field.setNumbers(1, "PointsList", trace_points)
        #gmsh.model.mesh.field.setNumbers(1, "CurvesList", trace_curves)
        gmsh.model.mesh.field.setNumbers(1, "SurfacesList", trace_total) ##traceAndGap
        #gmsh.model.mesh.field.setNumber(1, "Sampling", ceil(length_μm / l_trace))
        gmsh.model.mesh.field.setNumber(1, "Sampling", 100)

        gmsh.model.mesh.field.add("Constant", 2)
        gmsh.model.mesh.field.setNumber(2, "VIn", 1)
        gmsh.model.mesh.field.setNumber(2, "VOut", 0)
        gmsh.model.mesh.field.setNumbers(2, "VolumesList", metal_volumes)

        gmsh.model.mesh.field.add("MathEval", 3)
        gmsh.model.mesh.field.setString(3, "F",
            # Gmsh MathEval syntax: use F1 (distance) and F2 (volume)
            # F2*(l_trace+x*F1) refers to regions inside the volume,
            # (1-F2)*(l_trace+y*F1) to regions outside the volume
            # Se x = y non c'è differenza nello scaling della mesh tra interno e esterno.
            # Maggiore è il coefficiente moltiplicativo di F1, maggiori saranno le dimensioni delle mesh prodotte (e viceversa)
            string("F2*(", 0, "+1.4*F1) + (1-F2)*(", 0, "+1.6*F1)")
        )

        #Thr interna
        gmsh.model.mesh.field.add("Threshold", 4)
        gmsh.model.mesh.field.setNumber(4, "InField", 3)
        gmsh.model.mesh.field.setNumber(4, "SizeMin", l_trace) ##versione migliore: l_trace
        gmsh.model.mesh.field.setNumber(4, "SizeMax", 1.0 * metal_height_μm) ##versione migliore: metal_height_μm
        gmsh.model.mesh.field.setNumber(4, "DistMin", 1.0 * l_trace) ##versione migliore: 1.0 * l_trace
        gmsh.model.mesh.field.setNumber(4, "DistMax", 2.0 * metal_height_μm) ##versione migliore: 1.0 * metal_height_μm

        #Thr esterna
        gmsh.model.mesh.field.add("Threshold", 5)
        gmsh.model.mesh.field.setNumber(5, "InField", 3)
        gmsh.model.mesh.field.setNumber(5, "SizeMin", l_trace) ##versione migliore: l_trace
        gmsh.model.mesh.field.setNumber(5, "SizeMax", 572) ##versione migliore: l_farfield
        gmsh.model.mesh.field.setNumber(5, "DistMin", 80) ##versione migliore: 1.0 * trace_width_μm
        gmsh.model.mesh.field.setNumber(5, "DistMax", 200) ##versione migliore: 1.0 * sep_dy

        gmsh.model.mesh.field.add("Max", 101) ##versione migliore: Max
        gmsh.model.mesh.field.setNumbers(101, "FieldsList", [4, 5])
        gmsh.model.mesh.field.setAsBackgroundMesh(101)

        gmsh.option.setNumber("Mesh.Algorithm", 6)
        gmsh.option.setNumber("Mesh.Algorithm3D", 10)

    else

        #raw refinement

        # Generate mesh
        gmsh.option.setNumber("Mesh.MeshSizeMin", l_trace)
        gmsh.option.setNumber("Mesh.MeshSizeMax", l_farfield)
        gmsh.option.setNumber("Mesh.MeshSizeFromPoints", 0)
        gmsh.option.setNumber("Mesh.MeshSizeFromCurvature", 0)
        gmsh.option.setNumber("Mesh.MeshSizeExtendFromBoundary", 0)

        gmsh.model.mesh.field.add("Distance", 1)
        #gmsh.model.mesh.field.setNumbers(1, "PointsList", trace_points)
        #gmsh.model.mesh.field.setNumbers(1, "CurvesList", trace_curves)
        gmsh.model.mesh.field.setNumbers(1, "SurfacesList", trace_total)
        gmsh.model.mesh.field.setNumber(1, "Sampling", ceil(length_μm / l_trace))
        #gmsh.model.mesh.field.setNumber(1, "Sampling", 1000)

        #distMin = 1.0 * trace_width_μm
        #distMax = 3.0 * trace_width_μm
        #if distMax <= distMin
        #    distMin = 0.7 * distMin
        #    distMax = 1.3 * distMin
        #end
        #println("\ndistMin: ", distMin, " distMax: ", distMax)
        gmsh.model.mesh.field.add("Threshold", 2)
        gmsh.model.mesh.field.setNumber(2, "InField", 1)
        gmsh.model.mesh.field.setNumber(2, "SizeMin", l_trace) ##prima l_trace
        gmsh.model.mesh.field.setNumber(2, "SizeMax", l_farfield) ##prima l_farfield
        gmsh.model.mesh.field.setNumber(2, "DistMin", 1.0 * trace_width_μm) ##prima 1.0 * trace_width_μm
        gmsh.model.mesh.field.setNumber(2, "DistMax", 1.3 * trace_width_μm) ##prima 0.7 * sep_dy

        gmsh.model.mesh.field.add("MinAniso", 101)
        gmsh.model.mesh.field.setNumbers(101, "FieldsList", [2])
        gmsh.model.mesh.field.setAsBackgroundMesh(101)

        gmsh.option.setNumber("Mesh.Algorithm", 6)
        gmsh.option.setNumber("Mesh.Algorithm3D", 10)
    end

    gmsh.model.mesh.generate(3)
    gmsh.model.mesh.setOrder(order)

    # Save mesh
    gmsh.option.setNumber("Mesh.MshFileVersion", 2.2)
    gmsh.option.setNumber("Mesh.Binary", 0)
    meshFilename = replace(filename, ".msh2" => string("_lumped_mesh", l_trace, "um_port", length_μm*port_factor, "um_length", length_μm,"um.msh2"))
    gmsh.write(joinpath(@__DIR__, meshFilename))

    # Save geo file (geometry script)
    geo_file = replace(meshFilename, ".msh2" => ".geo_unrolled")
    gmsh.write(joinpath(@__DIR__, geo_file))

    # Print some information
    if verbose > 0
        println("\nFinished generating mesh. Physical group tags:")
        println("Air domain: ", air_domain_group)
        println("Si domain: ", si_domain_group)
        if length(metal_domains) > 0
            println("Metal domain: ", metal_domain_group)
        end
        println("Farfield boundaries: ", farfield_group)
        println("Metal boundaries: ", trace_group)
        if length(trace_top_group) > 0
            println("Extruded metal boundaries: ", trace_top_group)
        end
        println("Negative trace boundaries: ", gap_group)

        println("\nMultielement lumped ports:")
        println("Port 1: ", port1a_group, ", ", port1b_group)
        println("Port 2: ", port2a_group, ", ", port2b_group)
        println()
    end

    # Optionally launch GUI
    if gui
        gmsh.fltk.run()
    end

    return gmsh.finalize()
end


function calculate_Z0(;
    trace_width=90e-6,
    trace_thickness=20e-6,
    gap_width=100e-6,
    substrate_thickness=25e-6,
    dielectric_constant=3.3
)
    """
    Calculate Z0 from the predicted conductor loss using following paper IPC formula and comparing it with the usual conductor loss equation.
    "A new analytical, cad-oriented model for the ohmic and radiation losses of asymmetric coplanar waveguides in
    hybrid and monolithic mic's", G. Ghione, C.U. Naldi
    """

    a = trace_width / 2
    b = a + gap_width
    ks = a / b
    ksprime = sqrt(1 - ks*ks)
    k2 = sinh(pi * a / (2 * substrate_thickness)) / sinh(pi * b / (2 * substrate_thickness))
    k2prime = sqrt(1 - k2*k2)
    e_eff = 1 + ((dielectric_constant - 1) / 2) * (Elliptic.K(ksprime)/Elliptic.K(ks)) * (Elliptic.K(k2)/Elliptic.K(k2prime))
    #print(f'CPW impedance v1: {float((30 * pi) / sqrt(dielectric_constant) * (Elliptic.K(ksprime)/Elliptic.K(ks))):.2f} Ohm')
    Z0 = (1/trace_width) * (240 * pi * Elliptic.K(ks) * Elliptic.K(ksprime) * (1 - ks*ks)) / (sqrt(dielectric_constant) * (1/a * (pi + log(8 * pi * a * (1 - ks) / (trace_thickness * (1 + ks)))) + 1/b * (pi + log(8 * pi * b * (1 - ks) / (trace_thickness * (1 + ks))))))

    println("CPW impedance: $Z0 Ohm")

    return float(Z0)

end








#=
generate_cpw_wave_mesh(
    filename="cpw_port720x450.msh2",
    refinement=0,
    order=1,
    trace_width_μm=90.0,
    gap_width_μm=100.0,
    boundary_distance_um=300.0,
    ground_width_μm=90.0,
    substrate_height_μm=25.0,
    metal_height_μm=20.0,
    remove_metal_vol=false,
    length_μm=1000.0,
    air_distance_um=300.0,
    coax_ports=false,
    verbose=5,
    gui=true
)
=#

#=
generate_cpw_double_traces_lumped_mesh(
    filename="coplanar_waveguide_double_traces_lumped.msh2",
    refinement=1,
    order=1,
    trace_width_μm=90.0,
    gap_width_μm=100.0,
    boundary_distance_um=765.0,
    substrate_height_μm=25.0,
    metal_height_μm=20.0,
    remove_metal_vol=false,
    length_μm = 8000.0,
    verbose=5,
    gui=false
)
=#


generate_cpw_lumped_mesh(
    filename="cpw.msh2",
    refinement=0,
    order=1,
    trace_width_μm=90.0,
    gap_width_μm=100.0,
    boundary_distance_um=300.0,
    substrate_height_μm=25.0,
    metal_height_μm=20.0,
    remove_metal_vol=false,
    length_μm = 1000.0,
    air_distance_um=100.0,
    port_factor=0.075,
    verbose=5,
    gui=true
)


#=
generate_open_short_cpw_lumped_mesh(
    filename="short_cpw.msh2",
    refinement=0,
    order=1,
    trace_width_μm=116.0,
    gap_width_μm=65.0,
    ground_width_μm=300.0,
    boundary_distance_um=300.0,
    substrate_height_μm=525.0,
    metal_height_μm=0.5,
    remove_metal_vol=false,
    length_μm = 1350.0,
    air_distance_um=100.0,
    port_factor=0.075,
    open=false,
    verbose=5,
    gui=true
)
=#

#=
generate_gap_capacitor_cpw_lumped_mesh(
    filename="gap_cap_cpw.msh2",
    refinement=0,
    order=1,
    trace_width_μm=90.0,
    gap_width_μm=100.0,
    boundary_distance_um=300.0,
    substrate_height_μm=25.0,
    metal_height_μm=20.0,
    remove_metal_vol=false,
    length_μm = 500.0,
    air_distance_um=100.0,
    port_factor=0.075,
    gap_capacitor_length_um=20.0,
    verbose=5,
    gui=true
)
=#

#=
generate_diff_cpw_lumped_mesh(
    filename="diff_cpw.msh2",
    refinement=0,
    order=1,
    trace_width_μm=90.0,
    gap_width_μm=100.0,
    boundary_distance_um=300.0,
    substrate_height_μm=50.0,
    metal_height_μm=20.0,
    remove_metal_vol=false,
    length_μm = 1000.0,
    air_distance_um=100.0,
    port_factor=0.075,
    verbose=5,
    gui=true
)
=#

#=
generate_cpw_angle_mesh(
    filename="bent_cpw.msh2",
    refinement=0,
    order=1,
    trace_width_μm=90.0,
    gap_width_μm=100.0,
    boundary_distance_um=300.0,
    substrate_height_μm=25.0,
    metal_height_μm=20.0,
    remove_metal_vol=false,
    length_μm=100.0,
    bend_length_μm=100.0,  # <--- new: second segment length
    air_distance_um=100.0,
    port_factor=0.075,
    verbose=5,
    gui=true
)
=#

#=
generate_stepOut_cpw_mesh(
    filename="stepOut_cpw.msh2",
    refinement=0,
    order=1,
    trace_width_μm=90.0,
    gap_width_μm=100.0,
    boundary_distance_um=300.0,
    substrate_height_μm=25.0,
    metal_height_μm=20.0,
    remove_metal_vol=false,
    length_μm=1000.0,
    air_distance_um=100.0,
    port_factor=0.075,
    gnd_reduction=0.7,
    verbose=5,
    gui=true
)
=#

#=
generate_trapezoidal_cpw_lumped_mesh(
    filename="trapezoidal_cpw.msh2",
    refinement=0,
    order=1,
    trace_width_μm=90.0,
    gap_width_μm=100.0,
    boundary_distance_um=300.0,
    substrate_height_μm=25.0,
    metal_height_μm=20.0,
    remove_metal_vol=false,
    length_μm=1000.0,
    air_distance_um=100.0,
    port_factor=0.075,
    scale_factor=0.75,
    verbose=5,
    gui=true
)
=#