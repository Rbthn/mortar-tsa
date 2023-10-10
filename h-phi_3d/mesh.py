#!/bin/python3
import gmsh
import sys, os
import argparse
import numpy as np
from math import sqrt

def get_tags_from_dimtags(dimtags):
    return [dimtag[1] for dimtag in dimtags]


################################################################################
#################################### MORTAR ####################################
################################################################################

def mesh_mortar(refinement_count=0):
    gmsh.model.add("mag")

    ####
    #### TOP
    ####

    coil_top = gmsh.model.occ.addBox(size_x/2-w_x/2, size_y/2, 0, w_x, w_y/2, w_z)
    air_top = gmsh.model.occ.addBox(0, size_y/2, 0, size_x, size_y/2, size_z)

    ####
    #### BOTTOM
    ####

    coil_bottom_1 = gmsh.model.occ.addBox(size_x/2-w_x/2, size_y/2-w_y/2, 0, w_x, w_y/4, w_z)
    coil_bottom_2 = gmsh.model.occ.addBox(size_x/2-w_x/2, size_y/2-w_y/4, 0, w_x, w_y/4, w_z)
    air_bottom = gmsh.model.occ.addBox(0, 0, 0, size_x, size_y/2, size_z)
    split_tool_p1 = gmsh.model.occ.addPoint(size_x/2, 0, 0)
    split_tool_p2 = gmsh.model.occ.addPoint(size_x/2, size_y/2-w_y/2, 0)
    split_tool_l1 = gmsh.model.occ.addLine(split_tool_p1, split_tool_p2)
    split_tool = gmsh.model.occ.extrude([(1, split_tool_l1)], dx=0, dy=0, dz=size_z)
    split_tool = [dimtag for dimtag in split_tool if dimtag[0] == 2] # only dim==2

    ####
    #### Fragment, Translate
    ####

    [frag_vol, frag_map] = gmsh.model.occ.fragment([(3, air_top), (3, air_bottom)], [(3, coil_top), (3, coil_bottom_1), (3, coil_bottom_2)])
    air_top = 5
    air_bottom = 6
    gmsh.model.occ.synchronize()

    [frag_vol, frag_map] = gmsh.model.occ.fragment([(3, air_bottom), (3, coil_bottom_1)], split_tool)
    air_bottom = [6, 7]

    gmsh.model.occ.synchronize()
    top_tags = [(3, air_top), (3, coil_top)]
    gmsh.model.occ.translate(top_tags, dx=0, dy=1, dz=0)
    gmsh.model.occ.translate(top_tags, dx=0, dy=-1, dz=0)
    gmsh.model.occ.synchronize()

    air_bnd_top = gmsh.model.getBoundary([(3, air_top)], oriented=False)
    air_bnd_bottom = gmsh.model.getBoundary([(3, tag) for tag in air_bottom], oriented=False)
    coil_bnd_top = gmsh.model.getBoundary([(3, coil_top)], oriented=False)
    coil_bnd_bottom = gmsh.model.getBoundary([(3, coil_bottom_2)], oriented=False)
    gmsh.model.occ.synchronize()

    coil_interface_top = [(2,95)] #gmsh.model.occ.intersect(coil_bnd_top, [(3, coil_bottom)], removeObject=True, removeTool=False)[0]
    coil_interface_bottom = [(2, 52)] #gmsh.model.occ.intersect(coil_bnd_bottom, [(3, coil_top)], removeObject=True, removeTool=False)[0]
    air_interface_top = [(2, 86), (2, 92)] #gmsh.model.occ.intersect(air_bnd_top, [(3, air_bottom)], removeObject=True, removeTool=False)[0]
    air_interface_bottom = [(2, 33), (2, 39)] #gmsh.model.occ.intersect(air_bnd_bottom, [(3, air_top)], removeObject=True, removeTool=False)[0]

    ####
    #### Tags
    #####

    gmsh.model.addPhysicalGroup(dim=3, tags=[air_top], tag=phys_tag_air_top, name="Air top")
    gmsh.model.addPhysicalGroup(dim=3, tags=air_bottom, tag=phys_tag_air_bottom, name="Air bottom")
    gmsh.model.addPhysicalGroup(dim=3, tags=[coil_top], tag=phys_tag_coil_top, name="Coil top")
    gmsh.model.addPhysicalGroup(dim=3, tags=[coil_bottom_1, coil_bottom_2], tag=phys_tag_coil_bottom, name="Coil bottom")
    gmsh.model.addPhysicalGroup(dim=2, tags=get_tags_from_dimtags(coil_interface_top), tag=phys_tag_coil_interface_top, name="Coil interface - top")
    gmsh.model.addPhysicalGroup(dim=2, tags=get_tags_from_dimtags(coil_interface_bottom), tag=phys_tag_coil_interface_bottom, name="Coil interface - bottom")
    gmsh.model.addPhysicalGroup(dim=2, tags=get_tags_from_dimtags(air_interface_top), tag=phys_tag_air_interface_top, name="Air interface - top")
    gmsh.model.addPhysicalGroup(dim=2, tags=get_tags_from_dimtags(air_interface_bottom), tag=phys_tag_air_interface_bottom, name="Air interface - bottom")

    ####
    #### Mesh
    ####

    # setting different node counts depending on orientation of edge
    node_density = np.power(2, refinement_count) / INITIAL_MAX_MESH_SIZE

    transfinite_curves_top = get_tags_from_dimtags(gmsh.model.getBoundary(coil_interface_top, oriented=False))
    transfinite_curves_bottom = get_tags_from_dimtags(gmsh.model.getBoundary(coil_interface_bottom, oriented=False))

    for curve_list, factor in zip([transfinite_curves_bottom, transfinite_curves_top], [1., 1.2]):
        for c in curve_list:
            endpoints = get_tags_from_dimtags(gmsh.model.getBoundary([(1, c)]))
            coord_0 = gmsh.model.getValue(dim=0, tag=endpoints[0], parametricCoord=[])
            coord_1 = gmsh.model.getValue(dim=0, tag=endpoints[1], parametricCoord=[])
            diff = coord_0 - coord_1
            edge_length = sqrt(diff[0]**2+diff[1]**2+diff[2]**2)
            num_nodes = node_density * edge_length
            gmsh.model.mesh.setTransfiniteCurve(c, round(factor*num_nodes)+2)
    gmsh.model.mesh.setTransfiniteSurface(get_tags_from_dimtags(coil_interface_top)[0])
    gmsh.model.mesh.setTransfiniteSurface(get_tags_from_dimtags(coil_interface_bottom)[0])

    translation = [1, 0, 0, 0, 0, 1, 0, 0, 0, 0, 1, 0, 0, 0, 0, 1]
    gmsh.model.mesh.setPeriodic(2, get_tags_from_dimtags(air_interface_bottom), get_tags_from_dimtags(air_interface_top), translation)

    # recombine interface into quads
    for dimtag in coil_interface_top:
        gmsh.model.mesh.setRecombine(dimtag[0], dimtag[1])
    for dimtag in coil_interface_bottom:
        gmsh.model.mesh.setRecombine(dimtag[0], dimtag[1])

    gmsh.option.setNumber('Mesh.MeshSizeMax', INITIAL_MAX_MESH_SIZE)
    gmsh.model.mesh.generate()
    for i in range(refinement_count):
        gmsh.model.mesh.refine()

    ####
    #### Cohomology
    ####

    cuts_domain_helper = 50
    gmsh.model.addPhysicalGroup(dim=2, tags=[split_tool[0][1]], tag=phys_tag_cut_base, name="Location of cut")
    gmsh.model.addPhysicalGroup(dim=3, tags=[coil_bottom_2], tag=cuts_domain_helper, name="Helper for Cut")

    gmsh.model.mesh.addHomologyRequest(type="Cohomology", domainTags=[phys_tag_air_bottom, phys_tag_air_top, cuts_domain_helper], subdomainTags=[phys_tag_cut_base], dims=[1])
    cuts = gmsh.model.mesh.computeHomology()

    # rename cuts
    for i in range(len(cuts)):
        cut_dimtag = cuts[i]
        dim = cut_dimtag[0]
        tag = cut_dimtag[1]
        old_name = gmsh.model.getPhysicalName(dim, tag)
        gmsh.model.removePhysicalName(old_name)
        gmsh.model.setPhysicalName(dim=dim, tag=phys_tag_cut_base+i+1, name=("Cut_"+str(i+1)))

    ####
    #### Finalize
    ####

    gmsh.write(problem_dir + "/mag.msh")
    gmsh.finalize()



################################################################################
################################### REFERENCE ##################################
################################################################################

def mesh_reference(refinement_count=0):
    gmsh.model.add("mag_reference")

    ####
    #### Geometry
    ####

    coil = gmsh.model.occ.addBox(size_x/2-w_x/2, size_y/2-w_y/2, 0, w_x, w_y, w_z)
    air = gmsh.model.occ.addBox(0, 0, 0, size_x, size_y, size_z)
    split_tool_p1 = gmsh.model.occ.addPoint(size_x/2, 0, 0)
    split_tool_p2 = gmsh.model.occ.addPoint(size_x/2, size_y/2-w_y/2, 0)
    split_tool_l1 = gmsh.model.occ.addLine(split_tool_p1, split_tool_p2)
    split_tool = gmsh.model.occ.extrude([(1, split_tool_l1)], dx=0, dy=0, dz=size_z)
    split_tool = [dimtag for dimtag in split_tool if dimtag[0] == 2] # only dim==2

    [frag_vol, frag_map] = gmsh.model.occ.fragment([(3, air)], [(3, coil)])
    gmsh.model.occ.synchronize()

    ####
    #### Tags
    ####

    gmsh.model.addPhysicalGroup(dim=3, tags=[coil], tag=phys_tag_coil, name="Coil")
    gmsh.model.addPhysicalGroup(dim=3, tags=[air], tag=phys_tag_air, name="Air")

    ####
    #### Mesh
    ####

    gmsh.model.occ.synchronize()
    gmsh.option.setNumber('Mesh.MeshSizeMax', INITIAL_MAX_MESH_SIZE)
    gmsh.model.mesh.generate()
    for i in range(refinement_count):
        gmsh.model.mesh.refine()

    ####
    #### Cohomology
    ####

    gmsh.model.addPhysicalGroup(dim=2, tags=[split_tool[0][1]], tag=phys_tag_cut_base, name="Location of cut")

    gmsh.model.mesh.addHomologyRequest(type="Cohomology", domainTags=[phys_tag_air], subdomainTags=[phys_tag_cut_base], dims=[1])
    cuts = gmsh.model.mesh.computeHomology()

    # rename cuts
    for i in range(len(cuts)):
        cut_dimtag = cuts[i]
        dim = cut_dimtag[0]
        tag = cut_dimtag[1]
        old_name = gmsh.model.getPhysicalName(dim, tag)
        gmsh.model.removePhysicalName(old_name)
        gmsh.model.setPhysicalName(dim=dim, tag=phys_tag_cut_base+i+1, name=("Cut_"+str(i+1)))

    ####
    #### Finalize
    ####

    gmsh.write(problem_dir + "/mag_reference.msh")
    gmsh.finalize()


################################################################################
##################################### MAIN #####################################
################################################################################

### Get Parameters

problem_dir = os.path.dirname(__file__)

gmsh.initialize()
gmsh.parser.parse(problem_dir + "/geometry_parameters.pro")

size_x = gmsh.parser.getNumber("size_x")[0]
size_y = gmsh.parser.getNumber("size_y")[0]
size_z = gmsh.parser.getNumber("size_z")[0]
w_x = gmsh.parser.getNumber("w_x")[0]
w_y = gmsh.parser.getNumber("w_y")[0]
w_z = gmsh.parser.getNumber("w_z")[0]
#
phys_tag_air_top = int(gmsh.parser.getNumber("phys_tag_air_top")[0])
phys_tag_air_bottom = int(gmsh.parser.getNumber("phys_tag_air_bottom")[0])
phys_tag_coil_top = int(gmsh.parser.getNumber("phys_tag_coil_top")[0])
phys_tag_coil_bottom = int(gmsh.parser.getNumber("phys_tag_coil_bottom")[0])
phys_tag_coil_interface_top = int(gmsh.parser.getNumber("phys_tag_coil_interface_top")[0])
phys_tag_coil_interface_bottom = int(gmsh.parser.getNumber("phys_tag_coil_interface_bottom")[0])
phys_tag_air_interface_top = int(gmsh.parser.getNumber("phys_tag_air_interface_top")[0])
phys_tag_air_interface_bottom = int(gmsh.parser.getNumber("phys_tag_air_interface_bottom")[0])
phys_tag_cut_base = int(gmsh.parser.getNumber("phys_tag_cut_base")[0])
# reference
phys_tag_coil = int(gmsh.parser.getNumber("phys_tag_coil")[0])
phys_tag_air = int(gmsh.parser.getNumber("phys_tag_air")[0])

### Determine run conditions

USE_REFERENCE = False
N_REFINE = 0
INITIAL_MAX_MESH_SIZE = 0.1

parser = argparse.ArgumentParser(
    prog="Meshing script for H-Phi 3D Problem",
    description="",
)

parser.add_argument('--reference', action='store_true')
parser.add_argument('--refine', type=int)
parser.add_argument('--initsize', type=float)

args = parser.parse_args()
if args.reference:
    USE_REFERENCE = args.reference
if args.refine:
    N_REFINE = args.refine
if args.initsize:
    INITIAL_MESH_SIZE = args.initsize


### Run

if USE_REFERENCE:
    mesh_reference(N_REFINE)
else:
    mesh_mortar(N_REFINE)
