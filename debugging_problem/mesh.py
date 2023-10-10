#!/bin/python3

import sys, os
import gmsh
import argparse
import numpy as np

def mesh_reference(refinement_count=0):
    gmsh.model.add("thermal_reference")

    cl = min(cl_left, cl_right, cl_yoke, cl_air)
    cl_ref = cl / np.power(2, refinement_count)

    ### Left coil
    p101 = gmsh.model.geo.addPoint(0, 0, 0, cl_ref)
    p102 = gmsh.model.geo.addPoint(coil_xdim, 0, 0, cl_ref)
    p103 = gmsh.model.geo.addPoint(coil_xdim, coil_ydim, 0, cl_ref)
    p104 = gmsh.model.geo.addPoint(0, coil_ydim, 0, cl_ref)
    l101 = gmsh.model.geo.addLine(p101, p102)
    l102 = gmsh.model.geo.addLine(p102, p103)
    l103 = gmsh.model.geo.addLine(p103, p104)
    l104 = gmsh.model.geo.addLine(p104, p101)
    c101 = gmsh.model.geo.addCurveLoop([l101, l102, l103, l104])
    s101 = gmsh.model.geo.addPlaneSurface([c101])

    ### Right coil
    offset_r = coil_xdim + coil_wIns
    p201 = gmsh.model.geo.addPoint(offset_r, 0, 0, cl_ref)
    p202 = gmsh.model.geo.addPoint(offset_r+coil_xdim, 0, 0, cl_ref)
    p203 = gmsh.model.geo.addPoint(offset_r+coil_xdim, coil_ydim, 0, cl_ref)
    p204 = gmsh.model.geo.addPoint(offset_r, coil_ydim, 0, cl_ref)
    l201 = gmsh.model.geo.addLine(p201, p202)
    l202 = gmsh.model.geo.addLine(p202, p203)
    l203 = gmsh.model.geo.addLine(p203, p204)
    l204 = gmsh.model.geo.addLine(p204, p201)
    c201 = gmsh.model.geo.addCurveLoop([l201, l202, l203, l204])
    s201 = gmsh.model.geo.addPlaneSurface([c201])

    ### Yoke
    yoke_xoffset = coil_xdim + coil_wIns/2 - yoke_xdim/2
    yoke_yoffset = - yoke_wIns - yoke_ydim - yoke_dist
    p301 = gmsh.model.geo.addPoint(yoke_xoffset, yoke_yoffset, 0, cl_ref)
    p302 = gmsh.model.geo.addPoint(yoke_xoffset+yoke_xdim, yoke_yoffset, 0, cl_ref)
    p303 = gmsh.model.geo.addPoint(yoke_xoffset+yoke_xdim, yoke_yoffset+yoke_ydim, 0, cl_ref)
    p304 = gmsh.model.geo.addPoint(yoke_xoffset, yoke_yoffset+yoke_ydim, 0, cl_ref)
    l301 = gmsh.model.geo.addLine(p301, p302)
    l302 = gmsh.model.geo.addLine(p302, p303)
    l303 = gmsh.model.geo.addLine(p303, p304)
    l304 = gmsh.model.geo.addLine(p304, p301)
    c301 = gmsh.model.geo.addCurveLoop([l301, l302, l303, l304])
    s301 = gmsh.model.geo.addPlaneSurface([c301])

    # coil gap
    l501 = gmsh.model.geo.addLine(p102, p201)
    l502 = gmsh.model.geo.addLine(p103, l204)
    c501 = gmsh.model.geo.addCurveLoop([l501, -l204, -l502, -l102])
    s501 = gmsh.model.geo.addPlaneSurface([c501])

    # yoke ins
    p601 = gmsh.model.geo.addPoint(0, -yoke_wIns, 0, cl_ref)
    p602 = gmsh.model.geo.addPoint(2*coil_xdim+coil_wIns, -yoke_wIns, 0, cl_ref)
    l601 = gmsh.model.geo.addLine(p601, p602)
    l602 = gmsh.model.geo.addLine(p602, p202)
    l603 = gmsh.model.geo.addLine(p601, p101)
    c601 = gmsh.model.geo.addCurveLoop([l101, l501, l201, -l602, -l601, l603])
    s601 = gmsh.model.geo.addPlaneSurface([c601])

    # yoke gap
    l701 = gmsh.model.geo.addLine(p304, p601)
    l702 = gmsh.model.geo.addLine(p303, p602)
    c701 = gmsh.model.geo.addCurveLoop([-l303, l702, -l601, -l701])
    s701 = gmsh.model.geo.addPlaneSurface([c701])

    gmsh.model.geo.synchronize()
    gmsh.model.addPhysicalGroup(dim=2, tags=[s501], tag=phys_tag_coil_iso)
    gmsh.model.setPhysicalName(dim=2, tag=phys_tag_coil_iso, name="Insulation between coils")
    gmsh.model.addPhysicalGroup(dim=2, tags=[s601], tag=phys_tag_yoke_iso)
    gmsh.model.setPhysicalName(dim=2, tag=phys_tag_yoke_iso, name="Insulation between coils and yoke gap")
    gmsh.model.addPhysicalGroup(dim=2, tags=[s701], tag=phys_tag_yoke_gap)
    gmsh.model.setPhysicalName(dim=2, tag=phys_tag_yoke_gap, name="Yoke gap")

    ### Tags
    gmsh.model.addPhysicalGroup(dim=2, tags=[s101], tag=phys_tag_coil_left)
    gmsh.model.setPhysicalName(dim=2, tag=phys_tag_coil_left, name="Left coil")
    gmsh.model.addPhysicalGroup(dim=2, tags=[s201], tag=phys_tag_coil_right)
    gmsh.model.setPhysicalName(dim=2, tag=phys_tag_coil_right, name="Right coil")
    gmsh.model.addPhysicalGroup(dim=2, tags=[s301], tag=phys_tag_yoke)
    gmsh.model.setPhysicalName(dim=2, tag=phys_tag_yoke, name="Yoke")

    ### TMP
    gmsh.model.addPhysicalGroup(dim=1, tags=[l104], tag=101)
    gmsh.model.setPhysicalName(dim=1, tag=101, name="Dirichlet left")
    gmsh.model.addPhysicalGroup(dim=1, tags=[l202], tag=102)
    gmsh.model.setPhysicalName(dim=1, tag=102, name="Dirichlet right")

    gmsh.model.geo.synchronize()
    gmsh.model.mesh.generate()
    gmsh.write(problem_dir + "/thermal_reference.msh")
    gmsh.finalize()


def mesh_mortar(refinement_count=0):
    gmsh.model.add("thermal")

    cl_l = cl_left / np.power(2, refinement_count)
    cl_r = cl_right / np.power(2, refinement_count)
    cl_y = cl_yoke / np.power(2, refinement_count)
    cl_a = cl_y

    ### Left coil
    p101 = gmsh.model.geo.addPoint(0, 0, 0, cl_l)
    p102 = gmsh.model.geo.addPoint(coil_xdim, 0, 0, cl_l)
    p103 = gmsh.model.geo.addPoint(coil_xdim, coil_ydim, 0, cl_l)
    p104 = gmsh.model.geo.addPoint(0, coil_ydim, 0, cl_l)
    l101 = gmsh.model.geo.addLine(p101, p102)
    l102 = gmsh.model.geo.addLine(p102, p103)
    l103 = gmsh.model.geo.addLine(p103, p104)
    l104 = gmsh.model.geo.addLine(p104, p101)
    c101 = gmsh.model.geo.addCurveLoop([l101, l102, l103, l104])
    s101 = gmsh.model.geo.addPlaneSurface([c101])

    ### Right coil
    offset_r = coil_xdim + coil_wIns
    p201 = gmsh.model.geo.addPoint(offset_r, 0, 0, cl_r)
    p202 = gmsh.model.geo.addPoint(offset_r+coil_xdim, 0, 0, cl_r)
    p203 = gmsh.model.geo.addPoint(offset_r+coil_xdim, coil_ydim, 0, cl_r)
    p204 = gmsh.model.geo.addPoint(offset_r, coil_ydim, 0, cl_r)
    l201 = gmsh.model.geo.addLine(p201, p202)
    l202 = gmsh.model.geo.addLine(p202, p203)
    l203 = gmsh.model.geo.addLine(p203, p204)
    l204 = gmsh.model.geo.addLine(p204, p201)
    c201 = gmsh.model.geo.addCurveLoop([l201, l202, l203, l204])
    s201 = gmsh.model.geo.addPlaneSurface([c201])

    ### Yoke
    yoke_xoffset = coil_xdim + coil_wIns/2 - yoke_xdim/2
    yoke_yoffset = - yoke_wIns - yoke_ydim - yoke_dist
    p301 = gmsh.model.geo.addPoint(yoke_xoffset, yoke_yoffset, 0, cl_y)
    p302 = gmsh.model.geo.addPoint(yoke_xoffset+yoke_xdim, yoke_yoffset, 0, cl_y)
    p303 = gmsh.model.geo.addPoint(yoke_xoffset+yoke_xdim, yoke_yoffset+yoke_ydim, 0, cl_y)
    p304 = gmsh.model.geo.addPoint(yoke_xoffset, yoke_yoffset+yoke_ydim, 0, cl_y)
    l301 = gmsh.model.geo.addLine(p301, p302)
    l302 = gmsh.model.geo.addLine(p302, p303)
    l303 = gmsh.model.geo.addLine(p303, p304)
    l304 = gmsh.model.geo.addLine(p304, p301)
    c301 = gmsh.model.geo.addCurveLoop([l301, l302, l303, l304])
    s301 = gmsh.model.geo.addPlaneSurface([c301])

    # mesh air between yoke and yoke-coil TSA
    p401 = gmsh.model.geo.addPoint(0, -yoke_wIns, 0, cl_a)
    p402 = gmsh.model.geo.addPoint(coil_xdim, -yoke_wIns, 0, cl_a)
    p403 = gmsh.model.geo.addPoint(coil_xdim + coil_wIns, -yoke_wIns, 0, cl_a)
    p404 = gmsh.model.geo.addPoint(coil_xdim*2 + coil_wIns, -yoke_wIns, 0, cl_a)
    l401 = gmsh.model.geo.addLine(p303, p404)
    l402 = gmsh.model.geo.addLine(p404, p403)
    l403 = gmsh.model.geo.addLine(p403, p402)
    l404 = gmsh.model.geo.addLine(p402, p401)
    l405 = gmsh.model.geo.addLine(p401, p304)
    c401 = gmsh.model.geo.addCurveLoop([-l303, l401, l402, l403, l404, l405])
    s401 = gmsh.model.geo.addPlaneSurface([c401])

    ### Coil TSA
    coil_gap_x = coil_xdim+coil_wIns/2
    # aux. left
    p501 = gmsh.model.geo.addPoint(coil_gap_x, 0, 0, cl_l)
    p502 = gmsh.model.geo.addPoint(coil_gap_x, coil_ydim, 0, cl_l)
    l501 = gmsh.model.geo.addLine(p501, p502)

    # aux. right
    p511 = gmsh.model.geo.addPoint(coil_gap_x, 0, 0, cl_r)
    p512 = gmsh.model.geo.addPoint(coil_gap_x, coil_ydim, 0, cl_r)
    l511 = gmsh.model.geo.addLine(p511, p512)

    # aux. center
    p521 = gmsh.model.geo.addPoint(coil_gap_x, 0, 0, min(cl_l, cl_r))
    p522 = gmsh.model.geo.addPoint(coil_gap_x, coil_ydim, 0, min(cl_l, cl_r))
    l521 = gmsh.model.geo.addLine(p511, p512)

    ### Yoke TSA, left
    yoke_gap_y = - yoke_wIns / 2
    # aux. top
    p601 = gmsh.model.geo.addPoint(0, yoke_gap_y, 0, cl_l)
    p602 = gmsh.model.geo.addPoint(coil_xdim, yoke_gap_y, 0, cl_l)
    l601 = gmsh.model.geo.addLine(p601, p602)

    # aux. bottom
    p611 = gmsh.model.geo.addPoint(0, yoke_gap_y, 0, cl_a)
    p612 = gmsh.model.geo.addPoint(coil_xdim, yoke_gap_y, 0, cl_a)
    l611 = gmsh.model.geo.addLine(p611, p612)

    # aux. center
    p621 = gmsh.model.geo.addPoint(0, yoke_gap_y, 0, min(cl_l, cl_a))
    p622 = gmsh.model.geo.addPoint(coil_xdim, yoke_gap_y, 0, min(cl_l, cl_a))
    l621 = gmsh.model.geo.addLine(p621, p622)

    ### Yoke TSA, right
    # aux. top
    p701 = gmsh.model.geo.addPoint(coil_xdim+coil_wIns, yoke_gap_y, 0, cl_r)
    p702 = gmsh.model.geo.addPoint(coil_xdim*2+coil_wIns, yoke_gap_y, 0, cl_r)
    l701 = gmsh.model.geo.addLine(p701, p702)

    # aux. bottom
    p711 = gmsh.model.geo.addPoint(coil_xdim+coil_wIns, yoke_gap_y, 0, cl_a)
    p712 = gmsh.model.geo.addPoint(coil_xdim*2+coil_wIns, yoke_gap_y, 0, cl_a)
    l711 = gmsh.model.geo.addLine(p711, p712)

    # aux. center
    p721 = gmsh.model.geo.addPoint(coil_xdim+coil_wIns, yoke_gap_y, 0, min(cl_r, cl_a))
    p722 = gmsh.model.geo.addPoint(coil_xdim*2+coil_wIns, yoke_gap_y, 0, min(cl_r, cl_a))
    l721 = gmsh.model.geo.addLine(p721, p722)

    gmsh.model.geo.synchronize()
    gmsh.model.addPhysicalGroup(dim=2, tags=[s401], tag=phys_tag_yoke_gap)
    gmsh.model.setPhysicalName(dim=2, tag=phys_tag_yoke_gap, name="Yoke gap")

    gmsh.model.addPhysicalGroup(dim=1, tags=[l501], tag=phys_tag_coil_aux_left)
    gmsh.model.setPhysicalName(dim=1, tag=phys_tag_coil_aux_left, name="Aux. surface used for TSA between coils - left side")
    gmsh.model.addPhysicalGroup(dim=1, tags=[l511], tag=phys_tag_coil_aux_right)
    gmsh.model.setPhysicalName(dim=1, tag=phys_tag_coil_aux_right, name="Aux. surface used for TSA between coils - right side")
    gmsh.model.addPhysicalGroup(dim=1, tags=[l521], tag=phys_tag_coil_aux_center)
    gmsh.model.setPhysicalName(dim=1, tag=phys_tag_coil_aux_center, name="Aux. surface used for TSA between coils - center")
    gmsh.model.addPhysicalGroup(dim=1, tags=[l102], tag=phys_tag_coil_left_interface_right)
    gmsh.model.setPhysicalName(dim=1, tag=phys_tag_coil_left_interface_right, name="Left coil's interface with coil_aux_left")
    gmsh.model.addPhysicalGroup(dim=1, tags=[l204], tag=phys_tag_coil_right_interface_left)
    gmsh.model.setPhysicalName(dim=1, tag=phys_tag_coil_right_interface_left, name="Right coil's interface with coil_aux_right")

    gmsh.model.addPhysicalGroup(dim=1, tags=[l601], tag=phys_tag_yoke_left_aux_top)
    gmsh.model.setPhysicalName(dim=1, tag=phys_tag_yoke_left_aux_top, name="Aux. surface used for TSA between yoke and left coil - top side")
    gmsh.model.addPhysicalGroup(dim=1, tags=[l611], tag=phys_tag_yoke_left_aux_bottom)
    gmsh.model.setPhysicalName(dim=1, tag=phys_tag_yoke_left_aux_bottom, name="Aux. surface used for TSA between yoke and left coil - bottom side")
    gmsh.model.addPhysicalGroup(dim=1, tags=[l621], tag=phys_tag_yoke_left_aux_center)
    gmsh.model.setPhysicalName(dim=1, tag=phys_tag_yoke_left_aux_center, name="Aux. surface used for TSA between yoke and left coil - center")
    gmsh.model.addPhysicalGroup(dim=1, tags=[l101], tag=phys_tag_coil_left_interface_bottom)
    gmsh.model.setPhysicalName(dim=1, tag=phys_tag_coil_left_interface_bottom, name="Left coil's interface with yoke_left_aux_top")
    gmsh.model.addPhysicalGroup(dim=1, tags=[l404], tag=phys_tag_yoke_left_interface_top)
    gmsh.model.setPhysicalName(dim=1, tag=phys_tag_yoke_left_interface_top, name="Yoke air's interface with yoke_left_aux_bottom")

    gmsh.model.addPhysicalGroup(dim=1, tags=[l701], tag=phys_tag_yoke_right_aux_top)
    gmsh.model.setPhysicalName(dim=1, tag=phys_tag_yoke_right_aux_top, name="Aux. surface used for TSA between yoke and right coil - top side")
    gmsh.model.addPhysicalGroup(dim=1, tags=[l711], tag=phys_tag_yoke_right_aux_bottom)
    gmsh.model.setPhysicalName(dim=1, tag=phys_tag_yoke_right_aux_bottom, name="Aux. surface used for TSA between yoke and right coil - bottom side")
    gmsh.model.addPhysicalGroup(dim=1, tags=[l721], tag=phys_tag_yoke_right_aux_center)
    gmsh.model.setPhysicalName(dim=1, tag=phys_tag_yoke_right_aux_center, name="Aux. surface used for TSA between yoke and right coil - center")
    gmsh.model.addPhysicalGroup(dim=1, tags=[l201], tag=phys_tag_coil_right_interface_bottom)
    gmsh.model.setPhysicalName(dim=1, tag=phys_tag_coil_right_interface_bottom, name="Right coil's interface with yoke_right_aux_top")
    gmsh.model.addPhysicalGroup(dim=1, tags=[l402], tag=phys_tag_yoke_right_interface_top)
    gmsh.model.setPhysicalName(dim=1, tag=phys_tag_yoke_right_interface_top, name="Yoke air's interface with yoke_right_aux_bottom")

    ### Tags
    gmsh.model.addPhysicalGroup(dim=2, tags=[s101], tag=phys_tag_coil_left)
    gmsh.model.setPhysicalName(dim=2, tag=phys_tag_coil_left, name="Left coil")
    gmsh.model.addPhysicalGroup(dim=2, tags=[s201], tag=phys_tag_coil_right)
    gmsh.model.setPhysicalName(dim=2, tag=phys_tag_coil_right, name="Right coil")
    gmsh.model.addPhysicalGroup(dim=2, tags=[s301], tag=phys_tag_yoke)
    gmsh.model.setPhysicalName(dim=2, tag=phys_tag_yoke, name="Yoke")

    ### TMP
    gmsh.model.addPhysicalGroup(dim=1, tags=[l104], tag=101)
    gmsh.model.setPhysicalName(dim=1, tag=101, name="Dirichlet left")
    gmsh.model.addPhysicalGroup(dim=1, tags=[l202], tag=102)
    gmsh.model.setPhysicalName(dim=1, tag=102, name="Dirichlet right")

    gmsh.model.geo.synchronize()
    gmsh.model.mesh.generate()
    gmsh.write(problem_dir + "/thermal.msh")
    gmsh.finalize()


################################################################################
##################################### MAIN #####################################
################################################################################

problem_dir = os.path.dirname(__file__)

gmsh.initialize()

### Load parameters
gmsh.parser.parse(problem_dir + "/geometry_parameters.pro")

coil_xdim = gmsh.parser.getNumber("coil_xdim")[0]
coil_ydim = gmsh.parser.getNumber("coil_ydim")[0]
coil_wIns = gmsh.parser.getNumber("coil_wIns")[0]
cl_coil = gmsh.parser.getNumber("cl_coil")[0]
yoke_xdim = gmsh.parser.getNumber("yoke_xdim")[0]
yoke_ydim = gmsh.parser.getNumber("yoke_ydim")[0]
yoke_wIns = gmsh.parser.getNumber("yoke_wIns")[0]
yoke_dist = gmsh.parser.getNumber("yoke_dist")[0]
cl_yoke = gmsh.parser.getNumber("cl_yoke")[0]
cl_air = gmsh.parser.getNumber("cl_air")[0]


# Load tags
phys_tag_coil_left = int(gmsh.parser.getNumber("phys_tag_coil_left")[0])
phys_tag_coil_right = int(gmsh.parser.getNumber("phys_tag_coil_right")[0])
phys_tag_yoke = int(gmsh.parser.getNumber("phys_tag_yoke")[0])
phys_tag_yoke_gap = int(gmsh.parser.getNumber("phys_tag_yoke_gap")[0])
phys_tag_yoke_iso = int(gmsh.parser.getNumber("phys_tag_yoke_iso")[0])
phys_tag_coil_iso = int(gmsh.parser.getNumber("phys_tag_coil_iso")[0])

# mortar + TSA
phys_tag_coil_aux_left = int(gmsh.parser.getNumber("phys_tag_coil_aux_left")[0])
phys_tag_coil_aux_right = int(gmsh.parser.getNumber("phys_tag_coil_aux_right")[0])
phys_tag_coil_aux_center = int(gmsh.parser.getNumber("phys_tag_coil_aux_center")[0])
phys_tag_coil_left_interface_right = int(gmsh.parser.getNumber("phys_tag_coil_left_interface_right")[0])
phys_tag_coil_right_interface_left = int(gmsh.parser.getNumber("phys_tag_coil_right_interface_left")[0])

phys_tag_yoke_left_aux_top = int(gmsh.parser.getNumber("phys_tag_yoke_left_aux_top")[0])
phys_tag_coil_left_interface_bottom = int(gmsh.parser.getNumber("phys_tag_coil_left_interface_bottom")[0])
phys_tag_yoke_left_aux_bottom = int(gmsh.parser.getNumber("phys_tag_yoke_left_aux_bottom")[0])
phys_tag_yoke_left_interface_top = int(gmsh.parser.getNumber("phys_tag_yoke_left_interface_top")[0])
phys_tag_yoke_left_aux_center = int(gmsh.parser.getNumber("phys_tag_yoke_left_aux_center")[0])

phys_tag_yoke_right_aux_top = int(gmsh.parser.getNumber("phys_tag_yoke_right_aux_top")[0])
phys_tag_coil_right_interface_bottom = int(gmsh.parser.getNumber("phys_tag_coil_right_interface_bottom")[0])
phys_tag_yoke_right_aux_bottom = int(gmsh.parser.getNumber("phys_tag_yoke_right_aux_bottom")[0])
phys_tag_yoke_right_interface_top = int(gmsh.parser.getNumber("phys_tag_yoke_right_interface_top")[0])
phys_tag_yoke_right_aux_center = int(gmsh.parser.getNumber("phys_tag_yoke_right_aux_center")[0])

cl_left = cl_coil
cl_right = cl_coil/2.5

### Determine run conditions

USE_REFERENCE = False
N_REFINE = 0
INITIAL_MAX_MESH_SIZE = 0.1

parser = argparse.ArgumentParser(
    prog="Meshing script for Thermal Debugging Problem",
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
