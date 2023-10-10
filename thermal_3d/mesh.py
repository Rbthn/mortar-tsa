#!/bin/python3
import gmsh
import sys, os
import argparse

def mesh_mortar(refinement_count=0):
    gmsh.model.add("thermal")

    ####
    #### LEFT SUBDOMAIN
    ####

    # left port
    p101 = gmsh.model.geo.addPoint(0, port_yoffset, port_zoffset, cl1)
    p102 = gmsh.model.geo.addPoint(0, port_yoffset, port_zoffset+w_z, cl1)
    p103 = gmsh.model.geo.addPoint(0, port_yoffset+w_y, port_zoffset+w_z, cl1)
    p104 = gmsh.model.geo.addPoint(0, port_yoffset+w_y, port_zoffset, cl1)
    l101 = gmsh.model.geo.addLine(p101, p102)
    l102 = gmsh.model.geo.addLine(p102, p103)
    l103 = gmsh.model.geo.addLine(p103, p104)
    l104 = gmsh.model.geo.addLine(p104, p101)
    c101 = gmsh.model.geo.addCurveLoop([l101, l102, l103, l104])
    s101 = gmsh.model.geo.addPlaneSurface([c101])
    gmsh.model.addPhysicalGroup(dim=2, tags=[s101], tag=phys_tag_port_left)
    gmsh.model.setPhysicalName(dim=2, tag=phys_tag_port_left, name="Left port")

    # diagnoally cut face
    p105 = gmsh.model.geo.addPoint(L1-w_Ins/2, port_yoffset, port_zoffset, cl1)
    p106 = gmsh.model.geo.addPoint(L1-w_Ins/2, port_yoffset, port_zoffset+w_z, cl1)
    p107 = gmsh.model.geo.addPoint(L2-w_Ins/2, port_yoffset+w_y, port_zoffset+w_z, cl1)
    p108 = gmsh.model.geo.addPoint(L2-w_Ins/2, port_yoffset+w_y, port_zoffset, cl1)
    l105 = gmsh.model.geo.addLine(p105, p106)
    l106 = gmsh.model.geo.addLine(p106, p107)
    l107 = gmsh.model.geo.addLine(p107, p108)
    l108 = gmsh.model.geo.addLine(p108, p105)
    c102 = gmsh.model.geo.addCurveLoop([l105, l106, l107, l108])
    s102 = gmsh.model.geo.addPlaneSurface([c102])
    gmsh.model.addPhysicalGroup(dim=2, tags=[s102], tag=phys_tag_interface_left)
    gmsh.model.setPhysicalName(dim=2, tag=phys_tag_interface_left, name="Interface left")

    # connecting face, bottom
    l110 = gmsh.model.geo.addLine(p101, p105)
    l111 = gmsh.model.geo.addLine(p102, p106)
    c110 = gmsh.model.geo.addCurveLoop([l101, l111, -l105, -l110])
    s110 = gmsh.model.geo.addPlaneSurface([c110])
    # connecting face, top
    l112 = gmsh.model.geo.addLine(p104, p108)
    l113 = gmsh.model.geo.addLine(p103, p107)
    c111 = gmsh.model.geo.addCurveLoop([-l103, l113, l107, -l112])
    s111 = gmsh.model.geo.addPlaneSurface([c111])
    # connecting face, back
    c112 = gmsh.model.geo.addCurveLoop([l104, l110, -l108, -l112])
    s112 = gmsh.model.geo.addPlaneSurface([c112])
    # connecting face, front
    c113 = gmsh.model.geo.addCurveLoop([-l102, l111, l106, -l113])
    s113 = gmsh.model.geo.addPlaneSurface([c113])

    sl101 = gmsh.model.geo.addSurfaceLoop([s101, s102, s110, s111, s112, s113])
    v101 = gmsh.model.geo.addVolume([sl101])
    gmsh.model.addPhysicalGroup(dim=3, tags=[v101], tag=phys_tag_left)
    gmsh.model.setPhysicalName(dim=3, tag=phys_tag_left, name="Left")

    ####
    #### RIGHT SUBDOMAIN
    ####

    # right port
    p201 = gmsh.model.geo.addPoint(size_x, port_yoffset, port_zoffset, cl2)
    p202 = gmsh.model.geo.addPoint(size_x, port_yoffset, port_zoffset+w_z, cl2)
    p203 = gmsh.model.geo.addPoint(size_x, port_yoffset+w_y, port_zoffset+w_z, cl2)
    p204 = gmsh.model.geo.addPoint(size_x, port_yoffset+w_y, port_zoffset, cl2)
    l201 = gmsh.model.geo.addLine(p201, p202)
    l202 = gmsh.model.geo.addLine(p202, p203)
    l203 = gmsh.model.geo.addLine(p203, p204)
    l204 = gmsh.model.geo.addLine(p204, p201)
    c201 = gmsh.model.geo.addCurveLoop([l201, l202, l203, l204])
    s201 = gmsh.model.geo.addPlaneSurface([c201])
    gmsh.model.addPhysicalGroup(dim=2, tags=[s201], tag=phys_tag_port_right)
    gmsh.model.setPhysicalName(dim=2, tag=phys_tag_port_right, name="Right port")

    # diagnoally cut face
    p205 = gmsh.model.geo.addPoint(L1+w_Ins/2, port_yoffset, port_zoffset, cl2)
    p206 = gmsh.model.geo.addPoint(L1+w_Ins/2, port_yoffset, port_zoffset+w_z, cl2)
    p207 = gmsh.model.geo.addPoint(L2+w_Ins/2, port_yoffset+w_y, port_zoffset+w_z, cl2)
    p208 = gmsh.model.geo.addPoint(L2+w_Ins/2, port_yoffset+w_y, port_zoffset, cl2)
    l205 = gmsh.model.geo.addLine(p205, p206)
    l206 = gmsh.model.geo.addLine(p206, p207)
    l207 = gmsh.model.geo.addLine(p207, p208)
    l208 = gmsh.model.geo.addLine(p208, p205)
    c202 = gmsh.model.geo.addCurveLoop([l205, l206, l207, l208])
    s202 = gmsh.model.geo.addPlaneSurface([c202])
    gmsh.model.addPhysicalGroup(dim=2, tags=[s202], tag=phys_tag_interface_right)
    gmsh.model.setPhysicalName(dim=2, tag=phys_tag_interface_right, name="Interface right")

    # connecting face, bottom
    l210 = gmsh.model.geo.addLine(p201, p205)
    l211 = gmsh.model.geo.addLine(p202, p206)
    c210 = gmsh.model.geo.addCurveLoop([l201, l211, -l205, -l210])
    s210 = gmsh.model.geo.addPlaneSurface([c210])
    # connecting face, top
    l212 = gmsh.model.geo.addLine(p204, p208)
    l213 = gmsh.model.geo.addLine(p203, p207)
    c211 = gmsh.model.geo.addCurveLoop([-l203, l213, l207, -l212])
    s211 = gmsh.model.geo.addPlaneSurface([c211])
    # connecting face, back
    c212 = gmsh.model.geo.addCurveLoop([l204, l210, -l208, -l212])
    s212 = gmsh.model.geo.addPlaneSurface([c212])
    # connecting face, front
    c213 = gmsh.model.geo.addCurveLoop([-l202, l211, l206, -l213])
    s213 = gmsh.model.geo.addPlaneSurface([c213])

    sl201 = gmsh.model.geo.addSurfaceLoop([s201, s202, s210, s211, s212, s213])
    v201 = gmsh.model.geo.addVolume([sl201])
    gmsh.model.addPhysicalGroup(dim=3, tags=[v201], tag=phys_tag_right)
    gmsh.model.setPhysicalName(dim=3, tag=phys_tag_right, name="Right")


    ####
    #### Aux. TSA interfaces
    ####
    cl_tsa = min(cl1, cl2)

    # Shell_aux
    p1001 = gmsh.model.geo.addPoint(L1, port_yoffset, port_zoffset, cl_tsa)
    p1002 = gmsh.model.geo.addPoint(L1, port_yoffset, port_zoffset+w_z, cl_tsa)
    p1003 = gmsh.model.geo.addPoint(L2, port_yoffset+w_y, port_zoffset+w_z, cl_tsa)
    p1004 = gmsh.model.geo.addPoint(L2, port_yoffset+w_y, port_zoffset, cl_tsa)
    l1001 = gmsh.model.geo.addLine(p1001, p1002)
    l1002 = gmsh.model.geo.addLine(p1002, p1003)
    l1003 = gmsh.model.geo.addLine(p1003, p1004)
    l1004 = gmsh.model.geo.addLine(p1004, p1001)
    c1001 = gmsh.model.geo.addCurveLoop([l1001, l1002, l1003, l1004])
    s1001 = gmsh.model.geo.addPlaneSurface([c1001])
    gmsh.model.addPhysicalGroup(dim=2, tags=[s1001], tag=phys_tag_shell_aux)
    gmsh.model.setPhysicalName(dim=2, tag=phys_tag_shell_aux, name="Aux. shell for TSA")

    # left_aux
    p1011 = gmsh.model.geo.addPoint(L1, port_yoffset, port_zoffset, cl1)
    p1012 = gmsh.model.geo.addPoint(L1, port_yoffset, port_zoffset+w_z, cl1)
    p1013 = gmsh.model.geo.addPoint(L2, port_yoffset+w_y, port_zoffset+w_z, cl1)
    p1014 = gmsh.model.geo.addPoint(L2, port_yoffset+w_y, port_zoffset, cl1)
    l1011 = gmsh.model.geo.addLine(p1011, p1012)
    l1012 = gmsh.model.geo.addLine(p1012, p1013)
    l1013 = gmsh.model.geo.addLine(p1013, p1014)
    l1014 = gmsh.model.geo.addLine(p1014, p1011)
    c1011 = gmsh.model.geo.addCurveLoop([l1011, l1012, l1013, l1014])
    s1011 = gmsh.model.geo.addPlaneSurface([c1011])
    gmsh.model.addPhysicalGroup(dim=2, tags=[s1011], tag=phys_tag_left_aux)
    gmsh.model.setPhysicalName(dim=2, tag=phys_tag_left_aux, name="Left aux. shell for mortar")

    # right_aux
    p1021 = gmsh.model.geo.addPoint(L1, port_yoffset, port_zoffset, cl2)
    p1022 = gmsh.model.geo.addPoint(L1, port_yoffset, port_zoffset+w_z, cl2)
    p1023 = gmsh.model.geo.addPoint(L2, port_yoffset+w_y, port_zoffset+w_z, cl2)
    p1024 = gmsh.model.geo.addPoint(L2, port_yoffset+w_y, port_zoffset, cl2)
    l1021 = gmsh.model.geo.addLine(p1021, p1022)
    l1022 = gmsh.model.geo.addLine(p1022, p1023)
    l1023 = gmsh.model.geo.addLine(p1023, p1024)
    l1024 = gmsh.model.geo.addLine(p1024, p1021)
    c1021 = gmsh.model.geo.addCurveLoop([l1021, l1022, l1023, l1024])
    s1021 = gmsh.model.geo.addPlaneSurface([c1021])
    gmsh.model.addPhysicalGroup(dim=2, tags=[s1021], tag=phys_tag_right_aux)
    gmsh.model.setPhysicalName(dim=2, tag=phys_tag_right_aux, name="Right aux. shell for mortar")

    gmsh.model.geo.synchronize()
    gmsh.model.mesh.generate()
    for i in range(refinement_count):
        gmsh.model.mesh.refine()
    gmsh.write(problem_dir + "/thermal.msh")
    gmsh.finalize()

def mesh_reference(refinement_count=0):
    gmsh.model.add("thermal_reference")

    ####
    #### LEFT SUBDOMAIN
    ####

    # left port
    p101 = gmsh.model.geo.addPoint(0, port_yoffset, port_zoffset, cl1)
    p102 = gmsh.model.geo.addPoint(0, port_yoffset, port_zoffset+w_z, cl1)
    p103 = gmsh.model.geo.addPoint(0, port_yoffset+w_y, port_zoffset+w_z, cl1)
    p104 = gmsh.model.geo.addPoint(0, port_yoffset+w_y, port_zoffset, cl1)
    l101 = gmsh.model.geo.addLine(p101, p102)
    l102 = gmsh.model.geo.addLine(p102, p103)
    l103 = gmsh.model.geo.addLine(p103, p104)
    l104 = gmsh.model.geo.addLine(p104, p101)
    c101 = gmsh.model.geo.addCurveLoop([l101, l102, l103, l104])
    s101 = gmsh.model.geo.addPlaneSurface([c101])
    gmsh.model.addPhysicalGroup(dim=2, tags=[s101], tag=phys_tag_port_left)
    gmsh.model.setPhysicalName(dim=2, tag=phys_tag_port_left, name="Left port")

    # diagnoally cut face
    p105 = gmsh.model.geo.addPoint(L1-w_Ins/2, port_yoffset, port_zoffset, cl1)
    p106 = gmsh.model.geo.addPoint(L1-w_Ins/2, port_yoffset, port_zoffset+w_z, cl1)
    p107 = gmsh.model.geo.addPoint(L2-w_Ins/2, port_yoffset+w_y, port_zoffset+w_z, cl1)
    p108 = gmsh.model.geo.addPoint(L2-w_Ins/2, port_yoffset+w_y, port_zoffset, cl1)
    l105 = gmsh.model.geo.addLine(p105, p106)
    l106 = gmsh.model.geo.addLine(p106, p107)
    l107 = gmsh.model.geo.addLine(p107, p108)
    l108 = gmsh.model.geo.addLine(p108, p105)
    c102 = gmsh.model.geo.addCurveLoop([l105, l106, l107, l108])
    s102 = gmsh.model.geo.addPlaneSurface([c102])
    gmsh.model.addPhysicalGroup(dim=2, tags=[s102], tag=phys_tag_interface_left)
    gmsh.model.setPhysicalName(dim=2, tag=phys_tag_interface_left, name="Interface left")

    # connecting face, bottom
    l110 = gmsh.model.geo.addLine(p101, p105)
    l111 = gmsh.model.geo.addLine(p102, p106)
    c110 = gmsh.model.geo.addCurveLoop([l101, l111, -l105, -l110])
    s110 = gmsh.model.geo.addPlaneSurface([c110])
    # connecting face, top
    l112 = gmsh.model.geo.addLine(p104, p108)
    l113 = gmsh.model.geo.addLine(p103, p107)
    c111 = gmsh.model.geo.addCurveLoop([-l103, l113, l107, -l112])
    s111 = gmsh.model.geo.addPlaneSurface([c111])
    # connecting face, back
    c112 = gmsh.model.geo.addCurveLoop([l104, l110, -l108, -l112])
    s112 = gmsh.model.geo.addPlaneSurface([c112])
    # connecting face, front
    c113 = gmsh.model.geo.addCurveLoop([-l102, l111, l106, -l113])
    s113 = gmsh.model.geo.addPlaneSurface([c113])

    sl101 = gmsh.model.geo.addSurfaceLoop([s101, s102, s110, s111, s112, s113])
    v101 = gmsh.model.geo.addVolume([sl101])
    gmsh.model.addPhysicalGroup(dim=3, tags=[v101], tag=phys_tag_left)
    gmsh.model.setPhysicalName(dim=3, tag=phys_tag_left, name="Left")

    ####
    #### RIGHT SUBDOMAIN
    ####

    # right port
    p201 = gmsh.model.geo.addPoint(size_x, port_yoffset, port_zoffset, cl2)
    p202 = gmsh.model.geo.addPoint(size_x, port_yoffset, port_zoffset+w_z, cl2)
    p203 = gmsh.model.geo.addPoint(size_x, port_yoffset+w_y, port_zoffset+w_z, cl2)
    p204 = gmsh.model.geo.addPoint(size_x, port_yoffset+w_y, port_zoffset, cl2)
    l201 = gmsh.model.geo.addLine(p201, p202)
    l202 = gmsh.model.geo.addLine(p202, p203)
    l203 = gmsh.model.geo.addLine(p203, p204)
    l204 = gmsh.model.geo.addLine(p204, p201)
    c201 = gmsh.model.geo.addCurveLoop([l201, l202, l203, l204])
    s201 = gmsh.model.geo.addPlaneSurface([c201])
    gmsh.model.addPhysicalGroup(dim=2, tags=[s201], tag=phys_tag_port_right)
    gmsh.model.setPhysicalName(dim=2, tag=phys_tag_port_right, name="Right port")

    # diagnoally cut face
    p205 = gmsh.model.geo.addPoint(L1+w_Ins/2, port_yoffset, port_zoffset, cl2)
    p206 = gmsh.model.geo.addPoint(L1+w_Ins/2, port_yoffset, port_zoffset+w_z, cl2)
    p207 = gmsh.model.geo.addPoint(L2+w_Ins/2, port_yoffset+w_y, port_zoffset+w_z, cl2)
    p208 = gmsh.model.geo.addPoint(L2+w_Ins/2, port_yoffset+w_y, port_zoffset, cl2)
    l205 = gmsh.model.geo.addLine(p205, p206)
    l206 = gmsh.model.geo.addLine(p206, p207)
    l207 = gmsh.model.geo.addLine(p207, p208)
    l208 = gmsh.model.geo.addLine(p208, p205)
    c202 = gmsh.model.geo.addCurveLoop([l205, l206, l207, l208])
    s202 = gmsh.model.geo.addPlaneSurface([c202])
    gmsh.model.addPhysicalGroup(dim=2, tags=[s202], tag=phys_tag_interface_right)
    gmsh.model.setPhysicalName(dim=2, tag=phys_tag_interface_right, name="Interface right")

    # connecting face, bottom
    l210 = gmsh.model.geo.addLine(p201, p205)
    l211 = gmsh.model.geo.addLine(p202, p206)
    c210 = gmsh.model.geo.addCurveLoop([l201, l211, -l205, -l210])
    s210 = gmsh.model.geo.addPlaneSurface([c210])
    # connecting face, top
    l212 = gmsh.model.geo.addLine(p204, p208)
    l213 = gmsh.model.geo.addLine(p203, p207)
    c211 = gmsh.model.geo.addCurveLoop([-l203, l213, l207, -l212])
    s211 = gmsh.model.geo.addPlaneSurface([c211])
    # connecting face, back
    c212 = gmsh.model.geo.addCurveLoop([l204, l210, -l208, -l212])
    s212 = gmsh.model.geo.addPlaneSurface([c212])
    # connecting face, front
    c213 = gmsh.model.geo.addCurveLoop([-l202, l211, l206, -l213])
    s213 = gmsh.model.geo.addPlaneSurface([c213])

    sl201 = gmsh.model.geo.addSurfaceLoop([s201, s202, s210, s211, s212, s213])
    v201 = gmsh.model.geo.addVolume([sl201])
    gmsh.model.addPhysicalGroup(dim=3, tags=[v201], tag=phys_tag_right)
    gmsh.model.setPhysicalName(dim=3, tag=phys_tag_right, name="Right")

    # gap
    l301 = gmsh.model.geo.addLine(p106, p206)
    l302 = gmsh.model.geo.addLine(p105, p205)
    l303 = gmsh.model.geo.addLine(p107, p207)
    l304 = gmsh.model.geo.addLine(p108, p208)
    c301 = gmsh.model.geo.addCurveLoop([l301, l206, -l303, -l106])
    s301 = gmsh.model.geo.addPlaneSurface([c301])
    c302 = gmsh.model.geo.addCurveLoop([l302, -l208, -l304, l108])
    s302 = gmsh.model.geo.addPlaneSurface([c302])
    c303 = gmsh.model.geo.addCurveLoop([l301, -l205, -l302, l105])
    s303 = gmsh.model.geo.addPlaneSurface([c303])
    c304 = gmsh.model.geo.addCurveLoop([l303, l207, -l304, -l107])
    s304 = gmsh.model.geo.addPlaneSurface([c304])
    sl301 = gmsh.model.geo.addSurfaceLoop([s301, s302, s303, s304, s102, s202])
    v301 = gmsh.model.geo.addVolume([sl301])
    gmsh.model.addPhysicalGroup(dim=3, tags=[v301], tag=phys_tag_air, name="Air gap")


    gmsh.model.geo.synchronize()
    gmsh.model.mesh.generate()
    for i in range(refinement_count):
        gmsh.model.mesh.refine()
    gmsh.write(problem_dir + "/thermal_reference.msh")
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
L1 = gmsh.parser.getNumber("L1")[0]
L2 = gmsh.parser.getNumber("L2")[0]
w_y = gmsh.parser.getNumber("w_y")[0]
w_z = gmsh.parser.getNumber("w_z")[0]
w_Ins = gmsh.parser.getNumber("w_Ins")[0]
cl1 = gmsh.parser.getNumber("cl1")[0]
cl2 = gmsh.parser.getNumber("cl2")[0]
#
phys_tag_left = int(gmsh.parser.getNumber("phys_tag_left")[0])
phys_tag_right = int(gmsh.parser.getNumber("phys_tag_right")[0])
phys_tag_air = int(gmsh.parser.getNumber("phys_tag_air")[0])
phys_tag_port_left = int(gmsh.parser.getNumber("phys_tag_port_left")[0])
phys_tag_port_right = int(gmsh.parser.getNumber("phys_tag_port_right")[0])
phys_tag_interface_right = int(gmsh.parser.getNumber("phys_tag_interface_right")[0])
phys_tag_interface_left = int(gmsh.parser.getNumber("phys_tag_interface_left")[0])
phys_tag_shell_aux = int(gmsh.parser.getNumber("phys_tag_shell_aux")[0])
phys_tag_left_aux = int(gmsh.parser.getNumber("phys_tag_left_aux")[0])
phys_tag_right_aux = int(gmsh.parser.getNumber("phys_tag_right_aux")[0])

port_yoffset = (size_y - w_y)/2
port_zoffset = (size_z - w_z)/2

### Determine run conditions

USE_REFERENCE = False
N_REFINE = 0
INITIAL_MAX_MESH_SIZE = 0.1

parser = argparse.ArgumentParser(
    prog="Meshing script for Thermal 3D Problem",
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