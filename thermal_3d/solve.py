#!/bin/python3
import gmsh
import sys
import subprocess
import argparse

problem_name = "thermal.pro"
resolution_name = "heat_v"
postop_name = "Map"

### Determine run conditions
USE_REFERENCE = False

parser = argparse.ArgumentParser(
    prog="Meshing script for H-Phi 3D Problem",
    description="",
)

parser.add_argument('--reference', action='store_true')
parser.add_argument('--steady_state', action='store_true')

args = parser.parse_args()
USE_REFERENCE = args.reference
STEADY_STATE = args.steady_state

getdp_args = ["getdp", problem_name]
if USE_REFERENCE:
    getdp_args.append('-setnumber')
    getdp_args.append('USE_REFERENCE')
    getdp_args.append('1')
if STEADY_STATE:
    getdp_args.append('-setnumber')
    getdp_args.append('STEADY_STATE')
    getdp_args.append('1')
getdp_args.append("-solve")
getdp_args.append(resolution_name)
getdp_args.append("-pos")
getdp_args.append(postop_name)


process = subprocess.run(getdp_args)
exit(process.returncode)