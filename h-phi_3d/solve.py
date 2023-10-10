#!/bin/python3
import gmsh
import sys
import subprocess
import argparse

problem_name = "mag.pro"
resolution_name = "res"
postop_name = "map"

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
# set time params
getdp_args.append('-setnumber')
getdp_args.append('t_end')
getdp_args.append('7')
getdp_args.append('-setnumber')
getdp_args.append('t_step')
getdp_args.append('0.1')
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
