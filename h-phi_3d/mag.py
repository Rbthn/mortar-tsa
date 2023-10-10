#!/bin/python3
import sys
import gmsh
import subprocess

if "-nomesh" not in sys.argv:
    subprocess.run(["./mesh.py", "-nopopup"])

if "-nosolve" not in sys.argv:
    subprocess.run(["./solve.py", ""])

if "-nopopup" not in sys.argv:
    gmsh.initialize()
    gmsh.open("mag.brep")
    gmsh.merge("h.pos")
    gmsh.fltk.run()
    gmsh.finalize()