#!/usr/bin/python3
'''
setup_simulation.py

Reads setup.cfg file from current directory, creates necessary directories
and creates the solvers configuration files
'''
import os
import sys
from reader import read_input

if not os.path.exists("./setup.cfg"):
    print("Could not find 'setup.cfg' in current directory. Aborting")
    sys.exit(1)

for dir_name in ["./input","./output"]:
    if not os.path.exists(dir_name):
        os.mkdir(dir_name)
    if not os.path.isdir(dir_name):
        raise NotADirectoryError(dir_name)

input_file = open("./setup.cfg")
opts, geometry = read_input(input_file)

with open("./solver.cfg", "w") as out:
    out.write("# File generated automatically - do not change it\n"
              + "# To rebuild the configuration use the setup.cfg file\n"
              + "# and rerun the setup_simulation.py distributed\n"
              + "# with the package\n\n\n"
              )
    for key, value in sorted(opts.items()):
        out.write(key+" = "+value+"\n")

