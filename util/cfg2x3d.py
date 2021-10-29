import cfg
import sys
import math

x3d_atom_template = \
r"""<Transform translation='{position[0]:f} {position[1]:f} {position[2]:f}'>
<Transform rotation='0 0 1 {rotation[phi]:f}'> <!--phi-->
<Transform rotation='0 1 0 {rotation[theta]:f}'> <!--theta-->
<Transform scale='{scale:f} {scale:f} {scale:f}'>
<Transform rotation='1 0 0 1.570796325'>
<Transform scale='0.2 0.2 0.2'>
<Shape>
<Appearance>
<Material diffuseColor='0.7 0.7 0.7' ></Material>
</Appearance>
<Sphere></Sphere>
</Shape>
</Transform>
<Transform translation='0 -0.05 0'> <Transform scale='0.075 0.45 0.075'>
<Shape>
<Appearance>
<Material diffuseColor='{color[0]} {color[1]} {color[2]}' ></Material>
</Appearance>
<Cylinder></Cylinder>
</Shape>
</Transform> </Transform>
<Transform translation='0 0.4 0'> <Transform scale='0.15 0.1 0.15'>
<Shape>
<Appearance>
<Material diffuseColor='{color[0]} {color[1]} {color[2]}' ></Material>
</Appearance>
<Cone></Cone>
</Shape>
</Transform> </Transform>
</Transform>
</Transform>
</Transform>
</Transform>
</Transform>
"""

x3d_start_template = \
r"""<X3D version='3.2' profile='Immersive' xmlns:xsd='http://www.w3.org/2001/XMLSchema-instance' xsd:noNamespaceSchemaLocation='http://www.web3d.org/specifications/x3d-3.2.xsd'>
<Scene>
<Viewpoint bind='true' isActive='true' position='{position[0]:f} {position[1]:f} {position[2]:f}' centerOfRotation='{center[0]:f} {center[1]:f} {center[2]:f}'></Viewpoint>
"""

x3d_end_template = \
r"""</Scene>
</X3D>
"""

def norm(v):
    return (sum([x*x for x in v]))**0.5

def diff(arr):
    return [arr[i+1] - arr[i] for i in range(len(arr)-1)]

def main():
    # get parsers
    print("EXPERIMENTAL! USE WITH CAUTION! Email yifangu" + "@andrew.cmu.edu for questions.")
    parsers = {"atoms": cfg.load_atoms, "cells": cfg.load_cells}
    if len(sys.argv) != 4:
        print("Usage {:} [{:}] CELLS_OR_ATOMS_DATA_CFG CELLS_OR_ATOMS_COORDS_CFG".format(sys.argv[0], "|".join(parsers)))
        return -1
    mode = sys.argv[1]
    if mode not in parsers:
        print("{:} is not a valid mode")
        print("valid modes are {:}".format("|".join(parsers)))
        return -1
    filename = sys.argv[2]
    coords_filename = sys.argv[3]
    out_filename = filename + ".x3d"
    print("Processing {:}".format(filename))
    
    # load data
    try:
        data = []
        for d in parsers[mode](filename, coords_filename):
            if sum([abs(x) for x in d]) != 0:
                data.append(d)
    except Exception as e:
        print("Error while parsing file. Check selected mode and source files.")
        return -1

    # preprocess data
    data_t = list(zip(*data))
    xyz_sorted = [sorted(arr) for arr in data_t[:3]]
    mins = [K[0] for K in xyz_sorted]
    maxs = [K[-1] for K in xyz_sorted]
    centers = [(a+b)/2.0 for (a,b) in zip(maxs, mins)]
    dims = [a-b for (a,b) in zip(maxs, mins)]
    grid_sizes = [max(diff(arr)) for arr in xyz_sorted]

    x3d_start_data = {
        "center": centers,
        "position": [centers[0], centers[1], centers[2] + max(dims) * 2.0]
    }

    # output plt file
    with open(out_filename, "w") as f:
        f.write(x3d_start_template.format(**x3d_start_data))
        for d in data:
            (x,y,z,u,v,w) = d
            n = norm(d[3:])
            if n != 0:
                color = [(((k/n+1.0)/2.0)) for k in d[3:]]
                x3d_data = {
                    "position": (x, y, z),
                    "rotation": {
                        "theta":math.atan2((u**2.0+v**2.0)**0.5,w),
                        "phi":math.atan2(v,u)
                    },
                    "scale":0.8*min(grid_sizes),
                    "color": color
                }
                f.write(x3d_atom_template.format(**x3d_data))
            else:
                print("Invalid Magnetization! Skipping")
        f.write(x3d_end_template)
    print("Done, output {:}".format(out_filename))
    return 0

if __name__ == '__main__':
    main()
