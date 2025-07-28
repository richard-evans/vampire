import re

def load_array(lines, tuple_length=3, skip=0):
    regex = r"^([\-eE0-9\.]+)"
    for i in range(tuple_length+skip-1):
        regex += r"[ \t]+([\-eE0-9\.]+)"
    regex += r".*"
    compiled_regex = re.compile(regex)
    result = []
    for line in lines:
        m = compiled_regex.match(line)
        if m is not None:
            r = [float(m.group(i+skip+1)) for i in range(tuple_length)]
            result.append(tuple(r))
    return result

def load_cells(cells_filename, coords_filename = None):
    if coords_filename is None:
        coords_filename = 'cells-coords.cfg'
    cells_lines = []
    coords_lines = []
    with open(cells_filename) as f:
        cells_lines += f.readlines()
    with open(coords_filename) as f:
        coords_lines += f.readlines()
    cells_mags = load_array(cells_lines)
    cells_coords = load_array(coords_lines)
    assert(len(cells_mags) == len(cells_coords))
    result = []
    for i in range(len(cells_mags)):
        (x, y, z) = cells_coords[i]
        (u, v, w) = cells_mags[i]
        result.append((x,y,z,u,v,w))
    return result

def load_atoms(atoms_filename, coords_filename = None):
    if coords_filename is None:
        coords_filename = 'atoms-coords.cfg'
    atoms_lines = []
    coords_lines = []
    with open(atoms_filename) as f:
        atoms_lines += f.readlines()
    with open(coords_filename) as f:
        coords_lines += f.readlines()
    atoms_mags = load_array(atoms_lines)
    atoms_coords = load_array(coords_lines, skip=2)
    assert(len(atoms_mags) == len(atoms_coords))
    result = []
    for i in range(len(atoms_mags)):
        (x, y, z) = atoms_coords[i]
        (u, v, w) = atoms_mags[i]
        result.append((x,y,z,u,v,w))
    return result





