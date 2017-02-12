# -*- coding: utf-8 -*-
import argparse, sys, re

def read_xyz(fp):
    regexp_float = r'[+-]?(\d+\.)?\d+([deE][+-]?\d+)?'
    xyz = []
    is_first_line = True
    atom_number = 0
    for line in fp:
        if is_first_line:
            # Number of atoms only.
            m1 = re.search(r'(\d+)', line)
            if m1:
                num_atoms = int(m1.group(1))
                is_first_line = False
            else:
                print 'Header of xyz file not found'
                sys.exit(1)
            # Number of atoms and box size.
            m2 = re.search(r'(\d+)\s+(?P<box_x>%s)\s+(?P<box_y>%s)\s+(?P<box_z>%s)' %
                           (regexp_float, regexp_float, regexp_float), line)
            if m2:
                box = [float(m2.group('box_x')), float(m2.group('box_y')), float(m2.group('box_z'))]
            else:
                box = None
        else:
            m = re.search(r'(?P<element>\w+)\s+(?P<x>%s)\s+(?P<y>%s)\s+(?P<z>%s)' %
                         (regexp_float, regexp_float, regexp_float), line)
            if m:
                atom = (m.group('element'), float(m.group('x')), float(m.group('y')), float(m.group('z')))
                xyz.append(atom)
                atom_number += 1
        if len(xyz) >= num_atoms:
            break
    return xyz, box

def get_grid_division_num(xmin, xmax, max_bond_length):
    max_division = 10000
    div = int((xmax - xmin) / max_bond_length)
    return max(1, min(max_division, div))

def init_grid(xyz, box, max_bond_length):
    if box is None:  # Non-periodic mode.
        xs = map(lambda atom: atom[1], xyz)
        ys = map(lambda atom: atom[2], xyz)
        zs = map(lambda atom: atom[3], xyz)
        # eps is needed because get_grid_index_3arity(min(xs), max(xs), xdiv, max(xs)) == xdiv,
        # which means array index out of range.
        eps = 1.0e-8
        grid = {'xmin': min(xs), 'xmax': max(xs) + eps,
                'ymin': min(ys), 'ymax': max(ys) + eps,
                'zmin': min(zs), 'zmax': max(zs) + eps}
    else:  # Periodic mode.
        grid = {'xmin': 0., 'xmax': box[0],
                'ymin': 0., 'ymax': box[1],
                'zmin': 0., 'zmax': box[2]}
    grid['xdiv'] = get_grid_division_num(grid['xmin'], grid['xmax'], max_bond_length)
    grid['ydiv'] = get_grid_division_num(grid['ymin'], grid['ymax'], max_bond_length)
    grid['zdiv'] = get_grid_division_num(grid['zmin'], grid['zmax'], max_bond_length)
    # Initialize each division.
    grid['atoms'] = []
    for i in range(grid['xdiv'] * grid['ydiv'] * grid['zdiv']):
        grid['atoms'].append([])
    return grid

# Find grid index i which satisfies xmin + d * i <= x < xmin + d * (i + 1),
# where d is the grid width, (xmax - xmin) / xdiv.
# Note that the range is inclusve for the left and exclusive for the right.
def get_grid_index_3arity(xmin, xmax, xdiv, x):
    # Adjustment to the range xmin <= x < xmax. Needed only in the periodic mode.
    while x < xmin:
        x += xmax - xmin
    while xmax <= x:
        x -= xmax - xmin
    return int(xdiv * (x - xmin) / (xmax - xmin))

def grid_index_3to1(grid, ix, iy, iz):
    return ix * grid['ydiv'] * grid['zdiv'] + iy * grid['zdiv'] + iz

def grid_index_1to3(grid, i):
    ix = i / (grid['ydiv'] * grid['zdiv'])
    i -= ix * grid['ydiv'] * grid['zdiv']
    iy = i / grid['zdiv']
    iz = i - iy * grid['zdiv']
    return ix, iy, iz

def get_grid_index_1arity(grid, atom):
    ix = get_grid_index_3arity(grid['xmin'], grid['xmax'], grid['xdiv'], atom[1])
    iy = get_grid_index_3arity(grid['ymin'], grid['ymax'], grid['ydiv'], atom[2])
    iz = get_grid_index_3arity(grid['zmin'], grid['zmax'], grid['zdiv'], atom[3])
    return grid_index_3to1(grid, ix, iy, iz)

def allocate_atoms_to_grid(xyz, grid):
    for j, atom in enumerate(xyz):
        grid_index = get_grid_index_1arity(grid, atom)
        grid['atoms'][grid_index].append(j)

def get_neighbor_grid_indices_coord(cx, xdiv, is_non_periodic):
    if is_non_periodic:  # In the non-periodic mode, out-of-range grid indices are ignored.
        return range(max(0, cx - 1), min(xdiv, cx + 2))
    else:  # In the periodic mode, out-of-range grid indices should be wrapped.
        neighbors = []
        for i in range(cx - 1, cx + 2):
            if i == -1:
                neighbors.append(xdiv - 1)
            elif 0 <= i and i < xdiv:
                neighbors.append(i)
            elif i == xdiv:
                neighbors.append(0)
            else:
                assert(False)
        return neighbors

def get_neighbor_grid_indices(grid, index, is_non_periodic):
    cx, cy, cz = grid_index_1to3(grid, index)  # c means center.
    neighbors = []
    for ix in get_neighbor_grid_indices_coord(cx, grid['xdiv'], is_non_periodic):
        for iy in get_neighbor_grid_indices_coord(cy, grid['ydiv'], is_non_periodic):
            for iz in get_neighbor_grid_indices_coord(cz, grid['zdiv'], is_non_periodic):
                i = grid_index_3to1(grid, ix, iy, iz)
                neighbors.append(i)
    return neighbors

def read_bond_def(filename):
    rex_float = r'(\d+\.)?\d+([eE][+-]?\d+)?'
    rex_data = r'(?P<elem1>\w+)\s+(?P<elem2>\w+)\s+(?P<length>%s)' % rex_float
    bond_def = []
    with open(filename, 'r') as fp:
        for line in iter(fp):
            line = line.replace('d', 'e')
            m = re.match(rex_data, line)
            if m:
                data = {'elem1': m.group('elem1'),
                        'elem2': m.group('elem2'),
                        'length': float(m.group('length'))}
                bond_def.append(data)
    return bond_def

# Distance is treated in squared value to alleviate square root calculation,
def get_bond_length(bond_def, elem1, elem2):
    for bond in bond_def:
        if ((bond['elem1'] == elem1 and bond['elem2'] == elem2) or
            (bond['elem1'] == elem2 and bond['elem2'] == elem1)):
            return bond['length'] ** 2.0
    return 0.0

def get_distance(atom1, atom2):
    dist = (atom1[1] - atom2[1]) ** 2.0
    dist += (atom1[2] - atom2[2]) ** 2.0
    dist += (atom1[3] - atom2[3]) ** 2.0
    return dist

def is_close(bond_def, atom1, atom2):
    return get_distance(atom1, atom2) < get_bond_length(bond_def, atom1[0], atom2[0])

class UnionFind:
    def __init__(self, size):
        self.data = [-1] * size
    def union(self, x, y):
        xr = self.root(x)
        yr = self.root(y)
        if xr != yr:
            if self.data[xr] < self.data[yr]:
                xr, yr = yr, xr
            self.data[xr] = self.data[xr] + self.data[yr]
            self.data[yr] = xr
        return xr != yr
    def find(self, x, y):
        return self.root(x) == self.root(y)
    def root(self, x):
        visited = []
        px = self.data[x]
        while px >= 0:
            visited.append(x)
            x = px
            px = self.data[px]
        for v in visited:
            self.data[v] = x
        return x
    def size(self, x):
        return -self.data[self.root(x)]
    def print_groups(self, fp):
        num_atoms = len(self.data)
        roots = []
        for i in range(num_atoms):
            if self.data[i] < 0:
                roots.append(i)
        num_groups = len(roots)
        nums_group_atoms = filter(lambda n: n < 0, self.data)
        max_num_group_atoms = -min(nums_group_atoms)
        min_num_group_atoms = -max(nums_group_atoms)
        print '--- result ---'
        print 'number of groups:', num_groups
        print 'max number of atoms in group:', max_num_group_atoms
        print 'min number of atoms in group:', min_num_group_atoms
        fp.write('#\n')
        fp.write('# number of atoms, number of groups, maximum number of group atoms\n')
        fp.write('%d %d %d\n' % (num_atoms, num_groups, max_num_group_atoms))
        fp.write('#\n')
        for i in range(num_atoms):
            rx = roots.index(self.root(i))
            fp.write('%d %d\n' % (i + 1, rx + 1))
    def __str__(self):
        return str(self.data)

def get_groups(xyz, bond_def, box):
    size = len(xyz)
    groups = UnionFind(size)
    max_bond_length = max(map(lambda bond:bond['length'], bond_def))
    grid = init_grid(xyz, box, max_bond_length)
    print 'number of atoms: ', size
    print 'grid division for x y z total:', grid['xdiv'], grid['ydiv'], grid['zdiv'], \
        grid['xdiv'] * grid['ydiv'] * grid['zdiv']
    allocate_atoms_to_grid(xyz, grid)
    grid_sizes = map(len, grid['atoms'])
    print 'max number of atoms for grid divisions:', max(grid_sizes)
    assert(size == sum(grid_sizes))  # Check whether all atoms are allocated to one grid box.
    for atom1_i in range(size):
        atom1 = xyz[atom1_i]
        atom1_g = get_grid_index_1arity(grid, atom1)
        neighbor_grid_indices = get_neighbor_grid_indices(grid, atom1_g, box is None)
        for g in neighbor_grid_indices:
            for atom2_i in grid['atoms'][g]:
                atom2 = xyz[atom2_i]
                if is_close(bond_def, atom1, atom2):
                    groups.union(atom1_i, atom2_i)
    return groups

def print_bond_def(bond_def):
    print '# element1 - element2 bond_length [Angstrom]'
    for bond in bond_def:
        print '%s - %s %f' % (bond['elem1'], bond['elem2'], bond['length'])

# Unit of length is Angstrom.
default_bond_def = [{'elem1': 'C', 'elem2': 'H', 'length': 1.3},
                    {'elem1': 'C', 'elem2': 'C', 'length': 1.7},
                    {'elem1': 'C', 'elem2': 'O', 'length': 1.7},
                    {'elem1': 'O', 'elem2': 'H', 'length': 1.3}]

if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='')
    parser.add_argument('xyz_path', metavar='XYZ', type=str,
                        help='Input path of XYZ file')
    parser.add_argument('-b', metavar='BOND_DEF', dest='bond_def_path', type=str, default=None,
                        help='See the script source for default setting')
    parser.add_argument('-n', action='store_true', dest='is_non_periodic',
                        default=False, help='')
    parser.add_argument('-o', metavar='OUT', dest='out_path', type=str, default=None,
                        help='Output path of group id file')
    args = parser.parse_args()

    print 'set_group_id start for:', args.xyz_path
    # Read xyz.
    with open(args.xyz_path) as fp:
        xyz, box = read_xyz(fp)
    if box is None:
        print 'Non-periodic mode (no box info in xyz file)'
    elif args.is_non_periodic:
        box = None
        print 'Non-periodic mode (-n option specified)'
    else:
        print 'Periodic mode'

    # Read bond length definition.
    if args.bond_def_path is None:
        bond_def = default_bond_def
    else:
        bond_def = read_bond_def(args.bond_def_path)
    print_bond_def(bond_def)

    # Set output group id path.
    if args.out_path is None:
        out_path = re.sub('\.[^.]+$', '', args.xyz_path) + '_group_id.txt'
    else:
        out_path = args.out_path

    groups = get_groups(xyz, bond_def, box)
    with open(out_path, 'w') as fp:
        groups.print_groups(fp)
