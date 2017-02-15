# -*- coding: utf-8 -*-
import argparse
import math
import re
import sys

is_debug = False

# Relative atomic mass of elements.
# Retrieved from elses/src/elses-qm-geno-Huckel-atom-parameters.f90
element_to_mass = {'H': 1.0079,
                   'C': 12.0107,
                   'N': 14.0067,
                   'O': 15.9994}

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

def write_xyz(fp, xyz, box):
    if box is None:
        fp.write('%d\n#\n' % len(xyz))
    else:
        fp.write('%d %.16f %.16f %.16f\n#\n' % (len(xyz), box[0], box[1], box[2]))
    for e, x, y, z in xyz:
        fp.write('%s %.16f %.16f %.16f\n' % (e, x, y, z))

def read_group_id(fp):
    is_header_read = False
    dummy_group = -1
    for line in fp:
        if line[0] == '#':  # Comment line.
            continue
        if not is_header_read:
            m = re.search('(\d+)\s+(\d+)\s+(\d+)', line)
            if m:
                num_atoms = int(m.group(1))
                num_groups = int(m.group(2))
                num_max_group_size = int(m.group(3))
                atom_to_group = [dummy_group] * num_atoms
                is_header_read = True
        else:
            m = re.search('(\d+)\s+(\d+)', line)
            # Atom and group are 1-origin in the format, but 0-origin internally.
            if m:
                atom = int(m.group(1)) - 1
                group = int(m.group(2)) - 1
                atom_to_group[atom] = group
    if is_debug:
        print '[Debug] atom_to_group:', atom_to_group
    assert(dummy_group not in atom_to_group)  # All atoms must be belong to a group.
    return {'num_groups': num_groups,
            'atom_to_group': atom_to_group}

def get_distance(atom1, atom2):  # Get distance between atoms without considering box.
    acc = 0.
    for i in range(1, 4):
        acc += (atom1[i] - atom2[i]) ** 2.
    return math.sqrt(acc)

def get_mass_center_core(atoms):  # Get the center of mass for atoms without considering box.
    acc = [0., 0., 0.]  # x, y, z.
    total_mass = 0.
    for e, x, y, z in atoms:
        mass = element_to_mass[e]
        acc[0] += mass * x
        acc[1] += mass * y
        acc[2] += mass * z
        total_mass += mass
    return ('X', acc[0] / total_mass, acc[1] / total_mass, acc[2] / total_mass), total_mass

# Get width of the area where the atoms exist for each axis
# (in other words, the minimum size of covering cuboid).
def get_max_widths(atoms):
    mins = [1e100] * 3
    maxs = [-1e100] * 3
    for e, x, y, z in atoms:
        mins[0] = min(x, mins[0])
        mins[1] = min(y, mins[1])
        mins[2] = min(z, mins[2])
        maxs[0] = max(x, maxs[0])
        maxs[1] = max(y, maxs[1])
        maxs[2] = max(z, maxs[2])
    return [maxs[0] - mins[0], maxs[1] - mins[1], maxs[2] - mins[2]]

def get_mass_center(group, atoms, is_non_periodic, box):
    num_atoms = len(atoms)
    atom_represent = atoms[num_atoms / 2]  # Representation atom in the group.
    # In non-periodic system, shift for atoms means the same shift for center of mass.
    if is_non_periodic:
        def shift(atom, is_forward):
            return atom
    else:
        def shift(atom, is_forward):
            if is_forward:  # Forward shift to set the representation atom to the center of the box.
                x = atom[1] - atom_represent[1] + box[0] / 2.
                y = atom[2] - atom_represent[2] + box[1] / 2.
                z = atom[3] - atom_represent[3] + box[2] / 2.
                # Wrap over the box.
                x = x % box[0]
                y = y % box[1]
                z = z % box[2]
            else:  # In backward shift, wrapping is not needed.
                x = atom[1] + atom_represent[1] - box[0] / 2.
                y = atom[2] + atom_represent[2] - box[1] / 2.
                z = atom[3] + atom_represent[3] - box[2] / 2.
            return (atom[0], x, y, z)
    atoms_shifted = map(lambda a: shift(a, True), atoms)
    # Check whether the shift was appropriate.
    if not is_non_periodic:
        max_widths = get_max_widths(atoms_shifted)
        for axis, max_width, box_width in zip(['x', 'y', 'z'], max_widths, box):
            if max_width > box_width / 2.:
                print '[Warning] too distant atom from the representation atom in %s-axis' % axis
                print '[Warning]   group: %d' % (group + 1)
                print '[Warning]   representation atom: %s (%f, %f, %f)' % atom_represent
                print '[Warning]   half box width: %f' % (box_width / 2.)
                print '[Warning]   max width: %f' % max_width
    mass_center_shifted, total_mass = get_mass_center_core(atoms_shifted)
    mass_center = shift(mass_center_shifted, False)
    return mass_center, total_mass

if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='')
    parser.add_argument('xyz_path', metavar='XYZ', type=str,
                        help='Input path of XYZ file')
    parser.add_argument('-g', metavar='GROUP_ID', dest='group_id_path', type=str, default=None,
                        help='Group id file path. If not specified, all atoms belong to one group')
    parser.add_argument('-n', action='store_true', dest='is_non_periodic',
                        default=False, help='Force non-periodic mode')
    parser.add_argument('-o', metavar='OUTPUT', dest='output_path', type=str, default=None,
                        help='Output file path')
    parser.add_argument('--debug', action='store_true', dest='is_debug_mode',
                        default=False, help='Run in debug mode')
    args = parser.parse_args()

    if args.is_debug_mode:
        is_debug = True

    if is_debug:
        print 'Python mod specification check 1: -0.9 % 1. =', -0.9 % 1.
        print 'Python mod specification check 2: -1.1 % 1. =', -1.1 % 1.

    with open(args.xyz_path) as fp:
        xyz, box = read_xyz(fp)
    num_atoms_all = len(xyz)
    is_non_periodic = args.is_non_periodic or (box is None)
    if args.group_id_path is not None:
        with open(args.group_id_path) as fp:
            group_id = read_group_id(fp)
    else:
        group_id = {'num_groups': 1, 'atom_to_group': [0] * num_atoms_all}

    # Set output path.
    if args.output_path is None:
        output_path = re.sub('\.[^.]+$', '', args.xyz_path) + '_mass_center.txt'
    else:
        output_path = args.output_path

    # Print basic information.
    print 'number of atoms:', num_atoms_all
    print 'is non-periodic mode:', is_non_periodic
    if not is_non_periodic:
        print 'box:', box[0], box[1], box[2]
    if args.group_id_path is not None:
        print 'group_id path:', args.group_id_path
    print 'output path:', output_path

    # Group atoms by indices in group_id.
    group_to_atoms = []
    for group in range(group_id['num_groups']):
        group_to_atoms.append([])
    for group, atom in zip(group_id['atom_to_group'], xyz):
        group_to_atoms[group].append(atom)

    # Main procedure.
    mass_centers = []
    for group, atoms in enumerate(group_to_atoms):
        mass_centers.append(get_mass_center(group, atoms, is_non_periodic, box))

    # Output xyz representing result for debug.
    if is_debug:
        mass_centers_as_atoms = map(lambda m: m[0], mass_centers)
        with open('debug.xyz', 'w') as fp:
            write_xyz(fp, mass_centers_as_atoms, box)

    # Output result.
    with open(output_path, 'w') as fp:
        fp.write('# group x y z total_mass\n')
        for group, (mass_center, total_mass) in enumerate(mass_centers):
            x, y, z = mass_center[1], mass_center[2], mass_center[3]
            fp.write('%8d %14.6f %14.6f %14.6f %14.6f\n' % (group + 1, x, y, z, total_mass))
