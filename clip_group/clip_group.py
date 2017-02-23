# coding: utf-8
import argparse
import os.path
import re
import sys

NUM_ATOMS_DIGITS = 14
NUM_ATOMS_TEMPLATE = '%%%dd\n' % NUM_ATOMS_DIGITS

def check_atom_and_group_size(is_group_repeated, num_atoms, num_atoms_in_group):
    if is_group_repeated and num_atoms % num_atoms_in_group != 0:
        print '[Error] the number of atoms must be a multiple of the group size'
        print 'the number of atoms:', num_atoms
        print 'group size:', num_atoms_in_group
        sys.exit(1)
    elif not is_group_repeated and num_atoms != num_atoms_in_group:
        print '[Error] the number of atoms must be the same with the group size'
        print 'the number of atoms:', num_atoms
        print 'group size:', num_atoms_in_group
        sys.exit(1)

def main_gro(fp_in, fp_out, is_group_repeated, num_atoms_in_group, index_to_clipped):
    region = 'head'
    for line in fp_in:
        if region == 'head':
            fp_out.write(line)  # Comment.
            region = 'num_atoms'
        elif region == 'num_atoms':
            num_atoms_offset = fp_out.tell()
            fp_out.write(' ' * (NUM_ATOMS_DIGITS + 1))  # '+ 1' is for \n.
            num_atoms = int(line)
            check_atom_and_group_size(is_group_repeated, num_atoms, num_atoms_in_group)
            # Do not need conditional branch for is_group_repeated below (same for fp_out.write in body block).
            fp_out.write(NUM_ATOMS_TEMPLATE % (num_atoms / num_atoms_in_group * len(index_to_clipped)))
            atom_index = 0
            region = 'body'
        elif region == 'body':
            assert(atom_index < num_atoms)
            if atom_index % num_atoms_in_group in index_to_clipped:
                fp_out.write(line)
            atom_index += 1
            if atom_index == num_atoms:
                region = 'box'
        elif region == 'box':
            fp_out.write(line)
            region = 'head'  # Go to next step.

def main_xyz(fp_in, fp_out, is_group_repeated, num_atoms_in_group, index_to_clipped):
    region = 'head'
    for line in fp_in:
        if region == 'head':
            ss = line.split()
            num_atoms = int(ss[0])
            check_atom_and_group_size(is_group_repeated, num_atoms, num_atoms_in_group)
            # Do not need conditional branch for is_group_repeated below (same for fp_out.write in body block).
            num_atoms_clipped = num_atoms / num_atoms_in_group * len(index_to_clipped)
            if len(ss) == 1:
                fp_out.write('%d\n' % num_atoms_clipped)
            else:  # Extended xyz file (with box size).
                fp_out.write('%d %s\n' % (num_atoms_clipped, ' '.join(ss[1 :])))
            atom_index = 0
            region = 'comment'
        elif region == 'comment':
            fp_out.write(line)  # Copy comment.
            region = 'body'
        elif region == 'body':
            assert(atom_index < num_atoms)
            if atom_index % num_atoms_in_group in index_to_clipped:
                fp_out.write(line)
            atom_index += 1
            if atom_index == num_atoms:
                region = 'head'  # Go to next step.

def read_group_id(path, index):
    index_to_clipped = {}
    is_header_read = False
    with open(path) as fp:
        for line in fp:
            ss = line.split()
            if len(line) == 0 or line[0] == '#':  # Comment line.
                continue
            elif not is_header_read:
                num_atoms_in_group = int(ss[0])
                is_header_read = True
                assert(1 <= index and index <= int(ss[1]))  # Check whether group index is in collect range.
            elif is_header_read and int(ss[1]) == index:
                index_to_clipped[int(ss[0]) - 1] = True  # Convert 1-origin to 0-origin. True is a dummy value.
    return num_atoms_in_group, index_to_clipped

if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='')
    parser.add_argument('input_structure_path', metavar='XYZ/GRO', type=str,
                        help='Input structure file (format is xyz or gro)')
    parser.add_argument('group_id_path', metavar='GROUP_ID', type=str,
                        help='Input group ID file')
    parser.add_argument('group_index', metavar='GROUP', type=int,
                        help='Group index in group id file (1-origin)')
    parser.add_argument('-o', metavar='OUT', dest='output_path', type=str, default=None,
                        help='Output structure file')
    parser.add_argument('-r', action='store_true', dest='is_group_repeated', default=False,
                        help='Assume gro file is some number of repeat of group and copy group_id for them')
    parser.add_argument('-f', action='store_true', dest='to_force_overwrite', default=False,
                        help='Force overwrite to existing file')
    args = parser.parse_args()

    num_atoms_in_group, index_to_clipped = read_group_id(args.group_id_path, args.group_index)

    # Decide atomic structure file format. Default format is xyz and gro is optional.
    structure_format = 'xyz'
    if re.search('.gro\Z', args.input_structure_path):
        structure_format = 'gro'
    print 'in: ', args.input_structure_path

    if args.output_path is None:
        output_path = re.sub('\.[^.]+$', '', args.input_structure_path) + '_clip.' + structure_format
    else:
        output_path = args.output_path
    print 'out:  ', output_path
    if os.path.exists(output_path) and not args.to_force_overwrite:
        print 'output file "%s" already exists. Specify -f to force overwriting' % output_path
        sys.exit(1)

    with open(args.input_structure_path, 'r') as fp_in, open(output_path, 'w') as fp_out:
        if structure_format == 'xyz':
            main_xyz(fp_in, fp_out, args.is_group_repeated, num_atoms_in_group, index_to_clipped)
        elif structure_format == 'gro':
            main_gro(fp_in, fp_out, args.is_group_repeated, num_atoms_in_group, index_to_clipped)
        else:
            assert(False)  # Unknown format.
