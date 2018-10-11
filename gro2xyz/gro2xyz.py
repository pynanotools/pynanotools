# coding: utf-8
import argparse
import os
import os.path
import re
import sys

xyz_line_template = '%10s %8.5f %8.5f %8.5f\n'
XYZ_HEADER_WIDTH = 59
ANGSTROM_PER_NM = 10.0

def read_xyz(fp_xyz_check):
    check_xyz = []
    atom_info = fp_xyz_check.readline()
    comment = fp_xyz_check.readline()
    for line in fp_xyz_check:
        ss = line.split()
        element = ss[0]
        x = float(ss[1])
        y = float(ss[2])
        z = float(ss[3])
        coordinates = [x, y, z]
        check_xyz.append((element, coordinates))
    return atom_info, comment, check_xyz

def write_xyz(fp_write_xyz, atom_info, comment, check_xyz):
    fp_write_xyz.write("%s" % (atom_info))
    fp_write_xyz.write("%s" % (comment))
    for e, c in check_xyz:
        fp_write_xyz.write("%-3s %8.5f %8.5f %8.5f\n" % (e, c[0], c[1], c[2]))


def main(fp_gro, fp_xyz):
    gro_region = 'head'
    for line in fp_gro:
        if gro_region == 'head':
            gro_comment = line.strip()
            xyz_header_offset = fp_xyz.tell()
            fp_xyz.write(' ' * (XYZ_HEADER_WIDTH + len(gro_comment)))  # This region is overwritten after.
            gro_region = 'num_atoms'
            continue
        elif gro_region == 'num_atoms':
            num_atoms = int(line)
            atom_index = 0
            gro_region = 'body'
        elif gro_region == 'body':
            assert(atom_index < num_atoms)
            x, y, z = [float(line[20:28]), float(line[28:36]), float(line[36:44])]
            atom = line[8:15].strip()
            fp_xyz.write(xyz_line_template%(atom,
                                            x * ANGSTROM_PER_NM,
                                            y * ANGSTROM_PER_NM,
                                            z * ANGSTROM_PER_NM))
            atom_index += 1
            if atom_index == num_atoms:
                gro_region = 'box'
        elif gro_region == 'box':
            ss = line.split()
            cell_info = [float(ss[0]), float(ss[1]), float(ss[2])]
            xyz_tail_offset = fp_xyz.tell()
            fp_xyz.seek(xyz_header_offset)
            # Fixed character width (XYZ_HEADER_WIDTH).
            fp_xyz.write('%10d %14.5f %14.5f %14.5f\n# %s\n'%(num_atoms,
                                                              cell_info[0] * ANGSTROM_PER_NM,
                                                              cell_info[1] * ANGSTROM_PER_NM,
                                                              cell_info[2] * ANGSTROM_PER_NM,
                                                              gro_comment))
            fp_xyz.seek(xyz_tail_offset)
            gro_region = 'head'

if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='')
    parser.add_argument('input_path', metavar='GRO', type=str,
                        help='Input gro file path')
    parser.add_argument('-o', metavar='OUT', dest='output_path', type=str, default=None,
                        help='Output xyz file path')
    parser.add_argument('-f', action='store_true', dest='to_force_overwrite', default=False,
                        help='Force overwrite to existing file')
    args = parser.parse_args()

    if args.output_path is None:
        output_path = re.sub('\.[^.]+$', '', args.input_path) + '.xyz'
    else:
        output_path = args.output_path
    print 'in: ', args.input_path
    print 'out:  ', output_path
    if os.path.exists(output_path) and not args.to_force_overwrite:
        print 'output file "%s" already exists. Specify -f to force overwriting' % output_path
        sys.exit(1)
        
    with open(args.input_path) as fp_gro, open(output_path, 'w') as fp_xyz:
        main(fp_gro, fp_xyz)
    with open(output_path) as fp_xyz_check:
        atom_info, comment, check_xyz = read_xyz(fp_xyz_check)
    os.remove(output_path)
    with open(output_path, "w") as fp_write_xyz:
        write_xyz(fp_write_xyz, atom_info, comment, check_xyz)
