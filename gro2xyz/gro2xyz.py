# coding: utf-8
import os.path
import re
import sys

xyz_line_template = '%-3s %8.5f %8.5f %8.5f\n'
XYZ_HEADER_WIDTH = 59
ANGSTROM_PER_NM = 10.0

def main(input_gro_path):
    output_xyz_path = re.sub('\.[^.]+$', '', input_gro_path) + '.xyz'
    print 'in: ', input_gro_path
    print 'out:  ', output_xyz_path
    assert(not (os.path.exists(output_xyz_path)))
    with open(input_gro_path, 'r') as fp_gro:
        with open(output_xyz_path, 'w') as fp_xyz:
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
    if len(sys.argv) != 2:
        print 'usage: python %s <GRO>' % sys.argv[0]
        sys.exit(0)
    input_gro_path = sys.argv[1]
    main(input_gro_path)
