# coding: utf-8
import re
import sys

NM_PER_ANGSTROM = 0.1

if __name__ == '__main__':
    input_xyz_path = sys.argv[1]
    with open(input_xyz_path, 'r') as fp:
        line_num = 0
        atoms = []
        for line in fp:
            if line_num == 0:
                xyz_header = line.split()
                assert(len(xyz_header) == 1 or len(xyz_header) == 4)
                atom_num = int(xyz_header[0])
                if len(xyz_header) == 4:
                    box = map(lambda s: float(s), xyz_header[1 : 4])
                else:
                    box = None
            elif line_num == 1:
                comment = line
            elif line_num < atom_num + 2:
                xyz_body = line.split()
                element = xyz_body[0]
                coordinates = map(lambda s: float(s), xyz_body[1 : 4])
                atoms.append((element, coordinates))
            else:
                break
            line_num += 1

    output_gro_path = re.sub('\.[^.]+$', '.gro', input_xyz_path)
    with open(output_gro_path, 'w') as fp:
        template = "%8s%7s%5d%8.3f%8.3f%8.3f\n"
        if comment[0] == '#':
            comment = 'dummy\n'
        fp.write(comment)
        fp.write('    %d\n' % len(atoms))
        for i, (element, coordinates) in zip(range(len(atoms)), atoms):
            if True:
                tail = '  '
            else:
                tail = 'AX'  # Comes from Morino's bachelor thesis "A.1.10 convert xyz2gro.py".
            fp.write(template%("1PDB",
                               element + tail.strip(),
                               (i + 1) % 100000,
                               coordinates[0] * NM_PER_ANGSTROM,
                               coordinates[1] * NM_PER_ANGSTROM,
                               coordinates[2] * NM_PER_ANGSTROM))
        if box is None:
            xs = map(lambda a: a[1][0] * NM_PER_ANGSTROM, atoms)
            ys = map(lambda a: a[1][1] * NM_PER_ANGSTROM, atoms)
            zs = map(lambda a: a[1][2] * NM_PER_ANGSTROM, atoms)
            x_width = max(xs) - min(xs)
            y_width = max(ys) - min(ys)
            z_width = max(zs) - min(zs)
            fp.write('  %f  %f  %f\n' % (x_width, y_width, z_width))
        else:
            fp.write('  %f  %f  %f\n' % (box[0] * NM_PER_ANGSTROM,
                                         box[1] * NM_PER_ANGSTROM,
                                         box[2] * NM_PER_ANGSTROM))
