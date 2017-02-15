# get_mass_center.py
```
usage: python get_mass_center.py [-h] [-g GROUP_ID] [-n] [-o OUTPUT] [--debug] XYZ
```
Calculate the center of mass from input xyz file.

If a group ID file is given (`-g` option), the center of atom for each group is calculated.
For extended xyz file, periodic boundary condition is imposed.
Specify `-n` option to cancel bounday condition mode.

When `--debug` option is specified, `debug.xyz` is generated.
It contains imaginary atoms that represent the centers of mass for the atom groups.
