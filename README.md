# pynanotools

## set_group_id.py
```
usage: python set_group_id.py [-h] [-b BOND_DEF] [-o OUT] XYZ
```
Split an xyz file into some groups of atoms (for example, molecules).
The result of split is in the group_id format of ELSES.
For each pair of elements, a threshold for bond length is defined.
If there are a pair of atoms which is the distance between them is smaller than the threshold, the script infer that they are in the same group.
A default threshold table `default_bond_def` written in the script is used unless a user-defined one is given by `-b` option.
The following text is a sample bond definition file.
The unit of bond length is Angstrom.
```
C H 1.3
C C 1.6
C O 1.6
```
