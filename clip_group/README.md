# clip_group.py
```
usage: python clip_group.py [-h] [-o OUT] [-r] [-f] XYZ/GRO GROUP_ID GROUP
```
Clip atoms in a group from an atomic structure file.
`XYZ/GRO` is the path to the input structure file. xyz and gro file formats are available.
`GROUP_ID` is the path to group ID file describing atom-to-group mapping.
`GROUP` is the index of the clipped group.

When `-r` option is specified, atom-to-group mapping is repeated for M / N times
where M is the number of all atoms in .gro
and N is the number of atoms in the group.
For example, if the body of the group ID file is
```
1 1
2 2
3 2
```
and M = 6, atom indices in group 1 are {1, 4}
and those in group 2 are {2, 3, 5, 6}.
