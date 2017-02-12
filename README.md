# pynanotools

## About

Collection of 'nano' (tiny) Python scripts for nanomaterial simulations.
Developed on Python 2.7.

## List of tools
- set_group_id.py - Split xyz file into groups
- xyz2gro.py - Convert from xyz file format into gro file format.
- gro2xyz.py - Convert from gro file format into xyz file format.

## Note for the extended xyz file format
Among the tools, the xyz file format is extended for periodic system within orthorhombic cell.
The information of cell lengths, LX, LY, LZ can be written in the header (first line) of the file, 
such as '10 1.0 2.0 3.0'. 
The above header means that the system contains ten atoms and the cell lengths are defined as 
(LX, LY, LZ) = (1.0, 2.0, 3.0)
in the Angstrom unit.


