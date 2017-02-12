#!/bin/sh
python ../set_group_id.py benzene_dimer.xyz
# Output group_id file path can be changed by -o option.
python ../set_group_id.py -o group_id_test.txt F_Trip_03pvl_6mer_3x3x3.xyz
# Bond length definition can be given by -b option.
python ../set_group_id.py -b bond_def.txt pentacene_14112atom.xyz
