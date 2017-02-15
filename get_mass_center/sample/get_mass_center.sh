#!/bin/sh
# Run with printing debug information. debug.xyz is generated.
python ../get_mass_center.py --debug -g benzene_dimer_group_id.txt benzene_dimer.xyz

# No group ID file (all in one group) mode, and change output file path.
python ../get_mass_center.py -o benzene_dimer_mass_center_no_group.txt benzene_dimer.xyz

# This will generate warnings, but result is correct.
python ../get_mass_center.py -g benzene_dimer_group_id.txt benzene_dimer_warn.xyz
