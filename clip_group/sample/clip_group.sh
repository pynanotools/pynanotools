#!/bin/sh
# Clip only carbon atoms (group 2 in the group id file group_id_carbon_in_benzene.txt).
python ../clip_group.py benzene.xyz group_id_carbon_in_benzene.txt 2
# Clip only carbon atoms from a benzene dimer with the same group id file with above.
python ../clip_group.py -r benzene_dimer.xyz group_id_carbon_in_benzene.txt 2

# Do same procedure in gro format.
python ../clip_group.py benzene.gro group_id_carbon_in_benzene.txt 2
python ../clip_group.py -r benzene_dimer.gro group_id_carbon_in_benzene.txt 2
