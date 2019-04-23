# PDB_Utility

Contains modules to read pdb files and do many calculations and manipulations to the structure

pdb_read.py ----> reads pdb file and stores all the information in an object. This can then be used to do calculations or other manipulations

utilities.py ----> has functions to do calculations and manipulations.

List of functions

1. centerofmass: calculates center of mass

2. distance_between_two_points: calculates distance between two points

3. curate_pdb_segid: curates pdb with a list of given segment ids. Used for CHARMM pdb files where segment id is specified

4. rotate_structure: rotate the structure with certain angle and at certain direction (x, y, z)

5. translate: translate the coordinates in x, y, or z direction

6. mergecoords: merge a set of coordinates

7. calculate_helix_dipole: calculates dipole (ideally used for helix structure) for a segment of pdb giving a start and end residue number and the pdb object (from pdb_read.py)

