# distance between two atoms

import numpy as np
from math import sqrt, pi, sin, cos


def centerofmass(xcoords, ycoords, zcoords):
    x_com = np.mean(xcoords)
    y_com = np.mean(ycoords)
    z_com = np.mean(zcoords)
    return [x_com, y_com, z_com]

def distance_between_two_points(xyzcoords_1, xyzcoords_2):
    dist = sqrt((xyzcoords_1[0] - xyzcoords_2[0])**2 + (xyzcoords_1[1] - xyzcoords_2[1])**2 + (xyzcoords_1[2] - xyzcoords_2[2])**2)
    return dist

def curate_pdb_segid(readpdb_object, listofsegids):
    list_of_fields = []
    for ind, segid in enumerate(readpdb_object.atm_segid):
        for segid_from_list in listofsegids:
            if segid == segid_from_list:
                list_of_fields.append([readpdb_object.atm[ind], readpdb_object.atm_sernum[ind], readpdb_object.atm_name[ind],
                                       readpdb_object.atm_altloc[ind],
                                       readpdb_object.atm_resname[ind],
                                       readpdb_object.atm_chainid[ind], readpdb_object.atm_resseqnum[ind],
                                       readpdb_object.atm_cod[ind],
                                       readpdb_object.atm_xcoor[ind], readpdb_object.atm_ycoor[ind],
                                       readpdb_object.atm_zcoor[ind], readpdb_object.atm_occ[ind],
                                       readpdb_object.atm_tempf[ind], segid, readpdb_object.atm_elem[ind],
                                       readpdb_object.atm_charge[ind]])
    return list_of_fields


def curate_crd_segid(readcrd_object, listofsegids):
    list_of_fields = []
    for ind, segid in enumerate(readcrd_object.atm_segid):
        for segid_from_list in listofsegids:
            if segid == segid_from_list:
                list_of_fields.append([readcrd_object.atm_sernum[ind], readcrd_object.resnum_relative[ind],
                                       readcrd_object.atm_resname[ind], readcrd_object.atm_name[ind],
                                       readcrd_object.atm_xcoor[ind], readcrd_object.atm_ycoor[ind],
                                       readcrd_object.atm_zcoor[ind], segid, readcrd_object.atm_resseqnum[ind],
                                       readcrd_object.atm_charge[ind]])
    return list_of_fields



def rotation_matrix(angle=90, direction='x' or 'y' or 'z'):
    angle_rad = angle*pi/180
    cos_theta = cos(angle_rad)
    sin_theta = sin(angle_rad)
    if direction == 'x':
        rot_mat = np.array([[1, 0, 0], [0, cos_theta, sin_theta], [0, -1*sin_theta, cos_theta]], dtype='int')
    elif direction == 'y':
        rot_mat = np.array([[cos_theta, 0, -1*sin_theta], [0, 1, 0], [sin_theta, 0, cos_theta]], dtype='int')
    elif direction == 'z':
        rot_mat = np.array([[cos_theta, sin_theta, 0], [-1*sin_theta, cos_theta, 0], [0, 0, 1]], dtype='int')
    return rot_mat

def rotate_structure(coords, angle, direction):
    rot_mat = rotation_matrix(angle, direction)
    rot_coords = []
    for vector in coords:
        rot_coords.append(np.dot(vector, rot_mat))
    rot_coords = np.array(rot_coords)
    return rot_coords

def translate(coords, delta_x=6.0, delta_y=0.0, delta_z=0.0):
    trans_coords = []
    for vector in coords:
        trans_coords.append([vector[0]+delta_x, vector[1]+delta_y, vector[2]+delta_z])
    trans_coords = np.array(trans_coords)
    return trans_coords

def mergecoords(x, y, z):
    coor = []
    for i in range(len(x)):
        coor.append([x[i], y[i], z[i]])
    coor = np.array(coor, dtype='float')
    return coor

def calculate_helix_dipole(start_res_num, end_res_num, read_pdb_object):
    """
    calculate helix dipole
    :param read_pdb_object: read pdb object
    :param helix_start_res_num: helix start res num
    :param helix_end_res_num: helix end res num
    :return: dictionary of dipole. contains the total dipole along with start res, end res, r_ref for helix, dipole_arr
    """

    # define the backbone atom types
    # backbone atom names - N, HN, CA, HA, C, O

    dipole_dict = dict()

    backbone_atom_names = ['N', 'HN', 'CA', 'HA', 'C', 'O']
    # backbone_atom_names = ['N', 'CA', 'C', 'O']

    atom_names_arr = []
    x_coor_arr = []
    y_coor_arr = []
    z_coor_arr = []
    partial_charge_arr = []



    for ind, res_num in enumerate(read_pdb_object.atm_resseqnum):
        res_num = int(res_num)
        if res_num >= start_res_num:
            if res_num <= end_res_num:
                # print(res_num)
                if read_pdb_object.atm_name[ind] in backbone_atom_names:
                    atom_names_arr.append(read_pdb_object.atm_name[ind])
                    x_coor_arr.append(read_pdb_object.atm_xcoor[ind])
                    y_coor_arr.append(read_pdb_object.atm_ycoor[ind])
                    z_coor_arr.append(read_pdb_object.atm_zcoor[ind])
                    partial_charge_arr.append(read_pdb_object.atm_tempf[ind])


    xyz_ref = centerofmass(x_coor_arr, y_coor_arr, z_coor_arr)

    dipole_atom_arr = []

    for num, (partial_charge, xcoor, ycoor, zcoor) in enumerate(
            zip(partial_charge_arr, x_coor_arr, y_coor_arr, z_coor_arr)):
        atom_r = [xcoor, ycoor, zcoor]
        diff_coor = distance_between_two_points(atom_r, xyz_ref)
        diff_coor_meter = diff_coor * 1e-10
        partial_charge_e = partial_charge * 1.602e-19
        dipole_coul_meter = diff_coor_meter * partial_charge_e
        dipole = dipole_coul_meter * (1/3.33564e-30)
        dipole_atom_arr.append(dipole)

    dipole_ = np.sum(dipole_atom_arr)

    dipole_per_res = dipole_/(end_res_num- start_res_num)

    dipole_dict['ref_point'] = xyz_ref
    dipole_dict['start_res_num'] = start_res_num
    dipole_dict['end_res_num'] = end_res_num
    dipole_dict['dipole_atom_list'] = dipole_atom_arr
    dipole_dict['dipole'] = dipole_
    dipole_dict['dipole_per_res'] = dipole_per_res

    return dipole_dict