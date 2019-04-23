#pdb utility classes
#pdb read
#pdb write

import os
import datetime
import getpass
import numpy as np
from math import sqrt, pi, sin, cos
from fortranformat import FortranRecordReader, FortranRecordWriter

class ParsePDB(object):
    """
    Pdb parser object
    """

    def __init__(self, file, dirpath):
        """
        initialize the reader with file and directory of where the pdb file is
        and store the relative information
        :param file: input pdb file name
        :param dirpath: input directory path
        """
        self.pdbfile = file
        self.dirpath = dirpath

        # information to collect

        self.atm = []
        self.atm_sernum = []
        self.atm_space0 = []
        self.atm_name = []
        self.atm_altloc = []
        self.atm_resname = []
        self.atm_space1 = []
        self.atm_chainid = []
        self.atm_resseqnum = []
        self.atm_cod = []
        self.atm_space2 = []
        self.atm_xcoor = []
        self.atm_ycoor = []
        self.atm_zcoor = []
        self.atm_occ = []
        self.atm_tempf = []
        self.space3 = []
        self.atm_segid = []
        self.atm_elem = []
        self.atm_charge = []


    def read(self):
        with open(os.path.join(self.dirpath, self.pdbfile)) as pdbfile:
            pdbfile = pdbfile.read().splitlines(keepends=False)
            for line in pdbfile:
                if line.startswith('ATOM'):
                    self.atm.append(line[0:6].strip())
                    self.atm_sernum.append(line[6:11].strip())
                    self.atm_space0.append(line[12])
                    self.atm_name.append(line[12:16].strip())
                    self.atm_altloc.append(line[16])
                    self.atm_resname.append(line[17:21].strip())
                    self.atm_space1.append(line[20])
                    self.atm_chainid.append(line[21])
                    self.atm_resseqnum.append(line[22:26].strip())
                    self.atm_cod.append(line[26])
                    self.atm_space2.append(line[27:30])
                    self.atm_xcoor.append(float(line[30:38]))
                    self.atm_ycoor.append(float(line[38:46]))
                    self.atm_zcoor.append(float(line[46:54]))
                    self.atm_occ.append(float(line[54:60]))
                    self.atm_tempf.append(float(line[60:66]))
                    self.space3.append(line[66:72])
                    self.atm_segid.append(line[72:76].strip())
                    self.atm_elem.append(line[76:78].strip())
                    self.atm_charge.append(line[78:80])
        return self

class ParseCRD(object):
    """
    Parse CHARMM coordinate file .crd
    """

    def __init__(self, file, dirpath):
        """
        Initialize the reader with file name and directory
        :param file: filename
        :param dirpath: directorly absolute location
        """
        self.crdfile = file
        self.dirloc = dirpath

        # information to collect
        self.atm_sernum = []
        self.resnum_relative = []
        self.atm_resname = []
        self.atm_name = []
        self.atm_xcoor = []
        self.atm_ycoor = []
        self.atm_zcoor = []
        self.atm_segid = []
        self.atm_resseqnum = []
        self.atm_charge = []


    def read(self):
        with open(os.path.join(self.dirloc, self.crdfile)) as coorfile:
            coorfile = coorfile.read().splitlines(keepends=False)
            for line in coorfile[4:]:
                self.atm_sernum.append(int(line[:10].strip()))
                self.resnum_relative.append(int(line[10:20].strip()))
                self.atm_resname.append(str(line[22:27].strip()))
                self.atm_name.append(str(line[32:37].strip()))
                self.atm_xcoor.append(float(line[44:61].strip()))
                self.atm_ycoor.append(float(line[64:81].strip()))
                self.atm_zcoor.append(float(line[84:101].strip()))
                self.atm_segid.append(str(line[102:107].strip()))
                self.atm_resseqnum.append(int(line[110:116].strip()))
                self.atm_charge.append(float(line[124:141].strip()))

        return self

def write_crd(list_of_fields, fname, dirpath):
    crdwriter = FortranRecordWriter('(5X, I5, 5X, I5, 2X, A4, 6X, A4, 10X, F14.10, 6X, F14.10, 6X, F14.10, 2X, A4, 6X, A4, 11X, F13.10)')
    crdwriter_pre = FortranRecordWriter('(5X, I5, 1X, A4)')
    now = datetime.datetime.now()
    with open(os.path.join(dirpath, fname), 'w') as writecrd:
        writecrd.write('{} {}\n'.format('*', 'Suggie CRD writer'))
        writecrd.write('{} {} {} {} {} {}\n'.format('*', 'Date: ', now.strftime("%Y-%m-%d"),
                                                    now.strftime("%H:%M"), 'User: ', getpass.getuser()))
        writecrd.write('{}\n'.format('*'))

        numpy_list = np.asarray(list_of_fields)
        fortranformat_pre = crdwriter_pre.write([numpy_list[:, 0][-1], 'EXT'])
        writecrd.write(fortranformat_pre+'\n')

        for ind in range(len(list_of_fields)):
            atm_sernum = list_of_fields[ind][0]
            resnum_relative = list_of_fields[ind][1]
            atm_resname = list_of_fields[ind][2].ljust(4)
            atm_name = list_of_fields[ind][3].ljust(4)
            atm_xcoor = list_of_fields[ind][4]
            atm_ycoor = list_of_fields[ind][5]
            atm_zcoor = list_of_fields[ind][6]
            atm_segid = list_of_fields[ind][7]
            atm_resseqnum = list_of_fields[ind][8].ljust(4)
            atm_charge = list_of_fields[ind][9]

            fortranformat = crdwriter.write([atm_sernum, resnum_relative, atm_resname, atm_name, atm_xcoor, atm_ycoor, atm_zcoor,
                                             atm_segid, atm_resseqnum, atm_charge])

            writecrd.write(fortranformat+'\n')

    writecrd.close()

def convert_pdb_to_crd(readpdbobject, dirpath):
    list_of_fields = []
    for ind, atm_sernum in enumerate(readpdbobject.atm_sernum):
        list_of_fields.append([atm_sernum, readpdbobject.atm_resseqnum[ind],
                                       readpdbobject.atm_resname[ind], readpdbobject.atm_name[ind],
                                       readpdbobject.atm_xcoor[ind], readpdbobject.atm_ycoor[ind],
                                       readpdbobject.atm_zcoor[ind], readpdbobject.atm_segid[ind], readpdbobject.atm_resseqnum[ind],
                                       readpdbobject.atm_charge[ind]])

    outfilename = readpdbobject.pdbfile.strip('.')[0]+'_crd.crd'
    write_crd(list_of_fields, outfilename, dirpath)



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



def write_pdb(list_of_fields, fname, dirpath, elem_charge=True):
    now = datetime.datetime.now()
    with open(os.path.join(dirpath, fname), 'w') as writepdb:
        writepdb.write('{} {}\n'.format('REMARK', 'Suggie PDB writer'))
        writepdb.write('{} {} {} {} {} {}\n'.format('REMARK', 'Date: ', now.strftime("%Y-%m-%d"),
                                                           now.strftime("%H:%M"), 'User: ', getpass.getuser()))
        for ind in range(len(list_of_fields)):
            atm = list_of_fields[ind][0].ljust(6)
            atm_sernum = list_of_fields[ind][1].rjust(5)
            atm_space0 = ' '
            atm_name = list_of_fields[ind][2].ljust(4)
            atm_altloc = list_of_fields[ind][3]
            atm_resname = list_of_fields[ind][4].ljust(4)
            atm_space1 = ''
            atm_chainid = list_of_fields[ind][5]
            atm_resseqnum = list_of_fields[ind][6].rjust(4)
            atm_cod = list_of_fields[ind][7]
            atm_space2 = '   '
            atm_xcoor = str('%8.3f' % (float(list_of_fields[ind][8]))).rjust(8)
            atm_ycoor = str('%8.3f' % (float(list_of_fields[ind][9]))).rjust(8)
            atm_zcoor = str('%8.3f' % (float(list_of_fields[ind][10]))).rjust(8)
            atm_occ = str('%6.2f' % (float(list_of_fields[ind][11]))).rjust(6)
            atm_tempf = str('%6.2f' % (float(list_of_fields[ind][12]))).rjust(6)
            atm_space3 = '      '
            atm_segid = list_of_fields[ind][13].rjust(4)
            atm_elem = list_of_fields[ind][2][0].rjust(2)
            atm_charge = list_of_fields[ind][15].ljust(2)

            if elem_charge == True:
                writepdb.write(atm + atm_sernum + atm_space0 + atm_name + atm_altloc + atm_resname + atm_space1 +
                           atm_chainid + atm_resseqnum + atm_cod + atm_space2 + atm_xcoor + atm_ycoor + atm_zcoor +
                           atm_occ + atm_tempf + atm_space3 + atm_segid + atm_elem + atm_charge + '\n')
            else:
                writepdb.write(atm + atm_sernum + atm_space0 + atm_name + atm_altloc + atm_resname + atm_space1 +
                               atm_chainid + atm_resseqnum + atm_cod + atm_space2 + atm_xcoor + atm_ycoor + atm_zcoor +
                               atm_occ + atm_tempf + atm_space3 + atm_segid + '\n')

        ter = 'TER   ' + str(int(list_of_fields[-1][1]) + 1).rjust(5) + '      ' + str(list_of_fields[-1][4]).rjust(3) \
              + ' ' + str(list_of_fields[-1][5]) + list_of_fields[-1][6].rjust(4) + ' \n'
        writepdb.write(ter)
        end_model = 'ENDMDL\n'
        writepdb.write(end_model)
        writepdb.close()

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

def rotate_structure(coords, rot_mat):
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

    # print(dipole_)

    dipole_dict['ref_point'] = xyz_ref
    dipole_dict['start_res_num'] = start_res_num
    dipole_dict['end_res_num'] = end_res_num
    dipole_dict['dipole_atom_list'] = dipole_atom_arr
    dipole_dict['dipole'] = dipole_
    dipole_dict['dipole_per_res'] = dipole_per_res
    # print('heho')

    return dipole_dict




if __name__=='__main__':
    dirpath = r"C:\Users\sugyan\Documents\MembraneMDfiles\Serf_tmp_model\Replica_Exchange_Analysis_mod\step2_remd_coor_final"
    pdbfile = "step2_remd_tm_end_0.pdb"

    read_mod = ParsePDB(pdbfile, dirpath).read()