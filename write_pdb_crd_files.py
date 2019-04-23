import os
import datetime
import getpass
import numpy as np
from fortranformat import FortranRecordReader, FortranRecordWriter



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

            ter = 'TER   ' + str(int(list_of_fields[-1][1]) + 1).rjust(5) + '      ' + str(list_of_fields[-1][4]).rjust(
                3) \
                  + ' ' + str(list_of_fields[-1][5]) + list_of_fields[-1][6].rjust(4) + ' \n'
            writepdb.write(ter)
            end_model = 'ENDMDL\n'
            writepdb.write(end_model)
            writepdb.close()