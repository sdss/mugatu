import sys
import argparse
import os
import numpy as np
import glob

from astropy.io import fits
from astropy.table import Table

from sdssdb.peewee.sdss5db import targetdb
from mugatu.designs_to_targetdb import (make_design_field_targetdb,
                                        make_design_assignments_targetdb,
                                        TargetdbFieldIDs)


if __name__ == '__main__':
    parser = argparse.ArgumentParser(
        prog=os.path.basename(sys.argv[0]),
        description='In a batch, validate a set of designs')
    parser.add_argument('-l', '--loc', dest='loc',
                        type=str, help='local or utah',
                        choices=['local', 'utah'], required=True)
    parser.add_argument('-d', '--dir', dest='dir',
                        type=str, help='directory with design files',
                        required=True)

    args = parser.parse_args()
    loc = args.loc
    directory = args.dir

    if loc == 'local':
        targetdb.database.connect_from_parameters(user='sdss',
                                                  host='localhost',
                                                  port=7502)
    else:
        targetdb.database.connect_from_parameters(user='sdss',
                                                  host='operations.sdss.utah.edu',
                                                  port=5432)
    # get the files to ingest in the directory
    files = [file for file in glob.glob(directory + '*.fits')]

    # find the minimum field_id to start with
    fieldid = (TargetdbFieldIDs(fieldid_type='manual',
                                version_plan='manual').find_next_available())

    # get observatory insts
    obsDB = targetdb.Observatory()
    obs_inst = {}
    obs_inst['APO'] = obsDB.get(label='APO')
    obs_inst['LCO'] = obsDB.get(label='LCO')

    # get the instrument pks
    instr_pks = {}
    instr_pks['BOSS'] = targetdb.Instrument.get(label='BOSS').pk
    instr_pks['APOGEE'] = targetdb.Instrument.get(label='APOGEE').pk

    # create dict of fiber pks
    fiber_pks = {}
    fiber_pks['APO'] = {}
    fiber_pks['LCO'] = {}
    # get APO holes
    holes = (targetdb.Hole.select()
                          .where(targetdb.Hole.observatory ==
                                 obs_inst['APO'].pk))
    for hole in holes:
        fiber_pks['APO'][hole.holeid] = hole.pk
    # get LCO holes
    holes = (targetdb.Hole.select()
                          .where(targetdb.Hole.observatory ==
                                 obs_inst['LCO'].pk))
    for hole in holes:
        fiber_pks['LCO'][hole.holeid] = hole.pk

    # get plan pk
    ver_inst = targetdb.Version.get(plan='manual')

    # add designs in the files
    for file in files:
        slots_exposures = np.zeros((24, 2), dtype=int)

        # get header with field info
        head = fits.open(file)[0].header
        obs = head['obs'].strip().upper()
        # association between catalogid and instrument
        design_inst = fits.open(file)[1].data
        # catalogid assignment for each fiber
        design = fits.open(file)[2].data
        # get list of designmodes
        desmode_labels = head['DESMODE'].split(' ')

        # use mugatu function to create field in targetdb
        make_design_field_targetdb(cadence=head['FCADENCE'],
                                   fieldid=fieldid,
                                   plan=ver_inst,
                                   racen=head['RACEN'],
                                   deccen=head['DECCEN'],
                                   position_angle=head['PA'],
                                   observatory=obs_inst[obs],
                                   slots_exposures=slots_exposures)

        fieldid_inst = (targetdb.Field.select()
                                      .join(targetdb.Version)
                                      .where((targetdb.Field.field_id == fieldid) &
                                             (targetdb.Version.plan == plan)))

        # get number of exposures
        try:
            n_exp = len(design['robotID'][0, :])
        except IndexError:
            # some fields only have 1 design which is encoded as 1-D array apparently
            n_exp = 1

        for i in range(n_exp):
            # index correctly based on n_exp
            if n_exp == 1:
                roboIDs = design['robotID']
                holeIDs = design['holeID']
                desmode_label = desmode_labels[0]
            else:
                roboIDs = design['robotID'][:, i]
                holeIDs = design['holeID'][:, i]
                desmode_label = desmode_labels[i]
            # write exposure to targetdb
            design_id = make_design_assignments_targetdb(plan=ver_inst,
                                                         fieldid=fieldid_inst,
                                                         exposure=i,
                                                         desmode_label=desmode_label,
                                                         design_ids=design_inst['carton_to_target_pk'],
                                                         robotID=roboIDs,
                                                         holeID=holeIDs,
                                                         obsWavelength=design_inst['fiberType'],
                                                         carton=design_inst['carton'],
                                                         observatory=obs_inst[obs],
                                                         targetdb_ver=None,
                                                         instr_pks=instr_pks,
                                                         cart_pks=None,
                                                         fiber_pks=fiber_pks[obs],
                                                         idtype='carton_to_target')
            # create new entry in summary file
            dtype = np.dtype([('file_name', '<U50'),
                              ('exp', np.int32),
                              ('racen', np.float64),
                              ('deccen', np.float64),
                              ('designmode', '<U15'),
                              ('design_id', np.int32)])
            save_arr0 = np.zeros(1, dtype=dtype)
            save_arr0['file_name'][0] = os.path.split(file)[-1]
            save_arr0['exp'][0] = i + 1
            save_arr0['racen'][0] = head['RACEN']
            save_arr0['deccen'][0] = head['DECCEN']
            save_arr0['designmode'][0] = desmode_label
            save_arr0['design_id'][0] = design_id
            if 'save_arr' in locals():
                save_arr = np.append(save_arr,
                                     save_arr0)
            else:
                save_arr = save_arr0
        # add 1 for next fieldid
        fieldid += 1
    # write to fits file
    save_arr = Table(save_arr)
    save_arr.write(directory +
                   'design_ids_for_design_files.fits',
                   format='fits')
