import sys
import argparse
import os
import numpy as np
import glob

from astropy.io import fits

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
            make_design_assignments_targetdb(targetdb_ver=targetdb_ver,
                                             plan=ver_inst,
                                             fieldid=fieldid_inst,
                                             exposure=i,
                                             desmode_label=desmode_label,
                                             design_ids=design_inst['carton_to_target_pk'],
                                             robotID=roboIDs,
                                             holeID=holeIDs,
                                             obsWavelength=design_inst['fiberType'],
                                             carton=design_inst['carton'],
                                             observatory=obs_inst[obs],
                                             instr_pks=instr_pks,
                                             cart_pks=cart_pks,
                                             fiber_pks=fiber_pks,
                                             idtype='carton_to_target')
        # add 1 for next fieldid
        fieldid += 1
