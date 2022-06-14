# @Author: Ilija Medan and John Donor
# @Date: December 29, 2020
# @Filename: RS_to_targetdb.py
# @License: BSD 3-Clause
# @Copyright: Ilija Medan and John Donor

import sys
import argparse
import os
import warnings
import numpy as np

from astropy.io import fits

from sdssdb.peewee.sdss5db import targetdb
import sdss_access.path
from mugatu.designs_to_targetdb import (make_design_field_targetdb,
                                        make_design_assignments_targetdb,
                                        make_desigmmode_results_targetdb)
from mugatu.exceptions import MugatuWarning

sdss_path = sdss_access.path.Path(release='sdss5',
                                  preserve_envvars=True)

if __name__ == '__main__':

    # grabbed the parser from robostratgey code to keep consistent
    # to initiate code use
    # 'python RS_to_targetdb.py -p PLAN -o OBSERVATORY -v target_selection_version_pk'

    parser = argparse.ArgumentParser(
        prog=os.path.basename(sys.argv[0]),
        description='Write Robostratgey outputs to targetDB')

    parser.add_argument('-p', '--plan', dest='plan',
                        type=str, help='name of plan', required=True)
    parser.add_argument('-o', '--observatory', dest='observatory',
                        type=str, help='apo or lco',
                        choices=['apo', 'lco'], required=True)
    parser.add_argument('-t', '--tag', dest='tag',
                        type=str, help='tag for the plan', required=True)

    args = parser.parse_args()
    plan = args.plan
    observatory = args.observatory
    tag = args.tag

    # connect to targetdb
    targetdb.database.connect_from_parameters(user='sdss',
                                              host='operations.sdss.utah.edu',
                                              port=5432)

    # add new robostratgey version to targetDB if it doesnt exist
    ver = targetdb.Version.select().where(targetdb.Version.plan == plan)
    if ver.exists():
        flag = 'Robostrategy plan already in targetdb'
        warnings.warn(flag, MugatuWarning)
    else:
        targetdb.Version = targetdb.Version.create(plan=plan,
                                                   target_selection=False,
                                                   robostrategy=True,
                                                   tag=tag)

        targetdb.Version.save()

    # file with cadences for each field
    allocate_file = sdss_path.full('rsAllocationFinal', plan=plan,
                                   observatory=observatory)
    rsAllocation1 = fits.open(allocate_file)[1].data

    # get the instrument pks
    instr_pks = {}
    instr_pks['BOSS'] = targetdb.Instrument.get(label='BOSS').pk
    instr_pks['APOGEE'] = targetdb.Instrument.get(label='APOGEE').pk

    # get observatory pk
    obsDB = targetdb.Observatory()
    obs_inst = obsDB.get(label=observatory.upper())

    # create dict of fiber pks
    fiber_pks = {}
    holes = (targetdb.Hole.select()
                          .where(targetdb.Hole.observatory ==
                                 obs_inst.pk))
    for hole in holes:
        fiber_pks[hole.holeid] = hole.pk

    # get plan pk
    ver_inst = targetdb.Version.get(plan=plan)

    valid_results = fits.open(('/uufs/chpc.utah.edu/common/home/sdss50/sdsswork/'
                               'sandbox/mugatu/rs_plan_validations/{plan}/'
                               'rs_{plan}_{obs}_design_validation_results.fits'.format(
                                   plan=plan,
                                   obs=observatory)))[1].data

    for allo in rsAllocation1:
        # get fieldid and slots_exposure
        fieldid = allo["fieldid"]
        slots_exposures = [[int(i), int(j)]
                           for i, j in allo["slots_exposures"]]
        # now grab the assignment file for this field
        field_assigned_file = sdss_path.full('rsFieldAssignmentsFinal',
                                             plan=plan,
                                             observatory=observatory,
                                             fieldid=fieldid)

        # get header with field info
        head = fits.open(field_assigned_file)[0].header
        # association between catalogid and instrument
        design_inst = fits.open(field_assigned_file)[1].data
        # catalogid assignment for each fiber
        design = fits.open(field_assigned_file)[2].data
        # get list of designmodes
        desmode_labels = head['DESMODE'].split(' ')

        # use mugatu function to create field in targetdb
        make_design_field_targetdb(cadence=allo['cadence'],
                                   fieldid=fieldid,
                                   plan=ver_inst,
                                   racen=head['RACEN'],
                                   deccen=head['DECCEN'],
                                   position_angle=head['PA'],
                                   observatory=obs_inst,
                                   slots_exposures=slots_exposures)

        fieldid_inst = (targetdb.Field.select()
                                      .join(targetdb.Version)
                                      .switch(targetdb.Field)
                                      .join(targetdb.Cadence)
                                      .where((targetdb.Field.field_id == fieldid) &
                                             (targetdb.Version.plan == plan) &
                                             (targetdb.Cadence.label == allo['cadence'])))

        # get number of exposures
        n_exp = head['NEXP']

        # iterate over exposures for this field entry
        for i in range(allo['iexpst'], allo['iexpnd'] + 1):
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
            design_id = make_design_assignments_targetdb(
                plan=ver_inst,
                fieldid=fieldid_inst,
                exposure=i - allo['iexpst'],  # subtract off the min so it is 0-indexed
                field_exposure=i,
                desmode_label=desmode_label,
                design_ids=design_inst['carton_to_target_pk'],
                robotID=roboIDs,
                holeID=holeIDs,
                obsWavelength=design_inst['fiberType'],
                carton=design_inst['carton'],
                observatory=obs_inst,
                targetdb_ver=None,
                instr_pks=instr_pks,
                cart_pks=design_inst['carton_pk'],
                fiber_pks=fiber_pks,
                idtype='carton_to_target',
                return_design_id=True)

            file = sdss_path.full('rsFieldAssignmentsFinal',
                                  plan=plan,
                                  observatory=observatory,
                                  fieldid=fieldid).split('/')[-1]

            ind = np.where((valid_results['file_name'] == file) &
                           (valid_results['exp'] == i))[0][0]
            make_desigmmode_results_targetdb(
                design_id=design_id,
                design_pass=True,
                design_valid_file_row=valid_results[ind])
