# @Author: Ilija Medan and John Donor
# @Date: December 29, 2020
# @Filename: RS_to_targetdb.py
# @License: BSD 3-Clause
# @Copyright: Ilija Medan and John Donor

import sys
import argparse
import os
import numpy as np

import fitsio

from sdssdb.peewee.sdss5db import targetdb
import sdss_access.path

sdss_path = sdss_access.path.Path(release='sdss5')

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
    parser.add_argument('-v', '--targetdb_ver', dest='targetdb_ver',
                        type=str, help='target_selection version_pk',
                        required=True)

    args = parser.parse_args()
    plan = args.plan
    observatory = args.observatory
    targetdb_ver = args.targetdb_ver

    # COMMENT OUT FOR TEST
    # file with cadences for each field
    allocate_file = sdss_path.full('rsAllocation', plan=plan,
                                   observatory=observatory)

    # connect to targetdb
    targetdb.database.connect_from_parameters(user='sdss',
                                              host='operations.sdss.utah.edu',
                                              port=5432)

    # add new robostratgey version to targetDB if it doesnt exist
    try:
        versionDB = targetdb.Version()
        verpk = versionDB.get(plan=plan).pk
    except:
        versionDB = targetdb.Version.create(plan=plan,
                                            target_selection=False,
                                            robostrategy=True,
                                            tag="test")  # add test flag for now

        versionDB.save()

    # COMMENT OUT FOR TEST
    # data describing field (besides PA) stored here
    rsAllocation1 = fitsio.read(allocate_file, ext=1)
    # PAs are stored here
    # rsAllocation3 = fitsio.read(allocate_file, ext=3)

    # USING FOR TEST, CHANGE LATER
    # fieldids = [4373, 4673]
    fieldids = rsAllocation1["fieldid"]

    # need these databased for fks that will be stored in these tables
    cadenceDB = targetdb.Cadence()
    targetDB = targetdb.Target()
    carton_to_targetDB = targetdb.CartonToTarget()
    positionerDB = targetdb.Positioner()
    cartonDB = targetdb.Carton()

    # get the instrument pks
    instr_pks = {}
    instr_pks['BOSS'] = targetdb.Instrument.get(label='BOSS').pk
    instr_pks['APOGEE'] = targetdb.Instrument.get(label='APOGEE').pk

    # get observatory pk
    obsDB = targetdb.Observatory()
    obspk = obsDB.get(label=observatory.upper()).pk

    # get plan pk
    versionDB = targetdb.Version()
    verpk = versionDB.get(plan=plan).pk

    for fieldid in fieldids:
        # now grab the assignment file for this field
        field_assigned_file = sdss_path.full('rsFieldAssignments',
                                             plan=plan,
                                             observatory=observatory,
                                             fieldid=fieldid)
        # USING FOR TEST, CHANGE LATER
        # field_assigned_file = field_assigned_file.replace('/targets', '')
        # field_assigned_file = field_assigned_file[:-5] + '-example.fits'

        # get header with field info
        head = fitsio.read_header(field_assigned_file)
        # association between catalogid and instrument
        design_inst = fitsio.read(field_assigned_file, ext=1)
        # catalogid assignment for each fiber
        design = fitsio.read(field_assigned_file, ext=2)

        # print("loading field ", fieldid)
        dbCadence = cadenceDB.get(label=head['FCADENCE']).pk

        # grab all carton pks here
        cart_pks = {}
        for cart in np.unique(design_inst['carton']):
            # skip calibration from now
            if cart != 'CALIBRATION':
                cart_pks[cart] = (cartonDB.select(cartonDB.pk)
                                          .where((cartonDB.carton == cart) &
                                                 (cartonDB.version_pk == targetdb_ver))[0].pk)

        # try:
        #     dbCadence = cadenceDB.get(label=head['FCADENCE']).pk
        # except:  # dont know how to catch custom error 'CadenceDoesNotExist'
        #     print('No cadence found in db for field %d' % fieldid)
        #     # COMMENT OUT FOR TEST
        #     # continue
        #     # USING FOR TEST, CHANGE LATER
        #     dbCadence = -1

        # check if field exists
        field_test = (targetdb.Field
                      .select()
                      .where((targetdb.Field.field_id == fieldid) &
                             (targetdb.Field.version == verpk)))
        # creates new field in database if it doesnt exist
        if len(field_test) == 0:
            fieldDB = targetdb.Field.create(
                field_id=fieldid,
                racen=head['RACEN'],
                deccen=head['DECCEN'],
                position_angle=head['PA'],
                cadence=dbCadence,
                observatory=obspk,
                version=verpk)
            # save row in database
            fieldDB.save()

        # get number of exposures
        try:
            n_exp = len(design['robotID'][0, :])
        except IndexError:
            # some fields only have 1 design which is encoded as 1-D array apparently
            n_exp = 1

        for i in range(n_exp):
            # print("field {} pk {} design: {}".format(fieldid, fieldDB.pk, i))

            # creates new row in design database
            # NEED TO CHECK IF DESIGN ALREADY EXISTS?
            designDB = targetdb.Design.create(field=fieldDB.pk,
                                              exposure=i)
            # save row
            designDB.save()

            # add the assignments for the design to the assignment database
            rows = []
            for j in range(len(design['robotID'][:, i])):
                row_dict = {}

                # right now calibrations are fake, so need to skip
                if design['robotID'][j][i] != -1 and design_inst['program'][j] != 'CALIBRATION':
                    # get the pk for the positioner_info
                    # (where I assume the ID is just the
                    # row # in the fits file)

                    this_pos_DB = (positionerDB.get(
                        id=design['robotID'][j][i]).pk)

                    # get the instrument for fiber
                    inst_assign = design_inst['fiberType'][j]

                    # add db row info to dic
                    row_dict['design'] = designDB.pk
                    row_dict['instrument'] = instr_pks[inst_assign]
                    row_dict['positioner'] = this_pos_DB
                    cart_pk = cart_pks[design_inst['carton'][j]]
                    row_dict['carton_to_target'] = (carton_to_targetDB.select(
                        carton_to_targetDB.pk)
                        .join(targetDB,
                              on=(carton_to_targetDB.target_pk == targetDB.pk))
                        .where((targetDB.catalogid == design_inst['catalogid'][j]) &
                               (carton_to_targetDB.carton_pk == cart_pk))[0].pk)

                    rows.append(row_dict)

            # write all exposures for field to targetdb
            targetdb.Assignment.insert_many(rows).execute()
