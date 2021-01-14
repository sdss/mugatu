# @Author: Ilija Medan and John Donor
# @Date: December 29, 2020
# @Filename: RS_to_targetdb.py
# @License: BSD 3-Clause
# @Copyright: Ilija Medan and John Donor

import sys
import argparse
import os

from sdssdb.peewee.sdss5db import targetdb
import sdss_access.path
import fitsio

sdss_path = sdss_access.path.Path(release='sdss5')


if __name__ == '__main__':

    # grabbed the parser from robostratgey code to keep consistent
    # to initiate code use 'python RS_to_targetdb.py -p PLAN -o OBSERVATORY'

    parser = argparse.ArgumentParser(
        prog=os.path.basename(sys.argv[0]),
        description='Write Robostratgey outputs to targetDB')

    parser.add_argument('-p', '--plan', dest='plan',
                        type=str, help='name of plan', required=True)
    parser.add_argument('-o', '--observatory', dest='observatory',
                        type=str, help='apo or lco',
                        choices=['apo', 'lco'], required=True)

    args = parser.parse_args()
    plan = args.plan
    observatory = args.observatory

    # COMMENT OUT FOR TEST
    # file with cadences for each field
    # allocate_file = sdss_path.full('rsAllocation', plan=plan,
    #                               observatory=observatory)

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
                                            tag='test')  # add test flag for now

        versionDB.save()

    # COMMENT OUT FOR TEST
    # data describing field (besides PA) stored here
    # rsAllocation1 = fitsio.read(allocate_file, ext=1)
    # PAs are stored here
    # rsAllocation3 = fitsio.read(allocate_file, ext=3)

    # USING FOR TEST, CHANGE LATER
    fieldids = [4373, 4673]

    # need these databased for fks that will be stored in these tables
    cadenceDB = targetdb.Cadence()
    targetDB = targetdb.Target()
    positionerDB = targetdb.Positioner()

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
        field_assigned_file = field_assigned_file.replace('/targets', '')
        field_assigned_file = field_assigned_file[:-5] + '-example.fits'

        # get header with field info
        head = fitsio.read_header(field_assigned_file)
        # association between catalogid and instrument
        design_inst = fitsio.read(field_assigned_file, ext=1)
        # catalogid assignment for each fiber
        design = fitsio.read(field_assigned_file, ext=2)

        try:
            dbCadence = cadenceDB.get(label=head['FCADENCE']).pk
        except:  # dont know how to catch custom error 'CadenceDoesNotExist'
            print('No cadence found in db for field %d' % fieldid)
            # COMMENT OUT FOR TEST
            # continue
            # USING FOR TEST, CHANGE LATER
            dbCadence = -1

        # creates new field in database if it doesnt exist
        try:
            fieldDB = targetdb.Field.create(
                pk=fieldid,
                racen=head['RACEN'],
                deccen=head['DECCEN'],
                position_angle=head['PA'],
                cadence=dbCadence,
                observatory=obspk,
                version=verpk)
            # save row in database
            fieldDB.save()
        except:  # again, need to load correct error here
            pass

        # get number of exposures
        n_exp = len(design['robotID'][0, :])

        for i in range(n_exp):

            # creates new row in design database
            # NEED TO CHECK IF DESIGN ALREADY EXISTS?
            designDB = targetdb.Design.create(field=fieldid,
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
                    row_dict['target'] = (targetDB.get(
                        catalogid=design_inst['catalogid'][j]).pk)

                    rows.append(row_dict)

            # write all exposures for field to targetdb
            targetdb.Assignment.insert_many(rows).execute()
