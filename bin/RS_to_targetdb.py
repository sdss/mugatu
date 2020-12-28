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

    # file with cadences for each field
    allocate_file = sdss_path.full('rsAllocation', plan=plan,
                                   observatory=observatory)

    # connect to targetdb
    targetdb.database.connect_from_parameters(user='sdss_user',
                                              host='operations.sdss.utah.edu',
                                              port=5432)

    # add new robostratgey version to targetDB
    versionDB = targetdb.Version.create(plan=plan,
                                        target_selection=False,
                                        robostrategy=True,
                                        tag='')  # what is this supposed to be?

    versionDB.save()

    # data describing field (besides PA) stored here
    rsAllocation1 = fitsio.read(allocate_file, ext=1)
    # PAs are stored here
    rsAllocation3 = fitsio.read(allocate_file, ext=3)

    # need these databased for fks that will be stored in these tables
    cadenceDB = targetdb.Cadence()
    targetDB = targetdb.Target()
    positionerDB = targetdb.Positioner()
    postion_infoDB = targetdb.PositionerInfo()

    # get the instrument pks
    instr_pks = {}
    instr_pks['BOSS'] = targetdb.Instrument.get(label='BOSS').pk
    instr_pks['APOGEE'] = targetdb.Instrument.get(label='APOGEE').pk

    # get observatory pk
    obsDB = targetdb.Observatory()
    obspk = obsDB.get(label=observatory).pk

    # get plan pk
    versionDB = targetdb.Version()
    verpk = obsDB.get(plan=plan).pk

    for allo1, allo3 in zip(rsAllocation1, rsAllocation3):
        try:
            dbCadence = cadenceDB.get(label=allo1['cadence']).pk
        except:  # dont know how to catch custom error 'CadenceDoesNotExist'
            print('No cadence found in db for field %d' % allo1['fieldid'])
            continue

        # creates new row in database
        fieldDB = targetdb.Field.create(
            pk=allo1['fieldid'],
            racen=allo1['racen'],
            deccen=allo1['deccen'],
            position_angle=allo3['pa2'],
            cadence=dbCadence,
            observatory=obspk,
            version=verpk)
        # save row in database
        fieldDB.save()

        # now grab the assignment file for this field
        field_assigned_file = sdss_path.full('rsFieldAssignments',
                                             plan=plan,
                                             observatory=observatory,
                                             fieldid=allo1['fieldid'])

        # association between catalogid and instrument
        design_inst = fitsio.read(field_assigned_file, ext=1)
        # catalogid assignment for each fiber
        design = fitsio.read(field_assigned_file, ext=2)

        for i in range(len(design[0, :])):

            # creates new row in design database
            designDB = targetdb.Design.create(field=allo1['fieldid'],
                                              exposure=i)
            # save row
            designDB.save()

            # add the assignments for the design to the assignment database
            rows = []
            for j in range(len(design[:, i])):
                row_dict = {}

                # get the pk for the positioner_info
                # (where I assume the ID is just the row # in the fits file)
                this_pos_DB = positionerDB.get(id=j).pk

                # get the instrument for fiber
                inst_assign = design_inst['fiberType'][
                    design_inst['catalogid'] == design[j][i]][0]

                # add db row info to dic
                row_dict['design'] = designDB.pk
                row_dict['instrument'] = instr_pks[inst_assign]
                row_dict['positioner'] = this_pos_DB
                row_dict['target'] = targetDB.get(catalogid=design[j][i]).pk

                rows.append(row_dict)

            # write all exposures for field to targetdb
            targetdb.Assignment.insert_many(rows).execute()
