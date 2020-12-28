import sys
import argparse
import numpy as np

from sdssdb.peewee.sdss5db import targetdb
import sdss_access.path

sdss_path = sdss_access.path.Path(release='sdss5')


if __name__ == '__main__':

    # grabbed the parser from robostratgey code to keep consistent
    # this will fidn the correct files based on the robostrategy plan and observatory used as input
    # the way this code is writting right now, I think you would have to run apo first then lco

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

    allocate_file = sdss_path.full('rsAllocation', plan=plan,
                                   observatory=observatory)

    assignments_file = sdss_path.full('rsAssignments', plan=plan,
                                      observatory=observatory)

    cadences_file = sdss_path.full('rsCadences', plan=plan,
                                   observatory=observatory)

    targetdb.database.connect_from_parameters(user='sdss_user',
                                              host='operations.sdss.utah.edu',
                                              port=5432)

    # add new robostratgey version to targetDB
    versionDB = targetdb.Version.create(plan=plan,
                                        target_selection=False,
                                        robostrategy=True,
                                        tag='')  # what is this supposed to be?

    versionDB.save()

    # write the candeces to the cadenceDB

    # rsCadences=fitsio.read(cadences_file,ext=1)

    # need to pull out the pks for the instruments
    # instDB=targetdb.Instrument()

    # bosspk=instDB.get(instDB.label=='BOSS').pk
    # apogeepk=instDB.get(instDB.label=='APOGEE').pk

    # for cadence in rsCadences:
    #     #need to get instrument pk (so go from label to pk above)
    #     #setting blank indexes as -1, I dont know the right solution for this
    #     inst_pk=np.zeros(len(rsCadences['INSTRUMENT']),dtype=int)-1
    #     inst_pk[rsCadences['INSTRUMENT']=='BOSS']=bosspk
    #     inst_pk[rsCadences['INSTRUMENT']=='APOGEE']=apogeepk

    #     #creates new row in database
    #     cadenceDB=targetdb.Cadence.create(delta = rsCadences['DELTA'],
    #         delta_max = rsCadences['DELTA_MAX'],
    #         delta_min = rsCadences['DELTA_MIN'],
    #         instrument_pk = inst_pk,
    #         label = rsCadences['CADENCE'],
    #         nexposures = rsCadences['NEXPOSURES'],
    #         skybrightness = rsCadences['SKYBRIGHTNESS'])
    #     #save row in database
    #     cadenceDB.save()

    # write the fields to fieldDB and designs to designDB and
    # assignments all at once (I think this makes sense to do together)

    # most data needed is here
    rsAllocation1 = fitsio.read(allocate_file, ext=1)
    # PAs are stored here
    rsAllocation3 = fitsio.read(allocate_file, ext=3)

    # need these databased for fks that will be stored in these tables
    cadenceDB = targetdb.Cadence()
    targetDB = targetdb.Target()
    positionerDB = targetdb.Positioner()
    postion_infoDB = targetdb.PositionerInfo()

    # get observatory pk
    obsDB = targetdb.Observatory()
    obspk = obsDB.get(label=observatory).pk

    # get plan pk
    versionDB = targetdb.Version()
    verpk = obsDB.get(plan=plan).pk

    for allo1, allo3 in zip(rsAllocation1, rsAllocation3):
        try:
            dbCadence = cadenceDB.get(label=allo1['cadence']).pk
        except:
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

        # grab the second HDU to get the number of exposures needed
        # for the field (so number of designs in the field)
        # this will also have the catalogid assignments for each fiber
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
                this_pos_DB = positionerDB.get(id=j)

                row_dict['design'] = designDB.pk
                row_dict['instrument'] = ??
                row_dict['positioner'] = this_pos_DB
                row_dict['target'] = targetDB.get(catalogid=design[j][i])

                rows.append(row_dict)

                # check in the positioner_info DB is this fiber is an apogee or boss fiber
                # (and set key to whatever it coresponds to)

                # if postion_infoDB.get(postion_infoDB.pk=pos_infopk).apogee:
                #     inst_key=apogeepk
                # else:
                #     inst_key=bosspk

                # add the assignment
                # assignDB = targetdb.Assignment.create(
                #    design=designDB,  # can I call the row like this?
                #    instrument=??,
                #    positioner=this_pos_DB,
                #    target=targetDB.get(targetDB.catalogid=design[j][i]))

                # assignDB.save()

            targetdb.Assignment.insert_many(rows).execute()
