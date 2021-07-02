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
from mugatu.designs_to_targetdb import make_design_field_targetdb, make_design_assignments_targetdb

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

    args = parser.parse_args()
    plan = args.plan
    observatory = args.observatory

    # COMMENT OUT FOR TEST
    # file with cadences for each field
    allocate_file = sdss_path.full('rsAllocation', plan=plan,
                                   observatory=observatory)

    # get config file path
    # dont think is sdss_path
    stop = 0
    for i in range(len(allocate_file)):
        if allocate_file[i] == '/':
            stop = i + 1
    config_file = allocate_file[:stop] + 'robostrategy-' + plan + '.cfg'

    # connect to targetdb
    targetdb.database.connect_from_parameters(user='sdss',
                                              host='operations.sdss.utah.edu',
                                              port=5432)

    # get targetdb version pk by plan
    targetdb_ver_pk = {}
    ver = targetdb.Version.select()
    for v in ver:
        targetdb_ver_pk[v.plan] = v.pk

    # get the targetdb version pk for each carton from config file
    targetdb_ver = {}
    with open(config_file, 'r') as f:
        read = False
        for x in f:
            # remove \n from line
            if len(x) > 0:
                x = x[:-1]
            # stop reading at blank line
            if x == '':
                read = False
            # read if previously reached cartons line
            if read:
                # skip commented lines
                if x[0] != '#':
                    cart_line = x.split(' = ')
                    targetdb_ver[cart_line[0]] = targetdb_ver_pk[cart_line[1]]
            # start reading lines if reach cartons
            if x == '[Cartons]' or x == '[CartonsExtra]':
                read = True

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
    # cadenceDB = targetdb.Cadence()
    # targetDB = targetdb.Target()
    # carton_to_targetDB = targetdb.CartonToTarget()
    # positionerDB = targetdb.Positioner()
    cartonDB = targetdb.Carton()
    fieldDB = targetdb.Field()
    positionerDB = targetdb.Positioner()

    # create dict of fiber pks
    fiber_pks = {}
    for fp in range(500):
        fiber_pks[fp] = positionerDB.get(id=fp).pk

    # get the instrument pks
    instr_pks = {}
    instr_pks['BOSS'] = targetdb.Instrument.get(label='BOSS').pk
    instr_pks['APOGEE'] = targetdb.Instrument.get(label='APOGEE').pk

    # get observatory pk
    obsDB = targetdb.Observatory()
    obs_inst = obsDB.get(label=observatory.upper())

    # get plan pk
    versionDB = targetdb.Version()
    ver_inst = versionDB.get(plan=plan)

    for fieldid in fieldids:
        # now grab the assignment file for this field
        field_assigned_file = sdss_path.full('rsFieldAssignments',
                                             plan=plan,
                                             observatory=observatory,
                                             fieldid=fieldid)

        # get header with field info
        head = fitsio.read_header(field_assigned_file)
        # association between catalogid and instrument
        design_inst = fitsio.read(field_assigned_file, ext=1)
        # catalogid assignment for each fiber
        design = fitsio.read(field_assigned_file, ext=2)

        # use mugatu function to create field in targetdb
        make_design_field_targetdb(cadence=head['FCADENCE'],
                                   fieldid=fieldid,
                                   plan=ver_inst,
                                   racen=head['RACEN'],
                                   deccen=head['DECCEN'],
                                   position_angle=head['PA'],
                                   observatory=obs_inst)

        fieldid_inst = (fieldDB.select()
                               .join(versionDB)
                               .where((fieldDB.field_id=fieldid) &
                                      (versionDB.plan=plan)))

        # get number of exposures
        try:
            n_exp = len(design['robotID'][0, :])
        except IndexError:
            # some fields only have 1 design which is encoded as 1-D array apparently
            n_exp = 1

        # create the carton pk dict
        cart_pks = {}
        for cart in np.unique(design_inst['carton']):
            # skip calibration from now
            if cart != 'CALIBRATION':
                cart_pks[cart] = (cartonDB.select(cartonDB.pk)
                                          .where((cartonDB.carton == cart) &
                                                 (cartonDB.version_pk == targetdb_ver[cart]))[0].pk)

        for i in range(n_exp):
            # add design for ith exposure to targetdb
            make_design_assignments_targetdb(targetdb_ver=targetdb_ver,
                                             plan=ver_inst,
                                             fieldid=fieldid_inst,
                                             exposure=i,
                                             catalogID=design_inst['catalogid'],
                                             fiberID=design['robotID'][:, i],
                                             obsWavelength=design_inst['fiberType'],
                                             carton=design_inst['carton'],
                                             instr_pks=instr_pks,
                                             cart_pks=cart_pks,
                                             fiber_pks=fiber_pks)
