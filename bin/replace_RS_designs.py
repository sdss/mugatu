import sys
import argparse
import os
import numpy as np
import glob
import datetime

from astropy.io import fits
from astropy.table import Table

from sdssdb.peewee.sdss5db import targetdb
from mugatu.designs_to_targetdb import (make_design_field_targetdb,
                                        make_design_assignments_targetdb,
                                        TargetdbFieldIDs,
                                        make_desigmmode_results_targetdb,
                                        design_status_bitmask,
                                        make_designToField)
from mugatu.exceptions import MugatuError


if __name__ == '__main__':
    parser = argparse.ArgumentParser(
        prog=os.path.basename(sys.argv[0]),
        description='In a batch, validate a set of designs')
    parser.add_argument('-l', '--loc', dest='loc',
                        type=str, help='local or utah',
                        choices=['local', 'utah'], required=True)
    parser.add_argument('-p', '--plan', dest='plan',
                        type=str, help='name of RS plan',
                        required=True)
    parser.add_argument('-f', '--fieldids', dest='fieldids', nargs='+',
                        help='field_ids to replace)',
                        type=int, required=True)

    args = parser.parse_args()
    loc = args.loc
    plan = args.plan
    fieldids = args.fieldids

    if loc == 'local':
        targetdb.database.connect_from_parameters(user='sdss',
                                                  host='localhost',
                                                  port=7502)
    else:
        targetdb.database.connect_from_parameters(user='sdss',
                                                  host='operations.sdss.utah.edu',
                                                  port=5432)
    # get the files with designs
    replace_path = ('/uufs/chpc.utah.edu/common/home/sdss50/sdsswork/'
                    'target/robostrategy_replacement/{plan}/'.format(
                        plan=plan))
    file_save = replace_path + plan + '_field_replacement_log.fits'
    files = []
    for fid in fieldids:
        files += [file for file in glob.glob(replace_path +
                                             '{plan}_{fid}*.fits'.format(
                                                 plan=plan,
                                                 fid=fid))]
    for f in files:
        if 'validation' in f or 'status' in f:
            files.remove(f)
    for f in files:
        if 'validation' in f or 'status' in f:
            files.remove(f)
    files_valid = []
    for f in files:
        files_valid.append(f[:-5] + '_validation.fits')
    # get status files
    files_status = []
    for f in files:
        if os.path.exists(f[:-5] + '_designid_status.fits'):
            files_status.append(f[:-5] + '_designid_status.fits')
        else:
            raise MugatuError(message='No designid_status file for field')

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
    ver_inst = targetdb.Version.get(plan=plan)

    # add designs in the files
    for file, valid, status in zip(files, files_valid, files_status):
        # get header with field info
        head = fits.open(file)[0].header
        # association between catalogid and instrument
        design_inst = fits.open(file)[1].data
        # catalogid assignment for each fiber
        design = fits.open(file)[2].data
        # get list of designmodes
        desmode_labels = head['DESMODE'].split(' ')
        # get validation data
        valid_data = fits.open(valid)[1].data

        # get field info to be replaced
        fieldid_replace = int(file.split(replace_path)[1].split('_')[1])
        field_replace = (targetdb.Field.select().where(targetdb.Field.field_id == fieldid_replace,
                                                       targetdb.Field.version == ver_inst.pk)[0])
        # create replacement field
        # get new field_id only if racen, deccen and PA change
        if (head['RACEN'] == field_replace.racen and
            head['DECCEN'] == field_replace.deccen and
            head['PA'] == field_replace.position_angle):
            fieldid_new = fieldid_replace
        else:
            same_field = targetdb.Field.select().where(head['RACEN'] == targetdb.Field.racen,
                                                       head['DECCEN'] == targetdb.Field.deccen,
                                                       head['PA'] == targetdb.Field.position_angle,
                                                       targetdb.Field.field_id >= 100000)
            if len(same_field) > 0:
                fieldid_new = same_field[0].field_id
            else:
                fieldid_new = targetdb.FieldReservation.requestNext(
                     N=1, commit=True, commissioning=False)[0]
        racen_new = head['RACEN']
        deccen_new = head['DECCEN']
        pa_new = head['PA']
        # use mugatu function to create field in targetdb
        make_design_field_targetdb(cadence=field_replace.cadence.label,
                                   fieldid=fieldid_new,
                                   plan=ver_inst,
                                   racen=racen_new,
                                   deccen=deccen_new,
                                   position_angle=pa_new,
                                   observatory=obs_inst[field_replace.observatory.label],
                                   slots_exposures=field_replace.slots_exposures,
                                   replacement_field=True)

        # get info on field created
        fieldid_inst = (targetdb.Field.select()
                                      .join(targetdb.Version)
                                      .switch(targetdb.Field)
                                      .join(targetdb.Cadence)
                                      .where((targetdb.Field.field_id == fieldid_new) &
                                             (targetdb.Version.plan == plan) &
                                             (targetdb.Cadence.label == field_replace.cadence.label)))
        if len(fieldid_inst) > 1:
            pks = []
            for f in fieldid_inst:
                pks.append(f.pk)
            fieldid_inst = (targetdb.Field.select()
                                          .where(targetdb.Field.pk == max(pks)))

        # get number of exposures
        n_exp = head['NEXP']
        if os.path.exists(status):
            design_ids = fits.open(status)['STATUS'].data['designid']
        else:
            design_ids = np.zeros(n_exp, dtype=int) - 1

        # get number of exposures to be replaced
        des = targetdb.DesignToField.select().where(targetdb.DesignToField.field == field_replace.pk)
        n_exp_expc = len(des)
        # ensure exposures is correct
        if n_exp == n_exp_expc:
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
                if design_ids[i] == -1:
                    design_id = make_design_assignments_targetdb(
                        plan=ver_inst,
                        fieldid=fieldid_inst,
                        exposure=i,
                        field_exposure=i,  # safe to assume same for these?
                        desmode_label=desmode_label,
                        design_ids=design_inst['carton_to_target_pk'],
                        robotID=roboIDs,
                        holeID=holeIDs,
                        obsWavelength=design_inst['fiberType'],
                        carton=design_inst['carton'],
                        observatory=obs_inst[field_replace.observatory.label],
                        targetdb_ver=None,
                        instr_pks=instr_pks,
                        cart_pks=None,
                        fiber_pks=fiber_pks[field_replace.observatory.label],
                        idtype='carton_to_target',
                        return_design_id=True)
                    make_desigmmode_results_targetdb(
                        design_id=design_id,
                        design_pass=True,
                        design_valid_file_row=valid_data[i],
                        design_status=design_status_bitmask(replacement_design=True))
                else:
                    make_designToField(design=int(design_ids[i]),
                                       fieldid=fieldid_inst,
                                       exposure=i,
                                       field_exposure=i)
        elif n_exp == 1:
            for i in range(n_exp_expc):
                # repeat same design so always use only index
                roboIDs = design['robotID']
                holeIDs = design['holeID']
                desmode_label = desmode_labels[0]
                # write exposure to targetdb
                if design_ids[0] == -1:
                    design_id = make_design_assignments_targetdb(
                        plan=ver_inst,
                        fieldid=fieldid_inst,
                        exposure=i,
                        field_exposure=i,  # safe to assume same for these?
                        desmode_label=desmode_label,
                        design_ids=design_inst['carton_to_target_pk'],
                        robotID=roboIDs,
                        holeID=holeIDs,
                        obsWavelength=design_inst['fiberType'],
                        carton=design_inst['carton'],
                        observatory=obs_inst[field_replace.observatory.label],
                        targetdb_ver=None,
                        instr_pks=instr_pks,
                        cart_pks=None,
                        fiber_pks=fiber_pks[field_replace.observatory.label],
                        idtype='carton_to_target',
                        return_design_id=True)
                    make_desigmmode_results_targetdb(
                        design_id=design_id,
                        design_pass=True,
                        design_valid_file_row=valid_data[0],
                        design_status=design_status_bitmask(replacement_design=True))
                else:
                    make_designToField(design=int(design_ids[0]),
                                       fieldid=fieldid_inst,
                                       exposure=i,
                                       field_exposure=i)
        else:
            raise MugatuError(message='Incorrect number of exposures in design')
        # save field_id info
        dtype = np.dtype([('file_name', '<U50'),
                          ('field_id_old', np.int32),
                          ('field_pk_old', np.int32),
                          ('racen_old', np.float64),
                          ('deccen_old', np.float64),
                          ('pa_old', np.float64),
                          ('field_id_new', np.int32),
                          ('field_pk_new', np.int32),
                          ('racen_new', np.float64),
                          ('deccen_new', np.float64),
                          ('pa_new', np.float64),
                          ('date_replaced', '<U30')])
        field_arr = np.zeros(1, dtype=dtype)
        field_arr['file_name'][0] = os.path.split(file)[-1]
        # record old info
        field_arr['field_id_old'][0] = fieldid_replace
        field_arr['field_pk_old'][0] = field_replace.pk
        field_arr['racen_old'][0] = field_replace.racen
        field_arr['deccen_old'][0] = field_replace.deccen
        field_arr['pa_old'][0] = field_replace.position_angle
        # add new info
        field_arr['field_id_new'][0] = fieldid_new
        field_arr['field_pk_new'][0] = fieldid_inst[0].pk
        field_arr['racen_new'][0] = fieldid_inst[0].racen
        field_arr['deccen_new'][0] = fieldid_inst[0].deccen
        field_arr['pa_new'][0] = fieldid_inst[0].position_angle
        field_arr['date_replaced'][0] = str(datetime.datetime.now())
        if os.path.exists(file_save):
            curr_data = fits.open(file_save)[1].data
            curr_data_arr = np.zeros(len(curr_data), dtype=dtype)
            for i in range(len(curr_data)):
                for col in curr_data.columns:
                    curr_data_arr[col.name][i] = curr_data[col.name][i]
            field_arr = np.append(curr_data_arr, field_arr)
            field_arr = Table(field_arr)
            field_arr.write(file_save, format='fits', overwrite=True)
        else:
            field_arr = Table(field_arr)
            field_arr.write(file_save, format='fits')
