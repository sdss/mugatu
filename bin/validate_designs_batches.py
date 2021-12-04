import sys
import argparse
import os
import numpy as np
import glob
import time

from astropy.io import fits
from astropy.table import Table

from sdssdb.peewee.sdss5db import targetdb
from sdssdb.peewee.sdss5db import catalogdb
import sdss_access.path
import robostrategy.obstime as obstime
import coordio.time

from mugatu.fpsdesign import FPSDesign
from mugatu.exceptions import MugatuDesignError
from mugatu.designmode import DesignModeCheck
from mugatu.designmode import build_brigh_neigh_query, allDesignModes
from multiprocessing import Pool
from itertools import repeat


def validate_design(design_file, exp, obsTime,
                    db_query_results_boss, db_query_results_apogee,
                    desmode_manual):
    """
    Validate a design and record any errors or warnings
    from the validation
    """
    des = FPSDesign(design_pk=-1,
                    obsTime=obsTime,
                    design_file=design_file,
                    manual_design=True,
                    exp=exp)
    # set default
    decolide = True
    bright_safety = True
    # build design
    des.build_design_manual()
    # try to validate design and catch any design errors
    try:
        des.validate_design(db_query_results_boss=db_query_results_boss,
                            db_query_results_apogee=db_query_results_apogee,
                            desmode_manual=desmode_manual.todict())
    except MugatuDesignError as e:
        if 'Kaiju' in str(e):
            decolide = False
        if 'Bright' in str(e):
            bright_safety = False
    return des, decolide, bright_safety


def design_outputs_to_array(des, decolide,
                            bright_safety,
                            db_query_results_boss,
                            db_query_results_apogee,
                            desmode_manual):
    """
    Output validation parameters as a structured array
    """
    dtype = np.dtype([('file_name', '<U50'),
                      ('exp', np.int32),
                      ('racen', np.float64),
                      ('deccen', np.float64),
                      ('designmode', '<U15'),
                      ('decolide', bool),
                      ('bright_safety', bool),
                      ('bright_safety_pass', np.int32),
                      ('bright_safety_total', np.int32),
                      ('all_targets_assigned', bool),
                      ('all_targets_assigned_pass', np.int32),
                      ('all_targets_assigned_total', np.int32),
                      ('no_collisions', bool),
                      ('no_collisions_pass', np.int32),
                      ('no_collisions_total', np.int32),
                      ('boss_n_skies_min', bool),
                      ('boss_n_skies_min_value', np.int32),
                      ('apogee_n_skies_min', bool),
                      ('apogee_n_skies_min_value', np.int32),
                      ('boss_min_skies_fovmetric', bool),
                      ('boss_min_skies_fovmetric_value', np.float64),
                      ('apogee_min_skies_fovmetric', bool),
                      ('apogee_min_skies_fovmetric_value', np.float64),
                      ('boss_n_stds_min', bool),
                      ('boss_n_stds_min_value', np.int32),
                      ('apogee_n_stds_min', bool),
                      ('apogee_n_stds_min_value', np.int32),
                      ('boss_min_stds_fovmetric', bool),
                      ('boss_min_stds_fovmetric_value', np.float64),
                      ('apogee_min_stds_fovmetric', bool),
                      ('apogee_min_stds_fovmetric_value', np.float64),
                      ('boss_stds_mags', bool),
                      ('boss_stds_mags_pass', np.int32),
                      ('boss_stds_mags_total', np.int32),
                      ('apogee_stds_mags', bool),
                      ('apogee_stds_mags_pass', np.int32),
                      ('apogee_stds_mags_total', np.int32),
                      ('boss_bright_limit_targets', bool),
                      ('boss_bright_limit_targets_pass', np.int32),
                      ('boss_bright_limit_targets_total', np.int32),
                      ('apogee_bright_limit_targets', bool),
                      ('apogee_bright_limit_targets_pass', np.int32),
                      ('apogee_bright_limit_targets_total', np.int32),
                      ('boss_sky_neighbors_targets', bool),
                      ('boss_sky_neighbors_targets_pass', np.int32),
                      ('boss_sky_neighbors_targets_total', np.int32),
                      ('apogee_sky_neighbors_targets', bool),
                      ('apogee_sky_neighbors_targets_pass', np.int32),
                      ('apogee_sky_neighbors_targets_total', np.int32)])
    valid_arr = np.zeros(1, dtype=dtype)
    valid_arr['file_name'][0] = os.path.split(des.design_file)[-1]
    if des.exp == 0:
        valid_arr['exp'][0] = des.exp + 1
    else:
        valid_arr['exp'][0] = des.exp
    valid_arr['racen'][0] = des.racen
    valid_arr['deccen'][0] = des.deccen
    valid_arr['designmode'][0] = des.desmode_label
    valid_arr['decolide'][0] = decolide
    valid_arr['bright_safety'][0] = bright_safety
    # check the design warnings
    column_names = ['all_targets_assigned', 'no_collisions',
                    'boss_n_skies_min',
                    'apogee_n_skies_min', 'boss_min_skies_fovmetric',
                    'apogee_min_skies_fovmetric',
                    'boss_n_stds_min', 'apogee_n_stds_min',
                    'boss_min_stds_fovmetric', 'apogee_min_stds_fovmetric', 
                    'boss_stds_mags', 'apogee_stds_mags',
                    'boss_bright_limit_targets', 'apogee_bright_limit_targets',
                    'boss_sky_neighbors_targets', 'apogee_sky_neighbors_targets']
    warnings_order = ['all_targets_assigned', 'no_collisions', 'min_skies_boss',
                      'min_skies_apogee', 'fov_skies_boss', 'fov_skies_apogee',
                      'min_stds_boss', 'min_stds_apogee', 'fov_stds_boss',
                      'fov_stds_apogee', 'stds_mag_boss', 'stds_mag_apogee',
                      'sci_mag_boss', 'sci_mag_apogee', 'bright_neigh_boss',
                      'bright_neigh_apogee']
    for c, k in zip(column_names, warnings_order):
        valid_arr[c][0] = des.design_errors[k]
    # add in the metrics
    if not valid_arr['bright_safety'][0]:
        mode = DesignModeCheck(FPSDesign=des,
                               desmode_label=des.desmode_label,
                               db_query_results_boss=db_query_results_boss,
                               db_query_results_apogee=db_query_results_apogee,
                               desmode_manual=desmode_manual)
        bright_check_boss, hasFiber_boss = mode.bright_neighbors(instrument='BOSS',
                                                                 check_type='safety')
        check_tot = len(bright_check_boss[bright_check_boss &
                                          hasFiber_boss])
        design_tot = len(bright_check_boss[hasFiber_boss])
        valid_arr['bright_safety_pass'][0] = check_tot
        valid_arr['bright_safety_total'][0] = design_tot
        bright_check_apogee, hasFiber_apogee = mode.bright_neighbors(instrument='APOGEE',
                                                                     check_type='safety')
        check_tot = len(bright_check_apogee[bright_check_apogee &
                                            hasFiber_apogee])
        design_tot = len(bright_check_apogee[hasFiber_apogee])
        valid_arr['bright_safety_pass'][0] += check_tot
        valid_arr['bright_safety_total'][0] += design_tot
    # do assigned targets
    design_total = len(des.design['catalogID'][des.design['catalogID'] != -1])
    valid_arr['all_targets_assigned_pass'][0] = design_total - len(des.targets_unassigned)
    valid_arr['all_targets_assigned_total'][0] = design_total
    # do collisions
    valid_arr['no_collisions_pass'][0] = design_total - len(des.targets_collided)
    valid_arr['no_collisions_total'][0] = design_total
    # check the design warnings
    column_names = ['boss_n_skies_min',
                    'apogee_n_skies_min', 'boss_min_skies_fovmetric',
                    'apogee_min_skies_fovmetric',
                    'boss_n_stds_min', 'apogee_n_stds_min',
                    'boss_min_stds_fovmetric', 'apogee_min_stds_fovmetric', 
                    'boss_stds_mags', 'apogee_stds_mags',
                    'boss_bright_limit_targets', 'apogee_bright_limit_targets',
                    'boss_sky_neighbors_targets', 'apogee_sky_neighbors_targets']
    warnings_order = ['min_skies_boss',
                      'min_skies_apogee', 'fov_skies_boss', 'fov_skies_apogee',
                      'min_stds_boss', 'min_stds_apogee', 'fov_stds_boss',
                      'fov_stds_apogee', 'stds_mag_boss', 'stds_mag_apogee',
                      'sci_mag_boss', 'sci_mag_apogee', 'bright_neigh_boss',
                      'bright_neigh_apogee']
    for c, k in zip(column_names, warnings_order):
        if isinstance(des.design_errors[k + '_metric'], list):
            valid_arr[c + '_pass'][0] = des.design_errors[k + '_metric'][0]
            valid_arr[c + '_total'][0] = des.design_errors[k + '_metric'][1]
        else:
            valid_arr[c + '_value'][0] = des.design_errors[k + '_metric']
    return valid_arr


def valid_design_func(file, exp, obsTime, field_desmodes,
                      db_results_boss, db_results_apogee,
                      desmodes):
    dm = field_desmodes[exp]
    des, decolide, bright_safety = validate_design(file,
                                                   exp + 1,
                                                   obsTime,
                                                   db_results_boss[dm],
                                                   db_results_apogee[dm],
                                                   desmodes[dm])
    valid_arr_des = design_outputs_to_array(des, decolide,
                                            bright_safety,
                                            db_results_boss[dm],
                                            db_results_apogee[dm],
                                            desmodes[dm])
    return valid_arr_des


sdss_path = sdss_access.path.Path(release='sdss5',
                                  preserve_envvars=True)

if __name__ == '__main__':
    parser = argparse.ArgumentParser(
        prog=os.path.basename(sys.argv[0]),
        description='In a batch, validate a set of designs')
    parser.add_argument('-t', '--type', dest='type',
                        type=str, help='Validating files in directory (dir) or robostrategy (rs)', 
                        choices=['dir', 'rs'], required=True)
    parser.add_argument('-l', '--loc', dest='loc',
                        type=str, help='local or utah',
                        choices=['local', 'utah'], required=True)
    parser.add_argument('-d', '--dir', dest='dir',
                        type=str, help='directory with design files (for type=dir)', required=False)
    parser.add_argument('-p', '--plan', dest='plan',
                        type=str, help='name of plan (for type=rs)', required=False)
    parser.add_argument('-o', '--observatory', dest='observatory',
                        type=str, help='apo or lco (for type=rs)',
                        choices=['apo', 'lco'], required=False)
    parser.add_argument('-n', '--Ncores', dest='Ncores',
                        type=int, help='number of cores to use. If Ncores=1, then not run in parallal.',
                        default=1, nargs='?')

    args = parser.parse_args()
    vtype = args.type
    loc = args.loc
    directory = args.dir
    plan = args.plan
    observatory = args.observatory
    Ncores = args.Ncores

    if loc == 'local':
        targetdb.database.connect_from_parameters(user='sdss_user',
                                                  host='localhost',
                                                  port=7502)
        catalogdb.database.connect_from_parameters(user='sdss_user',
                                                  host='localhost',
                                                  port=7502)
    else:
        targetdb.database.connect_from_parameters(user='sdss_user',
                                                  host='operations.sdss.utah.edu',
                                                  port=5432)
        catalogdb.database.connect_from_parameters(user='sdss_user',
                                                  host='operations.sdss.utah.edu',
                                                  port=5432)
    if vtype == 'dir':
        files = [file for file in glob.glob(directory + '*.fits')]
    else:
        files = []

        allocate_file = sdss_path.full('rsAllocationFinal', plan=plan,
                                       observatory=observatory)
        if 'sas' in allocate_file:
            allocate_file = allocate_file[:32] + 'sdss50' + allocate_file[44:]

        # only grab unique fields for this
        rsAllocation1 = fits.open(allocate_file)[1].data

        fieldids = np.unique(rsAllocation1["fieldid"])

        for fieldid in fieldids:
            # now grab the assignment file for this field
            field_assigned_file = sdss_path.full('rsFieldAssignmentsFinal',
                                                 plan=plan,
                                                 observatory=observatory,
                                                 fieldid=fieldid)

            if 'sas' in field_assigned_file:
                field_assigned_file = field_assigned_file[:32] + 'sdss50' + field_assigned_file[44:]

            files.append(field_assigned_file)

    if vtype == 'dir':
        file_save = directory + 'design_validation_results.fits'
    else:
        file_save = 'rs_%s_%s_design_validation_results.fits' % (plan, observatory)
    start = time.time()
    # grab all designmodes
    desmodes = allDesignModes()
    # start validaitng designs
    for file in files:
        # get header info
        head = fits.open(file)[0].header
        racen = head['RACEN']
        deccen = head['DECCEN']
        ot = obstime.ObsTime(observatory=head['obs'].strip())
        obsTime = coordio.time.Time(ot.nominal(lst=racen)).jd
        n_exp = head['NEXP']
        field_desmodes = head['DESMODE'].split(' ')
        # do db query results for each desmode in field
        db_results_boss = {}
        db_results_apogee = {}
        for dm in np.unique(field_desmodes):
            db_results_boss[dm] = {}
            db_results_apogee[dm] = {}
            if 'bright' in dm:
                # no r_sdss for bright so do g band
                # this is hacky and needs to be fixed!!!
                mag_lim = desmodes[dm].bright_limit_targets['BOSS'][0][0]
            else:
                mag_lim = desmodes[dm].bright_limit_targets['BOSS'][1][0]
            db_results_boss[dm]['designmode'] = build_brigh_neigh_query('designmode',
                                                                        'BOSS',
                                                                        mag_lim,
                                                                        racen,
                                                                        deccen)
            db_results_boss[dm]['safety'] = build_brigh_neigh_query('safety',
                                                                    'BOSS',
                                                                    mag_lim,
                                                                    racen,
                                                                    deccen)
            mag_lim = desmodes[dm].bright_limit_targets['APOGEE'][-1][0]
            db_results_apogee[dm]['designmode'] = build_brigh_neigh_query('designmode',
                                                                          'APOGEE',
                                                                          mag_lim,
                                                                          racen,
                                                                          deccen)
            db_results_apogee[dm]['safety'] = build_brigh_neigh_query('safety',
                                                                      'APOGEE',
                                                                      mag_lim,
                                                                      racen,
                                                                      deccen)
        if n_exp == 1:
            exp = 0
            dm = field_desmodes[exp]
            des, decolide, bright_safety = validate_design(file,
                                                           exp,
                                                           obsTime,
                                                           db_results_boss[dm],
                                                           db_results_apogee[dm],
                                                           desmodes[dm])
            valid_arr_des = design_outputs_to_array(des, decolide,
                                                    bright_safety,
                                                    db_results_boss[dm],
                                                    db_results_apogee[dm],
                                                    desmodes[dm])
            if 'valid_arr' in locals():
                valid_arr = np.append(valid_arr,
                                      valid_arr_des)
            else:
                valid_arr = valid_arr_des
        else:
            if Ncores > 1:
                if n_exp > Ncores:
                    nPool = Ncores
                else:
                    nPool = n_exp
                with Pool(processes=nPool) as pool:
                    res = pool.starmap(valid_design_func, zip(repeat(file),
                                                              range(n_exp),
                                                              repeat(obsTime),
                                                              repeat(field_desmodes),
                                                              repeat(db_results_boss),
                                                              repeat(db_results_apogee),
                                                              repeat(desmodes)))
                for valid_arr_des in res:
                    if 'valid_arr' in locals():
                        valid_arr = np.append(valid_arr,
                                              valid_arr_des)
                    else:
                        valid_arr = valid_arr_des
            else:
                for exp in range(n_exp):
                    dm = field_desmodes[exp]
                    des, decolide, bright_safety = validate_design(file,
                                                                   exp + 1,
                                                                   obsTime,
                                                                   db_results_boss[dm],
                                                                   db_results_apogee[dm],
                                                                   desmodes[dm])
                    valid_arr_des = design_outputs_to_array(des, decolide,
                                                            bright_safety,
                                                            db_results_boss[dm],
                                                            db_results_apogee[dm],
                                                            desmodes[dm])
                    if 'valid_arr' in locals():
                        valid_arr = np.append(valid_arr,
                                              valid_arr_des)
                    else:
                        valid_arr = valid_arr_des
    # write to fits file
    valid_arr = Table(valid_arr)
    valid_arr.write(file_save, format='fits')
    print('Took %.3f minutes to validate designs' % ((time.time() - start) / 60))
