import sys
import argparse
import os
import numpy as np
import glob
import time

from astropy.io import fits
from astropy.table import Table

import sdss_access.path
import robostrategy.obstime as obstime
import coordio.time

from mugatu.exceptions import MugatuDesignError, MugatuError
from multiprocessing import Pool
from itertools import repeat
from functools import partial
from tqdm import tqdm

import mugatu
mugatu_ver = mugatu.__version__

import kaiju
kaiju_ver = kaiju.__version__

import coordio
coordio_ver = coordio.__version__

import fps_calibrations
fps_calib_ver = fps_calibrations.get_version()

primary_hdu = fits.PrimaryHDU()

primary_hdu.header['mugatu_version'] = mugatu_ver
primary_hdu.header['kaiju_version'] = kaiju_ver
primary_hdu.header['coordio_version'] = coordio_ver
primary_hdu.header['fps_calibrations_version'] = fps_calib_ver


def valid_field(all_files, offset_min_skybrightness, cache_bs):
    # need import here for create new connection
    from mugatu.fpsdesign import FPSDesign
    from mugatu.designmode import (build_brigh_neigh_query,
                                   DesignModeCheck,
                                   allDesignModes)

    def validate_design(design_file, exp, obsTime,
                        db_query_results_boss, db_query_results_apogee,
                        desmode_manual):
        """
        Validate a design and record any errors or warnings
        from the validation
        """
        if offset_min_skybrightness is None:
            des = FPSDesign(design_pk=-1,
                            obsTime=obsTime,
                            design_file=design_file,
                            manual_design=True,
                            exp=exp)
        else:
            des = FPSDesign(design_pk=-1,
                            obsTime=obsTime,
                            design_file=design_file,
                            manual_design=True,
                            exp=exp,
                            offset_min_skybrightness=offset_min_skybrightness)
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
        dtype = np.dtype([('file_name', '<U100'),
                          ('exp', np.int32),
                          ('racen', np.float64),
                          ('deccen', np.float64),
                          ('designmode', '<U15'),
                          ('designid_status', np.int32),
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
                          ('boss_sky_neighbors_targets_adj_mags', float,
                           (500,)),
                          ('apogee_sky_neighbors_targets', bool),
                          ('apogee_sky_neighbors_targets_pass', np.int32),
                          ('apogee_sky_neighbors_targets_total', np.int32),
                          ('apogee_sky_neighbors_targets_adj_mags', float,
                           (500,))])
        valid_arr = np.zeros(1, dtype=dtype)
        valid_arr['file_name'][0] = os.path.split(des.design_file)[-1]
        if des.exp == 0:
            valid_arr['exp'][0] = des.exp
        else:
            valid_arr['exp'][0] = des.exp - 1
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
                                   desmode_manual=desmode_manual.todict())
            bright_check_boss, hasFiber_boss, _, isassigned_boss = mode.bright_neighbors(
                instrument='BOSS', check_type='safety')
            check_tot = len(bright_check_boss[bright_check_boss &
                                              hasFiber_boss &
                                              isassigned_boss])
            design_tot = len(bright_check_boss[hasFiber_boss &
                                               isassigned_boss])
            valid_arr['bright_safety_pass'][0] = check_tot
            valid_arr['bright_safety_total'][0] = design_tot
            bright_check_apogee, hasFiber_apogee, _, isassigned_apogee = mode.bright_neighbors(
                instrument='APOGEE', check_type='safety')
            check_tot = len(bright_check_apogee[bright_check_apogee &
                                                hasFiber_apogee &
                                                isassigned_apogee])
            design_tot = len(bright_check_apogee[hasFiber_apogee &
                                                 isassigned_apogee])
            valid_arr['bright_safety_pass'][0] += check_tot
            valid_arr['bright_safety_total'][0] += design_tot
        else:
            valid_arr['bright_safety_pass'][0] = 500
            valid_arr['bright_safety_total'][0] = 500
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
                if len(des.design_errors[k + '_metric']) > 2:
                    valid_arr[c + '_adj_mags'][0] = list(des.design_errors[k + '_metric'][2])
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

    desmodes = allDesignModes()

    file = all_files[0]
    cache_file = all_files[1]
    if cache_bs:
        hdu_cache = fits.open(cache_file)

    head = fits.open(file)[0].header
    racen = head['RACEN']
    deccen = head['DECCEN']
    ot = obstime.ObsTime(observatory=head['obs'].strip())
    obsTime = coordio.time.Time(ot.nominal(lst=racen)).jd
    n_exp = head['NEXP']
    field_desmodes = head['DESMODE'].split(' ')
    try:
        design_ids = fits.open(file)['STATUS'].data['designid']
    except KeyError:
        design_ids = np.zeros(n_exp, dtype=np.int32) - 1
    # do db query results for each desmode in field
    db_results_boss = {}
    db_results_apogee = {}
    for dm in np.unique(field_desmodes):
        db_results_boss[dm] = {}
        db_results_apogee[dm] = {}
        if 'bright' in dm:
            mag_lim = desmodes[dm].bright_limit_targets['BOSS'][5][0]
        else:
            mag_lim = desmodes[dm].bright_limit_targets['BOSS'][1][0]
        if cache_bs:
            for idc in range(1, len(hdu_cache)):
                if hdu_cache[idc].header['FIBERTY'] == 'BOSS' and hdu_cache[idc].header['DESMODE'] == dm:
                    bs_field = hdu_cache[idc].data
                    db_results_boss[dm]['designmode'] = (bs_field['catalog_ra'],
                                                         bs_field['catalog_dec'],
                                                         bs_field['mag'],
                                                         bs_field['catalogid'],
                                                         bs_field['pmra'],
                                                         bs_field['pmdec'])
                    break
        else:
            db_results_boss[dm]['designmode'] = build_brigh_neigh_query('designmode',
                                                                        'BOSS',
                                                                        mag_lim,
                                                                        racen,
                                                                        deccen,
                                                                        head['obs'].strip().upper())
        db_results_boss[dm]['safety'] = build_brigh_neigh_query('safety',
                                                                'BOSS',
                                                                mag_lim,
                                                                racen,
                                                                deccen,
                                                                head['obs'].strip().upper())
        mag_lim = desmodes[dm].bright_limit_targets['APOGEE'][8][0]
        if cache_bs:
            for idc in range(1, len(hdu_cache)):
                if hdu_cache[idc].header['FIBERTY'] == 'APOGEE' and hdu_cache[idc].header['DESMODE'] == dm:
                    bs_field = hdu_cache[idc].data
                    db_results_apogee[dm]['designmode'] = (bs_field['catalog_ra'],
                                                           bs_field['catalog_dec'],
                                                           bs_field['mag'],
                                                           bs_field['catalogid'],
                                                           bs_field['pmra'],
                                                           bs_field['pmdec'])
                    break
        else:
            db_results_apogee[dm]['designmode'] = build_brigh_neigh_query('designmode',
                                                                          'APOGEE',
                                                                          mag_lim,
                                                                          racen,
                                                                          deccen,
                                                                          head['obs'].strip().upper())
        db_results_apogee[dm]['safety'] = build_brigh_neigh_query('safety',
                                                                  'APOGEE',
                                                                  mag_lim,
                                                                  racen,
                                                                  deccen,
                                                                  head['obs'].strip().upper())
    if n_exp == 1:
        exp = 0
        dm = field_desmodes[exp]
        des, decolide, bright_safety = validate_design(file,
                                                       exp,
                                                       obsTime,
                                                       db_results_boss[dm],
                                                       db_results_apogee[dm],
                                                       desmodes[dm])
        valid_arr = design_outputs_to_array(des, decolide,
                                            bright_safety,
                                            db_results_boss[dm],
                                            db_results_apogee[dm],
                                            desmodes[dm])

        valid_arr['designid_status'][0] = design_ids[0]
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
            valid_arr_des['designid_status'][0] = design_ids[exp]
            if exp == 0:
                valid_arr = valid_arr_des
            else:
                valid_arr = np.append(valid_arr,
                                      valid_arr_des)
    return valid_arr



sdss_path = sdss_access.path.Path(release='sdss5',
                                  preserve_envvars=True)

if __name__ == '__main__':
    parser = argparse.ArgumentParser(
        prog=os.path.basename(sys.argv[0]),
        description='In a batch, validate a set of designs')
    parser.add_argument('-t', '--type', dest='type',
                        type=str, help='Validating files in directory (dir) or robostrategy (rs)',
                        choices=['dir', 'rs', 'rs_replace', 'rs_catchup'], required=True)
    parser.add_argument('-d', '--dir', dest='dir',
                        type=str, help='directory with design files (for type=dir)',
                        required=False)
    parser.add_argument('-p', '--plan', dest='plan',
                        type=str, help='name of plan (for type=rs or type=rs_replace)',
                        required=False)
    parser.add_argument('-o', '--observatory', dest='observatory',
                        type=str, help='apo or lco (for type=rs)',
                        choices=['apo', 'lco'], required=False)
    parser.add_argument('-f', '--fieldids', dest='fieldids', nargs='+',
                        help='field_ids to validate (for type=rs_replace)',
                        type=int, required=False)
    parser.add_argument('-n', '--Ncores', dest='Ncores',
                        type=int, help='number of cores to use. If Ncores=1, then not run in parallal.',
                        default=1, nargs='?')
    parser.add_argument('-s','--offset_min_skybrightness', help='offset_min_skybrightness for design',
                        type=float, required=False)
    parser.add_argument('-c','--cache_bs', help='if want to use cache of bright star queries',
                        type=bool, required=False, default=False)
    parser.add_argument('-v', '--ver_catch', dest='ver_catch',
                        type=str, help='version of catchup (for type=rs_catchup)', required=False)

    args = parser.parse_args()
    vtype = args.type
    directory = args.dir
    plan = args.plan
    observatory = args.observatory
    fieldids = args.fieldids
    Ncores = args.Ncores
    offset_min_skybrightness = args.offset_min_skybrightness
    cache_bs = args.cache_bs
    ver_catch = args.ver_catch

    MUGATU_DATA = os.getenv('MUGATU_DATA')

    if vtype == 'dir':
        files = [file for file in glob.glob(directory + '*.fits')]
        cache_files = []
    elif vtype == 'rs':
        files = []

        cache_files = []
        cache_dir = os.popen('echo $ROBOSTRATEGY_DATA').read()[:-1] + '/allocations/eta-0-bs-cache/targets/'

        allocate_file = sdss_path.full('rsAllocationFinal', plan=plan,
                                       observatory=observatory)
        if 'sas' in allocate_file:
            allocate_file = allocate_file[:32] + 'sdss50' + allocate_file[44:]

        # only grab unique fields for this
        rsAllocation1 = fits.open(allocate_file)[1].data

        fieldids = np.unique(rsAllocation1["fieldid"])

        rs_fieldids = np.unique(rsAllocation1["rs_fieldid"])

        for fieldid, rs_fieldid in zip(fieldids, rs_fieldids):
            # now grab the assignment file for this field
            field_assigned_file = sdss_path.full('rsFieldAssignmentsFinal',
                                                 plan=plan,
                                                 observatory=observatory,
                                                 fieldid=fieldid)

            if 'sas' in field_assigned_file:
                field_assigned_file = field_assigned_file[:32] + 'sdss50' + field_assigned_file[44:]

            files.append(field_assigned_file)

            if cache_bs and 'eta' in plan:
                cache_file = cache_dir + 'rsBrightStars-eta-0-bs-cache-{obs}-{fid}.fits'.format(obs=observatory,
                                                                                               fid=rs_fieldid)
                cache_files.append(cache_file)
    elif vtype == 'rs_catchup':
        files = []

        cache_files = []
        cache_dir = os.popen('echo $ROBOSTRATEGY_DATA').read()[:-1] + '/allocations/eta-0-bs-cache/targets/'

        allocate_file = sdss_path.full('rsAllocationFinal', plan=plan,
                                       observatory=observatory)
        if 'sas' in allocate_file:
            allocate_file = allocate_file[:32] + 'sdss50' + allocate_file[44:]

        # get the catchup file
        allocate_file = allocate_file.replace('final', 'catchup').replace('Final', 'Catchup%s' % ver_catch)

        # only grab unique fields for this
        rsAllocation1 = fits.open(allocate_file)[1].data

        fieldids = np.unique(rsAllocation1["fieldid"])

        rs_fieldids = np.unique(rsAllocation1["rs_fieldid"])

        for fieldid, rs_fieldid in zip(fieldids, rs_fieldids):
            # now grab the assignment file for this field
            field_assigned_file = sdss_path.full('rsFieldAssignmentsFinal',
                                                 plan=plan,
                                                 observatory=observatory,
                                                 fieldid=fieldid)

            if 'sas' in field_assigned_file:
                field_assigned_file = field_assigned_file[:32] + 'sdss50' + field_assigned_file[44:]

            # get the catchup file
            field_assigned_file = field_assigned_file.replace('final', 'catchup').replace('Final', 'Catchup%s' % ver_catch)

            files.append(field_assigned_file)

            if cache_bs and 'eta' in plan:
                cache_file = cache_dir + 'rsBrightStars-eta-0-bs-cache-{obs}-{fid}.fits'.format(obs=observatory,
                                                                                               fid=rs_fieldid)
                cache_files.append(cache_file)
    elif vtype == 'rs_replace':
        replace_path = ('/uufs/chpc.utah.edu/common/home/sdss50/sdsswork/'
                        'target/robostrategy_replacement/{plan}/'.format(
                            plan=plan))
        files = []
        cache_files = []
        for fid in fieldids:
            files += [file for file in glob.glob(replace_path +
                                                 '{plan}_{fid}*.fits'.format(
                                                     plan=plan,
                                                     fid=fid))]
        for f in files:
            if 'validation' in f:
                files.remove(f)
    else:
        raise MugatuError(message='Improper Validation Type')

    if vtype == 'dir':
        file_save = directory + 'design_validation_results.fits'
    elif vtype == 'rs':
        if not os.path.isdir((MUGATU_DATA + '/rs_plan_validations/{plan}'
                              .format(plan=plan))):
            os.makedirs((MUGATU_DATA + '/rs_plan_validations/{plan}'
                         .format(plan=plan)))
        file_save = (MUGATU_DATA + '/rs_plan_validations/{plan}/'
                     'rs_{plan}_{obs}_design_validation_results.fits'.format(
                         plan=plan,
                         obs=observatory))
    elif vtype == 'rs_catchup':
        if not os.path.isdir((MUGATU_DATA + '/rs_plan_validations/{plan}'
                              .format(plan=plan))):
            os.makedirs((MUGATU_DATA + '/rs_plan_validations/{plan}'
                         .format(plan=plan)))
        file_save = (MUGATU_DATA + '/rs_plan_validations/{plan}/'
                     'rs_Catchup{ver_catch}_{plan}_{obs}_design_validation_results.fits'.format(
                         ver_catch=ver_catch,
                         plan=plan,
                         obs=observatory))
    elif vtype == 'rs_replace':
        replace_path = ('/uufs/chpc.utah.edu/common/home/sdss50/sdsswork/'
                        'target/robostrategy_replacement/{plan}/'.format(
                            plan=plan))
        file_save = []
        for f in files:
            file_save.append(f[:-5] + '_validation.fits')
    else:
        raise MugatuError(message='Improper Validation Type')

    start = time.time()
    # start validaitng designs
    with Pool(processes=Ncores) as pool:
        if len(cache_files) > 0:
            all_files = [(f, cf) for f, cf in zip(files, cache_files)]
        else:
            all_files = [(f, '') for f in files]
        res = tqdm(pool.imap(partial(valid_field, offset_min_skybrightness=offset_min_skybrightness,
                                     cache_bs=cache_bs),
                                     all_files),
                   total=len(files))
        res = [r for r in res]
    for i, r in enumerate(res):
        if vtype == 'rs_replace':
            valid_arr = Table(r)
            bin_fits = fits.BinTableHDU(valid_arr)
            hdu = fits.HDUList([primary_hdu, bin_fits])
            hdu.writeto(file_save[i])
        elif i == 0:
            valid_arr = r
        else:
            valid_arr = np.append(valid_arr,
                                  r)
    # write to fits file
    if vtype == 'dir' or vtype == 'rs' or vtype == 'rs_catchup':
        valid_arr = Table(valid_arr)
        bin_fits = fits.BinTableHDU(valid_arr)
        hdu = fits.HDUList([primary_hdu, bin_fits])
        hdu.writeto(file_save)
    print('Took %.3f minutes to validate designs' % ((time.time() - start) / 60))
