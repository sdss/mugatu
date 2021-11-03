import sys
import argparse
import os
import numpy as np
import glob

from astropy.io import fits
from astropy.table import Table

from sdssdb.peewee.sdss5db import targetdb
from sdssdb.peewee.sdss5db import catalogdb
import sdss_access.path
import robostrategy.obstime as obstime
import coordio.time

from mugatu.fpsdesign import FPSDesign
from mugatu.exceptions import MugatuDesignError


def validate_design(design_file, exp, obsTime):
    """
    Validate a design and record any errors or warnings
    from the validation
    """
    des = FPSDesign(design_pk=-1,
                    obsTime=obsTime,
                    design_file=design_file,
                    manual_design=True,
                    exp=exp,
                    idtype='catalogID')
    # set default
    decolide = True
    bright_safety = True
    # build design
    des.build_design_manual()
    # try to validate design and catch any design errors
    try:
        des.validate_design()
    except MugatuDesignError as e:
        if 'Kaiju' in str(e):
            decolide = False
        if 'Bright' in str(e):
            bright_safety = False
    return des, decolide, bright_safety


def design_outputs_to_array(des, decolide,
                            bright_safety):
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
                      ('all_targets_assigned', bool),
                      ('no_collisions', bool),
                      ('boss_n_skies_min', bool),
                      ('apogee_n_skies_min', bool),
                      ('boss_min_skies_fovmetric', bool),
                      ('apogee_min_skies_fovmetric', bool),
                      ('boss_n_stds_min', bool),
                      ('apogee_n_stds_min', bool),
                      ('boss_min_stds_fovmetric', bool),
                      ('apogee_min_stds_fovmetric', bool),
                      ('boss_stds_mags', bool),
                      ('apogee_stds_mags', bool),
                      ('boss_bright_limit_targets', bool),
                      ('apogee_bright_limit_targets', bool),
                      ('boss_sky_neighbors_targets', bool),
                      ('apogee_sky_neighbors_targets', bool)])
    valid_arr = np.zeros(1, dtype=dtype)
    valid_arr['file_name'][0] = os.path.split(des.design_file)[-1]
    if des.exp == 0:
        valid_arr['exp'][0] = des.exp + 1
    else:
        valid_arr['exp'][0] = des.exp
    valid_arr['racen'][0] = des.racen
    valid_arr['deccen'][0] = des.deccen
    valid_arr['designmode'][0] = des.desmode_label
    valid_arr['decolide'][0] = des.desmode_label
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
    return valid_arr


sdss_path = sdss_access.path.Path(release='sdss5')

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

    args = parser.parse_args()
    vtype = args.type
    loc = args.loc
    directory = args.dir
    plan = args.plan
    observatory = args.observatory

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

        allocate_file = sdss_path.full('rsAllocation', plan=plan,
                                       observatory=observatory)
        if 'sas' in allocate_file:
            allocate_file = allocate_file[:32] + 'sdss50' + allocate_file[44:]

        rsAllocation1 = fits.open(allocate_file)[1].data

        fieldids = rsAllocation1["fieldid"]

        for fieldid in fieldids:
            # now grab the assignment file for this field
            field_assigned_file = sdss_path.full('rsFieldAssignments',
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
    # start validaitng designs
    for file in files:
        # get header info
        head = fits.open(file)[0].header
        racen = head['RACEN']
        ot = obstime.ObsTime(observatory=head['obs'].strip())
        obsTime = coordio.time.Time(ot.nominal(lst=racen)).jd
        n_exp = head['NEXP']
        for exp in range(n_exp):
            if n_exp == 1:
                des, decolide, bright_safety = validate_design(file,
                                                               exp,
                                                               obsTime)
            else:
                des, decolide, bright_safety = validate_design(file,
                                                               exp + 1,
                                                               obsTime)
            valid_arr_des = design_outputs_to_array(des, decolide,
                                                    bright_safety)
            if 'valid_arr' in locals():
                valid_arr = np.append(valid_arr,
                                      valid_arr_des)
            else:
                valid_arr = valid_arr_des
    # write to fits file
    valid_arr = Table(valid_arr)
    valid_arr.write(file_save, format='fits')
