import os
import sys
import argparse
from astropy.io import fits

import sdss_access.path

from mugatu.fpsdesign import FPSDesign
from mugatu.designmode import DesignModeCheck
import robostrategy.obstime as obstime
import coordio.time

from sdssdb.peewee.sdss5db import targetdb

targetdb.database.connect_from_parameters(user='sdss_user',
                                          host='localhost',
                                          port=7502)


def designmode_manual(file, exp, obsTime, observatory,
                      designmode_label, desmode_manual):
    """
    test designmode check of manual design
    """
    des = FPSDesign(design_pk=-1,
                    obsTime=obsTime,
                    observatory=observatory.upper(),
                    design_file=file,
                    manual_design=True,
                    exp=exp,
                    idtype='carton_to_target')

    des.build_design_manual()
    des.validate_design()

    mode = DesignModeCheck(FPSDesign=des,
                           desmode_label=designmode_label,
                           desmode_manual=desmode_manual)

    mode.design_mode_check_all()


if __name__ == "__main__":
    sdss_path = sdss_access.path.Path(release='sdss5')

    parser = argparse.ArgumentParser(
        prog=os.path.basename(sys.argv[0]),
        description='Read Robostratgey design file to mugatu')

    parser.add_argument('-p', '--plan', dest='plan',
                        type=str, help='name of plan', required=True)
    parser.add_argument('-o', '--observatory', dest='observatory',
                        type=str, help='apo or lco',
                        choices=['apo', 'lco'], required=True)
    parser.add_argument('-f', '--fieldid', dest='fieldid',
                        type=int, help='fieldid',
                        required=True)
    parser.add_argument('-e', '--exp', dest='exp',
                        type=int, help='exposure of design',
                        required=True)
    parser.add_argument('-l', '--loc', dest='loc',
                        type=str, help='local or utah',
                        choices=['local', 'utah'], required=True)

    args = parser.parse_args()
    plan = args.plan
    observatory = args.observatory
    fieldid = args.fieldid
    exp = args.exp
    loc = args.loc

    if loc == 'local':
        file = 'rsFieldAssignments-%s-%s-%d.fits' % (plan, observatory, fieldid)
    else:
        file = sdss_path.full('rsFieldAssignments',
                              plan=plan,
                              observatory=observatory,
                              fieldid=fieldid)
        if 'sas' in file:
            file = file[:32] + 'sdss50' + file[44:]

    head = fits.open(file)[0].header
    racen = head['RACEN']
    designmode_label = head['DESMODE'].split(' ')[exp - 1]
    ot = obstime.ObsTime(observatory=observatory.lower())
    obsTime = coordio.time.Time(ot.nominal(lst=racen)).jd

    desmode_manual = {}
    data = fits.open(file)[3].data
    desmode_manual['boss_skies_min'] = data[data['label'] == designmode_label]['boss_n_skies_min'][0]
    desmode_manual['apogee_skies_min'] = data[data['label'] == designmode_label]['apogee_n_skies_min'][0]

    desmode_manual['boss_skies_fov'] = data[data['label'] == designmode_label]['boss_min_skies_fovmetric'][0]
    desmode_manual['apogee_skies_fov'] = data[data['label'] == designmode_label]['apogee_min_skies_fovmetric'][0]

    desmode_manual['boss_stds_min'] = data[data['label'] == designmode_label]['boss_n_stds_min'][0]
    desmode_manual['apogee_stds_min'] = data[data['label'] == designmode_label]['apogee_n_stds_min'][0]

    desmode_manual['boss_stds_fov'] = data[data['label'] == designmode_label]['boss_min_stds_fovmetric'][0]
    desmode_manual['apogee_stds_fov'] = data[data['label'] == designmode_label]['apogee_min_stds_fovmetric'][0]

    desmode_manual['boss_stds_mags'] = data[data['label'] == designmode_label]['boss_stds_mags'][0]
    desmode_manual['apogee_stds_mags'] = data[data['label'] == designmode_label]['apogee_stds_mags'][0]

    desmode_manual['boss_bright_limit_targets'] = data[data['label'] == designmode_label]['boss_bright_limit_targets'][0]
    desmode_manual['apogee_bright_limit_targets'] = data[data['label'] == designmode_label]['apogee_bright_limit_targets'][0]

    desmode_manual['boss_sky_neighbors_targets'] = data[data['label'] == designmode_label]['boss_sky_neighbors_targets'][0]
    desmode_manual['apogee_sky_neighbors_targets'] = data[data['label'] == designmode_label]['apogee_sky_neighbors_targets'][0]

    desmode_manual['apogee_trace_diff_targets'] = data[data['label'] == designmode_label]['apogee_trace_diff_targets'][0]

    designmode_manual(file=file, exp=exp, obsTime=obsTime,
                     observatory=observatory,
                     designmode_label=designmode_label,
                     desmode_manual=desmode_manual)
