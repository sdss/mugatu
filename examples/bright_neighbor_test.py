import os
import sys
import argparse
from astropy.io import fits

import sdss_access.path

from mugatu.fpsdesign import FPSDesign
from mugatu.designmode import DesignModeCheck, build_brigh_neigh_query
import robostrategy.obstime as obstime
import coordio.time
from sdssdb.peewee.sdss5db import targetdb
from sdssdb.peewee.sdss5db import catalogdb
import time


targetdb.database.connect_from_parameters(user='sdss_user',
                                          host='localhost',
                                          port=7502)

catalogdb.database.connect_from_parameters(user='sdss_user',
                                          host='localhost',
                                          port=7502)

def bright_neighbor_test(file, exp, obsTime,
                         designmode_label):
    des = FPSDesign(design_pk=-1,
                    obsTime=obsTime,
                    design_file=file,
                    manual_design=True,
                    exp=exp,
                    idtype='carton_to_target')

    start = time.time()
    des.build_design_manual()
    print('Took %.6f secs to build design' % (time.time() - start))
    start = time.time()
    des.validate_design(safety=False)
    print('Took %.6f secs to validate design' % (time.time() - start))

    mode = DesignModeCheck(FPSDesign=des,
                           desmode_label=designmode_label)

    if 'bright' in mode.desmode_label:
        # no r_sdss for bright so do g band
        # this is hacky and needs to be fixed!!!
        mag_lim = mode.bright_limit_targets['BOSS'][0][0]
    else:
        mag_lim = mode.bright_limit_targets['BOSS'][1][0]

    db_query_results_boss = build_brigh_neigh_query('designmode', 'BOSS', mag_lim,
                                                    des.racen, des.deccen)

    if 'bright' in mode.desmode_label:
        # no r_sdss for bright so do g band
        # this is hacky and needs to be fixed!!!
        mag_lim = mode.bright_limit_targets['APOGEE'][-1][0]
    else:
        mag_lim = mode.bright_limit_targets['APOGEE'][-1][0]

    db_query_results_apogee = build_brigh_neigh_query('designmode', 'APOGEE', mag_lim,
                                                    des.racen, des.deccen)

    des = FPSDesign(design_pk=-1,
                    obsTime=obsTime,
                    design_file=file,
                    manual_design=True,
                    exp=exp,
                    idtype='carton_to_target')

    start = time.time()
    des.build_design_manual()
    print('Took %.6f secs to build design' % (time.time() - start))
    start = time.time()
    des.validate_design(safety=False,
                        db_query_results_boss=db_query_results_boss,
                        db_query_results_apogee=db_query_results_apogee)
    print('Took %.6f secs to validate design' % (time.time() - start))

    mode = DesignModeCheck(FPSDesign=des,
                           desmode_label=designmode_label)

    print('Test without supplying queries')
    start = time.time()
    bright_check, hasFiber = mode.bright_neighbors(instrument='BOSS',
                                                   check_type='designmode')
    print('Took %.6f secs to do checkk' % (time.time() - start))
    print('%d fibers fail desmode bright check for BOSS' % len(bright_check[~bright_check & hasFiber]))
    print('')
    start = time.time()
    bright_check, hasFiber = mode.bright_neighbors(instrument='APOGEE',
                                                   check_type='designmode')
    print('Took %.6f secs to do checkk' % (time.time() - start))
    print('%d fibers fail desmode bright check for APOGEE' % len(bright_check[~bright_check & hasFiber]))
    print('')
    start = time.time()
    bright_check, hasFiber = mode.bright_neighbors(instrument='BOSS',
                                                   check_type='safety')
    print('Took %.6f secs to do checkk' % (time.time() - start))
    print('%d fibers fail safety bright check for BOSS' % len(bright_check[~bright_check & hasFiber]))
    print('')
    start = time.time()
    bright_check, hasFiber = mode.bright_neighbors(instrument='APOGEE',
                                                   check_type='safety')
    print('Took %.6f secs to do checkk' % (time.time() - start))
    print('%d fibers fail safety bright check for APOGEE' % len(bright_check[~bright_check & hasFiber]))
    print('')

    print('')
    print('Test with supplying queries')
    if 'bright' in mode.desmode_label:
        # no r_sdss for bright so do g band
        # this is hacky and needs to be fixed!!!
        mag_lim = mode.bright_limit_targets['BOSS'][0][0]
    else:
        mag_lim = mode.bright_limit_targets['BOSS'][1][0]

    db_query_results_boss = build_brigh_neigh_query('designmode', 'BOSS', mag_lim,
                                                    des.racen, des.deccen)

    start = time.time()
    mode = DesignModeCheck(FPSDesign=des,
                           desmode_label=designmode_label,
                           db_query_results_boss=db_query_results_boss)
    bright_check, hasFiber = mode.bright_neighbors(instrument='BOSS',
                                                   check_type='designmode')
    print('Took %.6f secs to do checkk' % (time.time() - start))
    print('%d fibers fail desmode bright check for BOSS' % len(bright_check[~bright_check & hasFiber]))
    print('')

    if 'bright' in mode.desmode_label:
        # no r_sdss for bright so do g band
        # this is hacky and needs to be fixed!!!
        mag_lim = mode.bright_limit_targets['APOGEE'][-1][0]
    else:
        mag_lim = mode.bright_limit_targets['APOGEE'][-1][0]

    db_query_results_apogee = build_brigh_neigh_query('designmode', 'APOGEE', mag_lim,
                                                    des.racen, des.deccen)

    start = time.time()
    mode = DesignModeCheck(FPSDesign=des,
                           desmode_label=designmode_label,
                           db_query_results_apogee=db_query_results_apogee)
    bright_check, hasFiber = mode.bright_neighbors(instrument='APOGEE',
                                                   check_type='designmode')
    print('Took %.6f secs to do checkk' % (time.time() - start))
    print('%d fibers fail desmode bright check for APOGEE' % len(bright_check[~bright_check & hasFiber]))
    print('')

    carts = ['ops_tycho2_brightneighbors',
             'ops_gaia_brightneighbors']
    mag_col = targetdb.Magnitude.gaia_g
    db_query = (targetdb.CartonToTarget.select(targetdb.Target.ra,
                                               targetdb.Target.dec,
                                               mag_col,
                                               targetdb.CartonToTarget.pk)
                                       .join(targetdb.Target)
                                       .switch(targetdb.CartonToTarget)
                                       .join(targetdb.Magnitude)
                                       .switch(targetdb.CartonToTarget)
                                       .join(targetdb.Carton)
                                       .where((targetdb.Target.cone_search(des.racen,
                                                                           des.deccen,
                                                                           1.5)) &
                                              (targetdb.Carton.carton.in_(carts))))
    db_query.execute()
    start = time.time()
    bright_check, hasFiber = mode.bright_neighbors(instrument='BOSS',
                                                   check_type='safety',
                                                   db_query=db_query)
    print('Took %.6f secs to do checkk' % (time.time() - start))
    print('%d fibers fail safety bright check for BOSS' % len(bright_check[~bright_check & hasFiber]))
    print('')

    carts = ['ops_2mass_psc_brightneighbors']
    mag_col = targetdb.Magnitude.h
    db_query = (targetdb.CartonToTarget.select(targetdb.Target.ra,
                                               targetdb.Target.dec,
                                               mag_col,
                                               targetdb.CartonToTarget.pk)
                                       .join(targetdb.Target)
                                       .switch(targetdb.CartonToTarget)
                                       .join(targetdb.Magnitude)
                                       .switch(targetdb.CartonToTarget)
                                       .join(targetdb.Carton)
                                       .where((targetdb.Target.cone_search(des.racen,
                                                                           des.deccen,
                                                                           1.5)) &
                                              (targetdb.Carton.carton.in_(carts))))
    db_query.execute()
    bright_check, hasFiber = mode.bright_neighbors(instrument='APOGEE',
                                                   check_type='safety',
                                                   db_query=db_query)
    print('Took %.6f secs to do checkk' % (time.time() - start))
    print('%d fibers fail safety bright check for APOGEE' % len(bright_check[~bright_check & hasFiber]))
    print('')


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

        # connect to targetdb
        targetdb.database.connect_from_parameters(user='sdss_user',
                                                  host='localhost',
                                                  port=7502)
        catalogdb.database.connect_from_parameters(user='sdss_user',
                                                  host='localhost',
                                                  port=7502)
    else:
        file = sdss_path.full('rsFieldAssignments',
                              plan=plan,
                              observatory=observatory,
                              fieldid=fieldid)
        if 'sas' in file:
            file = file[:32] + 'sdss50' + file[44:]
        targetdb.database.connect_from_parameters(user='sdss_user',
                                                  host='operations.sdss.utah.edu',
                                                  port=5432)
        catalogdb.database.connect_from_parameters(user='sdss_user',
                                                  host='operations.sdss.utah.edu',
                                                  port=5432)

    head = fits.open(file)[0].header
    racen = head['RACEN']
    designmode_label = head['DESMODE'].split(' ')[exp - 1]
    ot = obstime.ObsTime(observatory=observatory.lower())
    obsTime = coordio.time.Time(ot.nominal(lst=racen)).jd

    bright_neighbor_test(file=file, exp=exp, obsTime=obsTime,
                         designmode_label=designmode_label)
