import sys
import argparse
import os
import numpy as np
import pandas as pd

from astropy.io import fits

from sdssdb.peewee.sdss5db import targetdb
import sdss_access.path
import coordio
import datetime
from astropy.time import Time
from mugatu.fpsdesign import FPSDesign
from mugatu.designmode import DesignModeCheck
import time

sdss_path = sdss_access.path.Path(release='sdss5')


class ObsTime(object):
    """Class for finding appropriate observing times
    Parameters
    ----------
    observatory : str
        'apo' or 'lco'
    year : int
        nominal year to consider (default 2021)
    Attributes
    ----------
    observatory : str
        'apo' or 'lco'
    year : int
        nominal year to consider
    utcoff : int
        offset of local time from UTC
    transit_lst : ndarray of np.float64
        [365] LST (deg) transiting at each local standard midnight of year
    midnights : list of datetime.datetime objects
        [365] datetime format for each local standard midnight of year
    Methods
    -------
    nominal(lst=) : returns nominal observing time for a given RA
    Notes
    -----
    This class provides a way to assign a nominal observation time for
    a given LST.
    nominal() returns the local midnight at which the the LST is
    closest to transiting. It differs slightly from this at the 0/360
    deg boundary of LSTs.
    It uses SDSS's coordio for the astronomy calculation.
    This is taken from Robostrategy (not in current branch, so moved here)
    """
    def __init__(self, observatory='apo', year=2021):
        self.observatory = observatory
        self.year = year
        if(observatory == 'apo'):
            self.utcoff = - 7
        if(observatory == 'lco'):
            self.utcoff = - 4

        oneday = datetime.timedelta(days=1)
        onehour = datetime.timedelta(hours=1)

        site = coordio.site.Site(self.observatory.upper())

        self.transit_lst = np.zeros(365, dtype=np.float64)
        self.midnight = []

        day = datetime.datetime(year, 1, 1) - self.utcoff * onehour
        for n in range(365):
            midnight = day + oneday * n
            site.set_time(midnight)
            south = coordio.sky.Observed([[45., 180.]], site=site)
            self.transit_lst[n] = south.ra
            self.midnight.append(midnight)

        return

    def nominal(self, lst=None):
        """Return a nominal observation time for a given LST
        Parameters
        ----------
        lst : np.float64 or float
            LST desired for the observation (deg)
        Returns
        -------
        nominal_time : datetime object
            datetime object describing the midnight at which this LST
            is closest to transiting.
        Notes
        ------
        At 0/360 boundary picks the closest night to that boundary.
        This should be a very minor effect (few minutes).
"""
        imin = np.abs(self.transit_lst - lst).argmin()
        return(self.midnight[imin])


def arrayify(s):
    l = s.replace('[', '').replace(']', '').split(",")
    for i in range(len(l)):
        try:
            l[i] = float(l[i])
        except ValueError:
            l[i] = None
    return l


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
    field_file = sdss_path.full('rsFields', plan=plan,
                                observatory=observatory)

    # connect to targetdb
    targetdb.database.connect_from_parameters(user='sdss',
                                              host='operations.sdss.utah.edu',
                                              port=5432)

    rsField = fits.open(field_file)

    fieldids = rsField[1].data["fieldid"]

    designmodes = pd.read_csv('DesignMode_values.txt', sep=';')

    # setup dict
    desmode_manuals = {}
    for ind in range(len(designmodes)):
        desmode_manual = {}
        desmode_manual['boss_skies_min'] = designmodes.loc[ind,'BOSS_skies_min']
        desmode_manual['apogee_skies_min'] = designmodes.loc[ind,'APOGEE_skies_min']
        desmode_manual['boss_skies_fov'] = arrayify(designmodes.loc[ind,'BOSS_skies_FOV'])
        desmode_manual['apogee_skies_fov'] = arrayify(designmodes.loc[ind,'APOGEE_skies_FOV'])
        desmode_manual['boss_stds_min'] = designmodes.loc[ind,'BOSS_stds_min']
        desmode_manual['apogee_stds_min'] = designmodes.loc[ind,'APOGEE_stds_min']
        desmode_manual['boss_stds_fov'] = arrayify(designmodes.loc[ind,'BOSS_stds_FOV'])
        desmode_manual['apogee_stds_fov'] = arrayify(designmodes.loc[ind,'APOGEE_stds_FOV'])
        desmode_manual['boss_stds_mags'] = np.column_stack((arrayify(designmodes.loc[ind,'BOSS_stds_mags_min']),
                                                            arrayify(designmodes.loc[ind,'BOSS_stds_mags_max'])))
        desmode_manual['apogee_stds_mags'] = np.column_stack((arrayify(designmodes.loc[ind,'APOGEE_stds_mags_min']),
                                                              arrayify(designmodes.loc[ind,'APOGEE_stds_mags_max'])))
        desmode_manual['boss_bright_limit_targets'] = np.column_stack((arrayify(designmodes.loc[ind,'BOSS_bright_limit_targets_min']),
                                                                       arrayify(designmodes.loc[ind,'BOSS_bright_limit_targets_max'])))
        desmode_manual['apogee_bright_limit_targets'] = np.column_stack((arrayify(designmodes.loc[ind,'APOGEE_bright_limit_targets_min']),
                                                                         arrayify(designmodes.loc[ind,'APOGEE_bright_limit_targets_max'])))
        desmode_manual['boss_sky_neighbors_targets'] = np.array([5, 5, 5])
        desmode_manual['apogee_sky_neighbors_targets'] = np.array([5, 5, 5])
        desmode_manual['apogee_trace_diff_targets'] = None
        desmode_manuals[designmodes.loc[ind,'pk']] = desmode_manual

    start = time.time()
    bad = 0

    with open('rs_%s_%s_designmode_results_sdssdb_0_4_9.txt' % (plan, observatory), 'w') as f:
        header = ['fieldid',
                  'exp',
                  'cadence',
                  'designmode',
                  'BOSS_skies_min',
                  'BOSS_skies_FOV',
                  'APOGEE_skies_min',
                  'APOGEE_skies_FOV',
                  'BOSS_stds_min',
                  'BOSS_stds_mags',
                  'BOSS_stds_FOV',
                  'APOGEE_stds_min',
                  'APOGEE_stds_mags',
                  'APOGEE_stds_FOV',
                  'BOSS_bright_limit_targets',
                  'APOGEE_bright_limit_targets']
        f.write('\t'.join([str(x) for x in header]) + '\n')

        for fieldid in fieldids:
            # now grab the assignment file for this field
            field_assigned_file = sdss_path.full('rsFieldAssignments',
                                                 plan=plan,
                                                 observatory=observatory,
                                                 fieldid=fieldid)
            field_assign_1 = fits.open(field_assigned_file)[1].data
            field_assign_2 = fits.open(field_assigned_file)[2].data
            head = fits.open(field_assigned_file)[0].header

            carton_pks = np.zeros(len(field_assign_1['catalogid']), dtype=int)

            for cart in np.unique(field_assign_1['carton']):
                c = (targetdb.Carton.select()
                                    .join(targetdb.Version)
                                    .where((targetdb.Version.pk == 83) &
                                           (targetdb.Carton.carton == cart)))
                try:
                    carton_pks[field_assign_1['carton'] == cart] = c[0].pk
                except:
                    c = (targetdb.CartonToTarget.get(pk=field_assign_1['carton_to_target_pk'][field_assign_1['carton'] == cart][0])
                                    )
                    carton_pks[field_assign_1['carton'] == cart] = c.carton.pk

            # grab all magntiudes now
            mags = np.zeros((len(field_assign_1['catalogid']), 7)) - 9999.99
            for i in range(len(field_assign_1['catalogid'])):
                # do not query skies, they have no mag
                if 'sky' not in field_assign_1['category'][i] and field_assign_2['assigned'][i] == 1:
                    mag_query = (targetdb.Magnitude.select()
                                                   .where((targetdb.Magnitude.carton_to_target == field_assign_1['carton_to_target_pk'][i])))
                    if mag_query[0].optical_prov is None:
                        mags[i][0] = mag_query[0].g
                        mags[i][1] = mag_query[0].r
                        mags[i][2] = mag_query[0].i
                        mags[i][3] = mag_query[0].bp
                        mags[i][4] = mag_query[0].gaia_g
                        mags[i][5] = mag_query[0].rp
                        mags[i][6] = mag_query[0].h
                    elif 'psf' in mag_query[0].optical_prov:
                        mags[i][0] = mag_query[0].g
                        mags[i][1] = mag_query[0].r
                        mags[i][2] = mag_query[0].i
                        mags[i][3] = mag_query[0].bp
                        mags[i][4] = mag_query[0].gaia_g
                        mags[i][5] = mag_query[0].rp
                        mags[i][6] = mag_query[0].h
                    else:
                        mags[i][0] = None
                        mags[i][1] = None
                        mags[i][2] = None
                        mags[i][3] = None
                        mags[i][4] = None
                        mags[i][5] = None
                        mags[i][6] = None

            nd = len(field_assign_2['robotID'].shape)
            if nd == 1:
                ndim = 1
            else:
                ndim = len(field_assign_2['robotID'][0, :])

            for i in range(ndim):

                if nd == 1:
                    roboid = field_assign_2['robotID']
                else:
                    roboid = field_assign_2['robotID'][:, i]

                line = ['fieldid',
                        'exp',
                        'cadence',
                        'designmode',
                        'BOSS_skies_min',
                        'BOSS_skies_FOV',
                        'APOGEE_skies_min',
                        'APOGEE_skies_FOV',
                        'BOSS_stds_min',
                        'BOSS_stds_mags',
                        'BOSS_stds_FOV',
                        'APOGEE_stds_min',
                        'APOGEE_stds_mags',
                        'APOGEE_stds_FOV',
                        'BOSS_bright_limit_targets',
                        'APOGEE_bright_limit_targets']
                line[0] = fieldid
                line[1] = i
                line[2] = head['FCADENCE']
                if 'bright' in head['FCADENCE']:
                    line[3] = 'Bright Time'
                elif head['FCADENCE'] == 'dark_174x8' or head['FCADENCE'] == 'dark_100x8':
                    line[3] = 'Dark RM'
                elif head['FCADENCE'] == 'dark_10x4' or head['FCADENCE'] == 'dark_2x4':
                    line[3] = 'Dark Monitoring'
                elif head['FCADENCE'] == 'dark_2x2' or head['FCADENCE'] == 'dark_4x1':
                    line[3] = 'Dark Faint'
                else:
                    line[3] = 'Dark Plane'

                exp_ev = eval("roboid != -1")

                des = FPSDesign(design_pk=-1,
                                obsTime=Time(ObsTime().nominal(lst=head['RACEN'])).jd,
                                racen=head['RACEN'],
                                deccen=head['DECCEN'],
                                position_angle=head['PA'],
                                observatory=head['OBS'].upper(),
                                mode_pk=-1,
                                catalogids=field_assign_1['catalogid'][exp_ev],
                                ra=field_assign_1['ra'][exp_ev],
                                dec=field_assign_1['dec'][exp_ev],
                                fiberID=roboid[exp_ev],
                                obsWavelength=field_assign_1['fiberType'][exp_ev],
                                priority=field_assign_1['priority'][exp_ev],
                                carton_pk=carton_pks[exp_ev],
                                design_file=None,
                                manual_design=True)
                des.build_design_manual()
                des.design['x'] = field_assign_1['x'][exp_ev]
                des.design['y'] = field_assign_1['y'][exp_ev]

                des_check = DesignModeCheck(FPSDesign=des,
                                            desmode_label=line[3],
                                            desmode_manual=desmode_manuals[line[3]],
                                            mags=mags[exp_ev])

                # try:
                # do checks
                line[4] = int(des_check.skies_min('BOSS'))
                line[5] = int(des_check.skies_fov('BOSS'))
                if line[3] == 'Dark RM':
                    line[6] = int(True)
                    line[7] = int(True)
                else:
                    line[6] = int(des_check.skies_min('APOGEE'))
                    line[7] = int(des_check.skies_fov('APOGEE'))
                line[8] = int(des_check.stds_min('BOSS'))
                stds_mags_check = des_check.mag_limits(des_check.stds_mags['BOSS'],
                                                       'BOSS',
                                                       'std')
                check_tot = len(stds_mags_check[0][stds_mags_check[0]])
                design_tot = len(des_check.design['x'][(des_check.design['catalogID'] != -1) &
                                                       (np.isin(des_check.design['carton_pk'],
                                                                des_check.carton_classes['std'])) &
                                                       (des_check.design['obsWavelength'] == 'BOSS')])
                if design_tot == 0:
                    line[9] = -1
                else:
                    line[9] = check_tot / design_tot
                line[10] = int(des_check.stds_fov('BOSS'))
                if line[3] == 'Dark RM':
                    line[11] = int(True)
                    line[12] = 1.
                    line[13] = int(True)
                else:
                    line[11] = int(des_check.stds_min('APOGEE'))
                    stds_mags_check = des_check.mag_limits(des_check.stds_mags['APOGEE'],
                                                           'APOGEE',
                                                           'std')
                    check_tot = len(stds_mags_check[0][stds_mags_check[0]])
                    design_tot = len(des_check.design['x'][(des_check.design['catalogID'] != -1) &
                                                           (np.isin(des_check.design['carton_pk'],
                                                                    des_check.carton_classes['std'])) &
                                                           (des_check.design['obsWavelength'] == 'APOGEE')])
                    if design_tot == 0:
                        line[12] = -1
                    else:
                        line[12] = check_tot / design_tot
                    line[13] = int(des_check.stds_fov('APOGEE'))
                stds_mags_check = des_check.mag_limits(des_check.bright_limit_targets['BOSS'],
                                                       'BOSS',
                                                       'science')
                check_tot = len(stds_mags_check[0][stds_mags_check[0]])
                design_tot = len(des_check.design['x'][(des_check.design['catalogID'] != -1) &
                                                       (np.isin(des_check.design['carton_pk'],
                                                                des_check.carton_classes['science'])) &
                                                       (des_check.design['obsWavelength'] == 'BOSS')])
                if design_tot == 0:
                    line[14] = -1
                else:
                    line[14] = check_tot / design_tot
                if line[3] == 'Dark RM':
                    line[15] = 1.
                else:
                    stds_mags_check = des_check.mag_limits(des_check.bright_limit_targets['APOGEE'],
                                                           'APOGEE',
                                                           'science')
                    check_tot = len(stds_mags_check[0][stds_mags_check[0]])
                    design_tot = len(des_check.design['x'][(des_check.design['catalogID'] != -1) &
                                                           (np.isin(des_check.design['carton_pk'],
                                                                    des_check.carton_classes['science'])) &
                                                           (des_check.design['obsWavelength'] == 'APOGEE')])
                    if design_tot == 0:
                        line[15] = -1
                    else:
                        line[15] = check_tot / design_tot
                # except:
                #     line[4] = -1
                #     line[5] = -1
                #     line[6] = -1
                #     line[7] = -1
                #     line[8] = -1
                #     line[9] = -1
                #     line[10] = -1
                #     line[11] = -1
                #     line[12] = -1
                #     line[13] = -1
                #     line[14] = -1
                #     line[15] = -1
                #     bad += 1

                f.write('\t'.join([str(x) for x in line]) + '\n')

            print('Field: %d,Time: %.3f min' % (fieldid, (time.time() - start)/60), end='\n')
