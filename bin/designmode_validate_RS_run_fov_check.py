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
from sdssdb.peewee.sdss5db.targetdb import Carton, Category, Magnitude, CartonToTarget, Target
import time
from scipy.spatial import cKDTree

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


def stds_fov(instrument, design, carton_classes,
             min_stds_fovmetric):
    """
    Checks if design meets the FOV metric for the standards
    for some instrument. Returns True FOV metric met,
    False if not. Also, if no science targets in design
    returns True, while no standards returns False.

    Parameters
    ----------
    instrument: str
        Instrument to FOV metric.
        Must be 'BOSS' or 'APOGEE'.

    Returns
    -------
    : boolean
        True FOV metric met,
        False if not. Also, if no science targets in design
        returns True, while no standards returns False.
    """
    # get x,y of the standards
    x_std = design['x'][(design['catalogID'] != -1) &
                        (np.isin(design['carton_pk'],
                                 carton_classes['std'])) &
                        (design['obsWavelength'] == instrument)]
    y_std = design['y'][(design['catalogID'] != -1) &
                        (np.isin(design['carton_pk'],
                                 carton_classes['std'])) &
                        (design['obsWavelength'] == instrument)]

    x_sci = design['x'][(design['catalogID'] != -1) &
                        (np.isin(design['carton_pk'],
                                 carton_classes['science'])) &
                        (design['obsWavelength'] == instrument)]
    y_sci = design['y'][(design['catalogID'] != -1) &
                        (np.isin(design['carton_pk'],
                                 carton_classes['science'])) &
                        (design['obsWavelength'] == instrument)]

    # if no skies, dont do check
    if len(x_std) == 0:
        return -1.
    # if no science in band, doesnt matter?
    if len(x_sci) == 0:
        return -2.

    # create KDE tree
    tree = cKDTree(np.column_stack((x_std, y_std)))
    # get distances for nearest neighbors
    dd, ii = tree.query(np.column_stack((x_sci, y_sci)),
                        k=min_stds_fovmetric[instrument][0])
    # second column is the nth neighbor distance if k>1
    if min_stds_fovmetric[instrument][0] == 1:
        dists = dd
    else:
        dists = dd[:, 1]
    # this assumes percentile is on 0 to 100 scale
    perc_dist = np.percentile(dists,
                              min_stds_fovmetric[instrument][1])
    return perc_dist


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

    with open('rs_%s_%s_designmode_FOV_results_sdssdb_0_4_9_corr_fov_func.txt' % (plan, observatory), 'w') as f:
        header = ['fieldid',
                  'exp',
                  'cadence',
                  'designmode',
                  'n_APOGEE_stds',
                  'APOGGEE_std_FOV_r',
                  'n_BOSS_stds',
                  'BOSS_std_FOV_r']
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

            carton_classes = {}
            carton_classes['science'] = []
            carton_classes['sky'] = []
            carton_classes['std'] = []
            for pk in np.unique(carton_pks):
                carton = Category.select().join(Carton).where(Carton.pk == pk)[0]
                if 'sky' in carton.label:
                    carton_classes['sky'].append(pk)
                elif 'standard' in carton.label:
                    carton_classes['std'].append(pk)
                # I think else makes sense here as there are science
                # and open fiber labels?
                else:
                    carton_classes['science'].append(pk)

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
                        'n_APOGEE_stds',
                        'APOGGEE_std_FOV_r',
                        'n_BOSS_stds',
                        'BOSS_std_FOV_r']
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

                # carton_classes = {}
                # carton_classes['science'] = []
                # carton_classes['sky'] = []
                # carton_classes['std'] = []
                # for pk in np.unique(des.design['carton_pk'][des.design['catalogID'] != -1]):
                #     carton = Category.select().join(Carton).where(Carton.pk == pk)[0]
                #     if 'sky' in carton.label:
                #         carton_classes['sky'].append(pk)
                #     elif 'standard' in carton.label:
                #         carton_classes['std'].append(pk)
                #     # I think else makes sense here as there are science
                #     # and open fiber labels?
                #     else:
                #         carton_classes['science'].append(pk)

                min_stds_fovmetric = {}
                min_stds_fovmetric['BOSS'] = desmode_manuals[line[3]]['boss_stds_fov']
                min_stds_fovmetric['APOGEE'] = desmode_manuals[line[3]]['apogee_stds_fov']

                if line[3] != 'Dark RM':
                    line[5] = stds_fov('APOGEE',
                                       des.design,
                                       carton_classes,
                                       min_stds_fovmetric)
                    n_stds = len(des.design['catalogID'][(des.design['catalogID'] != -1) &
                                                         (np.isin(des.design['carton_pk'],
                                                                  carton_classes['std'])) &
                                                         (des.design['obsWavelength'] == 'APOGEE')])
                    line[4] = n_stds
                else:
                    line[5] = -1.
                    line[4] = -1.

                line[7] = stds_fov('BOSS',
                                   des.design,
                                   carton_classes,
                                   min_stds_fovmetric)
                n_stds = len(des.design['catalogID'][(des.design['catalogID'] != -1) &
                                                     (np.isin(des.design['carton_pk'],
                                                              carton_classes['std'])) &
                                                     (des.design['obsWavelength'] == 'BOSS')])
                line[6] = n_stds

                f.write('\t'.join([str(x) for x in line]) + '\n')

            print('Field: %d,Time: %.3f min' % (fieldid, (time.time() - start)/60), end='\n')
