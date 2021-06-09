# @Author: Ilija Medan
# @Date: May 18, 2021
# @Filename: obsmode.py
# @License: BSD 3-Clause
# @Copyright: Ilija Medan

import numpy as np
import warnings
from scipy.spatial import cKDTree

from mugatu.exceptions import MugatuError, MugatuWarning
# from sdssdb.peewee.sdss5db.targetdb import DesignMode
from sdssdb.peewee.sdss5db.targetdb import Carton, Category, Magnitude, CartonToTarget, Target
from sdssdb.peewee.sdss5db import catalogdb


def ang_sep(ra1, dec1, ra2, dec2):
    ra1 = np.radians(ra1)
    dec1 = np.radians(dec1)
    ra2 = np.radians(ra2)
    dec2 = np.radians(dec2)
    return (180 / np.pi) * np.arcos(np.sin(dec1) * np.sin(dec2) +
                                    np.cos(dec1) * np.cos(dec2) * np.cos(ra1 - ra2))


class DesignModeCheck(object):
    """
    Parameters
    ----------

    Attributes
    ----------

    Methods
    -------
    """

    def __init__(self, FPSDesign, desmode_label,
                 desmode_manual=None):
        self.design = FPSDesign.design
        self.desmode_label = desmode_label

        # grab the design mode params
        if desmode_manual is None:
            # uncomment once in database
            # desmode = DesignMode.get(DesignMode.label == desmode_label)
            desmode = desmode_manual
            self.n_skies_min = {}
            self.n_skies_min['BOSS'] = desmode.boss_skies_min
            self.n_skies_min['APOGEE'] = desmode.apogee_skies_min

            self.min_skies_fovmetric = {}
            self.min_skies_fovmetric['BOSS'] = desmode.boss_skies_fov
            self.min_skies_fovmetric['APOGEE'] = desmode.apogee_skies_fov

            self.n_stds_min = {}
            self.n_stds_min['BOSS'] = desmode.boss_stds_min
            self.n_stds_min['APOGEE'] = desmode.apogee_stds_min

            self.min_stds_fovmetric = {}
            self.min_stds_fovmetric['BOSS'] = desmode.boss_stds_fov
            self.min_stds_fovmetric['APOGEE'] = desmode.apogee_stds_fov

            self.stds_mags = {}
            self.stds_mags['BOSS'] = desmode.boss_stds_mags
            self.stds_mags['APOGEE'] = desmode.apogee_stds_mags

            self.bright_limit_targets = {}
            self.bright_limit_targets['BOSS'] = desmode.boss_bright_limit_targets
            self.bright_limit_targets['APOGEE'] = desmode.apogee_bright_limit_targets

            self.sky_neighbors_targets = {}
            self.sky_neighbors_targets['BOSS'] = desmode.boss_sky_neighbors_targets
            self.sky_neighbors_targets['APOGEE'] = desmode.apogee_sky_neighbors_targets

            self.trace_diff_targets = {}
            self.trace_diff_targets['APOGEE'] = desmode.apogee_trace_diff_targets
        else:
            self.n_skies_min = {}
            self.n_skies_min['BOSS'] = desmode_manual['boss_skies_min']
            self.n_skies_min['APOGEE'] = desmode_manual['apogee_skies_min']

            self.min_skies_fovmetric = {}
            self.min_skies_fovmetric['BOSS'] = desmode_manual['boss_skies_fov']
            self.min_skies_fovmetric['APOGEE'] = desmode_manual['apogee_skies_fov']

            self.n_stds_min = {}
            self.n_stds_min['BOSS'] = desmode_manual['boss_stds_min']
            self.n_stds_min['APOGEE'] = desmode_manual['apogee_stds_min']

            self.min_stds_fovmetric = {}
            self.min_stds_fovmetric['BOSS'] = desmode_manual['boss_stds_fov']
            self.min_stds_fovmetric['APOGEE'] = desmode_manual['apogee_stds_fov']

            self.stds_mags = {}
            self.stds_mags['BOSS'] = desmode_manual['boss_stds_mags']
            self.stds_mags['APOGEE'] = desmode_manual['apogee_stds_mags']

            self.bright_limit_targets = {}
            self.bright_limit_targets['BOSS'] = desmode_manual['boss_bright_limit_targets']
            self.bright_limit_targets['APOGEE'] = desmode_manual['apogee_bright_limit_targets']

            self.sky_neighbors_targets = {}
            self.sky_neighbors_targets['BOSS'] = desmode_manual['boss_sky_neighbors_targets']
            self.sky_neighbors_targets['APOGEE'] = desmode_manual['apogee_sky_neighbors_targets']

            self.trace_diff_targets = {}
            self.trace_diff_targets['APOGEE'] = desmode_manual['apogee_trace_diff_targets']

        # classify cartons as skies, standards or science
        self.carton_classes = {}
        self.carton_classes['science'] = []
        self.carton_classes['sky'] = []
        self.carton_classes['std'] = []
        for pk in np.unique(self.design['carton_pk'][self.design['catalogID'] != -1]):
            carton = Category.select().join(Carton).where(Carton.pk == pk)[0]
            if 'sky' in carton.label:
                self.carton_classes['sky'].append(pk)
            elif 'standard' in carton.label:
                self.carton_classes['std'].append(pk)
            # I think else makes sense here as there are science
            # and open fiber labels?
            else:
                self.carton_classes['science'].append(pk)

        # collect magntiudes of design
        # here I am doing g,r,i,BP,G,RP,H
        self.mags = np.zeros((len(self.design['catalogID']), 7)) - 9999.99
        for i in range(len(self.design['catalogID'])):
            if self.design['catalogID'][i] != -1:
                # do not query skies, they have no mag
                if self.design['carton_pk'][i] not in self.carton_classes['sky']:
                    mag_query = (Magnitude.select()
                                          .join(CartonToTarget)
                                          .join(Target)
                                          .where((Target.catalogid == self.design['catalogID'][i]) &
                                                 (CartonToTarget.carton_pk == self.design['carton_pk'][i])))
                    self.mags[i][0] = mag_query[0].g
                    self.mags[i][1] = mag_query[0].r
                    self.mags[i][2] = mag_query[0].i
                    self.mags[i][3] = mag_query[0].bp
                    self.mags[i][4] = mag_query[0].gaia_g
                    self.mags[i][5] = mag_query[0].rp
                    self.mags[i][6] = mag_query[0].h

    def skies_min(self, instrument):
        n_skies = len(self.design['catalogID'][(self.design['catalogID'] != -1) &
                                                     (np.isin(self.design['carton_pk'],
                                                              self.carton_classes['sky'])) &
                                                     (self.design['obsWavelength'] == instrument)])
        if n_skies >= self.n_skies_min[instrument]:
            return True
        else:
            return False

    def stds_min(self, instrument):
        n_stds = len(self.design['catalogID'][(self.design['catalogID'] != -1) &
                                                    (np.isin(self.design['carton_pk'],
                                                             self.carton_classes['std'])) &
                                                    (self.design['obsWavelength'] == instrument)])
        if n_stds >= self.n_stds_min[instrument]:
            return True
        else:
            return False

    def skies_fov(self, instrument):
        # get x,y of the skies
        x_sky = self.design['x'][(self.design['catalogID'] != -1) &
                                 (np.isin(self.design['carton_pk'],
                                          self.carton_classes['sky'])) &
                                 (self.design['obsWavelength'] == instrument)]
        y_sky = self.design['y'][(self.design['catalogID'] != -1) &
                                 (np.isin(self.design['carton_pk'],
                                          self.carton_classes['sky'])) &
                                 (self.design['obsWavelength'] == instrument)]

        x_sci = self.design['x'][(self.design['catalogID'] != -1) &
                                 (np.isin(self.design['carton_pk'],
                                          self.carton_classes['science'])) &
                                 (self.design['obsWavelength'] == instrument)]
        y_sci = self.design['y'][(self.design['catalogID'] != -1) &
                                 (np.isin(self.design['carton_pk'],
                                          self.carton_classes['science'])) &
                                 (self.design['obsWavelength'] == instrument)]

        # create KDE tree
        tree = cKDTree(np.column_stack((x_sky, y_sky)))
        # get distances for nearest neighbors
        dd, ii = tree.query(np.column_stack((x_sci, y_sci)),
                            k=self.min_skies_fovmetric[instrument][0])
        # second column is the nth neighbor distance if k>1
        if self.min_skies_fovmetric[instrument][0] == 1:
            dists = dd
        else:
            dists = dd[:, 1]
        # this assumes percentile is on 0 to 100 scale
        perc_dist = np.percentile(dists,
                                  self.min_stds_fovmetric[instrument][1])
        if perc_dist < self.min_stds_fovmetric[instrument][2]:
            return True
        else:
            return False

    def stds_fov(self, instrument):
        # get x,y of the standards
        x_std = self.design['x'][(self.design['catalogID'] != -1) &
                                 (np.isin(self.design['carton_pk'],
                                          self.carton_classes['std'])) &
                                 (self.design['obsWavelength'] == instrument)]
        y_std = self.design['y'][(self.design['catalogID'] != -1) &
                                 (np.isin(self.design['carton_pk'],
                                          self.carton_classes['std'])) &
                                 (self.design['obsWavelength'] == instrument)]

        x_sci = self.design['x'][(self.design['catalogID'] != -1) &
                                 (np.isin(self.design['carton_pk'],
                                          self.carton_classes['science'])) &
                                 (self.design['obsWavelength'] == instrument)]
        y_sci = self.design['y'][(self.design['catalogID'] != -1) &
                                 (np.isin(self.design['carton_pk'],
                                          self.carton_classes['science'])) &
                                 (self.design['obsWavelength'] == instrument)]

        # create KDE tree
        tree = cKDTree(np.column_stack((x_std, y_std)))
        # get distances for nearest neighbors
        dd, ii = tree.query(np.column_stack((x_sci, y_sci)),
                            k=self.min_skies_fovmetric[instrument][0])
        # second column is the nth neighbor distance if k>1
        if self.min_skies_fovmetric[instrument][0] == 1:
            dists = dd
        else:
            dists = dd[:, 1]
        # this assumes percentile is on 0 to 100 scale
        perc_dist = np.percentile(dists,
                                  self.min_skies_fovmetric[instrument][1])
        if perc_dist < self.min_skies_fovmetric[instrument][2]:
            return True
        else:
            return False

    def bright_limit(self, instrument):
        bright_checks = np.zeros(len(self.design['catalogID']),
                                 dtype=bool)
        # if complete, all bands for target present in check
        # if incomplete, then
        complete_check = np.array(['COMPLETE' for _ in range(len(self.design['catalogID']))],
                                  dtype='<U12')


        # check which limits are defined for mode
        check_inds = []
        for i in range(len(self.bright_limit[instrument])):
            if self.bright_limit[instrument][i] is not None:
                check_inds.append(i)

        # run checks
        for i in range(len(bright_checks)):
            if self.design['catalogID'][i] != -1:
                # check in each band that has check defined
                targ_check = np.zeros(len(check_inds), dtype=bool)
                for j, ind in enumerate(check_inds):
                    if self.mags[i][ind] is None:
                        complete_check[i] = 'INCOMPLETE'
                        # set True, no mag is not a fail
                        targ_check[j] = True
                    elif self.mags[i][ind] > self.bright_limit[instrument][ind]:
                        targ_check[j] = True
                # if all True, then passes
                if np.all(targ_check):
                    bright_checks[i] = True
            else:
                complete_check[i] = 'INCOMPLETE'

        return bright_checks, complete_check

    def faint_limit(self, instrument):
        faint_checks = np.zeros(len(self.design['catalogID']),
                                dtype=bool)
        # if complete, all bands for target present in check
        # if incomplete, then
        complete_check = np.array(['COMPLETE' for _ in range(len(self.design['catalogID']))],
                                  dtype='<U12')


        # check which limits are defined for mode
        check_inds = []
        for i in range(len(self.faint_limit[instrument])):
            if self.faint_limit[instrument][i] is not None:
                check_inds.append(i)

        # run checks
        for i in range(len(faint_checks)):
            if self.design['catalogID'][i] != -1:
                # check in each band that has check defined
                targ_check = np.zeros(len(check_inds), dtype=bool)
                for j, ind in enumerate(check_inds):
                    if self.mags[i][ind] is None:
                        complete_check[i] = 'INCOMPLETE'
                        # set True, no mag is not a fail
                        targ_check[j] = True
                    elif self.mags[i][ind] < self.faint_limit[instrument][ind]:
                        targ_check[j] = True
                # if all True, then passes
                if np.all(targ_check):
                    faint_checks[i] = True
            else:
                complete_check[i] = 'INCOMPLETE'

        return faint_checks, complete_check

    def mag_limits(self, mag_metric,
                   instrument, carton_class):
        mag_checks = np.zeros(len(self.design['catalogID']),
                              dtype=bool)
        # if complete, all bands for target present in check
        # if incomplete, then
        complete_check = np.array(['COMPLETE' for _ in range(len(self.design['catalogID']))],
                                  dtype='<U12')

        # check which limits are defined for mode
        check_inds = []
        for i in range(mag_metric.shape[0]):
            if mag_metric[i][0] is not None or mag_metric[i][1] is not None:
                check_inds.append(i)

        # run checks
        for i in range(len(mag_checks)):
            if (self.design['catalogID'][i] != -1 and
                self.design['carton_pk'][i] in self.carton_classes[carton_class] and
                self.design['obsWavelength'][i] == instrument):
                # check in each band that has check defined
                targ_check = np.zeros(len(check_inds), dtype=bool)
                for j, ind in enumerate(check_inds):
                    if self.mags[i][ind] is None:
                        complete_check[i] = 'INCOMPLETE'
                        # set True, no mag is not a fail
                        targ_check[j] = True
                    # check when greater than and less than
                    elif mag_metric[ind][0] is not None and mag_metric[ind][1] is not None:
                        if (mag_metric[ind][0] < self.mags[i][ind] < mag_metric[ind][1]):
                            targ_check[j] = True
                    # check when just greater than
                    elif mag_metric[ind][0] is not None:
                        if self.mags[i][ind] > mag_metric[ind][0]:
                            targ_check[j] = True
                    # check when less than
                    else:
                        if self.mags[i][ind] < mag_metric[ind][1]:
                            targ_check[j] = True
                # if all True, then passes
                if np.all(targ_check):
                    mag_checks[i] = True
            else:
                complete_check[i] = 'INCOMPLETE'
        return mag_checks, complete_check

    def sky_neighbors(self, instrument, catalogdb_ver):
        sky_checks = np.zeros(len(self.design['catalogID']),
                              dtype=bool)
        # set the columns/catalog based on instrument
        if instrument == 'BOSS':
            cat = catalogdb.Gaia_DR2
            ra_col = catalogdb.Gaia_DR2.ra
            dec_col = catalogdb.Gaia_DR2.dec
            mag_col = catalogdb.Gaia_DR2.phot_g_mean_mag
        else:
            cat = catalogdb.TwoMassPSC
            ra_col = catalogdb.TwoMassPSC.ra
            dec_col = catalogdb.TwoMassPSC.decl
            mag_col = catalogdb.TwoMassPSC.h_m

        # grab the params for check
        R_0 = self.sky_neighbors_targets[instrument][0]
        beta = self.sky_neighbors_targets[instrument][1]
        lim = self.sky_neighbors_targets[instrument][2]

        # go through all skies
        # grab magnitudes within 1' of skies
        # NOTE (is 1' big enough to catch everything?)
        for i in range(len(sky_checks)):
            if (self.design['catalogID'][i] != -1 and
                self.design['carton_pk'][i] in self.carton_classes['sky'] and
                self.design['obsWavelength'][i] == instrument):
                sky_neigh = (cat.select(ra_col,
                                        dec_col,
                                        mag_col,
                                        catalogdb.Catalog.catalogid)
                                .join(catalogdb.TIC_v8)
                                .join(catalogdb.CatalogToTIC_v8)
                                .join(catalogdb.Catalog)
                                .join(catalogdb.Version)
                                .where((cat.cone_search(self.design['ra'][i],
                                                        self.design['dec'][i],
                                                        1 / 60)) &
                                       (catalogdb.Version.id == catalogdb_ver)))
                # get result
                ras, decs, mags, catalogids = map(list, zip(*list(sky_neigh.tuples())))
                # remove the sky if in result
                ev_sky = np.isin(catalogids, [self.design['catalogID'][i]])
                ras = np.array(ras)[ev_sky]
                decs = np.array(decs)[ev_sky]
                mags = np.array(mags)[ev_sky]

                dists = 3600 * ang_sep(self.design['ra'][i],
                                       self.design['dec'][i],
                                       ras,
                                       decs)
                if np.all(dists > R_0 * (lim - mags) ** beta):
                    sky_checks[i] = True
        return sky_checks

    def design_mode_check_all(self, verbose=True):
        self.n_skies_min_check = {}
        self.n_skies_min_check['BOSS'] = self.skies_min(instrument='BOSS')
        self.n_skies_min_check['APOGEE'] = self.skies_min(instrument='APOGEE')

        self.min_skies_fovmetric_check = {}
        self.min_skies_fovmetric_check['BOSS'] = self.skies_fov(instrument='BOSS')
        self.min_skies_fovmetric_check['APOGEE'] = self.skies_fov(instrument='APOGEE')

        self.n_stds_min_check = {}
        self.n_stds_min_check['BOSS'] = self.stds_min(instrument='BOSS')
        self.n_stds_min_check['APOGEE'] = self.stds_min(instrument='APOGEE')

        self.min_stds_fovmetric_check = {}
        self.min_stds_fovmetric_check['BOSS'] = self.stds_fov(instrument='BOSS')
        self.min_stds_fovmetric_check['APOGEE'] = self.stds_fov(instrument='APOGEE')

        self.stds_mags_check = {}
        self.stds_mags_check['BOSS'] = self.mag_limits(self.stds_mags['BOSS'],
                                                       'BOSS',
                                                       'std')
        self.stds_mags_check['APOGEE'] = self.mag_limits(self.stds_mags['APOGEE'],
                                                         'APOGEE',
                                                         'std')

        self.bright_limit_targets_check = {}
        self.bright_limit_targets_check['BOSS'] = self.mag_limits(self.bright_limit_targets['BOSS'],
                                                                  'BOSS',
                                                                  'science')
        self.bright_limit_targets_check['APOGEE'] = self.mag_limits(self.bright_limit_targets['APOGEE'],
                                                                    'APOGEE',
                                                                    'science')

        if verbose:
            verbose_output = ''
            verbose_output += 'DesignMode Param                  | Pass Check?\n'
            verbose_output += '----------------------------------|-----------------\n'
            verbose_output += 'N Skies Min (BOSS):               | %s\n' % self.n_skies_min_check['BOSS']
            verbose_output += 'N Skies Min (APOGEE):             | %s\n' % self.n_skies_min_check['APOGEE']
            verbose_output += 'N Standards Min (BOSS):           | %s\n' % self.n_stds_min_check['BOSS']
            verbose_output += 'N Standards Min (APOGEE):         | %s\n' % self.n_stds_min_check['APOGEE']
            verbose_output += 'FOV Metric Skies (BOSS):          | %s\n' % self.min_skies_fovmetric_check['BOSS']
            verbose_output += 'FOV Metric Skies (APOGEE):        | %s\n' % self.min_skies_fovmetric_check['APOGEE']
            verbose_output += 'FOV Metric Standards (BOSS):      | %s\n' % self.min_stds_fovmetric_check['BOSS']
            verbose_output += 'FOV Metric Standards (APOGEE):    | %s\n' % self.min_stds_fovmetric_check['APOGEE']

            check_tot = len(self.stds_mags_check['BOSS'][0][self.stds_mags_check['BOSS'][0]])
            design_tot = len(self.design['x'][(self.design['catalogID'] != -1) &
                                              (np.isin(self.design['carton_pk'],
                                                       self.carton_classes['std'])) &
                                              (self.design['obsWavelength'] == 'BOSS')])
            verbose_output += 'Magnitude Limit Stds (BOSS):      | %d out of %d\n' % (check_tot, design_tot)

            check_tot = len(self.stds_mags_check['APOGEE'][0][self.stds_mags_check['APOGEE'][0]])
            design_tot = len(self.design['x'][(self.design['catalogID'] != -1) &
                                              (np.isin(self.design['carton_pk'],
                                                       self.carton_classes['std'])) &
                                              (self.design['obsWavelength'] == 'APOGEE')])
            verbose_output += 'Magnitude Limit Stds (APOGEE):    | %d out of %d\n' % (check_tot, design_tot)

            check_tot = len(self.bright_limit_targets_check['BOSS'][0][self.bright_limit_targets_check['BOSS'][0]])
            design_tot = len(self.design['x'][(self.design['catalogID'] != -1) &
                                              (np.isin(self.design['carton_pk'],
                                                       self.carton_classes['science'])) &
                                              (self.design['obsWavelength'] == 'BOSS')])
            verbose_output += 'Magnitude Limit Targets (BOSS):   | %d out of %d\n' % (check_tot, design_tot)

            check_tot = len(self.bright_limit_targets_check['APOGEE'][0][self.bright_limit_targets_check['APOGEE'][0]])
            design_tot = len(self.design['x'][(self.design['catalogID'] != -1) &
                                              (np.isin(self.design['carton_pk'],
                                                       self.carton_classes['science'])) &
                                              (self.design['obsWavelength'] == 'APOGEE')])
            verbose_output += 'Magnitude Limit Targets (APOGEE): | %d out of %d\n' % (check_tot, design_tot)

            print(verbose_output)
