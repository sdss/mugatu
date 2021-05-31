# @Author: Ilija Medan
# @Date: May 18, 2021
# @Filename: obsmode.py
# @License: BSD 3-Clause
# @Copyright: Ilija Medan

import numpy as np
import warnings
from scipy.spatial import cKDTree

from mugatu.exceptions import MugatuError, MugatuWarning
from sdssdb.peewee.sdss5db.targetdb import DesignMode, Carton, Category, Magnitude, CartonToTarget, Target


class ObsMode(object):
    """
    Parameters
    ----------

    Attributes
    ----------

    Methods
    -------
    """

    def __init__(self, FPSDesign, obsmode_label):
        self.design = FPSDesign.design
        self.obsmode_label = obsmode_label

        # grab the design mode params
        obsmode = DesignMode.get(DesignMode.label == obsmode_label)
        self.n_skies_min = {}
        self.n_skies_min['BOSS'] = obsmode.n_boss_skies_min
        self.n_skies_min['APOGEE'] = obsmode.n_apogee_skies_min

        self.min_skies_fovmetric = {}
        self.min_skies_fovmetric['BOSS'] = obsmode.min_boss_skies_fovmetric
        self.min_skies_fovmetric['APOGEE'] = obsmode.min_apogee_skies_fovmetric

        self.n_stds_min = {}
        self.n_stds_min['BOSS'] = obsmode.n_boss_stds_min
        self.n_stds_min['APOGEE'] = obsmode.n_apogee_stds_min

        self.min_stds_fovmetric = {}
        self.min_stds_fovmetric['BOSS'] = obsmode.min_boss_stds_fovmetric
        self.min_stds_fovmetric['APOGEE'] = obsmode.min_apogee_stds_fovmetric

        self.bright_limit = {}
        self.bright_limit['BOSS'] = obsmode.boss_bright_limit
        self.bright_limit['APOGEE'] = obsmode.apogee_bright_limit

        self.faint_limit = {}
        self.faint_limit['BOSS'] = obsmode.boss_faint_limit
        self.faint_limit['APOGEE'] = obsmode.apogee_faint_limit

        self.epoch_snrgoal = {}
        self.epoch_snrgoal['BOSS'] = obsmode.boss_epoch_snrgoal
        self.epoch_snrgoal['APOGEE'] = obsmode.apogee_epoch_snrgoal

        self.design_minsnr = {}
        self.design_minsnr['BOSS'] = obsmode.boss_design_minsnr
        self.design_minsnr['APOGEE'] = obsmode.apogee_design_minsnr

        # classify cartons as skies, standards or science
        self.carton_classes = {}
        self.carton_classes['science'] = []
        self.carton_classes['sky'] = []
        self.carton_classes['std'] = []
        for pk in np.unique(self.design['carton_pk'][self.design['catalogID'] != -1]):
            carton = Carton.select().join(Category).where(Carton.pk == pk)
            if 'sky' in carton.category.label:
                self.carton_classes['sky'].append(pk)
            elif 'standard' in carton.category.label:
                self.carton_classes['std'].append(pk)
            # I think else makes sense here as there are science
            # and open fiber labels?
            else:
                self.carton_classes['science'].append(pk)

        # collect magntiudes of design
        # here I am doing g,r,i,BP,G,RP
        self.mags = np.zeros((500, 6)) - 9999.99
        for i in range(500):
            if self.design['catalogID'][i] != -1:
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
        bright_checks = np.zeros(500, dtype=bool)
        # if complete, all bands for target present in check
        # if incomplete, then
        complete_check = np.array(['COMPLETE' for _ in range(500)],
                                  dtype='<U12')

        # check which limits are defined for mode
        check_inds = []
        for i in range(len(self.bright_limit[instrument])):
            if self.bright_limit[instrument][i] is not None:
                check_inds.append(i)

        # run checks
        for i in range(bright_checks):
            if self.design['catalogID'][i] != -1:
                # check in each band
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
        faint_checks = np.zeros(500, dtype=bool)
        # if complete, all bands for target present in check
        # if incomplete, then
        complete_check = np.array(['COMPLETE' for _ in range(500)],
                                  dtype='<U12')

        # check which limits are defined for mode
        check_inds = []
        for i in range(len(self.faint_limit[instrument])):
            if self.faint_limit[instrument][i] is not None:
                check_inds.append(i)

        # run checks
        for i in range(faint_checks):
            if self.design['catalogID'][i] != -1:
                # check in each band
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

        self.bright_limit_check = {}
        self.bright_limit_check['BOSS'] = self.bright_limit(instrument='BOSS')
        self.bright_limit_check['APOGEE'] = self.bright_limit(instrument='APOGEE')

        self.faint_limit_check = {}
        self.faint_limit_check['BOSS'] = self.faint_limit(instrument='BOSS')
        self.faint_limit_check['APOGEE'] = self.faint_limit(instrument='APOGEE')

        if verbose:
            design_tot = len(self.design['x'][self.design['catalogID'] != -1])
            verbose_output = ''
            verbose_output += 'DesignMode Param               | Pass Check?\n'
            verbose_output += '-------------------------------|-----------------\n'
            verbose_output += 'N Skies Min (BOSS):            | %s\n' % self.n_skies_min_check['BOSS']
            verbose_output += 'N Skies Min (APOGEE):          | %s\n' % self.n_skies_min_check['APOGEE']
            verbose_output += 'N Standards Min (BOSS):        | %s\n' % self.n_stds_min_check['BOSS']
            verbose_output += 'N Standards Min (APOGEE):      | %s\n' % self.n_stds_min_check['APOGEE']
            verbose_output += 'FOV Metric Skies (BOSS):       | %s\n' % self.min_skies_fovmetric_check['BOSS']
            verbose_output += 'FOV Metric Skies (APOGEE):     | %s\n' % self.min_skies_fovmetric_check['APOGEE']
            verbose_output += 'FOV Metric Standards (BOSS):   | %s\n' % self.min_stds_fovmetric_check['BOSS']
            verbose_output += 'FOV Metric Standards (APOGEE): | %s\n' % self.min_stds_fovmetric_check['APOGEE']
            check_tot = len(self.bright_limit_check['BOSS'][self.bright_limit_check['BOSS']])
            verbose_output += 'Bright Limit (BOSS):           | %d out of %d\n' % (check_tot, design_tot)
            check_tot = len(self.bright_limit_check['APOGEE'][self.bright_limit_check['APOGEE']])
            verbose_output += 'Bright Limit (APOGEE):         | %d out of %d\n' % (check_tot, design_tot)
            check_tot = len(self.faint_limit_check['BOSS'][self.faint_limit_check['BOSS']])
            verbose_output += 'Faint Limit (BOSS):            | %d out of %d\n' % (check_tot, design_tot)
            check_tot = len(self.faint_limit_check['APOGEE'][self.faint_limit_check['APOGEE']])
            verbose_output += 'Faint Limit (APOGEE):          | %d out of %d\n' % (check_tot, design_tot)
            print(verbose_output)
