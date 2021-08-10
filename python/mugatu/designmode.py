# @Author: Ilija Medan
# @Date: May 18, 2021
# @Filename: obsmode.py
# @License: BSD 3-Clause
# @Copyright: Ilija Medan

import numpy as np
import warnings
import fitsio
import collections
from scipy.spatial import cKDTree

from mugatu.exceptions import MugatuError, MugatuWarning
from sdssdb.peewee.sdss5db.targetdb import Carton, Category, Magnitude, CartonToTarget, Target
from sdssdb.peewee.sdss5db.targetdb import DesignMode as DesignModeDB
from sdssdb.peewee.sdss5db import catalogdb

try:
    from sdssdb.peewee.sdss5db import database
    database.set_profile('operations')
    _database = True
except:
    _database = False


def ang_sep(ra1, dec1, ra2, dec2):
    """
    Returns angular separation between objects.

    Parameters
    ----------
    ra1: float or np.array
        Right ascension(s) of first object(s).

    dec1: float or np.array
        Declination(s) of first object(s).

    ra2: float or np.array
        Right ascension(s) of second object(s).

    dec2: float or np.array
        Declination(s) of second object(s).

    Returns
    -------
    sep: float or np.array
        Angular seperation between objects in degrees.
    """
    ra1 = np.radians(ra1)
    dec1 = np.radians(dec1)
    ra2 = np.radians(ra2)
    dec2 = np.radians(dec2)
    sep = (180 / np.pi) * np.arcos(np.sin(dec1) * np.sin(dec2) +
                                   np.cos(dec1) * np.cos(dec2) * np.cos(ra1 - ra2))
    return sep


def allDesignModes(filename=None, ext=1):
    """Function to return a dictionary with all design modes

    Parameters
    ----------
    
    filename : str
        FITS file name to read array from (default None)

    ext : int or str
        FITS file extension to read array from (default 1)

    Comments
    --------

    Returns an OrderedDict
    
    If filename is not provided, it reads from that file name
    and the given extension. If not, it reads from the database.

    If neither the filename is given nor the db is available, it
    returns None.
"""
    if(filename is not None):
        dmas = fitsio.read(filename, ext=ext)
        dmd = collections.OrderedDict()
        for dma in dmas:
            dmd[dma['label']] = DesignMode()
            dmd[dma['label']].fromarray(dma)
        return(dmd)
            
    if(_database | (filename is not None)):
        desmodes = DesignModeDB.select()
        dmd = collections.OrderedDict()
        for desmode in desmodes:
            dmd[desmode.label] = DesignMode(label=desmode.label)
        return(dmd)

    print("No file name input or db access.")
    return(None)


class DesignMode(object):
    """Class to store parameters for a design mode

    Parameters
    ----------

    label : str
        label for design mode

    Attributes
    ----------

    n_stds_min: dict
        Dictonary with the minimum number of standards for a design
        for each instrument ('APOGEE' and 'BOSS').

    min_stds_fovmetric: dict
        Dictonary wih the FOV metric for the standards in a design
        for each instrument ('APOGEE' and 'BOSS').
        The FOV metric is described by three parameters, the
        nth neighbor to get distances to, the percentle distance
        to calculate and the distace to compare against for validation
        (in mm).

    stds_mags: dict
        Dictonary for the min/max magnitude for the standards in a
        design for each instrument ('APOGEE' and 'BOSS'). Indexes
        correspond to magntidues: [g, r, i, bp, gaia_g, rp, h].

    bright_limit_targets: dict
        Dictonary for the min/max magnitude for the science targets
        in adesign for each instrument ('APOGEE' and 'BOSS'). Indexes
        correspond to magntidues: [g, r, i, bp, gaia_g, rp, h].

    sky_neighbors_targets: dict
        Dictonary for the parameters used to check distance between
        skies and all possible sources in field for each instrument
        ('APOGEE' and 'BOSS'). Distances to targets (r, in arcseconds)
        must be r > R_0 * (lim - mags) ** beta, where mags is G band
        for BOSS and H band for APOGEE. Indexes correspond to:
        [R_0, beta, lim]

    trace_diff_targets: dict
        Dictonary for the maximum magnitude difference allowed between
        fibers next to each ther on the chip for each instrument
        ('APOGEE' and 'BOSS'). Here the magntidue difference is checked
        in the G band for BOSS and H band for APOGEE.
"""
    def __init__(self, label=None):
        if(label is not None):
            self.fromdb(label=label)
        return

    def fromdb(self, label=None):
        """Read in parameters for design mode from db

        Parameters
        ----------
        
        label : str
            name of design mode
"""
        if(_database is False):
            print("No database, cannot read.")
            return

        self.desmode_label = label

        desmode = DesignModeDB.select().where(DesignModeDB.label == self.desmode_label)[0]

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
        self.stds_mags['BOSS'] = np.zeros((len(desmode.boss_stds_mags_min), 2), dtype=np.float64)
        self.stds_mags['BOSS'][:, 0] = desmode.boss_stds_mags_min
        self.stds_mags['BOSS'][:, 1] = desmode.boss_stds_mags_max
        self.stds_mags['APOGEE'] = np.zeros((len(desmode.apogee_stds_mags_min), 2), dtype=np.float64)
        self.stds_mags['APOGEE'][:, 0] = desmode.apogee_stds_mags_min
        self.stds_mags['APOGEE'][:, 1] = desmode.apogee_stds_mags_max
        
        self.bright_limit_targets = {}
        self.bright_limit_targets['BOSS'] = np.zeros((len(desmode.boss_bright_limit_targets_min), 2), dtype=np.float64)
        self.bright_limit_targets['BOSS'][:, 0] = desmode.boss_bright_limit_targets_min
        self.bright_limit_targets['BOSS'][:, 1] = desmode.boss_bright_limit_targets_max
        self.bright_limit_targets['APOGEE'] = np.zeros((len(desmode.apogee_bright_limit_targets_min), 2), dtype=np.float64)
        self.bright_limit_targets['APOGEE'][:, 0] = desmode.apogee_bright_limit_targets_min
        self.bright_limit_targets['APOGEE'][:, 1] = desmode.apogee_bright_limit_targets_max
        
        self.sky_neighbors_targets = {}
        self.sky_neighbors_targets['BOSS'] = desmode.boss_sky_neighbors_targets
        self.sky_neighbors_targets['APOGEE'] = desmode.apogee_sky_neighbors_targets
        
        self.trace_diff_targets = {}
        self.trace_diff_targets['APOGEE'] = desmode.apogee_trace_diff_targets
        return

    def toarray(self):
        """Returns an ndarray form of the object
"""
        attrnames = ['n_skies_min', 'min_skies_fovmetric', 'n_stds_min',
                     'min_stds_fovmetric', 'stds_mags',
                     'bright_limit_targets', 'sky_neighbors_targets',
                     'trace_diff_targets']

        dtype = [('label', 'U30')]
        for attrname in attrnames:
            attr = getattr(self, attrname)
            for instrument in attr:
                name = instrument.lower() + '_' + attrname
                if(type(attr[instrument]) == np.ndarray):
                    dtype.append((name, np.float64,
                                  attr[instrument].shape))
                elif(type(attr[instrument]) == list):
                    dtype.append((name, np.float64,
                                  len(attr[instrument])))
                else:
                    dtype.append((name, np.float64))
        
        arr = np.zeros(1, dtype=dtype)
        arr['label'] = self.desmode_label
        for attrname in attrnames:
            attr = getattr(self, attrname)
            for instrument in attr:
                name = instrument.lower() + '_' + attrname
                arr[name] = attr[instrument]

        return(arr[0])

    def fromarray(self, arr=None):
        """Reads an ndarray form of the object created by toarray()
"""
        attrnames = ['n_skies_min', 'min_skies_fovmetric', 'n_stds_min',
                     'min_stds_fovmetric', 'stds_mags',
                     'bright_limit_targets', 'sky_neighbors_targets',
                     'trace_diff_targets']

        self.desmode_label = arr['label']
        for attrname in attrnames:
            setattr(self, attrname, dict())
            for instrument in ['BOSS', 'APOGEE']:
                name = instrument.lower() + '_' + attrname
                if(name in arr.dtype.names):
                    getattr(self, attrname)[instrument] = arr[name]
        
        return

    def todict(self):
        """Returns a dictionary form of the object
"""
        attrs = ['n_skies_min', 'min_skies_fovmetric', 'n_stds_min',
                 'min_stds_fovmetric', 'stds_mags', 'bright_limit_targets',
                 'sky_neighbors_targets', 'trace_diff_targets']

        dmd = dict()
        dmd['label'] = self.desmode_label
        for attr in attrs:
            dmda = getattr(self, attr)
            dmd[attr] = dict()
            for key in dmda:
                if(type(dmda[key]) == np.ndarray):
                    dmd[attr][key] = dmda[key].tolist()
                else:
                    dmd[attr][key] = dmda[key]

        return(dmd)

    def fromdict(self, designmode_dict=None):
        """Builds from a dictionary form of the object

        Parameters
        ----------
        
        designmode_dict : dict
            dictionary form of DesignMode object created by todict()
"""
        attrs = ['n_skies_min', 'min_skies_fovmetric', 'n_stds_min',
                 'min_stds_fovmetric', 'stds_mags', 'bright_limit_targets',
                 'sky_neighbors_targets', 'trace_diff_targets']

        self.desmod_label = designmode_dict['label']
        for param in attrs:
            setattr(self, param, dict())
            for instrument in designmode_dict[param]:
                if(type(designmode_dict[param][instrument]) == list):
                    getattr(self, param)[instrument] = np.array(designmode_dict[param][instrument])
                else:
                    getattr(self, param)[instrument] = designmode_dict[param][instrument]
        return

    def frommanual(self, label=None, desmode_manual=None):
        """Read in parameters for design mode from dictionary input

        Parameters
        ----------
        
        label : str
            name of design mode

        desmode_manual : dict
            Dictonary of DesignMode parameters to be used as manual
            inputs to validate the design, rather than from targetdb.
"""
        self.desmode_label = label

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
        return
        
        
def check_assign_mag_limit(mag_metric_min,
                           mag_metric_max,
                           assign_mag):
    """
    Checks the if magnitude of one assignment agrees with
    design mode for some instrument and carton class

    Parameters
    ----------
    mag_metric_min: float
        The minimum magitude for the specific DesignMode.

    mag_metric_min: float
        The maximum magitude for the specific DesignMode.

    assign_mag: float
        The magntiude of the assignment being checkaged against
        the DesignMode.

    Returns
    -------
    targ_check: boolean
        True of the assignment magntiude passed the DesignMode
        for the magnitude, False if not.

    complete_check: str
        'COMPLETE' if assign_mag is value, 'INCOMPLETE' if Null.
    """

    # set default values
    complete_check = 'COMPLETE'
    targ_check = False

    # do checks
    # define Null cases for targetdb.magnitude table
    cases = [-999, -9999, 999,
             0.0, np.nan, 99.9, None]
    if assign_mag in cases or np.isnan(assign_mag):
        complete_check = 'INCOMPLETE'
        # set True, no mag is not a fail
        targ_check = True
    # check when greater than and less than
    elif mag_metric_min != -999. and mag_metric_max != -999.:
        if (mag_metric_min < assign_mag < mag_metric_max):
            targ_check = True
    # check when just greater than
    elif mag_metric_min != -999.:
        if assign_mag > mag_metric_min:
            targ_check = True
    # check when less than
    else:
        if assign_mag < mag_metric_max:
            targ_check = True
    return targ_check, complete_check


class DesignModeCheck(DesignMode):
    """
    Parameters
    ----------
    FPSDesign: object
        mugatu.fpsdesign.FPSDesign object with a built design.

    desmode_label: str
        The DesignMode label from targetdb. Options for label are:
        bright_time, dark_plane, dark_monit, dark_rm and dark_faint.

    desmode_manual: dict
        Dictonary of DesignMode parameters to be used as manual
        inputs to validate the design, rather than from targetdb.

    mags: np.array
        Array of magntiudes of shape (N,M), where N is the length
        of mugatu.fpsdesign.FPSDesign.design and M=7. Columns of
        length M in array correspond to magntidues: [g, r, i, bp,
        gaia_g, rp, h].

    Attributes
    ----------
    n_skies_min: dict
        Dictonary with the minimum number of skies for a design
        for each instrument ('APOGEE' and 'BOSS').

    min_skies_fovmetric: dict
        Dictonary wih the FOV metric for the skies in a design
        for each instrument('APOGEE' and 'BOSS').
        The FOV metric is described by three parameters, the
        nth neighbor to get distances to, the percentle distance
        to calculate and the distace to compare against for validation
        (in mm).

    n_stds_min: dict
        Dictonary with the minimum number of standards for a design
        for each instrument ('APOGEE' and 'BOSS').

    min_stds_fovmetric: dict
        Dictonary wih the FOV metric for the standards in a design
        for each instrument ('APOGEE' and 'BOSS').
        The FOV metric is described by three parameters, the
        nth neighbor to get distances to, the percentle distance
        to calculate and the distace to compare against for validation
        (in mm).

    stds_mags: dict
        Dictonary for the min/max magnitude for the standards in a
        design for each instrument ('APOGEE' and 'BOSS'). Indexes
        correspond to magntidues: [g, r, i, bp, gaia_g, rp, h].

    bright_limit_targets: dict
        Dictonary for the min/max magnitude for the science targets
        in adesign for each instrument ('APOGEE' and 'BOSS'). Indexes
        correspond to magntidues: [g, r, i, bp, gaia_g, rp, h].

    sky_neighbors_targets: dict
        Dictonary for the parameters used to check distance between
        skies and all possible sources in field for each instrument
        ('APOGEE' and 'BOSS'). Distances to targets (r, in arcseconds)
        must be r > R_0 * (lim - mags) ** beta, where mags is G band
        for BOSS and H band for APOGEE. Indexes correspond to:
        [R_0, beta, lim]

    trace_diff_targets: dict
        Dictonary for the maximum magnitude difference allowed between
        fibers next to each ther on the chip for each instrument
        ('APOGEE' and 'BOSS'). Here the magntidue difference is checked
        in the G band for BOSS and H band for APOGEE.

    carton_classes: dict
        Dictonary of arrays for the carton pks in each category
        (i.e. science, standards and skies).
    """

    def __init__(self, FPSDesign, desmode_label,
                 desmode_manual=None,
                 mags=None):
        self.design = FPSDesign.design
        self.desmode_label = desmode_label

        # grab the design mode params
        if desmode_manual is None:
            self.fromdb(label=self.desmode_label)
        else:
            self.frommanual(label=self.desmode_label,
                            desmode_manual=desmode_manual)

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
        if mags is None:
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
        else:
            self.mags = mags

    def skies_min(self, instrument):
        """
        Checks if design has the required number of skies
        for some instrument. Returns True if number skies is
        greater than the minimum and False if not.

        Parameters
        ----------
        instrument: str
            Instrument to check number of sky fibers.
            Must be 'BOSS' or 'APOGEE'.

        Returns
        -------
        : boolean
            True if number skies is
            greater than the minimum and False if not.
        """
        n_skies = len(self.design['catalogID'][(self.design['catalogID'] != -1) &
                                                     (np.isin(self.design['carton_pk'],
                                                              self.carton_classes['sky'])) &
                                                     (self.design['obsWavelength'] == instrument)])
        if n_skies >= self.n_skies_min[instrument]:
            return True
        else:
            return False

    def stds_min(self, instrument):
        """
        Checks if design has the required number of standards
        for some instrument. Returns True if number standards is
        greater than the minimum and False if not.

        Parameters
        ----------
        instrument: str
            Instrument to check number of standard fibers.
            Must be 'BOSS' or 'APOGEE'.

        Returns
        -------
        : boolean
            True if number standards is
            greater than the minimum and False if not.
        """
        n_stds = len(self.design['catalogID'][(self.design['catalogID'] != -1) &
                                                    (np.isin(self.design['carton_pk'],
                                                             self.carton_classes['std'])) &
                                                    (self.design['obsWavelength'] == instrument)])
        if n_stds >= self.n_stds_min[instrument]:
            return True
        else:
            return False

    def skies_fov(self, instrument):
        """
        Checks if design meets the FOV metric for the skies
        for some instrument. Returns True FOV metric met,
        False if not. Also, if no science targets in design
        returns True, while no skies returns False.

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
            returns True, while no skies returns False.
        """
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

        # if no skies, dont do check
        if len(x_sky) == 0:
            return False
        # if no science in band, doesnt matter?
        if len(x_sci) == 0:
            return True

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
                                  self.min_skies_fovmetric[instrument][1])
        if perc_dist < self.min_skies_fovmetric[instrument][2]:
            return True
        else:
            return False

    def stds_fov(self, instrument):
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

        # if no stds, dont do check
        if len(x_std) == 0:
            return False
        # if no science in band, doesnt matter?
        if len(x_sci) == 0:
            return True

        # create KDE tree
        tree = cKDTree(np.column_stack((x_std, y_std)))
        # get distances for nearest neighbors
        dd, ii = tree.query(np.column_stack((x_sci, y_sci)),
                            k=self.min_stds_fovmetric[instrument][0])
        # second column is the nth neighbor distance if k>1
        if self.min_stds_fovmetric[instrument][0] == 1:
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

    def mag_limits(self, mag_metric,
                   instrument, carton_class):
        """
        Checks the if magnitude of assignments agree with
        design mode for some instrument and carton class

        Parameters
        ----------
        mag_metric: np.array
            Array of shape (N,M), where N=7 corresponds to
            magntiudes [g, r, i, bp, gaia_g, rp, h], and M=2
            where 0th column is minimum magnitude and 1st column
            is maximum magnitude. If no check in certain band,
            use None as value.

        instrument: str
            Instrument to check magnitudes.
            Must be 'BOSS' or 'APOGEE'.

        carton_class: str
            Carton class to check magnitude limits of.
            Must be 'std' or 'science'.

        Returns
        -------
        mag_checks: np.array
            Array of booleans equal to length of self.design.
            If True, assignment within magnitude limits, False
            if not. If assignment is not with instrument or in
            carton_class, will be False.

        complete_check: np.array
            Array of str equal to length of self.design.
            If 'INCOMPLETE', then assignment either did not have
            all required photometry for check, or assignment is
            not with instrument or in carton_class.
        """
        mag_checks = np.zeros(len(self.design['catalogID']),
                              dtype=bool)
        # if complete, all bands for target present in check
        # if incomplete, then
        complete_check = np.array(['COMPLETE' for _ in range(len(self.design['catalogID']))],
                                  dtype='<U12')

        # check which limits are defined for mode
        check_inds = []
        for i in range(mag_metric.shape[0]):
            if mag_metric[i][0] != -999. or mag_metric[i][1] != -999.:
                check_inds.append(i)

        # run checks
        for i in range(len(mag_checks)):
            if (self.design['catalogID'][i] != -1 and
                self.design['carton_pk'][i] in self.carton_classes[carton_class] and
                self.design['obsWavelength'][i] == instrument):
                # check in each band that has check defined
                targ_check = np.zeros(len(check_inds), dtype=bool)
                for j, ind in enumerate(check_inds):
                    # check the magntiude for this assignment
                    targ_check[j], complete_check[i] = check_assign_mag_limit(
                                                            mag_metric[ind][0],
                                                            mag_metric[ind][1],
                                                            self.mags[i][ind])
                # if all True, then passes
                if np.all(targ_check):
                    mag_checks[i] = True
            else:
                complete_check[i] = 'INCOMPLETE'
        return mag_checks, complete_check

    def sky_neighbors(self, instrument, catalogdb_ver):
        """
        Check if any skies are within some distance
        from another source where distances to targets
        (r, in arcseconds) must be r > R_0 * (lim - mags) ** beta,
        where mags is G band for BOSS and H band for APOGEE.

        Parameters
        ----------
        instrument: str
            Instrument to check sky distances.
            Must be 'BOSS' or 'APOGEE'.

        catalogdb_ver: int
            catalogdb.Version.id to check for targets
            in design field for check.

        Returns
        -------
        sky_checks: np.array
            Array of booleans equal to length of self.design.
            True if r > R_0 * (lim - mags) ** beta, False if not
            or assignment is not sky for specified instrument.
        """
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
        """
        Perform all DesignMode checks for the given design.

        Parameters
        ----------
        verbose: boolean
            True if want print statement summarizing results
            of all checks.

        Attributes
        ----------
        n_skies_min_check: dict
            Results of minumum sky check for each instrument.

        min_skies_fovmetric_check: dict
            Results of sky FOV metric check for each instrument.

        n_stds_min_check: dict
            Results of minumum standard check for each instrument.

        min_stds_fovmetric_check: dict
            Results of standard FOV metric check for each instrument.

        stds_mags_check: dict
            Results of standard magnitude limit check for
            each instrument.

        bright_limit_targets_check: dict
             Results of science magnitude limit check for
            each instrument.
        """
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
