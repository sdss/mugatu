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
from astropy.time import Time

from mugatu.exceptions import MugatuError, MugatuWarning
from coordio.utils import (radec2wokxy, wokxy2radec,
                           offset_definition, Moffat2dInterp)

fmagloss = Moffat2dInterp()

try:
    from sdssdb.peewee.sdss5db import database
    database.set_profile('operations')
    _database = True
except:
    _database = False

from sdssdb.peewee.sdss5db.targetdb import Carton, Category, Magnitude, CartonToTarget, Target
from sdssdb.peewee.sdss5db.targetdb import DesignMode as DesignModeDB
from sdssdb.peewee.sdss5db import catalogdb
from sdssdb.peewee.sdss5db import targetdb


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
    sep = (180 / np.pi) * np.arccos(np.sin(dec1) * np.sin(dec2) +
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
        correspond to magntidues: [g, r, i, z, bp, gaia_g, rp, J, H, K].

    bright_limit_targets: dict
        Dictonary for the min/max magnitude for the science targets
        in adesign for each instrument ('APOGEE' and 'BOSS'). Indexes
        correspond to magntidues: [g, r, i, z, bp, gaia_g, rp, J, H, K].

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


def bright_neigh_exclusion_r(mag_bs, mag_limit_r, lunation):
    """
    returns the exclusion radius for a fiber around a bright star
    in arcseconds for a given designmode. This is for the piecewise
    appromixation used based on:
    https://wiki.sdss.org/pages/viewpage.action?pageId=100173069

    Parameters
    ----------
    mag_bs: float or np.array
        The magniutde in the G band for the bright star(s)

    mag_limit_r: float
        Magnitude limit for the designmode in the r-SDSS band

    lunation: str:
        If the designmode is bright time ('bright') or dark
        time ('dark')

    Returns
    -------
    r_exclude: float or np.array
        exclusion radius in arcseconds around bright star(s)
    """
    flag = ('mugatu.designmode.bright_neigh_exclusion_r will be depriciated '
            'in future releases. For the same functionalitiy, please use '
            'coordio.utils.offset_definition')
    warnings.warn(flag, MugatuWarning)
    # linear portion in the wings
    r_wings = (mag_limit_r - mag_bs - 8.2) / 0.05
    # linear portion in transition area
    r_trans = (mag_limit_r - mag_bs - 4.5) / 0.25
    # core area
    if lunation == 'bright':
        r_core = 1.75 * (mag_limit_r - mag_bs) ** 0.6
    else:
        r_core = 1.5 * (mag_limit_r - mag_bs) ** 0.8
    # exlusion radius is the max of each section
    if isinstance(mag_bs, float):
        if mag_bs <= mag_limit_r:
            r_exclude = max(r_wings, r_trans, r_core)
        else:
            r_exclude = 0.
    else:
        r_exclude = np.nanmax(np.column_stack((r_wings,
                                               r_trans,
                                               r_core)),
                              axis=1)
        r_exclude[mag_bs > mag_limit_r] = 0.
    return r_exclude


def adjusted_brigh_neigh_mag(mag_bs, r, lunation):
    """
    returns the approximate adjusted magntidue of a
    bright source at the position of a nearby fiber.
    This is for the piecewise
    appromixation used based on:
    https://wiki.sdss.org/pages/viewpage.action?pageId=100173069

    Parameters
    ----------
    mag_bs: float
        The magniutde in the G band for the bright star

    r: float or np.array
        distance between fiber and bright star (in arcseconds)

    lunation: str:
        If the designmode is bright time ('bright') or dark
        time ('dark')

    Returns
    -------
    adjusted_mag_bs: float
        The adjusted magntiude of the bright star in SDSS r-band
    """
    mag_wings = 8.2 + 0.05 * r + mag_bs
    mag_trans = 4.5 + 0.25 * r + mag_bs
    if lunation == 'bright':
        mag_core = (r / 1.75) ** (1 / 0.6) + mag_bs
    else:
        mag_core = (r / 1.5) ** (1 / 0.8) + mag_bs
    if isinstance(r, float):
        adjusted_mag_bs = np.nanmin([mag_wings,
                                     mag_trans,
                                     mag_core])
    else:
        adjusted_mag_bs = np.nanmin(np.column_stack((mag_wings,
                                                     mag_trans,
                                                     mag_core)),
                                    axis=1)
    return adjusted_mag_bs


def build_brigh_neigh_query(check_type, instrument, mag_lim,
                            racen, deccen, observatory=None,
                            version_catdb='0.5.0'):
    """
    Builds the database query needed to run bright
    neighbor check

    Parameters
    ----------
    check_type: str
        Either 'designmode' for full designmode check with
        all stars down to mag_lim in catalogdb, or 'safety'
        for safety check using bright stars in targetdb.

    instrument: str
        Fiber instrument being checked. Either 'BOSS' or
        'APOGEE'.

    mag_lim: float
        Magnitude limit in the r_SDSS band for BOSS or
        H band for APOGEE.

    racen: float
        Field center in right ascension.

    deccen: float
        Feild center in declination

    observatory: str
        Observatory where observation is taking place, either
        'LCO' or 'APO'.

    version_catdb: str
        catalogdb.Version.plan to use for the query

    Outputs
    -------
    db_query_results: tuple
        Tuple of (ra, dec, mag, catalogid, pmra, pmdec) for the
        appropriate database query
    """
    # return empty tuple if no mag limit
    if mag_lim == -999.:
        return ()
    # change search radius based on observatory
    if observatory is None:
        r_search = 1.5
    elif observatory == 'APO':
        r_search = 1.5
    else:
        r_search = 1.0
    if check_type == 'designmode':
        if instrument == 'BOSS':
            cat = catalogdb.Gaia_DR2
            ra_col = catalogdb.Catalog.ra
            dec_col = catalogdb.Catalog.dec
            ra_col_str = 'ra'
            dec_col_str = 'dec'
            mag_col = catalogdb.Gaia_DR2.phot_g_mean_mag
            # run the query
            db_query_gaia = (catalogdb.Catalog.select(
                ra_col,
                dec_col,
                mag_col,
                catalogdb.Catalog.catalogid,
                catalogdb.Catalog.pmra,
                catalogdb.Catalog.pmdec)
                .join(catalogdb.Version)
                .switch(catalogdb.Catalog)
                .join(catalogdb.CatalogToTIC_v8)
                .join(catalogdb.TIC_v8)
                .join(catalogdb.Gaia_DR2)
                .where((cat.cone_search(racen,
                                        deccen,
                                        r_search,
                                        ra_col=ra_col_str,
                                        dec_col=dec_col_str)) &
                       (mag_col < mag_lim) &
                       (catalogdb.Version.plan == version_catdb)))
            rasg, decsg, magsg, catalogidsg, pmrasg, pmdecsg = map(list, zip(*list(db_query_gaia.tuples())))
            rasg = np.array(rasg, dtype=np.float64)
            decsg = np.array(decsg, dtype=np.float64)
            magsg = np.array(magsg, dtype=np.float64)
            catalogidsg = np.array(catalogidsg,
                                   dtype=int)
            pmrasg = np.array(pmrasg, dtype=np.float64)
            pmdecsg = np.array(pmdecsg, dtype=np.float64)

            cat = catalogdb.Tycho2
            ra_col = catalogdb.Catalog.ra
            dec_col = catalogdb.Catalog.dec
            ra_col_str = 'radeg'
            dec_col_str = 'dedeg'
            mag_colbt = catalogdb.Tycho2.btmag
            mag_colvt = catalogdb.Tycho2.vtmag
            # run the query
            db_query_tych = (catalogdb.Catalog.select(
                ra_col,
                dec_col,
                mag_colbt,
                mag_colvt,
                catalogdb.Catalog.catalogid,
                catalogdb.Catalog.pmra,
                catalogdb.Catalog.pmdec)
                .join(catalogdb.Version)
                .switch(catalogdb.Catalog)
                .join(catalogdb.CatalogToTycho2)
                .join(catalogdb.Tycho2)
                .where((cat.cone_search(racen,
                                        deccen,
                                        r_search,
                                        ra_col=ra_col_str,
                                        dec_col=dec_col_str)) &
                       (mag_colvt < mag_lim) &
                       (catalogdb.Version.plan == version_catdb)))
            rast, decst, magsbt, magsvt, catalogidst, pmrast, pmdecst = map(list, zip(*list(db_query_tych.tuples())))
            rast = np.array(rast, dtype=np.float64)
            decst = np.array(decst, dtype=np.float64)
            magsbt = np.array(magsbt, dtype=np.float64)
            magsvt = np.array(magsvt, dtype=np.float64)
            catalogidst = np.array(catalogidst,
                                   dtype=int)
            pmrast = np.array(pmrast, dtype=np.float64)
            pmdecst = np.array(pmdecst, dtype=np.float64)
            magsg_tych = np.zeros(len(magsbt))
            magsg_tych[magsbt != None] = (magsvt[magsbt != None] - 0.02051 -
                                         0.2706 * (magsbt[magsbt != None] -
                                                   magsvt[magsbt != None]) +
                                         0.03394 * (magsbt[magsbt != None] -
                                                    magsvt[magsbt != None]) ** 2 -
                                         0.05937 * (magsbt[magsbt != None] -
                                                    magsvt[magsbt != None]) ** 3)
            magsg_tych[magsbt == None] = magsvt[magsbt == None] - 1
            # set up the logic to only include tycho
            # stars not found in Gaia
            tych_in_gaia = np.isin(catalogidst, catalogidsg)
            # evals for each now
            tych_eval = (~tych_in_gaia)
            ras = np.append(rasg,
                            rast[tych_eval])
            decs = np.append(decsg,
                             decst[tych_eval])
            mags = np.append(magsg,
                             magsg_tych[tych_eval])
            catalogids = np.append(catalogidsg,
                                   catalogidst[tych_eval])
            pmras = np.append(pmrasg,
                              pmrast[tych_eval])
            pmdecs = np.append(pmdecsg,
                               pmdecst[tych_eval])
            db_query_results = (ras, decs, mags, catalogids, pmras, pmdecs)

        else:
            cat = catalogdb.TwoMassPSC
            ra_col = catalogdb.Catalog.ra
            dec_col = catalogdb.Catalog.dec
            ra_col_str = 'ra'
            dec_col_str = 'decl'
            mag_col = catalogdb.TwoMassPSC.h_m
            # run the query
            db_query = (catalogdb.Catalog.select(
                ra_col,
                dec_col,
                mag_col,
                catalogdb.Catalog.catalogid,
                catalogdb.Catalog.pmra,
                catalogdb.Catalog.pmdec)
                .join(catalogdb.Version)
                .switch(catalogdb.Catalog)
                .join(catalogdb.CatalogToTIC_v8)
                .join(catalogdb.TIC_v8)
                .join(cat)
                .where((cat.cone_search(racen,
                                        deccen,
                                        r_search,
                                        ra_col=ra_col_str,
                                        dec_col=dec_col_str)) &
                       (mag_col < mag_lim) &
                       (catalogdb.Version.plan == version_catdb)))
            ras, decs, mags, catalogids, pmras, pmdecs = map(list, zip(*list(db_query.tuples())))
            ras = np.array(ras, dtype=np.float64)
            decs = np.array(decs, dtype=np.float64)
            mags = np.array(mags, dtype=np.float64)
            catalogids = np.array(catalogids,
                                  dtype=int)
            pmras = np.array(pmras, dtype=np.float64)
            pmdecs = np.array(pmdecs, dtype=np.float64)
            db_query_results = (ras, decs, mags, catalogids, pmras, pmdecs)
    else:
        if instrument == 'BOSS':
            carts = ['ops_tycho2_brightneighbors',
                     'ops_gaia_brightneighbors']
            mag_col = targetdb.Magnitude.gaia_g
            # tycho query
            db_queryt = (targetdb.CartonToTarget.select(targetdb.Target.ra,
                                                        targetdb.Target.dec,
                                                        mag_col,
                                                        targetdb.Target.catalogid,
                                                        targetdb.Target.pmra,
                                                        targetdb.Target.pmdec)
                                                .join(targetdb.Target)
                                                .switch(targetdb.CartonToTarget)
                                                .join(targetdb.Magnitude)
                                                .switch(targetdb.CartonToTarget)
                                                .join(targetdb.Carton)
                                                .join(targetdb.Version)
                                                .where((targetdb.Target.cone_search(racen,
                                                                                    deccen,
                                                                                    r_search)) &
                                                       (targetdb.Carton.carton == carts[0]) &
                                                       (targetdb.Version.plan >= version_catdb) &
                                                       (~(targetdb.Version.plan % '%-test'))))
            rast, decst, magst, catalogidst, pmrast, pmdecst = map(list, zip(*list(db_queryt.tuples())))
            rast = np.array(rast, dtype=np.float64)
            decst = np.array(decst, dtype=np.float64)
            magst = np.array(magst, dtype=np.float64)
            catalogidst = np.array(catalogidst,
                                  dtype=int)
            pmrast = np.array(pmrast, dtype=np.float64)
            pmdecst = np.array(pmdecst, dtype=np.float64)

            # gaia query
            db_queryg = (targetdb.CartonToTarget.select(targetdb.Target.ra,
                                                        targetdb.Target.dec,
                                                        mag_col,
                                                        targetdb.Target.catalogid,
                                                        targetdb.Target.pmra,
                                                        targetdb.Target.pmdec)
                                                .join(targetdb.Target)
                                                .switch(targetdb.CartonToTarget)
                                                .join(targetdb.Magnitude)
                                                .switch(targetdb.CartonToTarget)
                                                .join(targetdb.Carton)
                                                .join(targetdb.Version)
                                                .where((targetdb.Target.cone_search(racen,
                                                                                    deccen,
                                                                                    r_search)) &
                                                       (targetdb.Carton.carton == carts[1]) &
                                                       (targetdb.Version.plan >= version_catdb) &
                                                       (~(targetdb.Version.plan % '%-test'))))
            rasg, decsg, magsg, catalogidsg, pmrasg, pmdecsg = map(list, zip(*list(db_queryg.tuples())))
            rasg = np.array(rasg, dtype=np.float64)
            decsg = np.array(decsg, dtype=np.float64)
            magsg = np.array(magsg, dtype=np.float64)
            catalogidsg = np.array(catalogidsg,
                                  dtype=int)
            pmrasg = np.array(pmrasg, dtype=np.float64)
            pmdecsg = np.array(pmdecsg, dtype=np.float64)

            # remove tycho objects in gaia
            tych_in_gaia = np.isin(catalogidst, catalogidsg)
            tych_eval = (~tych_in_gaia)
            ras = np.append(rasg,
                            rast[tych_eval])
            decs = np.append(decsg,
                             decst[tych_eval])
            mags = np.append(magsg,
                             magst[tych_eval])
            catalogids = np.append(catalogidsg,
                                   catalogidst[tych_eval])
            pmras = np.append(pmrasg,
                              pmrast[tych_eval])
            pmdecs = np.append(pmdecsg,
                               pmdecst[tych_eval])
        else:
            carts = ['ops_2mass_psc_brightneighbors']
            mag_col = targetdb.Magnitude.h
            # run the query
            db_query = (targetdb.CartonToTarget.select(targetdb.Target.ra,
                                                       targetdb.Target.dec,
                                                       mag_col,
                                                       targetdb.Target.catalogid,
                                                       targetdb.Target.pmra,
                                                       targetdb.Target.pmdec)
                                               .join(targetdb.Target)
                                               .switch(targetdb.CartonToTarget)
                                               .join(targetdb.Magnitude)
                                               .switch(targetdb.CartonToTarget)
                                               .join(targetdb.Carton)
                                               .join(targetdb.Version)
                                               .where((targetdb.Target.cone_search(racen,
                                                                                   deccen,
                                                                                   r_search)) &
                                                      (targetdb.Carton.carton == carts[0]) &
                                                      (targetdb.Version.plan >= version_catdb) &
                                                      (~(targetdb.Version.plan % '%-test'))))
            ras, decs, mags, catalogids, pmras, pmdecs = map(list, zip(*list(db_query.tuples())))
            ras = np.array(ras, dtype=np.float64)
            decs = np.array(decs, dtype=np.float64)
            mags = np.array(mags, dtype=np.float64)
            catalogids = np.array(catalogids,
                                  dtype=int)
            pmras = np.array(pmras, dtype=np.float64)
            pmdecs = np.array(pmdecs, dtype=np.float64)
        db_query_results = (ras, decs, mags, catalogids, pmras, pmdecs)
    return db_query_results


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

    db_query_results_boss: dict
        Database query results for BOSS bright neighbor check.
        Each index of dict is a tuple of (ras, decs, mags, catalogids)
        with one index for designmode and the other safety.

    db_query_results_apogee: dict
        Database query results for APOGEE bright neighbor check.
        Each index of dict is a tuple of (ras, decs, mags, catalogids)
        with one index for designmode and the other safety.

    Attributes
    ----------
    design: dict
        FPSDesign dictonary.

    racen: float
        RA center of the field (degrees)

    deccen: float
        DEC center of the field (degrees)

    position_angle: float
        Position angle of the field E of N in degrees.

    observatory: str
        Observatory where observation is taking place, either
        'LCO' or 'APO'.

    obsTime: float
        Julian date of the observation.

    rg: kaiju.robotGrid
        Kaiju robotGrid object for the design.

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
                 db_query_results_boss=None,
                 db_query_results_apogee=None):
        # grab needed info from FPSDesign object
        self.design = FPSDesign.design
        self.racen = FPSDesign.racen
        self.deccen = FPSDesign.deccen
        self.position_angle = FPSDesign.position_angle
        self.observatory = FPSDesign.observatory
        self.obsTime = FPSDesign.obsTime
        self.rg = FPSDesign.rg
        self.desmode_label = desmode_label
        self.db_query_results_boss = db_query_results_boss
        self.db_query_results_apogee = db_query_results_apogee

        # grab the design mode params
        if desmode_manual is None:
            self.fromdb(label=self.desmode_label)
        else:
            self.fromdict(designmode_dict=desmode_manual)

        # classify cartons as skies, standards or science
        self.carton_classes = {}
        self.carton_classes['science'] = []
        self.carton_classes['sky'] = []
        self.carton_classes['std'] = []
        for cat in np.unique(self.design['category'][self.design['catalogID'] != -1]):
            if 'sky' in cat:
                self.carton_classes['sky'] += list(self.design['category'][(self.design['catalogID'] != -1) &
                                                                            (self.design['category'] == cat)])
            elif 'standard' in cat:
                self.carton_classes['std'] += list(self.design['category'][(self.design['catalogID'] != -1) &
                                                                            (self.design['category'] == cat)])
            # I think else makes sense here as there are science
            # and open fiber labels?
            else:
                self.carton_classes['science'] += list(self.design['category'][(self.design['catalogID'] != -1) &
                                                                                (self.design['category'] == cat)])

        # collect magntiudes of design
        # here I am doing g,r,i,BP,G,RP,H
        self.mags = self.design['magnitudes']

    def skies_min(self, instrument, return_metric=False):
        """
        Checks if design has the required number of skies
        for some instrument. Returns True if number skies is
        greater than the minimum and False if not.

        Parameters
        ----------
        instrument: str
            Instrument to check number of sky fibers.
            Must be 'BOSS' or 'APOGEE'.

        return_metric: boolean
            If True, will return the number of skies
            in the design.

        Returns
        -------
        : boolean
            True if number skies is
            greater than the minimum and False if not.

        n_skies: int
            Number of skies in design, only returned if
            return_metric=True.
        """
        n_skies = len(self.design['catalogID'][(self.design['catalogID'] != -1) &
                                                     (np.isin(self.design['category'],
                                                              self.carton_classes['sky'])) &
                                                     (self.design['obsWavelength'] == instrument)])
        if return_metric:
            if n_skies >= self.n_skies_min[instrument]:
                return True, n_skies
            else:
                return False, n_skies
        else:
            if n_skies >= self.n_skies_min[instrument]:
                return True
            else:
                return False

    def stds_min(self, instrument, return_metric=False):
        """
        Checks if design has the required number of standards
        for some instrument. Returns True if number standards is
        greater than the minimum and False if not.

        Parameters
        ----------
        instrument: str
            Instrument to check number of standard fibers.
            Must be 'BOSS' or 'APOGEE'.

        return_metric: boolean
            If True, will return the number of standards
            in the design.

        Returns
        -------
        : boolean
            True if number standards is
            greater than the minimum and False if not.

        n_stds: int
            Number of standards in design, only returned if
            return_metric=True.
        """
        n_stds = len(self.design['catalogID'][(self.design['catalogID'] != -1) &
                                                    (np.isin(self.design['category'],
                                                             self.carton_classes['std'])) &
                                                    (self.design['obsWavelength'] == instrument)])
        if return_metric:
            if n_stds >= self.n_stds_min[instrument]:
                return True, n_stds
            else:
                return False, n_stds
        else:
            if n_stds >= self.n_stds_min[instrument]:
                return True
            else:
                return False

    def skies_fov(self, instrument, return_metric=False):
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

        return_metric: boolean
            If True, will return the FOV metric
            for the design.

        Returns
        -------
        : boolean
            True FOV metric met,
            False if not. Also, if no science targets in design
            returns True, while no skies returns False.

        perc_dist: float
            Percentile distance between science assignments and
            skies (FOV metric) for the design. Only returned if
            return_metric=True.
        """
        # for now, return True (passing) if no FOV metric supplied
        if np.all(self.min_skies_fovmetric[instrument] == -999.):
            if return_metric:
                return True, -1.
            else:
                return True
        # get x,y of the skies
        x_sky = self.design['x'][(self.design['catalogID'] != -1) &
                                 (np.isin(self.design['category'],
                                          self.carton_classes['sky'])) &
                                 (self.design['obsWavelength'] == instrument)]
        y_sky = self.design['y'][(self.design['catalogID'] != -1) &
                                 (np.isin(self.design['category'],
                                          self.carton_classes['sky'])) &
                                 (self.design['obsWavelength'] == instrument)]

        x_sci = self.design['x'][(self.design['catalogID'] != -1) &
                                 (~np.isin(self.design['category'],
                                           self.carton_classes['sky'])) &
                                 (self.design['obsWavelength'] == instrument)]
        y_sci = self.design['y'][(self.design['catalogID'] != -1) &
                                 (~np.isin(self.design['category'],
                                           self.carton_classes['sky'])) &
                                 (self.design['obsWavelength'] == instrument)]

        # if no science in band, doesnt matter?
        if len(x_sci) == 0:
            if return_metric:
                return True, -1.
            else:
                return True
        # if no skies required, dont do check
        if self.n_skies_min[instrument] == 0:
            if return_metric:
                return True, -1.
            else:
                return True
        # if no skies, dont do check
        if len(x_sky) == 0:
            if return_metric:
                return False, -1.
            else:
                return False

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
        if return_metric:
            if perc_dist < self.min_skies_fovmetric[instrument][2]:
                return True, perc_dist
            else:
                return False, perc_dist
        else:
            if perc_dist < self.min_skies_fovmetric[instrument][2]:
                return True
            else:
                return False

    def stds_fov(self, instrument, return_metric=False):
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

        return_metric: boolean
            If True, will return the FOV metric
            for the design.

        Returns
        -------
        : boolean
            True FOV metric met,
            False if not. Also, if no science targets in design
            returns True, while no standards returns False.

        perc_dist: float
            Percentile distance between science assignments and
            standards (FOV metric) for the design. Only returned if
            return_metric=True.
        """
        # for now, return True (passing) if no FOV metric supplied
        if np.all(self.min_stds_fovmetric[instrument] == -999.):
            if return_metric:
                return True, -1.
            else:
                return True
        # get x,y of the standards
        x_std = self.design['x'][(self.design['catalogID'] != -1) &
                                 (np.isin(self.design['category'],
                                          self.carton_classes['std'])) &
                                 (self.design['obsWavelength'] == instrument)]
        y_std = self.design['y'][(self.design['catalogID'] != -1) &
                                 (np.isin(self.design['category'],
                                          self.carton_classes['std'])) &
                                 (self.design['obsWavelength'] == instrument)]

        x_sci = self.design['x'][(self.design['catalogID'] != -1) &
                                 (~np.isin(self.design['category'],
                                           self.carton_classes['std'])) &
                                 (self.design['obsWavelength'] == instrument)]
        y_sci = self.design['y'][(self.design['catalogID'] != -1) &
                                 (~np.isin(self.design['category'],
                                           self.carton_classes['std'])) &
                                 (self.design['obsWavelength'] == instrument)]

        # if no science in band, doesnt matter?
        if len(x_sci) == 0:
            if return_metric:
                return True, -1.
            else:
                return True
        # if no stds required, dont do check
        if self.n_stds_min[instrument] == 0:
            if return_metric:
                return True, -1.
            else:
                return True
        # if no stds, dont do check
        if len(x_std) == 0:
            if return_metric:
                return False, -1.
            else:
                return False

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
        if return_metric:
            if perc_dist < self.min_stds_fovmetric[instrument][2]:
                return True, perc_dist
            else:
                return False, perc_dist
        else:
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
            Array of shape (N,M), where N=10 corresponds to
            magntiudes [g, r, i, z, bp, gaia_g, rp, J, H, K], and M=2
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
                self.design['category'][i] in self.carton_classes[carton_class] and
                self.design['obsWavelength'][i] == instrument):
                # don't do check and make true if offset target
                if self.design['offset'][i] and self.design['offset_flag'][i] == 0:
                    mag_checks[i] = True
                else:
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

    def bright_neighbors(self, instrument, check_type='designmode',
                         db_query=None):
        """
        Check if any fibers are placed near a bright star

        Parameters
        ----------
        instrument: str
            Instrument to check sky distances.
            Must be 'BOSS' or 'APOGEE'.

        check_type: str
            Either 'designmode' for checking bright
            stars in catalogdb down to designmode magnitude
            limit, or 'safety' for checking bright stars in
            targetdb in bright star cartons.

        db_query: peewee query
            Optionally pass the pre-queried bright stars
            from catalogdb or targetdb. Must include ra, dec,
            magnitude and id. Other option is to provide tuple
            of (ra,dec,magntiude,id).

        Returns
        -------
        neigh_checks_des: np.array
            Array of booleans equal to length of 500, where order
            is robotID 1 to 500 robotID.
            True assignment valid (i.e. not near bright star),
            False if not valid.

        hasFiber: np.array
            Array of booleans equal to length of 500, where order
            is robotID 1 to 500 robotID. True if robotID has intrument
            on fiber, False if not. If False, do not consider result from
            neigh_checks_des for this robotID

        mag_adj_robo: np.array
            If neigh_checks_des == False, reports the magntiude of nearby
            bright source at the fiber position. If neigh_checks_des == True,
            then no value reported (where NULL = -9999.).

        isassigned: np.array
            Array of booleans equal to length of 500, where order
            is robotID 1 to 500 robotID. True if robotID is assigned in grid
            and False if it is not.
        """
        # get xPos and yPos from robotGrid
        xrobo = np.zeros(500)
        yrobo = np.zeros(500)
        hasFiber = np.zeros(500, dtype=bool) + True
        isassigned = np.zeros(500, dtype=bool) + True
        mag_adj_robo = np.zeros(500) - 9999.
        for i, robotID in enumerate(self.rg.robotDict):
            if instrument == 'BOSS':
                xrobo[i] = self.rg.robotDict[robotID].bossWokXYZ[0]
                yrobo[i] = self.rg.robotDict[robotID].bossWokXYZ[1]
            else:
                hasFiber[i] = self.rg.robotDict[robotID].hasApogee
                xrobo[i] = self.rg.robotDict[robotID].apWokXYZ[0]
                yrobo[i] = self.rg.robotDict[robotID].apWokXYZ[1]
            isassigned[i] = self.rg.robotDict[robotID].isAssigned()

        ra_robo, dec_robo, fieldWarn = wokxy2radec(xWok=xrobo,
                                                   yWok=yrobo,
                                                   waveName=instrument.title(),
                                                   raCen=self.racen,
                                                   decCen=self.deccen,
                                                   obsAngle=self.position_angle,
                                                   obsSite=self.observatory,
                                                   obsTime=self.obsTime)

        neigh_checks = np.zeros(500, dtype=bool) + True
        mag_limits = self.bright_limit_targets[instrument][:, 0]
        # set magntiude limit for instrument and lunation
        if instrument == 'Apogee':
            # 2MASS H
            mag_lim = mag_limits[8]
        elif 'bright' in self.desmode_label:
            # Gaia G
            mag_lim = mag_limits[5]
        else:
            # SDSS r
            mag_lim = mag_limits[1]
        # run query for field if not supplied
        if (self.db_query_results_boss is not None and
           instrument == 'BOSS'):
            db_query = self.db_query_results_boss[check_type]
        elif (self.db_query_results_apogee is not None and
           instrument == 'APOGEE'):
            db_query = self.db_query_results_apogee[check_type]
        elif db_query is None:
            db_query = build_brigh_neigh_query(check_type, instrument,
                                               mag_lim, self.racen,
                                               self.deccen, self.observatory)

        # only do check if any stars returned
        if len(db_query) > 0:
            if isinstance(db_query, tuple):
                ras, decs, mags, catalogids, pmras, pmdecs = db_query
            else:
                ras, decs, mags, catalogids, pmras, pmdecs = map(list, zip(*list(db_query.tuples())))

            # set nan pms to 0
            pmras[np.isnan(pmras)] = 0.
            pmdecs[np.isnan(pmdecs)] = 0.
            # convert to x,y
            radVel = (np.zeros(len(ras),
                               dtype=np.float64) + 1.e-4)
            parallax = (np.zeros(len(ras),
                                 dtype=np.float64) + 1.e-4)
            res = radec2wokxy(
                ra=ras,
                dec=decs,
                coordEpoch=Time(np.array([2015.5] * len(ras)),
                                format='decimalyear').jd,
                waveName=instrument.title(),
                raCen=self.racen,
                decCen=self.deccen,
                obsAngle=self.position_angle,
                obsSite=self.observatory,
                obsTime=self.obsTime,
                pmra=pmras,
                pmdec=pmdecs,
                parallax=parallax,
                radVel=radVel)
            # now convert back to ra,dec
            ras, decs, fieldWarn = wokxy2radec(xWok=res[0],
                                               yWok=res[1],
                                               waveName=instrument.title(),
                                               raCen=self.racen,
                                               decCen=self.deccen,
                                               obsAngle=self.position_angle,
                                               obsSite=self.observatory,
                                               obsTime=self.obsTime)


            if 'bright' in self.desmode_label:
                r_exclude, _ = offset_definition(mags,
                                                 mag_limits,
                                                 lunation='bright',
                                                 waveName=instrument.title(),
                                                 fmagloss=fmagloss)
            else:
                r_exclude, _ = offset_definition(mags,
                                                 mag_limits,
                                                 lunation='dark',
                                                 waveName=instrument.title(),
                                                 fmagloss=fmagloss)

            # check if fibers too close to bright neighbors
            for i in range(len(r_exclude)):
                # only do check if exclusion radius larger
                # than 0" (otherwise star not bright enough)
                if r_exclude[i] > 0.:
                    dist = ang_sep(ras[i], decs[i],
                                   ra_robo, dec_robo) * 3600.
                    if 'bright' in self.desmode_label:
                        mag_adj = adjusted_brigh_neigh_mag(mags[i],
                                                           dist,
                                                           lunation='bright')
                    else:
                        mag_adj = adjusted_brigh_neigh_mag(mags[i],
                                                           dist,
                                                           lunation='dark')
                    neigh_checks[dist < r_exclude[i]] = False
                    mag_adj_robo[dist < r_exclude[i]] = mag_adj[dist < r_exclude[i]]
        # dont accoutn for places where no fiber
        mag_adj_robo[~hasFiber] = -9999.

        return neigh_checks, hasFiber, mag_adj_robo, isassigned

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
        result = self.skies_min(instrument='BOSS',
                                return_metric=True)
        self.n_skies_min_check['BOSS'] = result[0]
        self.n_skies_min_check['BOSS_metric'] = result[1]
        result = self.skies_min(instrument='APOGEE',
                                return_metric=True)
        self.n_skies_min_check['APOGEE'] = result[0]
        self.n_skies_min_check['APOGEE_metric'] = result[1]

        self.min_skies_fovmetric_check = {}
        result = self.skies_fov(instrument='BOSS',
                                return_metric=True)
        self.min_skies_fovmetric_check['BOSS'] = result[0]
        self.min_skies_fovmetric_check['BOSS_metric'] = result[1]
        result = self.skies_fov(instrument='APOGEE',
                                return_metric=True)
        self.min_skies_fovmetric_check['APOGEE'] = result[0]
        self.min_skies_fovmetric_check['APOGEE_metric'] = result[1]

        self.n_stds_min_check = {}
        result = self.stds_min(instrument='BOSS',
                               return_metric=True)
        self.n_stds_min_check['BOSS'] = result[0]
        self.n_stds_min_check['BOSS_metric'] = result[1]
        result = self.stds_min(instrument='APOGEE',
                               return_metric=True)
        self.n_stds_min_check['APOGEE'] = result[0]
        self.n_stds_min_check['APOGEE_metric'] = result[1]

        self.min_stds_fovmetric_check = {}
        result = self.stds_fov(instrument='BOSS',
                               return_metric=True)
        self.min_stds_fovmetric_check['BOSS'] = result[0]
        self.min_stds_fovmetric_check['BOSS_metric'] = result[1]
        result = self.stds_fov(instrument='APOGEE',
                               return_metric=True)
        self.min_stds_fovmetric_check['APOGEE'] = result[0]
        self.min_stds_fovmetric_check['APOGEE_metric'] = result[1]

        self.stds_mags_check = {}
        self.stds_mags_check['BOSS'] = self.mag_limits(self.stds_mags['BOSS'],
                                                       'BOSS',
                                                       'std')
        check_tot = len(self.stds_mags_check['BOSS'][0][self.stds_mags_check['BOSS'][0]])
        design_tot = len(self.design['x'][(self.design['catalogID'] != -1) &
                                          (np.isin(self.design['category'],
                                                   self.carton_classes['std'])) &
                                          (self.design['obsWavelength'] == 'BOSS')])
        self.stds_mags_check['BOSS_metric'] = [check_tot, design_tot]
        self.stds_mags_check['APOGEE'] = self.mag_limits(self.stds_mags['APOGEE'],
                                                         'APOGEE',
                                                         'std')
        check_tot = len(self.stds_mags_check['APOGEE'][0][self.stds_mags_check['APOGEE'][0]])
        design_tot = len(self.design['x'][(self.design['catalogID'] != -1) &
                                          (np.isin(self.design['category'],
                                                   self.carton_classes['std'])) &
                                          (self.design['obsWavelength'] == 'APOGEE')])
        self.stds_mags_check['APOGEE_metric'] = [check_tot, design_tot]

        self.bright_limit_targets_check = {}
        self.bright_limit_targets_check['BOSS'] = self.mag_limits(self.bright_limit_targets['BOSS'],
                                                                  'BOSS',
                                                                  'science')
        check_tot = len(self.bright_limit_targets_check['BOSS'][0][self.bright_limit_targets_check['BOSS'][0]])
        design_tot = len(self.design['x'][(self.design['catalogID'] != -1) &
                                          (np.isin(self.design['category'],
                                                   self.carton_classes['science'])) &
                                          (self.design['obsWavelength'] == 'BOSS')])
        self.bright_limit_targets_check['BOSS_metric'] = [check_tot, design_tot]
        self.bright_limit_targets_check['APOGEE'] = self.mag_limits(self.bright_limit_targets['APOGEE'],
                                                                    'APOGEE',
                                                                    'science')
        check_tot = len(self.bright_limit_targets_check['APOGEE'][0][self.bright_limit_targets_check['APOGEE'][0]])
        design_tot = len(self.design['x'][(self.design['catalogID'] != -1) &
                                          (np.isin(self.design['category'],
                                                   self.carton_classes['science'])) &
                                          (self.design['obsWavelength'] == 'APOGEE')])
        self.bright_limit_targets_check['APOGEE_metric'] = [check_tot, design_tot]

        self.bright_neighbor_check = {}
        self.bright_neighbor_check['BOSS'] = self.bright_neighbors(instrument='BOSS',
                                                                   check_type='designmode')
        check_tot = len(self.bright_neighbor_check['BOSS'][0][self.bright_neighbor_check['BOSS'][0] &
                                                              self.bright_neighbor_check['BOSS'][1] &
                                                              self.bright_neighbor_check['BOSS'][3]])
        design_tot = len(self.bright_neighbor_check['BOSS'][0][self.bright_neighbor_check['BOSS'][1] &
                                                               self.bright_neighbor_check['BOSS'][3]])
        mag_adj_near_bs = self.bright_neighbor_check['BOSS'][2]
        self.bright_neighbor_check['BOSS_metric'] = [check_tot, design_tot, mag_adj_near_bs]
        self.bright_neighbor_check['APOGEE'] = self.bright_neighbors(instrument='APOGEE',
                                                                     check_type='designmode')
        check_tot = len(self.bright_neighbor_check['APOGEE'][0][self.bright_neighbor_check['APOGEE'][0] &
                                                                self.bright_neighbor_check['APOGEE'][1] &
                                                                self.bright_neighbor_check['APOGEE'][3]])
        design_tot = len(self.bright_neighbor_check['APOGEE'][0][self.bright_neighbor_check['APOGEE'][1] &
                                                                 self.bright_neighbor_check['APOGEE'][3]])
        mag_adj_near_bs = self.bright_neighbor_check['APOGEE'][2]
        self.bright_neighbor_check['APOGEE_metric'] = [check_tot, design_tot, mag_adj_near_bs]

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

            check_tot = self.stds_mags_check['BOSS_metric'][0]
            design_tot = self.stds_mags_check['BOSS_metric'][1]
            verbose_output += 'Magnitude Limit Stds (BOSS):      | %d out of %d\n' % (check_tot, design_tot)

            check_tot = self.stds_mags_check['APOGEE_metric'][0]
            design_tot = self.stds_mags_check['APOGEE_metric'][1]
            verbose_output += 'Magnitude Limit Stds (APOGEE):    | %d out of %d\n' % (check_tot, design_tot)

            check_tot = self.bright_limit_targets_check['BOSS_metric'][0]
            design_tot = self.bright_limit_targets_check['BOSS_metric'][1]
            verbose_output += 'Magnitude Limit Targets (BOSS):   | %d out of %d\n' % (check_tot, design_tot)

            check_tot = self.bright_limit_targets_check['APOGEE_metric'][0]
            design_tot = self.bright_limit_targets_check['APOGEE_metric'][1]
            verbose_output += 'Magnitude Limit Targets (APOGEE): | %d out of %d\n' % (check_tot, design_tot)

            check_tot = self.bright_neighbor_check['BOSS_metric'][0]
            design_tot = self.bright_neighbor_check['BOSS_metric'][1]
            verbose_output += 'Bright Neighbor Check (BOSS):     | %d out of %d\n' % (check_tot, design_tot)

            check_tot = self.bright_neighbor_check['APOGEE_metric'][0]
            design_tot = self.bright_neighbor_check['APOGEE_metric'][1]
            verbose_output += 'Bright Neighbor Check (APOGEE):   | %d out of %d\n' % (check_tot, design_tot)

            print(verbose_output)
