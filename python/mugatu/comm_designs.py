# @Author: Ilija Medan
# @Date: February 4, 2021
# @Filename: comm_designs.py
# @License: BSD 3-Clause
# @Copyright: Ilija Medan

import numpy as np
import warnings

from sdssdb.peewee.sdss5db.targetdb import Carton, Target, CartonToTarget
from mugatu.exceptions import MugatuError, MugatuWarning
from mugatu.fpsdesign import FPSDesign
from coordio.utils import radec2wokxy
import robostrategy.field as field
import coordio
import datetime


def all_sky_design_RS(racen, deccen, position_angle, observatory,
                      obsTime, n_sky_apogee, n_sky_boss, cadence):
    """
    creates a design that consists of all sky targets

    Parameters
    ----------

    racen: np.float64
        Right Ascension of center of the field

    deccen: np.float64
        Declination of center of the field

    position_angle: np.float64
        Position angle of the field

    observatory: str
        Observatory where observation is taking place

    obsTime: np.float64
        Julian date of the observation

    n_sky_apogee: int
        number of fibers that should be apogee skies

    n_sky_boss: int
        number of fibers that should be boss skies

    cadence: str
        Name of the cadence for the design

    Returns
    -------

    f: robostrategy.Field()
        robostrategy.Field object with the assignmends
        for the design

    fps_design: mugatu.FPSDesign()
        mugatu.FPSDesign() object created using the assingments
        from robostrategy
    """

    # one not for the below is it assumes the field has
    # apogee and boss. There are some fields with no
    # boss skies, will these ever be targeted?

    # make sure n skies is = 500
    if n_sky_apogee + n_sky_boss != 500:
        raise MugatuError(message='n_sky_apogee + n_sky_boss != 500')

    # base search radius on field size
    if observatory == 'APO':
    	r_search = 1.49
    else:
    	r_serach = 0.95

    # grab the skies from targetdb
    apogee_sky = (Target.select(Target.catalogid,
                                Target.ra,
                                Target.dec,
                                Target.pk)
                        .join(CartonToTarget)
                        .join(Carton)
                        .where((Carton.carton == 'ops_sky_apogee') & 
                               (Target.cone_search(racen, deccen, r_search))))

    if len(apogee_sky) < n_sky_apogee:
        raise MugatuError(message='Not enough apogee skies in field')

    boss_sky = (Target.select(Target.catalogid,
                              Target.ra,
                              Target.dec,
                              Target.pk)
                      .join(CartonToTarget)
                      .join(Carton)
                      .where((Carton.carton == 'ops_sky_boss') & 
                             (Target.cone_search(racen, deccen, r_search))))

    if len(boss_sky) < n_sky_boss:
        raise MugatuError(message='Not enough boss skies in field')

    catalogid_apogee, ra_apogee, dec_apogee, target_pk_apogee = map(list, zip(*list(apogee_sky.tuples())))
    catalogid_apogee = np.array(catalogid_apogee, dtype=np.int64)
    ra_apogee = np.array(ra_apogee)
    dec_apogee = np.array(dec_apogee)
    target_pk_apogee = np.array(target_pk_apogee, dtype=np.int64)
    x_apogee, y_apogee, fieldWarn, HA, PA_coordio = radec2wokxy(ra=ra_apogee,
                                                                dec=dec_apogee,
                                                                coordEpoch=np.array([2457174] * len(ra_apogee)),
                                                                waveName=np.array(list(map(lambda x:x.title(), ['APOGEE'] * len(ra_apogee)))),
                                                                raCen=racen,
                                                                decCen=deccen,
                                                                obsAngle=position_angle,
                                                                obsSite=observatory,
                                                                obsTime=obsTime)


    catalogid_boss, ra_boss, dec_boss, target_pk_boss = map(list, zip(*list(boss_sky.tuples())))
    catalogid_boss = np.array(catalogid_boss, dtype=np.int64)
    ra_boss = np.array(ra_boss)
    dec_boss = np.array(dec_boss)
    target_pk_boss = np.array(target_pk_boss, dtype=np.int64)
    x_boss, y_boss, fieldWarn, HA, PA_coordio = radec2wokxy(ra=ra_boss,
                                                            dec=dec_boss,
                                                            coordEpoch=np.array([2457174] * len(ra_boss)),
                                                            waveName=np.array(list(map(lambda x:x.title(), ['BOSS'] * len(ra_boss)))),
                                                            raCen=racen,
                                                            decCen=deccen,
                                                            obsAngle=position_angle,
                                                            obsSite=observatory,
                                                            obsTime=obsTime)

    # remove apogee skies that are also boss ids
    # dont know if this is right or not
    idx_dup = np.isin(catalogid_apogee, catalogid_boss)
    catalogid_apogee = catalogid_apogee[~idx_dup]
    ra_apogee = ra_apogee[~idx_dup]
    dec_apogee = dec_apogee[~idx_dup]
    target_pk_apogee = target_pk_apogee[~idx_dup]
    x_apogee = x_apogee[~idx_dup]
    y_apogee = y_apogee[~idx_dup]

    # create field object
    f = field.Field(racen=racen, deccen=deccen, pa=position_angle,
                    field_cadence=cadence, observatory=observatory.lower())

    # set the required skies
    f.required_calibrations['sky_boss'] = n_sky_boss
    f.required_calibrations['sky_apogee'] = n_sky_apogee

    # create array for RS field
    N = len(ra_apogee) + len(ra_boss)
    # these are datatypes from robostrategy.Field
    targets_dtype = np.dtype([('ra', np.float64),
                              ('dec', np.float64),
                              ('x', np.float64),
                              ('y', np.float64),
                              ('within', np.int32),
                              ('incadence', np.int32),
                              ('priority', np.int32),
                              ('value', np.float32),
                              ('program', np.unicode_, 30),
                              ('carton', np.unicode_, 30),
                              ('category', np.unicode_, 30),
                              ('cadence', np.unicode_, 30),
                              ('fiberType', np.unicode_, 10),
                              ('catalogid', np.int64),
                              ('rsid', np.int64),
                              ('target_pk', np.int64)])
    targs = np.zeros(N, dtype=targets_dtype)
    targs['ra'] = np.append(ra_apogee, ra_boss)
    targs['dec'] = np.append(dec_apogee, dec_boss)
    targs['x'] = np.append(x_apogee, x_boss)
    targs['y'] = np.append(y_apogee, y_boss)
    targs['within'] = np.zeros(len(targs), dtype=np.int32) + 1
    targs['incadence'] = np.zeros(len(targs), dtype=np.int32) + 1
    targs['priority'] = np.zeros(len(targs), dtype=np.int32) + 100000
    targs['value'] = np.zeros(len(targs), dtype=np.int32) + 1
    targs['program'] = np.array(['CALIBRATION'] * len(targs),
                                dtype='<U30')
    targs['carton'] = np.array(['CALIBRATION'] * len(targs),
                               dtype='<U30')
    targs['category'] = np.append(np.array(['sky_apogee'] * len(ra_apogee),
                                           dtype='<U30'),
                                  np.array(['sky_boss'] * len(ra_boss),
                                           dtype='<U30'))
    targs['cadence'] = np.array([cadence] * len(targs),
                                dtype='<U30')
    targs['fiberType'] = np.append(np.array(['APOGEE'] * len(ra_apogee),
                                            dtype='<U10'),
                                   np.array(['BOSS'] * len(ra_boss),
                                            dtype='<U10'))
    targs['catalogid'] = np.append(catalogid_apogee, catalogid_boss)
    targs['rsid'] = np.arange(len(targs), dtype=np.int64) + 1
    targs['target_pk'] = np.append(target_pk_apogee, target_pk_boss)

    # assign targets
    f.targets_fromarray(targs)

    f.assign()
    # create mugatu FPSDesign object
    catalogids = np.zeros(500, dtype=np.int64) - 1
    ra = np.zeros(500, dtype=float) - 9999.99
    dec = np.zeros(500, dtype=float) - 9999.99
    fiberID = np.zeros(500, dtype=np.int64) - 1
    obsWavelength = np.zeros(500, dtype='<U10')
    priority = np.zeros(500, dtype=int) - 1

    for i in range(len(f.targets)):
        if f.assignments[i][2][0] != -1:
            rid = f.assignments[i][2][0]
            catalogids[rid] = f.targets[i]['catalogid']
            ra[rid] = f.targets[i]['ra']
            dec[rid] = f.targets[i]['dec']
            fiberID[rid] = rid
            obsWavelength[rid] = f.targets[i]['fiberType']
            priority[rid] = f.targets[i]['priority']

    fps_design = FPSDesign(design_pk=-1,
                           obsTime=obsTime,
                           racen=racen,
                           deccen=deccen,
                           position_angle=position_angle,
                           observatory=observatory,
                           mode_pk=None,
                           catalogids=catalogids,
                           ra=ra,
                           dec=dec,
                           fiberID=fiberID,
                           obsWavelength=obsWavelength,
                           priority=priority,
                           design_file=None,
                           manual_design=True)

    return f, fps_design


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
