# @Author: Ilija Medan
# @Date: December 29, 2020
# @Filename: fpsdesign.py
# @License: BSD 3-Clause
# @Copyright: Ilija Medan

import numpy as np
import warnings
from astropy.io import fits
from astropy.time import Time
from peewee import JOIN

import kaiju
import kaiju.robotGrid
from mugatu.exceptions import (MugatuError, MugatuWarning,
                               MugatuDesignError, MugatuDesignWarning,
                               MugatuDesignModeWarning)
from coordio.utils import radec2wokxy
from mugatu.designs_to_targetdb import (make_design_assignments_targetdb,
                                        make_design_field_targetdb)
from mugatu.designmode import DesignModeCheck

try:
    from sdssdb.peewee.sdss5db import database
    database.set_profile('operations')
    _database = True
except:
    _database = False

from sdssdb.peewee.sdss5db.targetdb import (Design, Field, Observatory,
                                            Assignment, Instrument, Target,
                                            Hole, CartonToTarget, Carton,
                                            Magnitude, Category)


class FPSDesign(object):
    """
    Class to load, validate and export a design.

    Parameters
    ----------
    design_pk: int
        The pk of the design as it appears in targetdb.

    obsTime: float
        Julian date of the observation.

    racen: np.float64
        Right Ascension of center of the field in degrees.

    deccen: np.float64
        Declination of center of the field in degrees.

    position_angle: np.float64
        Position angle of the field E of N in degrees.

    observatory: str
        Observatory where observation is taking place, either
        'LCO' or 'APO'.

    desmode_label: int
        The pk in targetdb for the observing mode for the design.

    idtype: str
        ID type used for catalogids in design. Must be 'catalogID'
        or 'carton_to_target'.

    catalogids: np.array
        List of ids for a manual design in db. Default from idtype
        is carton_to_target.

    ra: np.array
        List of right ascensions that correspond to catalogids.
        If list not provided, but catalogids for a manual design
        are provided, right ascensions will be pulled from
        targetdb.

    dec: np.array
        List of declinations that correspond to catalogids.
        If list not provided, but catalogids for a manual design
        are provided, declinations will be pulled from
        targetdb.

    pmra: np.array
        Array of proper motions in the RA axis for the targets, in
        milliarcsec/yr. Must be a true angle, i.e, it must include
        the cos(dec) term. If no proper motion, set index to zero.

    pmdec: np.array
        Array of proper motions in the DEC axis for the targets, in
        milliarcsec/yr. If no proper motion, set index to zero.

    delta_ra: np.array
        List of offsets in right ascensions that correspond to catalogids.
        If list not provided, but catalogids for a manual design
        are provided, right ascensions will be pulled from
        targetdb.

    delta_dec: np.array
        List of offsets in declinations that correspond to catalogids.
        If list not provided, but catalogids for a manual design
        are provided, declinations will be pulled from
        targetdb.

    epoch: np.array
        Array of epochs for the coordinates in decimal years.

    robotID: np.array
        Fiber assignement for each catalogid target in the
        manual design.

    obsWavelength: np.array
        Wavelength of observation for each fiber
        that correspond to catalogids.
        Must be either 'BOSS' or 'APOGEE'.

    priority: np.array
        Priorties for targets that correspond to catalogids.

    carton_pk: np.array
        The carton_pks (as in targetdb) for the targets
        that correspond to catalogids.

    category: np.array
        The category for each target. Can be 'science', 'sky_INSTRUMENT'
        or 'standard_INSTRUMENT'

    magnitudes: np.array
        Magnitudes of the targets. Should be of size (N, 7), where
        columns correspond to g, r, i, bp, gaia_g, rp and h band
        magnitudes.

    design_file: str
        FITS file with a manual design. The file must be the same
        format as the rsFieldAssignments file.

    manual_design: boolean
        Boolean if the design being validated is manual
        (manual_design=True) or in targetdb (manual_design=False)

    exp: int
        Exposure number for design in file. If exp=0, assumes
        only 1 exposure in the design file. if exp>0, then will
        choose exposure exp = 1 to N.

    collisionBuffer: float
        collisionBuffer parmameter for Kaiju RobotGrid.

    Attributes
    ----------
    design: dict
        contains all of the design inputs for Kaiju needed to validate
        a design.

    rg: Kaiju RobotGrid object
        Kaiju RobotGrid with the assignments from design.

    targets_unassigned: list
        catalogid of targets that could not be assigned due to assigned
        robot not being able to reach the assigned target.

    holeID_mapping: np.array
        Mapping between robotID (index + 1) and holeID

    targets_collided: list
        catalogid of targets that could not be assigned due to assigned
        robot to assigned target resulting in a collision.

    valid_design: dictonary
        Same format as design dictonary, except this is validated
        design by Kaiju with collision/delock assignments removed.

    design_built: booleen
        True if design dictonary has been populated, False if not.

    hourAngle: float
        Hour angle of field center from coordio.utils.radec2wokxy.

    positionAngle_coordio: float
        position angle of field center from coordio.utils.radec2wokxy.

    design_errors: dict
        dictonary with errors for design. Created and filled after validation.
    """

    def __init__(self, design_pk, obsTime, racen=None, deccen=None,
                 position_angle=None, observatory=None, desmode_label=None,
                 idtype='carton_to_target', catalogids=None, ra=None, dec=None,
                 pmra=None, pmdec=None, delta_ra=None, delta_dec=None,
                 epoch=None, robotID=None, obsWavelength=None,
                 priority=None, carton_pk=None, category=None, magnitudes=None,
                 design_file=None, manual_design=False, exp=0,
                 collisionBuffer=2.):
        if idtype != 'catalogID' and idtype != 'carton_to_target':
            message = 'idtype must be catalogID or carton_to_target'
            raise MugatuError(message=message)
        self.design_pk = design_pk
        self.obsTime = obsTime
        self.design = {}
        # either set field params or pull from db is design
        # is in targetdb
        if manual_design:
            if design_file is None:
                self.racen = racen
                self.deccen = deccen
                self.position_angle = position_angle
                self.observatory = observatory
                self.desmode_label = desmode_label
            else:
                head = fits.open(design_file)[0].header
                desmode_labels = head['DESMODE'].split(' ')
                self.racen = head['RACEN']
                self.deccen = head['DECCEN']
                self.position_angle = head['PA']
                self.observatory = head['obs'].strip().upper()
                self.desmode_label = desmode_labels[exp - 1]
        else:
            design_field_db = (
                Design.select(Design.design_id,
                              Design.design_mode,
                              Field.racen,
                              Field.deccen,
                              Field.position_angle,
                              Observatory.label)
                      .join(Field,
                            on=(Design.field_pk == Field.pk))
                      .join(Observatory,
                            on=(Field.observatory == Observatory.pk))
                      .where(Design.design_id == self.design_pk))

            self.racen = design_field_db.objects()[0].racen
            self.deccen = design_field_db.objects()[0].deccen
            self.position_angle = design_field_db.objects()[0].position_angle
            self.observatory = design_field_db.objects()[0].label
            self.desmode_label = design_field_db.objects()[0].design_mode.label

        # should these be catalogids or carton_to_target?
        self.idtype = idtype
        self.catalogids = catalogids
        self.ra = ra
        self.dec = dec
        self.pmra = pmra
        self.pmdec = pmdec
        self.delta_ra = delta_ra
        self.delta_dec = delta_dec
        self.epoch = epoch
        self.robotID = robotID
        self.obsWavelength = obsWavelength
        self.priority = priority
        self.carton_pk = carton_pk
        self.category = category
        self.magnitudes = magnitudes
        self.design_file = design_file
        self.manual_design = manual_design
        self.exp = exp
        self.collisionBuffer = collisionBuffer

        # set dummy value for collision for now
        # this may want to be a input, not sure the standard here
        # initialize robotGrid
        if self.observatory == 'APO':
            self.rg = kaiju.robotGrid.RobotGridAPO(
                collisionBuffer=self.collisionBuffer,
                stepSize=0.05)
        else:
            self.rg = kaiju.robotGrid.RobotGridLCO(
                collisionBuffer=self.collisionBuffer,
                stepSize=0.05)
        self.holeID_mapping = np.zeros(500, dtype='<U10')
        for i, robotID in enumerate(range(1, 501)):
            self.holeID_mapping[i] = self.rg.robotDict[robotID].holeID
        # this is in Conor's test, I'm not quite sure what it does
        # but without paths wont generate
        for k in self.rg.robotDict.keys():
            self.rg.homeRobot(k)
        # for rID in self.rg.robotDict:
            # robot = self.rg.getRobot(rID)
            # robot.setXYUniform()
        self.targets_unassigned = []
        self.targets_collided = []

        self.valid_design = {}

        self.design_built = False

    def _offset_radec(self, ra=None, dec=None, delta_ra=0., delta_dec=0.):
        """Offsets ra and dec according to specified amount. From Mike's
        robostrategy.Field object

        Parameters
        ----------
        ra : np.float64 or ndarray of np.float64
        right ascension, deg
        dec : np.float64 or ndarray of np.float64
            declination, deg
        delta_ra : np.float64 or ndarray of np.float64
            right ascension direction offset, arcsec
        delta_dec : np.float64 or ndarray of np.float64
            declination direction offset, arcsec

        Returns
        -------
        offset_ra : np.float64 or ndarray of np.float64
            offset right ascension, deg
        offset_dec : np.float64 or ndarray of np.float64
            offset declination, deg

        Notes
        -----
        Assumes that delta_ra, delta_dec are in proper coordinates; i.e.
        an offset of delta_ra=1 arcsec represents the same angular separation
        on the sky at any declination.
        Carefully offsets in the local directions of ra, dec based on
        the local tangent plane (i.e. does not just scale delta_ra by
        1/cos(dec))
        """
        deg2rad = np.pi / 180.
        arcsec2rad = np.pi / 180. / 3600.
        x = np.cos(dec * deg2rad) * np.cos(ra * deg2rad)
        y = np.cos(dec * deg2rad) * np.sin(ra * deg2rad)
        z = np.sin(dec * deg2rad)
        ra_x = - np.sin(ra * deg2rad)
        ra_y = np.cos(ra * deg2rad)
        ra_z = 0.
        dec_x = - np.sin(dec * deg2rad) * np.cos(ra * deg2rad)
        dec_y = - np.sin(dec * deg2rad) * np.sin(ra * deg2rad)
        dec_z = np.cos(dec * deg2rad)
        xoff = x + (ra_x * delta_ra + dec_x * delta_dec) * arcsec2rad
        yoff = y + (ra_y * delta_ra + dec_y * delta_dec) * arcsec2rad
        zoff = z + (ra_z * delta_ra + dec_z * delta_dec) * arcsec2rad
        offnorm = np.sqrt(xoff**2 + yoff**2 + zoff**2)
        xoff = xoff / offnorm
        yoff = yoff / offnorm
        zoff = zoff / offnorm
        decoff = np.arcsin(zoff) / deg2rad
        raoff = ((np.arctan2(yoff, xoff) / deg2rad) + 360.) % 360.
        return(raoff, decoff)

    def build_design_db(self):
        """
        Populate the design dictonary for design in targetdb.

        Notes
        -----
        This method will build a design from the db if manual_design=False

        """

        if self.manual_design:
            flag = 'Building database design with manual_design=True flag'
            warnings.warn(flag, MugatuWarning)

        # initialize dict for the design
        # not using None or nan for no assignments
        # using -1 (for int) and -9999.99 (for float) for None assignment
        self.design['design_pk'] = self.design_pk
        self.design['desmode_label'] = self.desmode_label
        self.design['catalogID'] = np.zeros(500, dtype=np.int64) - 1
        self.design['robotID'] = np.zeros(500, dtype=np.int64) - 1
        self.design['holeID'] = np.zeros(500, dtype='<U10')
        self.design['obsWavelength'] = np.zeros(500, dtype='<U6')
        self.design['priority'] = np.zeros(500, dtype=int) - 1
        self.design['carton_pk'] = np.zeros(500, dtype=int) - 1
        self.design['category'] = np.zeros(500, dtype='<U10')
        self.design['ra'] = np.zeros(500, dtype=float) - 9999.99
        self.design['dec'] = np.zeros(500, dtype=float) - 9999.99
        self.design['delta_ra'] = np.zeros(500, dtype=float) - 9999.99
        self.design['delta_dec'] = np.zeros(500, dtype=float) - 9999.99
        self.design['epoch'] = np.zeros(500, dtype=float) - 9999.99
        self.design['ra_off'] = np.zeros(500, dtype=float) - 9999.99
        self.design['dec_off'] = np.zeros(500, dtype=float) - 9999.99
        self.design['pmra'] = np.zeros(500, dtype=float) - 9999.99
        self.design['pmdec'] = np.zeros(500, dtype=float) - 9999.99
        self.design['x'] = np.zeros(500, dtype=float) - 9999.99
        self.design['y'] = np.zeros(500, dtype=float) - 9999.99
        self.design['magnitudes'] = np.zeros((500, 7), dtype=float) - 9999.99

        # need to add wokHole to query when in db (not there now)
        # I need to test this when v05 is up, im unsure about Joins
        design_targ_db = (
            Assignment.select(Target.catalogid,
                              Hole.holeid,
                              Instrument.label,
                              CartonToTarget.priority,
                              Target.ra,
                              Target.dec,
                              CartonToTarget.carton,
                              Carton.pk.alias('carton_pk'),
                              CartonToTarget.pk,
                              Magnitude.g,
                              Magnitude.r,
                              Magnitude.i,
                              Magnitude.bp,
                              Magnitude.gaia_g,
                              Magnitude.rp,
                              Magnitude.h,
                              CartonToTarget.delta_ra,
                              CartonToTarget.delta_dec,
                              Target.pmra,
                              Target.pmdec,
                              Target.epoch,
                              Category.label.alias('cat_lab'))
                      .join(Hole)
                      .switch(Assignment)
                      .join(Instrument)
                      .switch(Assignment)
                      .join(CartonToTarget)
                      .join(Target)
                      .switch(CartonToTarget)
                      .join(Magnitude, JOIN.LEFT_OUTER)
                      .switch(CartonToTarget)
                      .join(Carton, JOIN.LEFT_OUTER)
                      .join(Category, JOIN.LEFT_OUTER)
                      .where(Assignment.design_id == self.design_pk))

        for d in design_targ_db.objects():
            # assign to index that corresponds to fiber assignment
            # index should match length of arrays
            pos_id = np.where(self.holeID_mapping == d.holeid)[0][0]
            if self.idtype == 'carton_to_target':
                self.design['catalogID'][pos_id] = d.pk
            else:
                self.design['catalogID'][pos_id] = d.catalogid

            self.design['robotID'][pos_id] = pos_id + 1
            self.design['holeID'][pos_id] = d.holeid
            # design['wokHoleID'][i] = design_targ_db[i]
            self.design['obsWavelength'][pos_id] = d.label
            # catch targets with no assigned priority
            try:
                self.design['priority'][pos_id] = d.priority
            except AttributeError:
                self.design['priority'][pos_id] = -1
            self.design['carton_pk'][pos_id] = d.carton_pk
            self.design['category'][pos_id] = d.cat_lab
            self.design['ra'][pos_id] = d.ra
            self.design['dec'][pos_id] = d.dec
            self.design['delta_ra'][pos_id] = d.delta_ra
            self.design['delta_dec'][pos_id] = d.delta_dec
            self.design['pmra'][pos_id] = d.pmra
            self.design['pmdec'][pos_id] = d.pmdec
            self.design['epoch'][pos_id] = d.epoch
            self.design['magnitudes'][pos_id][0] = d.g
            self.design['magnitudes'][pos_id][1] = d.r
            self.design['magnitudes'][pos_id][2] = d.i
            self.design['magnitudes'][pos_id][3] = d.bp
            self.design['magnitudes'][pos_id][4] = d.gaia_g
            self.design['magnitudes'][pos_id][5] = d.rp
            self.design['magnitudes'][pos_id][6] = d.h

        # set nan pm tp zero
        self.design['pmra'][np.isnan(self.design['pmra'])] = 0.
        self.design['pmdec'][np.isnan(self.design['pmdec'])] = 0.
        # here convert ra/dec to x/y based on field/time of observation
        # I think I need to add inertial in here at some point,
        # dont see this in targetdb though
        ev = eval("(self.design['ra'] != -9999.99)")
        res = self._offset_radec(ra=self.design['ra'][ev],
                                 dec=self.design['dec'][ev],
                                 delta_ra=self.design['delta_ra'][ev],
                                 delta_dec=self.design['delta_dec'][ev])
        self.design['ra_off'][ev], self.design['dec_off'][ev] = res
        with warnings.catch_warnings(record=True) as w:
            res = radec2wokxy(
                ra=self.design['ra_off'][ev],
                dec=self.design['dec_off'][ev],
                coordEpoch=Time(self.design['epoch'][ev],
                                format='decimalyear').jd,
                waveName=np.array(list(map(lambda x: x.title(),
                                           self.design['obsWavelength'][ev]))),
                raCen=self.racen,
                decCen=self.deccen,
                obsAngle=self.position_angle,
                obsSite=self.observatory,
                obsTime=self.obsTime)  #,
                # pmra=self.design['pmra'],
                # pmdec=self.design['pmdec'])
            self.design['x'][ev] = res[0]
            self.design['y'][ev] = res[1]
            fieldWarn = res[2]
            self.hourAngle = res[3]
            self.positionAngle_coordio = res[4]
            if len(w) > 0:
                for wm in w:
                    if 'iauPmsafe return' in wm.message.args[0]:
                        flag = ('Coordio xy coordinates converted '
                                'should be eyed with suspicion.')
                        warnings.warn(flag, MugatuDesignWarning)

        if np.any(fieldWarn):
            flag = ('Coordio xy coordinates converted '
                    'should be eyed with suspicion.')
            warnings.warn(flag, MugatuDesignWarning)

        self.design_built = True

        return

    def build_design_manual(self):
        """
        Populate the design dictonary for manual design.

        Notes
        -----
        This function creates a manual design whether it is from
        user inputted catalogids, or if it is a FITS file
        if manual_design=True.

        """

        if not self.manual_design:
            flag = 'Building manual design with manual_design=False flag'
            warnings.warn(flag, MugatuWarning)

        # initialize dict for the design
        # for manual design, dont make arrays since length
        # will be defined by user
        self.design['design_pk'] = self.design_pk
        self.design['desmode_label'] = self.desmode_label

        if self.design_file is None:
            # manual design with targets in targetdb
            self.design['catalogID'] = self.catalogids
            self.design['obsWavelength'] = self.obsWavelength
            self.design['priority'] = self.priority
            self.design['carton_pk'] = self.carton_pk
            self.design['category'] = self.category

            if self.ra is None:
                for i in range(len(self.design['catalogID'])):
                    targ_db = Target.get(catalogid=self.design['catalogID'][i])
                    self.design['ra'][i] = targ_db.ra
                    self.design['dec'][i] = targ_db.dec
                    self.design['pmra'][i] = targ_db.pmra
                    self.design['pmdec'][i] = targ_db.pmdec
                    self.design['delta_ra'][i] = targ_db.delta_ra
                    self.design['delta_dec'][i] = targ_db.delta_dec
            else:
                self.design['ra'] = self.ra
                self.design['dec'] = self.dec
                self.design['pmra'] = self.pmra
                self.design['pmdec'] = self.pmdec
                self.design['delta_ra'] = self.delta_ra
                self.design['delta_dec'] = self.delta_dec
                self.design['epoch'] = self.epoch

            self.design['robotID'] = self.robotID
            self.design['holeID'] = np.zeros(len(self.design['robotID']),
                                             dtype='<U10')
            for i in range(len(self.design['robotID'])):
                self.design['holeID'][i] = self.holeID_mapping[self.design['robotID'][i] - 1]
            self.design['magnitudes'] = self.magnitudes
        else:
            # manual design from flat file
            # get header with field info
            head = fits.open(self.design_file)[0].header
            desmode_labels = head['DESMODE'].split(' ')
            # association between catalogid and instrument
            design_inst = fits.open(self.design_file)[1].data
            # catalogid assignment for each fiber
            design = fits.open(self.design_file)[2].data

            # grab obs info
            self.racen = head['RACEN']
            self.deccen = head['DECCEN']
            self.position_angle = head['PA']
            self.observatory = head['obs'].strip().upper()
            self.desmode_label = desmode_labels[self.exp - 1]

            # grab assignment info
            if self.exp == 0:
                roboIDs = design['robotID']
            else:
                roboIDs = design['robotID'][:, self.exp - 1]
            if self.idtype == 'catalogID':
                self.design['catalogID'] = (design_inst['catalogid']
                                            [roboIDs != -1])
            else:
                self.design['catalogID'] = (design_inst['carton_to_target_pk']
                                            [roboIDs != -1])
            self.design['ra'] = design_inst['ra'][roboIDs != -1]
            self.design['dec'] = design_inst['dec'][roboIDs != -1]
            self.design['delta_ra'] = design_inst['delta_ra'][roboIDs != -1]
            self.design['delta_dec'] = design_inst['delta_dec'][roboIDs != -1]
            self.design['pmra'] = design_inst['pmra'][roboIDs != -1]
            self.design['pmdec'] = design_inst['pmdec'][roboIDs != -1]
            self.design['epoch'] = design_inst['epoch'][roboIDs != -1]
            self.design['robotID'] = roboIDs[roboIDs != -1]
            self.design['holeID'] = np.zeros(len(self.design['robotID']),
                                             dtype='<U10')
            for i in range(len(self.design['robotID'])):
                self.design['holeID'][i] = self.holeID_mapping[self.design['robotID'][i] - 1]
            self.design['obsWavelength'] = (design_inst['fiberType']
                                            [roboIDs != -1])
            self.design['priority'] = design_inst['priority'][roboIDs != -1]
            # need to change this
            self.design['carton_pk'] = np.arange(0,
                                                 len(self.design['catalogID']),
                                                 1,
                                                 dtype=int)
            self.design['category'] = design_inst['category'][roboIDs != -1]
            self.design['magnitudes'] = design_inst['magnitude'][roboIDs != -1]

        # set nan pm tp zero
        self.design['pmra'][np.isnan(self.design['pmra'])] = 0.
        self.design['pmdec'][np.isnan(self.design['pmdec'])] = 0.
        # make empty x,y arrays
        self.design['x'] = (np.zeros(len(self.design['catalogID']),
                            dtype=float) -
                            9999.99)
        self.design['y'] = (np.zeros(len(self.design['catalogID']),
                            dtype=float) -
                            9999.99)
        self.design['ra_off'] = (np.zeros(len(self.design['catalogID']),
                                 dtype=float) -
                                 9999.99)
        self.design['dec_off'] = (np.zeros(len(self.design['catalogID']),
                                  dtype=float) -
                                  9999.99)

        # here convert ra/dec to x/y based on field/time of observation
        ev = eval("(self.design['ra'] != -9999.99)")
        res = self._offset_radec(ra=self.design['ra'][ev],
                                 dec=self.design['dec'][ev],
                                 delta_ra=self.design['delta_ra'][ev],
                                 delta_dec=self.design['delta_dec'][ev])
        self.design['ra_off'][ev], self.design['dec_off'][ev] = res
        with warnings.catch_warnings(record=True) as w:
            res = radec2wokxy(
                ra=self.design['ra_off'][ev],
                dec=self.design['dec_off'][ev],
                coordEpoch=Time(self.design['epoch'][ev],
                                format='decimalyear').jd,
                waveName=np.array(list(map(lambda x: x.title(),
                                           self.design['obsWavelength'][ev]))),
                raCen=self.racen,
                decCen=self.deccen,
                obsAngle=self.position_angle,
                obsSite=self.observatory,
                obsTime=self.obsTime)  #,
                # pmra=self.design['pmra'],
                # pmdec=self.design['pmdec'])
            self.design['x'][ev] = res[0]
            self.design['y'][ev] = res[1]
            fieldWarn = res[2]
            self.hourAngle = res[3]
            self.positionAngle_coordio = res[4]
            if len(w) > 0:
                for wm in w:
                    if 'iauPmsafe return' in wm.message.args[0]:
                        flag = ('Coordio xy coordinates converted '
                                'should be eyed with suspicion.')
                        warnings.warn(flag, MugatuDesignWarning)

        if np.any(fieldWarn):
            flag = ('Coordio xy coordinates converted '
                    'should be eyed with suspicion.')
            warnings.warn(flag, MugatuDesignWarning)

        self.design_built = True

        return

    def design_to_RobotGrid(self):
        """
        Add assignments to Kaiju RobotGrid.

        Notes
        -----
        Adds targets to the kaiju.robotGrid.RobotGridFilledHex and
        assigns them to fibers based on the design dictonary.

       """

        is_unassigned = False

        for i in range(len(self.design['x'])):
            if self.design['robotID'][i] != -1:
                if self.design['obsWavelength'][i] == 'BOSS':
                    self.rg.addTarget(targetID=self.design['catalogID'][i],
                                      x=self.design['x'][i],
                                      y=self.design['y'][i],
                                      priority=self.design['priority'][i],
                                      fiberType=kaiju.cKaiju.BossFiber)
                else:
                    self.rg.addTarget(targetID=self.design['catalogID'][i],
                                      x=self.design['x'][i],
                                      y=self.design['y'][i],
                                      priority=self.design['priority'][i],
                                      fiberType=kaiju.cKaiju.ApogeeFiber)
        for i in range(len(self.design['x'])):
            if self.design['robotID'][i] > 0:
                try:
                    self.rg.assignRobot2Target(self.design['robotID'][i],
                                               self.design['catalogID'][i])
                except RuntimeError:
                    # this catches the fact that robot cant be
                    # assigned to given fiber
                    self.targets_unassigned.append(self.design['catalogID'][i])
                    is_unassigned = True

        if is_unassigned:
            flag = 'Some targets could not be assigned to fiber'
            warnings.warn(flag, MugatuDesignWarning)

        return

    def decollide_grid(self):
        """
        Check to see if any collisions in grid, and if so,
        decollide grid and record assignments that were
        removed
        """
        # decollide the unassigned robots
        for robotID in self.rg.robotDict:
            if(self.rg.robotDict[robotID].isAssigned() is False):
                self.rg.decollideRobot(robotID)
                # need to still set alpha/beta
                # robot = self.rg.getRobot(robotID)
                # robot.setDestinationAlphaBeta(0, 180)
        if self.rg.getNCollisions() > 0:
            self.design_errors['no_collisions'] = False
            self.rg.decollideGrid()

            # check if de-collision was successful
            if not self.rg.getNCollisions() == 0:
                raise MugatuDesignError(message='Kaiju decollideGrid failed')

            flag = 'Some targets removed from design due to collisions'
            warnings.warn(flag, MugatuDesignWarning)

            # grab all of the targets removed due to collisions
            for i in self.rg.robotDict:
                fiber_idx = np.where(self.design['robotID'] == i)[0]
                if (self.rg.robotDict[i].assignedTargetID == -1 and
                   len(fiber_idx) > 0):
                    targ_remove = self.design['catalogID'][fiber_idx[0]]
                    if targ_remove not in self.targets_unassigned:
                        self.targets_collided.append(targ_remove)
        else:
            self.design_errors['no_collisions'] = True
        return

    def designmode_validate(self, db_query_results_boss=None,
                            db_query_results_apogee=None,
                            desmode_manual=None):
        """
        Check all designmode parameters for a design
        """
        mode = DesignModeCheck(FPSDesign=self,
                               desmode_label=self.desmode_label,
                               db_query_results_boss=db_query_results_boss,
                               db_query_results_apogee=db_query_results_apogee,
                               desmode_manual=desmode_manual)
        mode.design_mode_check_all(verbose=False)
        self.design_errors['min_skies_boss'] = mode.n_skies_min_check['BOSS']
        self.design_errors['min_skies_boss_metric'] = (mode
                                                       .n_skies_min_check
                                                       ['BOSS_metric'])
        if self.design_errors['min_skies_boss'] is False:
            flag = 'Design does not meet minimum BOSS skies for DesignMode'
            warnings.warn(flag, MugatuDesignModeWarning)
        self.design_errors['min_skies_apogee'] = (mode
                                                  .n_skies_min_check['APOGEE'])
        self.design_errors['min_skies_apogee_metric'] = (mode
                                                         .n_skies_min_check
                                                         ['APOGEE_metric'])
        if self.design_errors['min_skies_apogee'] is False:
            flag = 'Design does not meet minimum APOGEE skies for DesignMode'
            warnings.warn(flag, MugatuDesignModeWarning)

        self.design_errors['fov_skies_boss'] = (mode
                                                .min_skies_fovmetric_check
                                                ['BOSS'])
        self.design_errors['fov_skies_boss_metric'] = (mode
                                                       .min_skies_fovmetric_check
                                                       ['BOSS_metric'])
        if self.design_errors['fov_skies_boss'] is False:
            flag = ('Design does not meet FOV criteria '
                    'for BOSS skies for DesignMode')
            warnings.warn(flag, MugatuDesignModeWarning)
        self.design_errors['fov_skies_apogee'] = (mode
                                                  .min_skies_fovmetric_check
                                                  ['APOGEE'])
        self.design_errors['fov_skies_apogee_metric'] = (mode
                                                         .min_skies_fovmetric_check
                                                         ['APOGEE_metric'])
        if self.design_errors['fov_skies_apogee'] is False:
            flag = ('Design does not meet FOV criteria '
                    'for APOGEE skies for DesignMode')
            warnings.warn(flag, MugatuDesignModeWarning)

        self.design_errors['min_stds_boss'] = mode.n_stds_min_check['BOSS']
        self.design_errors['min_stds_boss_metric'] = (mode
                                                      .n_stds_min_check
                                                      ['BOSS_metric'])
        if self.design_errors['min_stds_boss'] is False:
            flag = 'Design does not meet minimum BOSS standards for DesignMode'
            warnings.warn(flag, MugatuDesignModeWarning)
        self.design_errors['min_stds_apogee'] = mode.n_stds_min_check['APOGEE']
        self.design_errors['min_stds_apogee_metric'] = (mode
                                                        .n_stds_min_check
                                                        ['APOGEE_metric'])
        if self.design_errors['min_stds_apogee'] is False:
            flag = ('Design does not meet minimum '
                    'APOGEE standards for DesignMode')
            warnings.warn(flag, MugatuDesignModeWarning)

        self.design_errors['fov_stds_boss'] = (mode
                                               .min_stds_fovmetric_check
                                               ['BOSS'])
        self.design_errors['fov_stds_boss_metric'] = (mode
                                                      .min_stds_fovmetric_check
                                                      ['BOSS_metric'])
        if self.design_errors['fov_stds_boss'] is False:
            flag = ('Design does not meet FOV criteria '
                    'for BOSS standards for DesignMode')
            warnings.warn(flag, MugatuDesignModeWarning)
        self.design_errors['fov_stds_apogee'] = (mode
                                                 .min_stds_fovmetric_check
                                                 ['APOGEE'])
        self.design_errors['fov_stds_apogee_metric'] = (mode
                                                        .min_stds_fovmetric_check
                                                        ['APOGEE_metric'])
        if self.design_errors['fov_stds_apogee'] is False:
            flag = ('Design does not meet FOV criteria '
                    'for APOGEE standards for DesignMode')
            warnings.warn(flag, MugatuDesignModeWarning)

        self.design_errors['stds_mag_boss'] = np.all(
            mode.stds_mags_check['BOSS'][0][(self.design['catalogID'] != -1) &
                                            (self.design['category'] == 'standard_boss')])
        self.design_errors['stds_mag_boss_metric'] = mode.stds_mags_check['BOSS_metric']
        if self.design_errors['stds_mag_boss'] is False:
            flag = ('Design has BOSS standard '
                    'assignments too bright for DesignMode')
            warnings.warn(flag, MugatuDesignModeWarning)
        self.design_errors['stds_mag_apogee'] = np.all(
            mode.stds_mags_check['APOGEE'][0][(self.design['catalogID'] != -1) &
                                              (self.design['category'] == 'standard_apogee')])
        self.design_errors['stds_mag_apogee_metric'] = mode.stds_mags_check['APOGEE_metric']
        if self.design_errors['stds_mag_apogee'] is False:
            flag = ('Design has APOGEE standard '
                    'assignments too bright for DesignMode')
            warnings.warn(flag, MugatuDesignModeWarning)

        self.design_errors['sci_mag_boss'] = np.all(
            mode.bright_limit_targets_check['BOSS'][0][(self.design['catalogID'] != -1) &
                                                       (self.design['category'] == 'science') &
                                                       (self.design['obsWavelength'] == 'BOSS')])
        self.design_errors['sci_mag_boss_metric'] = mode.bright_limit_targets_check['BOSS_metric']
        if self.design_errors['sci_mag_boss'] is False:
            flag = ('Design has BOSS science assignments '
                    ' too bright for DesignMode')
            warnings.warn(flag, MugatuDesignModeWarning)
        self.design_errors['sci_mag_apogee'] = np.all(
            mode.bright_limit_targets_check['APOGEE'][0][(self.design['catalogID'] != -1) &
                                                         (self.design['category'] == 'science') &
                                                         (self.design['obsWavelength'] ==
                                                          'APOGEE')])
        self.design_errors['sci_mag_apogee_metric'] = (mode
                                                       .bright_limit_targets_check
                                                       ['APOGEE_metric'])
        if self.design_errors['sci_mag_apogee'] is False:
            flag = ('Design has APOGEE science assignments '
                    'too bright for DesignMode')
            warnings.warn(flag, MugatuDesignModeWarning)

        self.design_errors['bright_neigh_boss'] = np.all(
            mode.bright_neighbor_check['BOSS'][0][mode.bright_neighbor_check['BOSS'][1]])
        self.design_errors['bright_neigh_boss_metric'] = mode.bright_neighbor_check['BOSS_metric']
        if self.design_errors['bright_neigh_boss'] == False:
            flag = ('Design has BOSS fibers too near '
                    'bright source for DesignMode')
            warnings.warn(flag, MugatuDesignModeWarning)

        self.design_errors['bright_neigh_apogee'] = np.all(
            mode.bright_neighbor_check['APOGEE'][0][mode.bright_neighbor_check['APOGEE'][1]])
        self.design_errors['bright_neigh_apogee_metric'] = (mode
                                                            .bright_neighbor_check
                                                            ['APOGEE_metric'])
        if self.design_errors['bright_neigh_apogee'] == False:
            flag = ('Design has APOGEE fibers too '
                    'near bright source for DesignMode')
            warnings.warn(flag, MugatuDesignModeWarning)
        return

    def bright_neigh_safety(self, db_query_results_boss=None,
                            db_query_results_apogee=None,
                            desmode_manual=None):
        """
        Perform safety version of bright neighbor check
        on the design
        """
        mode = DesignModeCheck(FPSDesign=self,
                               desmode_label=self.desmode_label,
                               db_query_results_boss=db_query_results_boss,
                               db_query_results_apogee=db_query_results_apogee,
                               desmode_manual=desmode_manual)
        bright_check_boss, hasFiber_boss = mode.bright_neighbors(instrument='BOSS',
                                                                 check_type='safety')
        bright_check_apogee, hasFiber_apogee = mode.bright_neighbors(instrument='APOGEE',
                                                                     check_type='safety')
        if (len(bright_check_boss[~bright_check_boss & hasFiber_boss]) > 0 or
           len(bright_check_apogee[~bright_check_apogee & hasFiber_apogee]) > 0):
            message = 'Bright Neighbor Safety Checked Failed,'
            message += (' %d BOSS and %d APOGEE fibers near bright sources' %
                        (len(bright_check_boss[~bright_check_boss & hasFiber_boss]),
                         len(bright_check_apogee[~bright_check_apogee & hasFiber_apogee])))
            raise MugatuDesignError(message=message)
        return

    def RobotGrid_to_valid_design(self):
        """
        Construct valid design from Kaiju Robotgrid
       """

        # initialize dict for the validated design
        # not using None or nan for no assignments
        # using -1 (for int) and -9999.99 (for float) for None assignment
        self.valid_design['design_pk'] = self.design_pk
        self.valid_design['desmode_label'] = self.desmode_label
        self.valid_design['catalogID'] = np.zeros(500, dtype=np.int64) - 1
        self.valid_design['robotID'] = np.zeros(500, dtype=np.int64) - 1
        self.valid_design['holeID'] = np.zeros(500, dtype='<U10')
        self.valid_design['obsWavelength'] = np.zeros(500, dtype='<U6')
        self.valid_design['priority'] = np.zeros(500, dtype=int) - 1
        self.valid_design['carton_pk'] = np.zeros(500, dtype=int) - 1
        self.valid_design['category'] = np.zeros(500, dtype='<U10')
        self.valid_design['ra'] = np.zeros(500, dtype=float) - 9999.99
        self.valid_design['dec'] = np.zeros(500, dtype=float) - 9999.99
        self.valid_design['pmra'] = np.zeros(500, dtype=float) - 9999.99
        self.valid_design['pmdec'] = np.zeros(500, dtype=float) - 9999.99
        self.valid_design['delta_ra'] = np.zeros(500, dtype=float) - 9999.99
        self.valid_design['delta_dec'] = np.zeros(500, dtype=float) - 9999.99
        self.valid_design['ra_off'] = np.zeros(500, dtype=float) - 9999.99
        self.valid_design['dec_off'] = np.zeros(500, dtype=float) - 9999.99
        self.valid_design['epoch'] = np.zeros(500, dtype=float) - 9999.99
        self.valid_design['x'] = np.zeros(500, dtype=float) - 9999.99
        self.valid_design['y'] = np.zeros(500, dtype=float) - 9999.99

        for i, rid in enumerate(self.rg.robotDict):
            self.valid_design['catalogID'][i] = (self.rg.robotDict[rid]
                                                 .assignedTargetID)
            if self.valid_design['catalogID'][i] != -1:
                self.valid_design['robotID'][i] = rid
                self.valid_design['holeID'][i] = self.holeID_mapping[rid - 1]
                # is below necessary? i dont know if decollide ever reassigns
                # or just removes
                cond = eval("self.design['catalogID'] == self.valid_design['catalogID'][i]")
                self.valid_design['obsWavelength'][i] = self.design['obsWavelength'][cond][0]
                self.valid_design['priority'][i] = self.design['priority'][cond][0]
                self.valid_design['carton_pk'][i] = self.design['carton_pk'][cond][0]
                self.valid_design['category'][i] = self.design['category'][cond][0]
                self.valid_design['ra'][i] = self.design['ra'][cond][0]
                self.valid_design['dec'][i] = self.design['dec'][cond][0]
                self.valid_design['pmra'][i] = self.design['pmra'][cond][0]
                self.valid_design['pmdec'][i] = self.design['pmdec'][cond][0]
                self.valid_design['delta_ra'][i] = self.design['delta_ra'][cond][0]
                self.valid_design['delta_dec'][i] = self.design['delta_dec'][cond][0]
                self.valid_design['ra_off'][i] = self.design['ra_off'][cond][0]
                self.valid_design['dec_off'][i] = self.design['dec_off'][cond][0]
                self.valid_design['epoch'][i] = self.design['epoch'][cond][0]
            self.valid_design['x'][i] = self.rg.robotDict[rid].xPos
            self.valid_design['y'][i] = self.rg.robotDict[rid].yPos

        return

    def validate_design(self, designmode=True, safety=True,
                        db_query_results_boss=None,
                        db_query_results_apogee=None,
                        desmode_manual=None):
        """
        Validate design for deadlocks and collisions using Kaiju.

        Parameters
        ----------
        designmode: bool
            Check designmodes for the design.

        safety: bool
            Check for bright neighbors using carton of brightest Gaia/2MASS
            stars in targetdb.

        db_query_results_boss: dict
            Database query results for BOSS bright neighbor check.
            Each index of dict is a tuple of (ras, decs, mags, catalogids)
            with one index for designmode and the other safety.

        db_query_results_apogee: dict
            Database query results for APOGEE bright neighbor check.
            Each index of dict is a tuple of (ras, decs, mags, catalogids)
            with one index for designmode and the other safety.
        """
        # make dict to store design errors that may come up
        self.design_errors = {}
        # build the design with all needed parameters if it has not been
        # built already (save time by doing this check)
        if not self.design_built:
            if self.manual_design:
                self.build_design_manual()
            else:
                self.build_design_db()

        # construct the Kaiju robotGrid
        self.design_to_RobotGrid()
        if len(self.targets_unassigned) > 0:
            self.design_errors['all_targets_assigned'] = False
        else:
            self.design_errors['all_targets_assigned'] = True

        # validate the design

        # de-collide the grid if collisions exist
        # and check for targets removed
        self.decollide_grid()

        # generate paths
        # self.rg.pathGen()
        # if self.rg.didFail:
        #     raise MugatuError(message='Kaiju pathGen failed')

        if designmode:
            self.designmode_validate(db_query_results_boss=db_query_results_boss,
                                     db_query_results_apogee=db_query_results_apogee,
                                     desmode_manual=desmode_manual)

        # do the safety check for design
        if safety:
            self.bright_neigh_safety(db_query_results_boss=db_query_results_boss,
                                     db_query_results_apogee=db_query_results_apogee,
                                     desmode_manual=desmode_manual)

        # I imagine that the above step would manipulate the robogrid based on
        # collisions and deadlocks, so the below would take these Kaiju results
        # and put them back as a design object which will be the output of the
        # function
        self.RobotGrid_to_valid_design()

        # another note to add to this above design dictonary, is how will this
        # include paths? This seems important for manual designs and to send
        # to Jaeger

        return

    def design_to_targetdb(self, cadence, fieldid, exposure):
        """
        Write a validated esign to targetdb

        Parameters
        ----------

        cadence: str
            Cadence label for the field.

        fieldid: int
            fieldid for the field.

        exposure: int
            The exposure of this set of designs, 0th indexed.

        Notes
        -----
        Version above should allow for manual designs to be added under
        version = 'manual' to seperate them from robostrategy ingested designs
        """
        if not self.manual_design:
            raise MugatuError(message='Can only write manual designs to targetdb')

        if len(self.valid_design) == 0:
            raise MugatuError(message='Need to validate design first')

        # create the manual field for the design
        make_design_field_targetdb(cadence=cadence,
                                   fieldid=fieldid,
                                   plan='manual',
                                   racen=self.racen,
                                   deccen=self.deccen,
                                   position_angle=self.position_angle,
                                   observatory=self.observatory,
                                   slots_exposures=None)

        # create dictonary for unique carton pks for idtype is catalogID
        if self.idtype == 'catalogID':
            cart_pks = {}
            targetdb_ver = {}
            for pk in np.unique(self.valid_design['carton_pk'][self.valid_design['catalogID'] != -1]):
                cart_pks[pk] = pk
                targetdb_ver[pk] = Carton.get(pk).version_pk
        else:
            cart_pks = None
            targetdb_ver = None

        # add the design to targetdb
        make_design_assignments_targetdb(
            plan='manual',
            fieldid=fieldid,
            exposure=exposure,
            desmode_label=self.desmode_label,
            catalogID=self.valid_design['catalogID'][self.valid_design['catalogID'] != -1],
            robotID=self.valid_design['robotID'][self.valid_design['catalogID'] != -1],
            holeID=self.valid_design['holeID'][self.valid_design['catalogID'] != -1],
            obsWavelength=self.valid_design['obsWavelength'][self.valid_design['catalogID'] != -1],
            carton=self.valid_design['carton_pk'][self.valid_design['catalogID'] != -1],
            observatory=self.observatory,
            targetdb_ver=targetdb_ver,
            instr_pks=None,
            cart_pks=cart_pks,
            fiber_pks=None,
            idtype=self.idtype)
        return

    def design_to_opsdb(self, design):
        """
        Write a validated design to opsdb

        Notes
        -----
        I dont know what needs to be included in opsdb yet, need to look
        at the schema
        """
        return
