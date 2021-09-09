# @Author: Ilija Medan
# @Date: December 29, 2020
# @Filename: fpsdesign.py
# @License: BSD 3-Clause
# @Copyright: Ilija Medan

import numpy as np
import warnings
from astropy.io import fits

import kaiju
import kaiju.robotGrid
# import coordio
from sdssdb.peewee.sdss5db.targetdb import Design, Field, Observatory, Assignment, Instrument, Target, Positioner, CartonToTarget, Carton
import fitsio
from mugatu.exceptions import MugatuError, MugatuWarning
from coordio.utils import radec2wokxy, wokxy2radec
from mugatu.designs_to_targetdb import make_design_assignments_targetdb, make_design_field_targetdb


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

    mode_pk: int
        The pk in targetdb for the observing mode for the design.

    idtype: str
        ID type used for catalogids in design. Must be 'catalogID'
        or 'carton_to_target'.

    catalogids: np.array
        List of catalogids for a manual design in db.

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

    fiberID: np.array
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
    """

    def __init__(self, design_pk, obsTime, racen=None, deccen=None,
                 position_angle=None, observatory=None, mode_pk=None,
                 idtype='carton_to_target', catalogids=None, ra=None, dec=None,
                 pmra=None, pmdec=None, fiberID=None, obsWavelength=None,
                 priority=None, carton_pk=None, category=None, design_file=None,
                 manual_design=False, exp=0,
                 collisionBuffer=2.):
        if idtype != 'catalogID' and idtype != 'carton_to_target':
            raise MugatuError(message='idtype must be catalogID or carton_to_target')
        self.design_pk = design_pk
        self.obsTime = obsTime
        self.design = {}
        # either set field params or pull from db is design
        # is in targetdb
        if manual_design:
            self.racen = racen
            self.deccen = deccen
            self.position_angle = position_angle
            self.observatory = observatory
            self.mode_pk = mode_pk
        else:
            design_field_db = (
                Design.select(Design.pk,
                              # Design.mode_pk,
                              Field.racen,
                              Field.deccen,
                              Field.position_angle,
                              Observatory.label)
                      .join(Field,
                            on=(Design.field == Field.pk))
                      .join(Observatory,
                            on=(Field.observatory == Observatory.pk))
                      .where(Design.pk == self.design_pk))

            self.racen = design_field_db[0].field.racen
            self.deccen = design_field_db[0].field.deccen
            self.position_angle = design_field_db[0].field.position_angle
            self.observatory = design_field_db[0].field.observatory.label
            # no mode right now in targetdb
            # self.mode_pk = design_field_db[0].mode_pk
            self.mode_pk = None

        # should these be catalogids or carton_to_target?
        self.idtype = idtype
        self.catalogids = catalogids
        self.ra = ra
        self.dec = dec
        self.pmra = pmra
        self.pmdec = pmdec
        self.fiberID = fiberID
        self.obsWavelength = obsWavelength
        self.priority = priority
        self.carton_pk = carton_pk
        self.category = category
        self.design_file = design_file
        self.manual_design = manual_design
        self.exp = exp
        self.collisionBuffer = collisionBuffer

        # set dummy value for collision for now
        # this may want to be a input, not sure the standard here
        # initialize robotGrid
        if self.observatory == 'APO':
            self.rg = kaiju.robotGrid.RobotGridAPO(collisionBuffer=self.collisionBuffer,
                                                   stepSize=0.05)
        else:
            self.rg = kaiju.robotGrid.RobotGridLCO(collisionBuffer=self.collisionBuffer,
                                                   stepSize=0.05)
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
        self.design['mode_pk'] = self.mode_pk
        self.design['catalogID'] = np.zeros(500, dtype=np.int64) - 1
        self.design['fiberID'] = np.zeros(500, dtype=np.int64) - 1
        # self.design['wokHoleID'] = np.zeros(500, dtype=np.int64) - 1
        self.design['obsWavelength'] = np.zeros(500, dtype='<U6')
        self.design['priority'] = np.zeros(500, dtype=int) - 1
        self.design['carton_pk'] = np.zeros(500, dtype=int) - 1
        self.design['category'] = np.zeros(500, dtype='<U10')
        self.design['ra'] = np.zeros(500, dtype=float) - 9999.99
        self.design['dec'] = np.zeros(500, dtype=float) - 9999.99
        self.design['pmra'] = np.zeros(500, dtype=float) - 9999.99
        self.design['pmdec'] = np.zeros(500, dtype=float) - 9999.99
        self.design['x'] = np.zeros(500, dtype=float) - 9999.99
        self.design['y'] = np.zeros(500, dtype=float) - 9999.99

        # need to add wokHole to query when in db (not there now)
        # I need to test this when v05 is up, im unsure about Joins
        design_targ_db = (
            Assignment.select(Target.catalogid,
                              Positioner.id,
                              Positioner.pk,
                              Instrument.label,
                              CartonToTarget.priority,
                              Target.ra,
                              Target.dec,
                              CartonToTarget.carton)
                      .join(Positioner,
                            on=(Assignment.positioner_pk == Positioner.pk))
                      .switch(Assignment)
                      .join(Instrument,
                            on=(Assignment.instrument_pk == Instrument.pk))
                      .switch(Assignment)
                      .join(CartonToTarget,
                            on=(Assignment.carton_to_target_pk == CartonToTarget.pk))
                      .join(Target,
                            on=(CartonToTarget.target_pk == Target.pk))
                      .where(Assignment.design_pk == self.design_pk))

        for i in range(len(design_targ_db)):
            # assign to index that corresponds to fiber assignment
            # index should match length of arrays
            pos_id = design_targ_db[i].positioner.id
            if self.idtype == 'carton_to_target':
                self.design['catalogID'][pos_id] = (design_targ_db[i]
                                                    .carton_to_target
                                                    .target.catalogid.catalogid)
            else:
                self.design['catalogID'][pos_id] = (design_targ_db[i]
                                                    .carton_to_target
                                                    .pk)
            self.design['fiberID'][pos_id] = pos_id
            # design['wokHoleID'][i] = design_targ_db[i]
            self.design['obsWavelength'][pos_id] = (design_targ_db[i]
                                                    .instrument.label)
            # catch targets with no assigned priority
            try:
                self.design['priority'][pos_id] = (design_targ_db[i]
                                                   .carton_to_target.priority)
            except AttributeError:
                self.design['priority'][pos_id] = -1
            self.design['carton_pk'][pos_id] = (design_targ_db[i]
                                                .carton_to_target.carton.pk)
            self.design['category'][pos_id] = (design_targ_db[i]
                                                .carton_to_target.carton.category.label)
            self.design['ra'][pos_id] = (design_targ_db[i].carton_to_target
                                                          .target.ra)
            self.design['dec'][pos_id] = (design_targ_db[i].carton_to_target
                                                           .target.dec)
            self.design['pmra'][pos_id] = (design_targ_db[i].carton_to_target
                                                            .target.pmra)
            self.design['pmdec'][pos_id] = (design_targ_db[i].carton_to_target
                                                             .target.pmdec)

        # here convert ra/dec to x/y based on field/time of observation
        # I think I need to add inertial in here at some point, dont see this in targetdb though
        ev = eval("(self.design['ra'] != -9999.99)")
        self.design['x'][ev], self.design['y'][ev], fieldWarn, self.hourAngle, self.positionAngle_coordio = radec2wokxy(ra=self.design['ra'][ev],
                                                                                                                        dec=self.design['dec'][ev],
                                                                                                                        coordEpoch=2457205.9999942128,
                                                                                                                        waveName=np.array(list(map(lambda x:x.title(), self.design['obsWavelength'][ev]))),
                                                                                                                        raCen=self.racen,
                                                                                                                        decCen=self.deccen,
                                                                                                                        obsAngle=self.position_angle,
                                                                                                                        obsSite=self.observatory,
                                                                                                                        obsTime=self.obsTime)

        if np.any(fieldWarn):
            flag = 'Coordio xy coordinates converted should be eyed with suspicion.'
            warnings.warn(flag, MugatuWarning)

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
        self.design['mode_pk'] = self.mode_pk

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
            else:
                self.design['ra'] = self.ra
                self.design['dec'] = self.dec
                self.design['pmra'] = self.pmra
                self.design['pmdec'] = self.pmdec

            # here somehow assign these
            self.design['fiberID'] = self.fiberID
        else:
            # manual design from flat file
            # get header with field info
            head = fits.open(self.design_file)[0].header
            # association between catalogid and instrument
            design_inst = fits.open(self.design_file)[1].data
            # catalogid assignment for each fiber
            design = fits.open(self.design_file)[2].data

            # grab obs info
            self.racen = head['RACEN']
            self.deccen = head['DECCEN']
            self.position_angle = head['PA']
            self.observatory = head['obs'].strip().upper()

            # grab assignment info
            if self.exp == 0:
                roboIDs = design['robotID']
            else:
                roboIDs = design['robotID'][:, self.exp - 1]
            if self.idtype == 'catalogID':
                self.design['catalogID'] = design_inst['catalogid'][roboIDs != -1]
            else:
                self.design['catalogID'] = design_inst['carton_to_target_pk'][roboIDs != -1]
            self.design['ra'] = design_inst['ra'][roboIDs != -1]
            self.design['dec'] = design_inst['dec'][roboIDs != -1]
            self.design['pmra'] = design_inst['pmra'][roboIDs != -1]
            self.design['pmdec'] = design_inst['pmdec'][roboIDs != -1]
            self.design['fiberID'] = roboIDs[roboIDs != -1]
            self.design['obsWavelength'] = design_inst['fiberType'][roboIDs != -1]
            self.design['priority'] = design_inst['priority'][roboIDs != -1]
            # need to change this
            self.design['carton_pk'] = np.arange(0, len(self.design['catalogID']), 1, dtype=int)
            self.design['category'] = design_inst['category'][roboIDs != -1]

        # make empty x,y arrays
        self.design['x'] = np.zeros(len(self.design['catalogID']), dtype=float) - 9999.99
        self.design['y'] = np.zeros(len(self.design['catalogID']), dtype=float) - 9999.99

        # here convert ra/dec to x/y based on field/time of observation
        ev = eval("(self.design['ra'] != -9999.99)")
        self.design['x'][ev], self.design['y'][ev], fieldWarn, self.hourAngle, self.positionAngle_coordio = radec2wokxy(ra=self.design['ra'][ev],
                                                                                                                        dec=self.design['dec'][ev],
                                                                                                                        coordEpoch=2457205.9999942128,  # this is roughly 2015.5, need to ask about this and change it
                                                                                                                        waveName=np.array(list(map(lambda x:x.title(), self.design['obsWavelength'][ev]))),
                                                                                                                        raCen=self.racen,
                                                                                                                        decCen=self.deccen,
                                                                                                                        obsAngle=self.position_angle,
                                                                                                                        obsSite=self.observatory,
                                                                                                                        obsTime=self.obsTime)

        if np.any(fieldWarn):
            flag = 'Coordio xy coordinates converted should be eyed with suspicion.'
            warnings.warn(flag, MugatuWarning)

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
            if self.design['fiberID'][i] != -1:
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
            if self.design['fiberID'][i] > 0:
                try:
                    self.rg.assignRobot2Target(self.design['fiberID'][i],
                                               self.design['catalogID'][i])
                except RuntimeError:
                    # this catches the fact that robot cant be
                    # assigned to given fiber
                    self.targets_unassigned.append(self.design['catalogID'][i])
                    is_unassigned = True

        if is_unassigned:
            flag = 'Some targets could not be assigned to fiber'
            warnings.warn(flag, MugatuWarning)

        return

    def RobotGrid_to_valid_design(self):
        """
        Construct valid design from Kaiju Robotgrid
       """

        # initialize dict for the validated design
        # not using None or nan for no assignments
        # using -1 (for int) and -9999.99 (for float) for None assignment
        self.valid_design['design_pk'] = self.design_pk
        self.valid_design['mode_pk'] = self.mode_pk
        self.valid_design['catalogID'] = np.zeros(500, dtype=np.int64) - 1
        self.valid_design['fiberID'] = np.zeros(500, dtype=np.int64) - 1
        # self.valid_design['wokHoleID'] = np.zeros(500, dtype=np.int64) - 1
        self.valid_design['obsWavelength'] = np.zeros(500, dtype='<U6')
        self.valid_design['priority'] = np.zeros(500, dtype=int) - 1
        self.valid_design['carton_pk'] = np.zeros(500, dtype=int) - 1
        self.valid_design['category'] = np.zeros(500, dtype='<U10')
        self.valid_design['ra'] = np.zeros(500, dtype=float) - 9999.99
        self.valid_design['dec'] = np.zeros(500, dtype=float) - 9999.99
        self.valid_design['x'] = np.zeros(500, dtype=float) - 9999.99
        self.valid_design['y'] = np.zeros(500, dtype=float) - 9999.99

        for i, rid in enumerate(self.rg.robotDict):
            self.valid_design['catalogID'][i] = (self.rg.robotDict[rid]
                                                 .assignedTargetID)
            if self.valid_design['catalogID'][i] != -1:
                self.valid_design['fiberID'][i] = rid
                # is below necessary? i dont know if decollide ever reassigns
                # or just removes
                cond = eval("self.design['catalogID'] == self.valid_design['catalogID'][i]")
                # self.valid_design['wokHoleID'][i] = self.design['wokHoleID'][cond][0]
                self.valid_design['obsWavelength'][i] = self.design['obsWavelength'][cond][0]
                self.valid_design['priority'][i] = self.design['priority'][cond][0]
                self.valid_design['carton_pk'][i] = self.design['carton_pk'][cond][0]
                self.valid_design['category'][i] = self.design['category'][cond][0]
                self.valid_design['ra'][i] = self.design['ra'][cond][0]
                self.valid_design['dec'][i] = self.design['dec'][cond][0]
                self.valid_design['x'][i] = self.design['x'][cond][0]
                self.valid_design['y'][i] = self.design['y'][cond][0]

        return

    def validate_design(self):
        """
        Validate design for deadlocks and collisions using Kaiju.
        """

        # build the design with all needed parameters if it has not been
        # built already (save time by doing this check)
        if not self.design_built:
            if self.manual_design:
                self.build_design_manual()
            else:
                self.build_design_db()

        # construct the Kaiju robotGrid
        self.design_to_RobotGrid()

        # validate the design

        # decollide the unassigned robots
        for robotID in self.rg.robotDict:
            if(self.rg.robotDict[robotID].isAssigned() == False):
                self.rg.decollideRobot(robotID)
                # need to still set alpha/beta
                # robot = self.rg.getRobot(robotID)
                # robot.setDestinationAlphaBeta(0, 180)

        # de-collide the grid if collisions exist
        # and check for targets removed
        if self.rg.getNCollisions() > 0:
            self.rg.decollideGrid()

            # check if de-collision was successful
            if not self.rg.getNCollisions() == 0:
                raise MugatuError(message='Kaiju decollideGrid failed')

            flag = 'Some targets removed from design due to collisions'
            warnings.warn(flag, MugatuWarning)

            # grab all of the targets removed due to collisions
            for i in self.rg.robotDict:
                fiber_idx = np.where(self.design['fiberID'] == i)[0]
                if (self.rg.robotDict[i].assignedTargetID == -1
                    and len(fiber_idx) > 0):
                    targ_remove = self.design['catalogID'][fiber_idx[0]]
                    if targ_remove not in self.targets_unassigned:
                        self.targets_collided.append(targ_remove)

        # generate paths
        # self.rg.pathGen()
        # if self.rg.didFail:
        #     raise MugatuError(message='Kaiju pathGen failed')

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
                                   observatory=self.observatory)

        # create dictonary for unique carton pks
        cart_pks = {}
        targetdb_ver = {}
        for pk in np.unique(self.valid_design['carton_pk'][self.valid_design['catalogID'] != -1]):
            cart_pks[pk] = pk
            targetdb_ver[pk] = Carton.get(pk).version_pk

        # add the design to targetdb
        make_design_assignments_targetdb(targetdb_ver=targetdb_ver,
                                         plan='manual',
                                         fieldid=fieldid,
                                         exposure=exposure,
                                         catalogID=self.valid_design['catalogID'][self.valid_design['catalogID'] != -1],
                                         fiberID=self.valid_design['fiberID'][self.valid_design['catalogID'] != -1],
                                         obsWavelength=self.valid_design['obsWavelength'][self.valid_design['catalogID'] != -1],
                                         carton=self.valid_design['carton_pk'][self.valid_design['catalogID'] != -1],
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
