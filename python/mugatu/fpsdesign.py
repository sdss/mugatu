# @Author: Ilija Medan
# @Date: December 29, 2020
# @Filename: fpsdesign.py
# @License: BSD 3-Clause
# @Copyright: Ilija Medan

import numpy as np
import warnings

import kaiju
import kaiju.robotGrid
# import coordio
from sdssdb.peewee.sdss5db.targetdb import Design, Field, Observatory, Assignment, Instrument, Target, Positioner, CartonToTarget
import fitsio
from mugatu.exceptions import MugatuError, MugatuWarning
from coordio.utils import radec2wokxy, wokxy2radec


class FPSDesign(object):
    """
    Parameters
    ----------

    design_pk: int
        The pk of the design as it appears in targetdb

    obsTime: np.float64
        Julian date of the observation

    design: dict
        Dictonary that contains all parameters related to the design.
        This is initiallized as empty and is populated using either
        build_design_db() or build_design_manual().

    racen: np.float64
        Right Ascension of center of the field

    deccen: np.float64
        Declination of center of the field

    position_angle: np.float64
        Position angle of the field

    observatory: str
        Observatory where observation is taking place

    mode_pk: int
        The pk in targetdb for the observing mode for the design

    catalogids: np.array
        List of catalogids for a manual design in db.
        Length of array must be n=500.

    ra: np.array
        List of right ascensions that correspond to catalogids.
        If list not provided, but catalogids for a manual design
        are provided, right ascensions will be pulled from
        targetdb. Length of array must be n=500.

    dec: np.array
        List of declinations that correspond to catalogids.
        If list not provided, but catalogids for a manual design
        are provided, declinations will be pulled from
        targetdb. Length of array must be n=500.

    fiberID: np.array
        Fiber assignement for each catalogid target in the
        manual design. Length of array must be n=500.

    obsWavelength: np.array
        Wavelength of observation for each fiber in the
        manual design in db. Length of array must be n=500.

    priority: np.array
        Priorties for targets in a manual design in db.
        Length of array must be n=500.

    design_file: str
        Flate file with a manual design not in the db. The
        file should have the columns: RA, DEC, obsWavelength,
        priority, ???

    manual_design: boolean
        Boolean if the design being validated is manual
        (manual_design=True) or in targetdb (manual_design=False)

    Attributes
    ----------

    design: dictonary
        contains all of the design inputs for Kaiju needed to validate
        a design

    rg: Kaiju RobotGrid object
        MIKE NOTE: I imay just want to inherit this? need to think about
        this when I get to writing Kaiju part of code.

    targets_unassigned: list
        catalogid of targets that could not be assigned due to assigned
        robot not being able to reach the assigned target

    targets_collided: list
        catalogid of targets that could not be assigned due to assigned
        robot to assigned target resulting in a collision

    valid_design: dictonary
        same form as design dictonary, except this design has been
        validated by Kaiju and collision/delock assignments removed

    design_built: booleen
        keeps track if design from db or file has already been built. If so,
        this will save time on validating the design by not requiring it
        to be rebuilt.

    hourAngle: float
        Hour angle of field center from coordio radec2wokxy

    positionAngle_coordio: float
        position angle of field center from coordio radec2wokxy

    Methods
    -------

    build_design_db(): compile the parameters for a design that exists in
        targetdb. This method will generally be triggered when
        manual_design=False

    build_design_manual(): compile the parameters for a design that has
        been manually created. This method with work for manual designs that
        are flat files (not in targetdb) and those targets in targetdb

    design_to_RobotGrid(): construct the Kaiju RobotGrid object from a design

    RobotGrid_to_valid_design(): construct valid design from validated
        RobotGrid object

    validate_design(): calls appropriate functions to carry out creating
        the design object (manual or from db), converting coordinates,
        creating Kaiju.RobotGrid, assigning targets, checking for
        collisions and generating paths for design

    design_to_targetdb(): write a design to targetdb

    design_to_opsdb(): write a validated design to opsdb. NOTE: Need to
        address here that "designs" do not go to opsdb, "configurations"
        do.

    """

    def __init__(self, design_pk, obsTime, racen=None, deccen=None,
                 position_angle=None, observatory=None, mode_pk=None,
                 catalogids=None, ra=None, dec=None, fiberID=None,
                 obsWavelength=None, priority=None, design_file=None,
                 manual_design=False):
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
                              Design.mode_pk,
                              Field.racen,
                              Field.deccen,
                              Field.position_angle,
                              Observatory.label)
                      .join(Field,
                            on=(Design.field_pk == Field.pk))
                      .join(Observatory,
                            on=(Field.observatory_pk == Observatory.pk))
                      .where(Design.pk == self.design_pk))

            self.racen = design_field_db[0].field.racen
            self.deccen = design_field_db[0].field.deccen
            self.position_angle = design_field_db[0].field.position_angle
            self.observatory = design_field_db[0].field.observatory.label
            self.mode_pk = design_field_db[0].mode_pk
        # should these be catalogids or carton_to_target?
        # either way, I think I need to specify carton info
        # for checking modes here
        self.catalogids = catalogids
        self.ra = ra
        self.dec = dec
        self.fiberID = fiberID
        self.obsWavelength = obsWavelength
        self.priority = priority
        self.design_file = design_file
        self.manual_design = manual_design

        # set dummy value for collision for now
        # this may want to be a input, not sure the standard here
        # initialize robotGrid
        self.rg = kaiju.robotGrid.RobotGridFilledHex()
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
        compile the parameters for a design that exists in targetdb

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
        self.design['ra'] = np.zeros(500, dtype=float) - 9999.99
        self.design['dec'] = np.zeros(500, dtype=float) - 9999.99
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
                              Target.dec)
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
            self.design['catalogID'][pos_id] = (design_targ_db[i]
                                                .cartontotarget
                                                .target.catalogid.catalogid)
            self.design['fiberID'][pos_id] = pos_id
            # design['wokHoleID'][i] = design_targ_db[i]
            self.design['obsWavelength'][pos_id] = (design_targ_db[i]
                                                    .instrument.label)
            # catch targets with no assigned priority
            try:
                self.design['priority'][pos_id] = (design_targ_db[i]
                                                   .cartontotarget.priority)
            except AttributeError:
                self.design['priority'][pos_id] = -1
            self.design['ra'][pos_id] = (design_targ_db[i].cartontotarget
                                                          .target.ra)
            self.design['dec'][pos_id] = (design_targ_db[i].cartontotarget
                                                           .target.dec)

        # here convert ra/dec to x/y based on field/time of observation
        ev = eval("(self.design['ra'] != -9999.99)")
        self.design['x'][ev], self.design['y'][ev], fieldWarn, self.hourAngle, self.positionAngle_coordio = radec2wokxy(ra=self.design['ra'][ev],
                                                                                                                        dec=self.design['dec'][ev],
                                                                                                                        coordEpoch=np.array([2457174] * len(self.design['ra'][ev])),
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
        compile the parameters for a manual design

        Notes
        -----
        This function creates a manual design whether it is from
        user inputted catalogids, or if it is a flat file list
        of coordinates (and fiber assignments?), if manual_design=True

        """

        if not self.manual_design:
            flag = 'Building manual design with manual_design=False flag'
            warnings.warn(flag, MugatuWarning)

        # initialize dict for the design
        # not using None or nan for no assignments
        # using -1 (for int) and -9999.99 (for float) for None assignment
        self.design['design_pk'] = self.design_pk
        self.design['mode_pk'] = self.mode_pk
        self.design['catalogID'] = np.zeros(500, dtype=np.int64) - 1
        self.design['fiberID'] = np.zeros(500, dtype=np.int64) - 1
        # self.design['wokHoleID'] = np.zeros(500, dtype=np.int64) - 1
        self.design['obsWavelength'] = np.zeros(500, dtype=str)
        self.design['priority'] = np.zeros(500, dtype=int) - 1
        self.design['ra'] = np.zeros(500, dtype=float) - 9999.99
        self.design['dec'] = np.zeros(500, dtype=float) - 9999.99
        self.design['x'] = np.zeros(500, dtype=float) - 9999.99
        self.design['y'] = np.zeros(500, dtype=float) - 9999.99

        if self.design_file is None:
            # manual design with targets in targetdb
            self.design['catalogID'] = self.catalogids
            self.design['obsWavelength'] = self.obsWavelength
            self.design['priority'] = self.priority

            if self.ra is None:
                for i in range(len(self.design['catalogID'])):
                    targ_db = Target.get(catalogid=self.design['catalogID'][i])
                    self.design['ra'][i] = targ_db.ra
                    self.design['dec'][i] = targ_db.dec
            else:
                self.design['ra'] = self.ra
                self.design['dec'] = self.dec

            # here somehow assign these
            self.design['fiberID'] = self.fiberID
        else:
            # manual design from flat file
            man_des = fitsio.read(self.design_file)

            self.design['catalogID'] = man_des['catalogID']
            self.design['fiberID'] = man_des['fiberID']
            self.design['obsWavelength'] = man_des['obsWavelength']
            self.design['priority'] = man_des['priority']
            self.design['ra'] = man_des['ra']
            self.design['dec'] = man_des['dec']

        # here convert ra/dec to x/y based on field/time of observation
        ev = eval("(self.design['ra'] != -9999.99)")
        self.design['x'][ev], self.design['y'][ev], fieldWarn, self.hourAngle, self.positionAngle_coordio = radec2wokxy(ra=self.design['ra'][ev],
                                                                                                                        dec=self.design['dec'][ev],
                                                                                                                        coordEpoch=np.array([2457174] * len(self.design['ra'][ev])),  # this is roughly 2015.5, need to ask about this and change it
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
        contruct a Kaiju RobotGrid object

        Notes
        -----
        Adds targets to the kaiju.robotGrid.RobotGridFilledHex and
        assigns them to fibers based on the design dict parameters

       """

        is_unassigned = False

        for i in range(len(self.design['x'])):
            if self.design['fiberID'][i] != -1:
                if self.design['obsWavelength'][i] == 'BOSS':
                    self.rg.addTarget(targetID=self.design['catalogID'][i],
                                      x=self.design['x'][i],
                                      y=self.design['y'][i],
                                      priority=self.design['priority'][i],
                                      fiberType=kaiju.BossFiber)
                else:
                    self.rg.addTarget(targetID=self.design['catalogID'][i],
                                      x=self.design['x'][i],
                                      y=self.design['y'][i],
                                      priority=self.design['priority'][i],
                                      fiberType=kaiju.ApogeeFiber)
        for i in range(len(self.design['x'])):
            if self.design['fiberID'][i] != -1:
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
        construct design from RobotGrid object

        Notes
        -----
        I don't really know how this will work yet

       """

        # initialize dict for the validated design
        # not using None or nan for no assignments
        # using -1 (for int) and -9999.99 (for float) for None assignment
        self.valid_design['design_pk'] = self.design_pk
        self.valid_design['mode_pk'] = self.mode_pk
        self.valid_design['catalogID'] = np.zeros(500, dtype=np.int64) - 1
        self.valid_design['fiberID'] = np.zeros(500, dtype=np.int64) - 1
        # self.valid_design['wokHoleID'] = np.zeros(500, dtype=np.int64) - 1
        self.valid_design['obsWavelength'] = np.zeros(500, dtype=str)
        self.valid_design['priority'] = np.zeros(500, dtype=int) - 1
        self.valid_design['ra'] = np.zeros(500, dtype=float) - 9999.99
        self.valid_design['dec'] = np.zeros(500, dtype=float) - 9999.99
        self.valid_design['x'] = np.zeros(500, dtype=float) - 9999.99
        self.valid_design['y'] = np.zeros(500, dtype=float) - 9999.99

        for i in self.rg.robotDict:
            self.valid_design['catalogID'][i] = (self.rg.robotDict[i]
                                                 .assignedTargetID)
            if self.valid_design['catalogID'][i] != -1:
                self.valid_design['fiberID'][i] = i
                # is below necessary? i dont know if decollide ever reassigns
                # or just removes
                cond = eval("self.design['catalogID'] == self.valid_design['catalogID'][i]")
                # self.valid_design['wokHoleID'][i] = self.design['wokHoleID'][cond][0]
                self.valid_design['obsWavelength'][i] = self.design['obsWavelength'][cond][0]
                self.valid_design['priority'][i] = self.design['priority'][cond][0]
                self.valid_design['ra'][i] = self.design['ra'][cond][0]
                self.valid_design['dec'][i] = self.design['dec'][cond][0]
                self.valid_design['x'][i] = self.design['x'][cond][0]
                self.valid_design['y'][i] = self.design['y'][cond][0]

        return

    def validate_design(self):
        """
        validate the design

        Notes
        -----
        This function will call all the necessary steps to create the design
        file for Kaiju and validate the design (following steps 3-7 in Conor's
        outline). This will utlize multiple functions that will be above.

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
        self.rg.pathGen()
        if self.rg.didFail:
            raise MugatuError(message='Kaiju pathGen failed')

        # I imagine that the above step would manipulate the robogrid based on
        # collisions and deadlocks, so the below would take these Kaiju results
        # and put them back as a design object which will be the output of the
        # function
        self.RobotGrid_to_valid_design()

        # another note to add to this above design dictonary, is how will this
        # include paths? This seems important for manual designs and to send
        # to Jaeger

        return

    def design_to_targetdb(self, design, version):
        """
        write a design to targetdb

        Notes
        -----
        version above should allow for manual designs to be added under
        version = 'manual' to seperate them from robostrategy ingested designs
        """
        return

    def design_to_opsdb(self, design):
        """
        write a validated design to opsdb

        Notes
        -----
        I dont know what needs to be included in opsdb yet, need to look
        at the schema
        """
        return
