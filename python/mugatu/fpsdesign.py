# @Author: Ilija Medan
# @Date: December 22, 2020
# @Filename: fpsdesign.py
# @License: BSD 3-Clause
# @Copyright: Ilija Medan

import numpy as np

import kaiju
import coordio
from sdssdb.peewee.sdss5db.targetdb import Design, Field, Observatory, Assignment, Instrument, Target, Positioner, CartonToTarget


class FPSDesign(object):
    """
    Parameters:
    -----------

    design_pk: int
        The pk of the design as it appears in targetdb

    hour_angle: np.float64
        Hour angle of the observation

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

    catalogids: np.array
        List of catalogids for a manual design in db

    obsWavelength: np.array
        Wavelength of observation for each fiber in the
        manual design in db

    priority: np.array
        Priorties for targets in a manual design in db

    design_file: str
        Flate file with a manual design not in the db. The
        file should have the columns: RA, DEC, obsWavelength,
        priority, ???

    manual_design: boolean
        Boolean if the design being validated is manual
        (manual_design=True) or in targetdb (manual_design=False)

    Attributes:
    -----------

    design: dictonary
        contains all of the design inputs for Kaiju needed to validate
        a design

    robogrid: Kaiju RobotGrid object
        MIKE NOTE: I imay just want to inherit this? need to think about
        this when I get to writing Kaiju part of code.

    valid_design: dictonary
        same form as design dictonary, except this design has been
        validated by Kaiju and collision/delock assignments removed

    Methods:
    --------

    radec_to_xy(): convert ra/dec to x/y using field information and
        coordio

    build_design_db(): compile the parameters for a design that exists in
        targetdb. This method will generally be triggered when
        manual_design=False

    build_design_manual(): compile the parameters for a design that has
        been manually created. This method with work for manual designs that
        are flat files (not in targetdb) and those targets in targetdb

    design_to_RobotGrid(): construct the Kaiju RobotGrid object from a design

    RobotGrid_to_design(): construct design from RobotGrid object

    validate_design(): calls appropriate functions to carry out creating
        the design object (manual or from db), converting coordinates,
        creating Kaiju.RobotGrid, assigning targets, checking for
        collisions and generating paths for design

    design_to_targetdb(): write a design to targetdb

    design_to_opsdb(): write a validated design to opsdb. NOTE: Need to
        address here that "designs" do not go to opsdb, "configurations"
        do.

    """

    def _init_(self, design_pk, hour_angle, racen=None, deccen=None,
               position_angle=None, observatory=None, catalogids=None,
               obsWavelength=None, priority=None, design_file=None,
               manual_design=False):
        self.design_pk = design_pk
        self.hour_angle = hour_angle
        self.design = {}
        if manual_design:
            self.racen = racen
            self.deccen = deccen
            self.position_angle = position_angle
            self.observatory = observatory
        else:
            design_field_db = (Design.select(Design.pk,
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
            self.observatory = design_field_db[0].observatory.label
        self.catalogids = catalogids
        self.obsWavelength = obsWavelength
        self.priority = priority
        self.design_file = design_file
        self.manual_design = manual_design

    def radec_to_xy(self, ra, dec):
        """
        convert ra/dec to x/y using field information and
        coordio

        Notes:
        ------
        This function will mainly rely on calling coordio for the
        conversions

        """
        return x, y

    def build_design_db(self):
        """
        compile the parameters for a design that exists in targetdb

        Notes:
        ------
        This method will build a design from the db if manual_design=False

        """

        self.design['design_pk'] = self.design_pk
        self.design['catalogID'] = np.array(500, dtype=np.int64)
        self.design['fiberID'] = np.array(500, dtype=np.int64)
        self.design['wokHoleID'] = np.array(500, dtype=np.int64)
        self.design['obsWavelength'] = np.array(500, dtype=str)
        self.design['priority'] = np.array(500, dtype=int)
        self.design['ra'] = np.array(500, dtype=float)
        self.design['dec'] = np.array(500, dtype=float)
        self.design['x'] = np.array(500, dtype=float)
        self.design['y'] = np.array(500, dtype=float)

        # need to add wokHole to query when in db (not there now)
        design_targ_db = (Assignment.select(Target.catalogid,
                                            Positioner.id,
                                            Instrument.label,
                                            CartonToTarget.priority,
                                            Target.ra,
                                            Target.dec)
                                    .join(Target,
                                          on=(Assignment.target_pk == Target.pk))
                                    .join(Positioner,
                                          on=(Assignment.positioner_pk == Positioner.pk))
                                    .join(Instrument,
                                          on=(Assignment.instrument_pk == Instrument.pk))
                                    .join(CartonToTarget,
                                          on=(Target.pj == CartonToTarget.target_pk))
                                    .where(Assignment.design_pk == self.design_pk))

        for i in range(len(design_targ_db)):
            self.design['catalogID'][i] = design_targ_db[i].target.catalogid
            self.design['fiberID'][i] = design_targ_db[i].positioner.id
            # design['wokHoleID'][i] = design_targ_db[i]
            self.design['obsWavelength'][i] = design_targ_db[i].instrument.label
            self.design['priority'][i] = design_targ_db[i].cartontotarget.priority
            self.design['ra'][i] = design_targ_db[i].target.ra
            self.design['dec'][i] = design_targ_db[i].target.dec

        # here convert ra/dec to x/y based on field/HA observation
        self.design['x'], self.design['y'] = self.radec_to_xy()

        return

    def build_design_manual(self):
        """
        compile the parameters for a manual design

        Notes:
        ------
        This function creates a manual design whether it is from
        user inputted catalogids, or if it is a flat file list
        of coordinates (and fiber assignments?), if manual_design=True

        """

        self.design['design_pk'] = self.design_pk
        self.design['catalogID'] = np.array(500, dtype=np.int64)
        self.design['fiberID'] = np.array(500, dtype=np.int64)
        self.design['wokHoleID'] = np.array(500, dtype=np.int64)
        self.design['obsWavelength'] = np.array(500, dtype=str)
        self.design['priority'] = np.array(500, dtype=int)
        self.design['ra'] = np.array(500, dtype=float)
        self.design['dec'] = np.array(500, dtype=float)
        self.design['x'] = np.array(500, dtype=float)
        self.design['y'] = np.array(500, dtype=float)

        if self.design_file is None:
            # manual design with targets in targetdv
            self.design['catalogID'] = self.catalogids
            self.design['obsWavelength'] = self.obsWavelength
            self.design['priority'] = self.priority

            for i in range(len(self.design['catalogID'])):
                targ_db = Target.get(catalogid=self.design['catalogID'][i])
                self.design['ra'][i] = targ_db.ra
                self.design['dec'][i] = targ_db.dec

            # here somehow assign these
            self.design['fiberID'], self.design['wokHoleID'] = something()
        else:
            # manual design from flat file
            man_des = load(self.design_file)

            # now need to load all the params to dictonary

        # here convert ra/dec to x/y based on field/HA observation
        self.design['x'], self.design['y'] = self.radec_to_xy()

        return

    def design_to_RobotGrid(self, design):
        """
        contruct a Kaiju RobotGrid object

        Notes:
        -----
        I don't really know how this will work yet

       """

        return kaiju.RobotGrid

    def RobotGrid_to_design(self, robogrid):
        """
        construct design from RobotGrid object

        Notes:
        -----
        I don't really know how this will work yet

       """

        return design

    def validate_design(self):
        """
        validate the design

        Notes:
        ------
        This function will call all the necessary steps to create the design
        file for Kaiju and validate the design (following steps 3-7 in Conor's
        outline). This will utlize multiple functions that will be above.

        """

        # build the design with all needed parameters
        design = self.build_design()

        # construct the Kaiju robotGrid
        robogrid = self.design_to_RobotGrid(design=design)

        # validate the design
        # I don't know what would go here yet, but I imagine it would be calls
        # to Kaiju that would look for collisions and generate the paths for
        # the design that was built into the above robogrid

        # I imagine that the above step would manipulate the robogrid based on
        # collisions and deadlocks, so the below would take these Kaiju results
        # and put them back as a design object which will be the output of the
        # function
        valid_design = self.RobotGrid_to_design(robogrid=robogrid)

        # another note to add to this above design dictonary, is how will this
        # include paths? This seems important for manual designs and to send
        # to Jaeger

        return valid_design

    def design_to_targetdb(self, design, version):
        """
        write a design to targetdb

        Notes:
        -----
        version above should allow for manual designs to be added under
        version = 'manual' to seperate them from robostrategy ingested designs
        """
        return

    def design_to_opsdb(self, design):
        """
        write a validated design to opsdb

        Notes:
        -----
        I dont know what needs to be included in opsdb yet, need to look
        at the schema
        """
        return
