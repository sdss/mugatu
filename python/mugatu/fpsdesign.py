# @Author: Ilija Medan
# @Date: December 22, 2020
# @Filename: fpsdesign.py
# @License: BSD 3-Clause
# @Copyright: Ilija Medan

import numpy as np

import kaiju
import coordio
import sdssdb


class FPSDesign(object):
    """
    Parameters:
    -----------

    racen: np.float64
        Right Ascension of center of the field

    deccen: np.float64
        Declination of center of the field

    PA: np.float64
        Position angle of the field

    HA: np.float64
        Hour angle of the observation

    observatory: str
        Observatory where observation is taking place

    design_pk: int
        The pk of the design as it appears in targetdb

    target_ids: np.array
        List of target ids (or catalog ids?) for a manual
        design in db

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

    design_to_opsdb(): write a validated design to opsdb

    """

    def _init_(self, racen, deccen, PA, HA, observatory,
               design_pk, target_ids=None,
               obsWavelength=None, priority=None, design_file=None,
               manual_design=False):
        self.racen = racen
        self.deccen = deccen
        self.PA = PA
        self.HA = HA
        self.observatory = observatory
        self.design_pk = design_pk
        self.target_ids = target_ids
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
        This method will generally be triggered when manual_design=False

        """

        design = {}
        design['design_pk'] = self.design_pk
        design['targetID'] = np.array(500, dtype=np.int64)
        design['fiberID'] = np.array(500, dtype=np.int64)
        design['wokHoleID'] = np.array(500, dtype=np.int64)
        design['obsWavelength'] = np.array(500, dtype=str)
        design['priority'] = np.array(500, dtype=int)
        design['x'] = np.array(500, dtype=float)
        design['y'] = np.array(500, dtype=float)

        # something here that queries the db to pull targetID, fiberID,
        # wokHoleID and obsWavelength based on design_pk

        # here grab the ra/decs from targetdb based on targetID
        ra, dec = something(design['targetID'])

        # here convert ra/dec to x/y based on field/HA observation
        design['x'], design['y'] = self.radec_to_xy(ra, dec)

        return design

    def build_design_manual(self):
        """
        compile the parameters for a design

        Notes:
        ------
        This function will be able to do two things:

        (1) Build a manual design by collecting ra/decs of
        input targets and randomally (or algorthmically?) assign
        fiberIDs (and WokHoleIDs) based on priority

        (2) Build a design from a design_pk in targetdb. This will be just
        collecting information that already exists in the db.

        """

        design = {}
        design['design_pk'] = self.design_pk
        design['targetID'] = np.array(500, dtype=np.int64)
        design['fiberID'] = np.array(500, dtype=np.int64)
        design['wokHoleID'] = np.array(500, dtype=np.int64)
        design['obsWavelength'] = np.array(500, dtype=str)
        design['priority'] = np.array(500, dtype=int)
        design['x'] = np.array(500, dtype=float)
        design['y'] = np.array(500, dtype=float)

        if self.design_file is None:
            # manual design with targets in targetdv
            design['targetID'] = self.target_ids
            design['obsWavelength'] = self.obsWavelength
            design['priority'] = self.priority

            # here somehow assign these
            design['fiberID'], design['wokHoleID'] = something()

            # here grab the ra/decs from targetdb based on targetID
            ra, dec = something(design['targetID'])
        else:
            # manual design from flat file
            man_des = load(self.design_file)

            # now need to load all the params to dictonary

        # here convert ra/dec to x/y based on field/HA observation
        design['x'], design['y'] = self.radec_to_xy(ra, dec)

        return design

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
