# @Author: Ilija Medan
# @Date: February 4, 2021
# @Filename: fpsdesign.py
# @License: BSD 3-Clause
# @Copyright: Ilija Medan

import numpy as np
import warnings
from scipy.spatial import cKDTree

from sdssdb.peewee.sdss5db.targetdb import Carton, CartonToTarget, Magnitude, Mode, Assignment, Instrument, Positioner

from mugatu.exceptions import MugatuError, MugatuWarning


def mode_validation(design):
    """
    Validate that a design meets the criteria of
    an assigned observing mode

    Parameters
    ----------
    design: dict
        A design dictonary made from a FPSDesign object

    Outputs
    -------
    didPass: boolean
        Did the design pass the mode validation
    """

    # grab the design mode info
    des_mode = Mode.get(design['mode_pk'])

    design_db = (
        Assignment.select(Positioner.id,
                          Instrument.label,
                          Magnitude.gaia_g,
                          Carton.carton,
                          Carton.program)
                  .join(Positioner,
                        on=(Assignment.positioner_pk == Positioner.pk))
                  .switch(Assignment)
                  .join(Instrument,
                        on=(Assignment.instrument_pk == Instrument.pk))
                  .switch(Assignment)
                  .join(CartonToTarget,
                        on=(Assignment.carton_to_target_pk == CartonToTarget.pk))
                  .join(Magnitude,
                        on=(CartonToTarget.pk == Magnitiude.carton_to_target_pk))
                  .switch(CartonToTarget)
                  .join(Carton,
                        on=(CartonToTarget.carton_pk == Carton.pkk))
                  .where(Assignment.design_pk == design['design_pk']))

    pos_id, instr_label, gaia_g, carton, program = map(np.array, zip(*np.array(design_db.tuples())))

    # setup booleans for if design passes each mode param
    # default to false, which is fail
    skybrightness = False
    boss_skies = False
    boss_stds = False
    boss_bright_limit = False
    apogee_skies = False
    apogee_stds = False
    apogee_bright_limit = False

    # check the boss modes
    if len(program[(program == 'ops_sky') & (instr_label == 'BOSS')]) >= des_mode.boss_skies:
        boss_skies = True

    if len(program[(program == 'ops_std') & (instr_label == 'BOSS')]) >= des_mode.boss_stds:
        boss_stds = True

    gaia_g_boss = gaia_g[instr_label == 'BOSS']
    if np.any(gaia_g_boss[gaia_g_boss <= des_mode.boss_bright_limit]):
        boss_bright_limit = True

    # check the apogee modes
    if len(program[(program == 'ops_sky') & (instr_label == 'APOGEE')]) >= des_mode.apogee_skies:
        apogee_skies = True

    if len(program[(program == 'ops_std') & (instr_label == 'BOSS')]) >= des_mode.apogee_stds:
        apogee_stds = True

    gaia_g_apogee = gaia_g[instr_label == 'APOGEE']
    if np.any(gaia_g_apogee[gaia_g_apogee <= des_mode.apogee_bright_limit]):
        apogee_bright_limit = True

    # spit out an error if any of the above failed
    mode_checks = np.array([skybrightness, boss_skies,
                            boss_stds, boss_bright_limit,
                            apogee_skies, apogee_stds,
                            apogee_bright_limit])
    if np.any(~mode_checks):
        message = 'Design not compliant with observing mode. Design fails: '
        if not skybrightness:
            message += 'skybrightness '
        if not boss_skies:
            message += 'boss_skies '
        if not boss_stds:
            message += 'boss_stds '
        if not boss_bright_limit:
            message += 'boss_bright_limit '
        if not apogee_skies:
            message += 'apogee_skies '
        if not apogee_stds:
            message += 'apogee_stds '
        if not apogee_bright_limit:
            message += 'apogee_bright_limit '
        warnings.warn(message, MugatuWarning)
        return False
    else:
        return True


def bright_neighbors(fps_design, bright_limit, bright_neighbor_limit):
    """
    looks for bright neighbors some distance away from all of
    the assigned targets in a design

    Parameters
    ----------
    fps_design: object
        FPSDesign object

    bright_limit: float
        Bright limit for the design in mags. This will search
        for Gaia(?) mags brighter than this

    bright_neighbor_limit: float
        the max distance a fiber should be away from a bright target.
        this is provided in degrees.

    Outputs
    -------
    fiber_bright: np.array
        Array of booleens to signify if fiber close to a bright neighbor.
        True means too close to bright neighbor.
    """

    # first thing here will be to query bright objects
    # dont really know where this will come from, so leave
    # blank for now
    query_bright = something()

    ra_bright, dec_bright, mag_bright = map(np.array, zip(*np.array(query_bright.tuples())))

    # make trees for bright targets and designs
    tree_bright = cKDTree(np.column_stack((ra_bright, dec_bright)))
    tree_design = cKDTree(np.column_stack((fps_design.design['ra'],
                                           fps_design.design['dec'])))

    # get idnexes of fibers close to right targets
    index_bright = tree_design.query_ball_tree(tree_bright,
                                               r=bright_neighbor_limit)

    # set fibers close enough to bright objects to True
    fiber_bright = np.array([False] * len(index_bright))
    for i in range(len(index_bright)):
        if len(index_bright[i]) > 0.:
            fiber_bright[i] = True

    # don't count any fibers unassigned
    fiber_bright[fps_design.design['fiberID'] == -1] = False

    return fiber_bright
