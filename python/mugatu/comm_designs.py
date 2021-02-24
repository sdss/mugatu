# @Author: Ilija Medan
# @Date: February 4, 2021
# @Filename: comm_designs.py
# @License: BSD 3-Clause
# @Copyright: Ilija Medan

import numpy as np
import warnings
import random

import kaiju
import kaiju.robotGrid
from sdssdb.peewee.sdss5db.targetdb import Carton, Target, CartonToTarget
from mugatu.exceptions import MugatuError, MugatuWarning
from mugatu.fpsdesign import FPSDesign
from coordio.utils import radec2wokxy


def all_sky_design(racen, deccen, position_angle, observatory,
                   obsTime, n_sky_apogee, n_sky_boss):
    """
    creates a design that consists of all sky targets

    Parameters:
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
                                Target.dec)
                        .join(CartonToTarget)
                        .join(Carton)
                        .where((Carton.carton == 'ops_sky_apogee') & 
                               (Target.cone_search(racen, deccen, r_search))))

    if len(apogee_sky) < n_sky_apogee:
        raise MugatuError(message='Not enough apogee skies in field')

    boss_sky = (Target.select(Target.catalogid,
                              Target.ra,
                              Target.dec)
                      .join(CartonToTarget)
                      .join(Carton)
                      .where((Carton.carton == 'ops_sky_boss') & 
                             (Target.cone_search(racen, deccen, r_search))))

    if len(boss_sky) < n_sky_boss:
        raise MugatuError(message='Not enough boss skies in field')

    catalogid_apogee, ra_apogee, dec_apogee = map(list, zip(*list(apogee_sky.tuples())))
    catalogid_apogee = np.array(catalogid_apogee, dtype=int)
    ra_apogee = np.array(ra_apogee)
    dec_apogee = np.array(dec_apogee)
    x_apogee, y_apogee, fieldWarn, HA, PA_coordio = radec2wokxy(ra=ra_apogee,
                                                                dec=dec_apogee,
                                                                coordEpoch=np.array([2457174] * len(ra_apogee)),
                                                                waveName=np.array(list(map(lambda x:x.title(), ['APOGEE'] * len(ra_apogee)))),
                                                                raCen=racen,
                                                                decCen=deccen,
                                                                obsAngle=position_angle,
                                                                obsSite=observatory,
                                                                obsTime=obsTime)


    catalogid_boss, ra_boss, dec_boss = map(list, zip(*list(boss_sky.tuples())))
    catalogid_boss = np.array(catalogid_boss, dtype=int)
    ra_boss = np.array(ra_boss)
    dec_boss = np.array(dec_boss)
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

    # make a robot grid
    rg = kaiju.robotGrid.RobotGridFilledHex()
    for rID in rg.robotDict:
        robot = rg.getRobot(rID)
        robot.setXYUniform()

    # add all targets to grid
    for i in range(len(ra_apogee)):
        rg.addTarget(targetID=catalogid_apogee[i],
                     x=x_apogee[i],
                     y=y_apogee[i],
                     priority=1,
                     fiberType=kaiju.ApogeeFiber)

    for i in range(len(ra_boss)):
        rg.addTarget(targetID=catalogid_boss[i],
                     x=x_boss[i],
                     y=y_boss[i],
                     priority=1,
                     fiberType=kaiju.BossFiber)

    # this assignment below basically works
    # issues is idk if assigning boss and then apogee makes sense
    # also, a few tests show that a lot of collisions can still happen
    # may need to be a reassign step after this?

    # assign BOSS first as there are less of them
    fiberID_boss = np.zeros(len(catalogid_boss), dtype=int) - 9999
    n_boss_assigned = 0
    for i in range(len(rg.robotDict)):
        if n_boss_assigned == n_sky_boss:
            break
        valid_targs = np.isin(rg.robotDict[i].validTargetIDs,
                              catalogid_boss[fiberID_boss == -9999])
        valid_targs_arr = np.array(rg.robotDict[i].validTargetIDs, dtype=int)
        if np.any(valid_targs):
            fiberID_boss[catalogid_boss == valid_targs_arr[valid_targs][0]] = i
            n_boss_assigned += 1

    if n_boss_assigned < n_sky_boss:
        flag = 'Not enough boss skies could be assigned'
        warnings.warn(flag, MugatuWarning)

    # now assign the apogee fibers
    fiberID_apogee = np.zeros(len(catalogid_apogee), dtype=int) - 9999
    n_apogee_assigned = 0
    for i in range(len(rg.robotDict)):
        if n_apogee_assigned == n_sky_apogee:
            break
        elif i in fiberID_boss:
            pass
        else:
            valid_targs = np.isin(rg.robotDict[i].validTargetIDs,
                                  catalogid_apogee[fiberID_apogee == -9999])
            valid_targs_arr = np.array(rg.robotDict[i].validTargetIDs, dtype=int)
            if np.any(valid_targs):
                fiberID_apogee[catalogid_apogee == valid_targs_arr[valid_targs][0]] = i
                n_apogee_assigned += 1

    if n_apogee_assigned < n_sky_apogee:
        flag = 'Not enough apogee skies could be assigned'
        warnings.warn(flag, MugatuWarning)

    # get indexes for assigments
    idx_apogee = np.where(fiberID_apogee != -9999)

    idx_boss = np.where(fiberID_boss != -9999)

    catalogids = np.append(catalogid_apogee[idx_apogee], catalogid_boss[idx_boss])
    ras = np.append(ra_apogee[idx_apogee], ra_boss[idx_boss])
    decs = np.append(dec_apogee[idx_apogee], dec_boss[idx_boss])
    fiberIDs = np.append(fiberID_apogee[idx_apogee], fiberID_boss[idx_boss])
    obsWaves = np.array(['APOGEE'] * n_apogee_assigned + ['BOSS'] * n_boss_assigned)

    if len(catalogids) < 500:
        size = len(catalogids)
        catalogids = np.append(catalogids, [-1] * (500 - size))
        ras = np.append(ras, [-9999.99] * (500 - size))
        decs = np.append(decs, [-9999.99] * (500 - size))
        fiberIDs = np.append(fiberIDs, [-1] * (500 - size))
        obsWaves = np.append(obsWaves, [''] * (500 - size))

    # create the initial design
    init_design = FPSDesign(-1, obsTime, racen=racen, deccen=deccen,
                            position_angle=position_angle,
                            observatory=observatory, mode_pk=None,
                            catalogids=catalogids,
                            ra=ras,
                            dec=decs,
                            fiberID=fiberIDs,
                            obsWavelength=obsWaves,
                            priority=np.array([1] * 500), design_file=None,
                            manual_design=True)

    # validate the design
    init_design.validate_design()

    return init_design
