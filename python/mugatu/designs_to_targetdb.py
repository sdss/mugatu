import warnings
import numpy as np

from sdssdb.peewee.sdss5db import targetdb
from mugatu.exceptions import MugatuError, MugatuWarning


def make_design_field_targetdb(cadence, fieldid, plan,
                               racen, deccen, position_angle,
                               observatory):
    """
    Create a new field in targetdb. Will return warning
    if the field already exists in targetdb

    Parameters
    ----------
    cadence: str or targetdb.Cadence instance
        Either label of the cadence for the field (str) or
        a targetdb.Cadence.get instance for the label that can
        be used to get the cadence pk

    fieldid: int
        The fieldid for the field

    plan: str or targetdb.Version instance
        Either robostratgegy plan as a str or a targetdb.Version.get
        instance for the plan that can be used to get the version pk

    racen: float
        Right Ascension center of the field (degrees)

    deccen: float
        Declination center of the field (degrees)

    position_angle: float
        Position angle of the field, East of North (degrees)

    observatory: str or targetdb.Observatory instance
        Either label of the observatory for the field (str; either
        'apo' or 'lco') or a targetdb.Observatory.get instance
        for the observatory label that can be used to get the
        observatory pk
    """

    # get the field cadence pk
    if isinstance(cadence, targetdb.Cadence):
        dbCadence = cadence.pk
    else:
        cadenceDB = targetdb.Cadence()
        dbCadence = cadenceDB.get(label=cadence).pk

    # get the observatory pk
    if isinstance(observatory, targetdb.Observatory):
        obspk = observatory.pk
    else:
        obsDB = targetdb.Observatory()
        obspk = obsDB.get(label=observatory.upper()).pk

    # get the version pk based on the plan
    if isinstance(plan, targetdb.Version):
        verpk = plan.pk
    else:
        versionDB = targetdb.Version()
        verpk = versionDB.get(plan=plan).pk

    # check if field exists
    field_test = (targetdb.Field
                  .select()
                  .where((targetdb.Field.field_id == fieldid) &
                         (targetdb.Field.version == verpk)))
    # creates new field in database if it doesnt exist
    if field_test.exists():
        flag = 'Field already exists in targetdb'
        warnings.warn(flag, MugatuWarning)
    else:
        fieldDB = targetdb.Field.create(
            field_id=fieldid,
            racen=racen,
            deccen=deccen,
            position_angle=position_angle,
            cadence=dbCadence,
            observatory=obspk,
            version=verpk)
        # save row in database
        fieldDB.save()


def make_design_assignments_targetdb(targetdb_ver, plan, fieldid, exposure,
                                     catalogID, fiberID, obsWavelength,
                                     carton, instr_pks=None, cart_pks=None,
                                     fiber_pks=None):
    """
    Add assignments for a design to targetdb.

    Parameters
    ----------
    targetdb_ver: dict
        dictonary of pks for the targetdb version of each carton
        used in this design

    plan: str or targetdb.Version instance
        Either robostratgegy plan as a str or a targetdb.Version.get
        instance for the plan that can be used to get the version pk

    fieldid: int or targetdb.Field instance
        The fieldid for the field (int) or a targetdb.Field
        instance for the field that can be used to get the field pk

    exposure: int
        The exposure of this set of designs. 0th indexed

    catalogID: np.array
        Array of catalogids for the design of length N

    fiberID: np.array
        Array of the fiberIDs (robotIDs in robostrategy)
        for the design of length N

    obsWavelength: np.array
        Array of obsWavelength for the design (choice of
        'BOSS' or 'APOGEE') for the design of legnth N

    carton: np.array
        Array of cartons for the design of length N

    instr_pks: dict
        Optional dictonary with the isntrument pks from
        targetdb

    cart_pks: dict
        Optional dictonary with the possible carton pks
        for the design

    fiber_pks: dict
        Optional dictonary with the fiber pks
    """

    # grab the targetdb tables
    carton_to_targetDB = targetdb.CartonToTarget()
    positionerDB = targetdb.Positioner()

   # get the version pk based on the plan
    if isinstance(plan, targetdb.Version):
        verpk = plan.pk
    else:
        versionDB = targetdb.Version()
        verpk = versionDB.get(plan=plan).pk

    # get the instrument pks
    if instr_pks is None:
        instr_pks = {}
        instr_pks['BOSS'] = targetdb.Instrument.get(label='BOSS').pk
        instr_pks['APOGEE'] = targetdb.Instrument.get(label='APOGEE').pk

    # grab all carton pks here
    if cart_pks is None:
        cart_pks = {}
        for cart in np.unique(carton):
            # skip calibration from now
            if cart != 'CALIBRATION':
                cart_pks[cart] = (targetdb.Carton.select(targetdb.Carton.pk)
                                                 .where((targetdb.Carton.carton == cart) &
                                                        (targetdb.Carton.version_pk == targetdb_ver[cart]))[0].pk)

    # get the fieldpk
    if isinstance(fieldid, targetdb.Field):
        fieldpk = fieldid[0].pk
    else:
        fieldDB = targetdb.Field()
        field = (targetdb.Field.select()
                               .join(targetdb.Version)
                               .where((targetdb.Field.field_id == fieldid) &
                                      (targetdb.Version.plan == plan)))
        fieldpk = field[0].pk

    designDB = targetdb.Design.create(field=fieldpk,
                                      exposure=exposure)
    # save row
    designDB.save()

    # add the assignments for the design to the assignment database
    rows = []
    for j in range(len(fiberID)):
        row_dict = {}

        # right now calibrations are fake, so need to skip
        if fiberID[j] != -1 and carton[j] != 'CALIBRATION':
            # get the pk for the positioner_info
            # (where I assume the ID is just the
            # row # in the fits file)

            if fiber_pks is None:
                this_pos_DB = (targetdb.Positioner.get(
                    id=fiberID[j]).pk)
            else:
                this_pos_DB = fiber_pks[fiberID[j]]

            # get the instrument for fiber
            inst_assign = obsWavelength[j]

            # add db row info to dic
            row_dict['design'] = designDB.pk
            row_dict['instrument'] = instr_pks[inst_assign]
            row_dict['positioner'] = this_pos_DB
            cart_pk = cart_pks[carton[j]]
            row_dict['carton_to_target'] = (targetdb.CartonToTarget.select(
                targetdb.CartonToTarget.pk)
                .join(targetdb.Target,
                      on=(targetdb.CartonToTarget.target_pk == targetdb.Target.pk))
                .where((targetdb.Target.catalogid == catalogID[j]) &
                       (targetdb.CartonToTarget.carton_pk == cart_pk))[0].pk)

            rows.append(row_dict)

    # write all exposures for field to targetdb
    targetdb.Assignment.insert_many(rows).execute()
