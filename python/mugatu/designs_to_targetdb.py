import warnings

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
    cadence: str
        The label for the field cadence

    fieldid: int
        The fieldid for the field

    plan: str
        the robostrategy plan that this field corresponds to

    racen: float
        Right Ascension center of the field (degrees)

    deccen: float
        Declination center of the field (degrees)

    position_angle: float
        Position angle of the field, East of North (degrees)

    observatory: str
        Observatory for the field, 'apo' or 'lco'
    """

    # grab the targetdb tables
    obsDB = targetdb.Observatory()
    cadenceDB = targetdb.Cadence()
    versionDB = targetdb.Version()

    # get the field cadence pk
    dbCadence = cadenceDB.get(label=cadence).pk

    # get the observatory pk
    obspk = obsDB.get(label=observatory.upper()).pk

    # get the version pk based on the plan
    verpk = versionDB.get(plan=plan).pk

    # check if field exists
    field_test = (targetdb.Field
                  .select()
                  .where((targetdb.Field.field_id == fieldid) &
                         (targetdb.Field.version == verpk)))
    # creates new field in database if it doesnt exist
    if len(field_test) == 0:
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
    else:
        flag = 'Field already exists in targetdb'
        warnings.warn(flag, MugatuWarning)


def make_design_assignments_targetdb(targetdb_ver, plan, fieldid, exposure,
                                     catalogID, fiberID, obsWavelength,
                                     carton):
    """
    targetdb_ver: int
        pk for the targetdb version of this design

    plan: str
        the robostrategy plan that this field corresponds to

    fieldid: int
        The fieldid for the field

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
    """

    # grab the targetdb tables
    versionDB = targetdb.Version()
    fieldDB = targetdb.Field()
    cartonDB = targetdb.Carton()
    carton_to_targetDB = targetdb.CartonToTarget()
    positionerDB = targetdb.Positioner()

    # get the version pk based on the plan
    verpk = versionDB.get(plan=plan).pk

    # get the instrument pks
    instr_pks = {}
    instr_pks['BOSS'] = targetdb.Instrument.get(label='BOSS').pk
    instr_pks['APOGEE'] = targetdb.Instrument.get(label='APOGEE').pk

    # grab all carton pks here
    cart_pks = {}
    for cart in np.unique(carton):
        # skip calibration from now
        if cart != 'CALIBRATION':
            cart_pks[cart] = (cartonDB.select(cartonDB.pk)
                                      .where((cartonDB.carton == cart) &
                                             (cartonDB.version_pk == targetdb_ver))[0].pk)

    # get the fieldpk
    field = (fieldDB.select()
                    .join(versionDB)
                    .where((fieldDB.field_id=fieldid) &
                           (versionDB.plan=plan)))
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

            this_pos_DB = (positionerDB.get(
                id=fiberID[j]).pk)

            # get the instrument for fiber
            inst_assign = obsWavelength[j]

            # add db row info to dic
            row_dict['design'] = designDB.pk
            row_dict['instrument'] = instr_pks[inst_assign]
            row_dict['positioner'] = this_pos_DB
            cart_pk = cart_pks[carton[j]]
            row_dict['carton_to_target'] = (carton_to_targetDB.select(
                carton_to_targetDB.pk)
                .join(targetDB,
                      on=(carton_to_targetDB.target_pk == targetDB.pk))
                .where((targetDB.catalogid == catalogID[j]) &
                       (carton_to_targetDB.carton_pk == cart_pk))[0].pk)

            rows.append(row_dict)

    # write all exposures for field to targetdb
    targetdb.Assignment.insert_many(rows).execute()
