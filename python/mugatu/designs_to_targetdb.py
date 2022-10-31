import warnings
import numpy as np
import datetime
import json
import hashlib
from mugatu.exceptions import MugatuError, MugatuWarning

try:
    from sdssdb.peewee.sdss5db import database
    database.set_profile('operations')
    _database = True
except:
    _database = False

from sdssdb.peewee.sdss5db import targetdb
from peewee import fn
import mugatu


mugatu_version = mugatu.__version__


def make_design_field_targetdb(cadence, fieldid, plan,
                               racen, deccen, position_angle,
                               observatory, slots_exposures,
                               replacement_field=False):
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

    replacement_field: boolean
        If the field is a replacement field. If True then will
        ignore check for if field exists.
    """

    # get the field cadence pk
    if isinstance(cadence, str):
        cadenceDB = targetdb.Cadence()
        dbCadence = cadenceDB.get(label=cadence).pk
    else:
        dbCadence = cadence.pk

    # get the observatory pk
    if isinstance(observatory, str):
        obsDB = targetdb.Observatory()
        obspk = obsDB.get(label=observatory.upper()).pk
    else:
        obspk = observatory.pk

    # get the version pk based on the plan
    if isinstance(plan, str):
        versionDB = targetdb.Version()
        verpk = versionDB.get(plan=plan).pk
    else:
        verpk = plan.pk

    if replacement_field:
        fieldDB = targetdb.Field.create(
                field_id=fieldid,
                racen=racen,
                deccen=deccen,
                position_angle=position_angle,
                slots_exposures=slots_exposures,
                cadence=dbCadence,
                observatory=obspk,
                version=verpk)
        # save row in database
        fieldDB.save()
    else:
        # check if field exists
        field_test = (targetdb.Field
                      .select()
                      .where((targetdb.Field.field_id == fieldid) &
                             (targetdb.Field.version == verpk) &
                             (targetdb.Field.cadence == dbCadence)))
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
                slots_exposures=slots_exposures,
                cadence=dbCadence,
                observatory=obspk,
                version=verpk)
            # save row in database
            fieldDB.save()


def design_status_bitmask(revalidated_design=False,
                          replacement_design=False):
    """
    function to create bitmask for design_status in
    targetdb.DesignModeCheckResults

    Parameters
    ----------
    revalidated_design: boolean
        Has this design been revalidated

    replacement_design: boolean
        Is this design a replacement design
        for the robostrategy run

    Returns
    -------
    bitmask: int
        bitmask for the desin
    """
    bitmask = 0
    if revalidated_design:
        bitmask += int(2 ** 0)
    if replacement_design:
        bitmask += int(2 ** 1)
    return bitmask


def make_desigmmode_results_targetdb(design_id, design_pass,
                                     design_valid_file_row=None,
                                     design_status=None,
                                     boss_skies_min_pass=None,
                                     boss_skies_min_value=None,
                                     boss_skies_fov_pass=None,
                                     boss_skies_fov_value=None,
                                     apogee_skies_min_pass=None,
                                     apogee_skies_min_value=None,
                                     apogee_skies_fov_pass=None,
                                     apogee_skies_fov_value=None,
                                     boss_stds_min_pass=None,
                                     boss_stds_min_value=None,
                                     boss_stds_fov_pass=None,
                                     boss_stds_fov_value=None,
                                     apogee_stds_min_pass=None,
                                     apogee_stds_min_value=None,
                                     apogee_stds_fov_pass=None,
                                     apogee_stds_fov_value=None,
                                     boss_stds_mags_pass=None,
                                     apogee_stds_mags_pass=None,
                                     boss_bright_limit_targets_pass=None,
                                     apogee_bright_limit_targets_pass=None,
                                     boss_sky_neighbors_targets_pass=None,
                                     apogee_sky_neighbors_targets_pass=None,
                                     apogee_trace_diff_targets_pass=None):
    """
    add designmode check results to targetdb
    """
    if design_valid_file_row is not None:
        boss_skies_min_pass = design_valid_file_row['boss_n_skies_min']
        boss_skies_min_value = design_valid_file_row['boss_n_skies_min_value']
        boss_skies_fov_pass = design_valid_file_row['boss_min_skies_fovmetric']
        boss_skies_fov_value = design_valid_file_row['boss_min_skies_fovmetric_value']
        apogee_skies_min_pass = design_valid_file_row['apogee_n_skies_min']
        apogee_skies_min_value = design_valid_file_row['apogee_n_skies_min_value']
        apogee_skies_fov_pass = design_valid_file_row['apogee_min_skies_fovmetric']
        apogee_skies_fov_value = design_valid_file_row['apogee_min_skies_fovmetric_value']
        boss_stds_min_pass = design_valid_file_row['boss_n_stds_min']
        boss_stds_min_value = design_valid_file_row['boss_n_stds_min_value']
        boss_stds_fov_pass = design_valid_file_row['boss_min_stds_fovmetric']
        boss_stds_fov_value = design_valid_file_row['boss_min_stds_fovmetric_value']
        apogee_stds_min_pass = design_valid_file_row['apogee_n_stds_min']
        apogee_stds_min_value = design_valid_file_row['apogee_n_stds_min_value']
        apogee_stds_fov_pass = design_valid_file_row['apogee_min_stds_fovmetric']
        apogee_stds_fov_value = design_valid_file_row['apogee_min_stds_fovmetric_value']
        boss_stds_mags_pass = design_valid_file_row['boss_stds_mags']
        apogee_stds_mags_pass = design_valid_file_row['apogee_stds_mags']
        boss_bright_limit_targets_pass = design_valid_file_row['boss_bright_limit_targets']
        apogee_bright_limit_targets_pass = design_valid_file_row['apogee_bright_limit_targets']
        boss_sky_neighbors_targets_pass = design_valid_file_row['boss_sky_neighbors_targets']
        apogee_sky_neighbors_targets_pass = design_valid_file_row['apogee_sky_neighbors_targets']
        apogee_trace_diff_targets_pass = None  # not checking currently
    # write to database
    design_checkDB = targetdb.DesignModeCheckResults.create(
        design=design_id,
        design_pass=design_pass,
        design_status=design_status,
        boss_skies_min_pass=boss_skies_min_pass,
        boss_skies_min_value=boss_skies_min_value,
        boss_skies_fov_pass=boss_skies_fov_pass,
        boss_skies_fov_value=boss_skies_fov_value,
        apogee_skies_min_pass=apogee_skies_min_pass,
        apogee_skies_min_value=apogee_skies_min_value,
        apogee_skies_fov_pass=apogee_skies_fov_pass,
        apogee_skies_fov_value=apogee_skies_fov_value,
        boss_stds_min_pass=boss_stds_min_pass,
        boss_stds_min_value=boss_stds_min_value,
        boss_stds_fov_pass=boss_stds_fov_pass,
        boss_stds_fov_value=boss_stds_fov_value,
        apogee_stds_min_pass=apogee_stds_min_pass,
        apogee_stds_min_value=apogee_stds_min_value,
        apogee_stds_fov_pass=apogee_stds_fov_pass,
        apogee_stds_fov_value=apogee_stds_fov_value,
        boss_stds_mags_pass=boss_stds_mags_pass,
        apogee_stds_mags_pass=apogee_stds_mags_pass,
        boss_bright_limit_targets_pass=boss_bright_limit_targets_pass,
        apogee_bright_limit_targets_pass=apogee_bright_limit_targets_pass,
        boss_sky_neighbors_targets_pass=boss_sky_neighbors_targets_pass,
        apogee_sky_neighbors_targets_pass=apogee_sky_neighbors_targets_pass,
        apogee_trace_diff_targets_pass=apogee_trace_diff_targets_pass)
    design_checkDB.save()


def assignment_hash(ids, holeIDs):
    """
    Create a unique hash for the design

    Parameters
    ----------
    ids: np.array
        carton_to_target_pks for the assignments.

    holeIDs: np.array
        corresponding holeIDs for the assignments.

    Returns
    -------
    assign_hash: str
        The assignment_hash for the design
    """
    assign_dict = {}
    idx = np.argsort(holeIDs)
    assign_dict['carton_to_target'] = ids[idx].tolist()
    assign_dict['hole_id'] = holeIDs[idx].tolist()
    assign_hash = (hashlib.md5(json.dumps(assign_dict)
                               .encode('utf-8')).hexdigest())
    return assign_hash


def designToField_exists(design_id, field_id, plan):
    """
    Check if designToField entry already exists to
    to ensure unique entries (i.e. one unqiue
    combo of design_id, field_id and version)

    Parameters
    ----------
    design_id: int
        design_id for the design

    field_id: int
        field_id for the design

    plan: str or targetdb.Version instance
        Either robostratgegy plan as a str or a targetdb.Version.get
        instance for the plan that can be used to get the version pk

    Returns
    -------
    exists: boolean
        whether or not the combo of design_id, field_id and version
        already exists in targetdb
    """
    # get the version pk based on the plan
    if isinstance(plan, str):
        versionDB = targetdb.Version()
        verpk = versionDB.get(plan=plan).pk
    else:
        verpk = plan.pk

    # create the query
    check_designToField = (targetdb.DesignToField.select()
                                                 .join(targetdb.Field)
                                                 .where((targetdb.DesignToField.design == design_id) &
                                                        (targetdb.Field.field_id == field_id) &
                                                        (targetdb.Field.version == verpk)))
    exists = check_designToField.exists()
    return exists


def make_designToField(design, fieldid, exposure, field_exposure):
    """
    Create a new designToField entry

    Parameters
    ----------
    design: int or targetdb.Design instance
        The design_id for the design (can be int or peewee
        instance)

    fieldid: int or targetdb.Field instance
        The fieldid for the field (int) or a targetdb.Field
        instance for the field that can be used to get the field pk

    exposure: int
        The exposure of this set of designs. 0th indexed

    field_exposure: int
        The exposure of this set of designs as listed in the
        robostrategy design file. 0th indexed
    """
    # get the fieldpk
    if isinstance(design, int):
        design_id = design
    else:
        design_id = design.design_id

    if isinstance(fieldid, int):
        field = (targetdb.Field.select()
                               .where((targetdb.Field.field_id == fieldid) &
                                      (targetdb.Field.version == verpk)))
        fieldpk = field[0].pk
    else:
        fieldpk = fieldid[0].pk

    designToFieldDB = targetdb.DesignToField.create(design=design_id,
                                                    field=fieldpk,
                                                    exposure=exposure,
                                                    field_exposure=field_exposure)
    designToFieldDB.save()


def make_design_assignments_targetdb(plan, fieldid, exposure, field_exposure,
                                     desmode_label,
                                     design_ids, robotID, holeID,
                                     obsWavelength,
                                     carton, observatory, targetdb_ver=None,
                                     instr_pks=None,
                                     cart_pks=None, fiber_pks=None,
                                     idtype='carton_to_target',
                                     return_design_id=False):
    """
    Add assignments for a design to targetdb.

    Parameters
    ----------
    plan: str or targetdb.Version instance
        Either robostratgegy plan as a str or a targetdb.Version.get
        instance for the plan that can be used to get the version pk

    fieldid: int or targetdb.Field instance
        The fieldid for the field (int) or a targetdb.Field
        instance for the field that can be used to get the field pk

    exposure: int
        The exposure of this set of designs. 0th indexed

    field_exposure: int
        The exposure of this set of designs as listed in the
        robostrategy design file. 0th indexed

    desmode_label: str
        DesignMode labe for the design.

    design_ids: np.array
        Array of catalogids or carton_to_target_pks for the
        design of length N

    robotID: np.array
        Array of the robotIDs (robotIDs in robostrategy)
        for the design of length N

    holeID: np.array
        Array of holeIDs for the design of length N

    carton: np.array
        Array of cartons for the design of length N

    obsWavelength: np.array
        Array of obsWavelength for the design (choice of
        'BOSS' or 'APOGEE') for the design of legnth N

    targetdb_ver: dict
        Optional dictonary of pks for the targetdb version of each carton
        used in this design. Dict is indexed by carton names.
        Only needed if idtype='catalogID'.

    instr_pks: dict
        Optional dictonary with the isntrument pks from
        targetdb. Dict is indexed by instrument names.

    cart_pks: dict or array
        Optional dictonary with the possible carton pks
        for the design. If dict, then indexed by carton name.
        Optionally can be array of carton pks
        same length as design entries. Only needed if idtype='catalogID'.

    fiber_pks: dict
        Optional dictonary with the holeID pks. Dict is indexed by holeID.

    idtype: str
        Defines the id type used in defining the design_ids.
        Must be 'catalogID' or 'carton_to_target'.

    return_design_id: boolean
        Optionally return the design_id for the new entry in the database.
    """
    # make sure idtype is catalogID or carton_to_target
    if idtype != 'catalogID' and idtype != 'carton_to_target':
        raise MugatuError(message='idtype must be catalogID or carton_to_target')

    # get the version pk based on the plan
    if isinstance(plan, str):
        versionDB = targetdb.Version()
        verpk = versionDB.get(plan=plan).pk
    else:
        verpk = plan.pk

    # get the observatory pk
    if isinstance(observatory, str):
        obsDB = targetdb.Observatory()
        obspk = obsDB.get(label=observatory.upper()).pk
    else:
        obspk = observatory.pk

    # get the instrument pks
    if instr_pks is None:
        instr_pks = {}
        instr_pks['BOSS'] = targetdb.Instrument.get(label='BOSS').pk
        instr_pks['APOGEE'] = targetdb.Instrument.get(label='APOGEE').pk

    # grab all carton pks here
    if cart_pks is None and idtype == 'catalogID':
        cart_pks = {}
        for cart in np.unique(carton):
            # skip calibration from now
            if cart != 'CALIBRATION':
                cart_pks[cart] = (targetdb.Carton.select(targetdb.Carton.pk)
                                                 .where((targetdb.Carton.carton ==
                                                         cart) &
                                                        (targetdb.Carton.version_pk ==
                                                         targetdb_ver[cart]))[0].pk)

    # get the fieldpk
    if isinstance(fieldid, int):
        field = (targetdb.Field.select()
                               .where((targetdb.Field.field_id == fieldid) &
                                      (targetdb.Field.version == verpk)))
        fieldpk = field[0].pk
    else:
        fieldpk = fieldid[0].pk

    assign_hash = assignment_hash(design_ids[robotID != -1],
                                  holeID[robotID != -1])

    designDB = targetdb.Design.create(design_mode=desmode_label,
                                      mugatu_version=mugatu_version,
                                      run_on=datetime.datetime.now(),
                                      assignment_hash=assign_hash,
                                      version=verpk)
    # save row
    designDB.save()

    # add the designToField entry
    make_designToField(designDB.design_id, fieldid, exposure, field_exposure)

    # add the assignments for the design to the assignment database
    rows = []
    for j in range(len(robotID)):
        row_dict = {}

        # right now calibrations are fake, so need to skip
        if robotID[j] != -1:
            # get the pk for the positioner_info
            # (where I assume the ID is just the
            # row # in the fits file)

            if fiber_pks is None:
                this_pos_DB = (targetdb.Hole.get(
                    (targetdb.Hole.holeid == holeID[j]) &
                    (targetdb.Hole.observatory == obspk)).pk)
            else:
                this_pos_DB = fiber_pks[holeID[j]]

            # get the instrument for fiber
            inst_assign = obsWavelength[j]

            # add db row info to dic
            row_dict['design'] = designDB.design_id
            row_dict['instrument'] = instr_pks[inst_assign]
            row_dict['hole'] = this_pos_DB
            if idtype == 'catalogID':
                if isinstance(cart_pks, dict):
                    cart_pk = cart_pks[carton[j]]
                else:
                    cart_pk = cart_pks[j]
                row_dict['carton_to_target'] = (targetdb.CartonToTarget.select(
                    targetdb.CartonToTarget.pk)
                    .join(targetdb.Target,
                          on=(targetdb.CartonToTarget.target_pk == targetdb.Target.pk))
                    .where((targetdb.Target.catalogid == design_ids[j]) &
                           (targetdb.CartonToTarget.carton_pk == cart_pk))[0].pk)
            if idtype == 'carton_to_target':
                row_dict['carton_to_target'] = design_ids[j]

            rows.append(row_dict)

    # write all exposures for field to targetdb
    targetdb.Assignment.insert_many(rows).execute()
    # return design_id if requested
    if return_design_id:
        return designDB.design_id


class TargetdbFieldIDs(object):
    """
    Class to find available fieldids in targetdb

    Parameters
    ----------
    fieldid_type: str
        Which type of fields you are considering.
        Either 'manual' for manual commissioning
        fields or 'survey' for fields used in survey
        operations.

    version_plan: str
        The targetdb.Version.plan for the fields you
        want to consider.
    """
    def __init__(self, fieldid_type=None,
                 version_plan=None):
        self.fieldid_type = fieldid_type
        self.version_plan = version_plan
        # get the reserved field_ids
        field_reserve = targetdb.FieldReservation.select()
        self.field_reserve = [f.field_id for f in field_reserve]

    def find_next_available(self):
        """
        Find the next availble fieldid for the fieldid_type.
        Optionally can be some version plan.

        Returns
        -------
        fieldid: int
            Next available fieldid for the fieldid_type.
            This accounts for any current gaps that may be
            in range of fieldids.
        """
        if self.fieldid_type == 'manual':
            fieldid_bounds = [16000, 100000]
        else:
            fieldid_bounds = [100000, -999]
        # find the minimum field_id to start with
        if fieldid_bounds[1] == -999:
            arg = (targetdb.Field.field_id >= fieldid_bounds[0])
        else:
            arg = ((targetdb.Field.field_id >= fieldid_bounds[0]) &
                   (targetdb.Field.field_id < fieldid_bounds[1]))
        if self.version_plan is not None:
            arg = arg & (targetdb.Version.plan == self.version_plan)
        field_start = (targetdb.Field.select(
            targetdb.Field.field_id)
            .join(targetdb.Version)
            .where(arg))
        if len(field_start) > 0:
            fieldids = [f.field_id for f in field_start]
            # add reserved fieldids to list
            other_reserved = list(set(self.field_reserve) - set(fieldids))
            for oth in other_reserved:
                if fieldid_bounds[1] == -999:
                    if oth >= fieldid_bounds[0]:
                        fieldids.append(oth)
                else:
                    if oth >= fieldid_bounds[0] and oth < fieldid_bounds[1]:
                        fieldids.append(oth)
            fieldids = list(set(fieldids))
            # sort in order
            fieldids.sort()
            # look for breaks
            stop = False
            for i in range(1, len(fieldids)):
                if fieldids[i] != fieldids[i - 1] + 1:
                    fieldid = fieldids[i - 1] + 1
                    stop = True
                    break
            # if no stop, then max + 1
            if not stop:
                fieldid = max(fieldids) + 1
        else:
            fieldid = fieldid_bounds[0]
        return fieldid

    def check_availability(self, fieldid):
        """
        check if fieldid(s) are currently available in
        targetdb.

        Parameters
        ----------
        fieldid: int or list
            Single or multiple fieldids to check

        Returns
        -------
        field_avail: bool or np.array
            Single or array of booleans if the fieldid(s)
            are available (True) or not (False).
        """
        if isinstance(fieldid, int):
            arg = (targetdb.Field.field_id == fieldid)
        else:
            arg = (targetdb.Field.field_id.in_(fieldid))
        if self.version_plan is not None:
            arg = arg & (targetdb.Version.plan == self.version_plan)
        field = (targetdb.Field.select(
            targetdb.Field.field_id)
            .join(targetdb.Version)
            .where(arg))
        fieldids = [f[0] for f in field.tuples()]
        # add reserved fieldids to list
        other_reserved = list(set(self.field_reserve) - set(fieldids))
        fieldids += other_reserved
        if isinstance(fieldid, int):
            field_avail = fieldid not in fieldids
        else:
            field_avail = ~np.isin(fieldid, fieldids)
        return field_avail
