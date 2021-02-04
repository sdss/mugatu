# @Author: Ilija Medan
# @Date: February 4, 2021
# @Filename: fpsdesign.py
# @License: BSD 3-Clause
# @Copyright: Ilija Medan

import numpy as np

from sdssdb.peewee.sdss5db.targetdb import Design, Carton, CartonToTarget, Magnitude, Mode, Assignment, Instrument, Target, Positioner

def mode_validation(design):
    """
    Validate that a design meets the criteria of
    an assigned observing mode

    Parameters
    ----------

    design: dict
        A design dictonary made from a FPSDesign object
    """

    design_targ_db = (
        Assignment.select(Carton.carton,
                          Carton.program,
                          Instrument.label,
                          Magnitude.gaia_g,
                          Target.ra,
                          Target.dec,
                          Positioner.id)
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
                  .switch(CartonToTarget)
                  .join(Magnitude,
                        on=(CartonToTarget.pk == Magnitiude.carton_to_target_pk))
                  .switch(CartonToTarget)
                  .join(Carton,
                        on=(CartonToTarget.carton_pk == Carton.pkk))
                  .where(Assignment.design_pk == design['design_pk']))
