import sys
import argparse
import os
import warnings

from sdssdb.peewee.sdss5db import targetdb
from robostrategy.params import RobostrategyParams
from mugatu.exceptions import MugatuError, MugatuWarning
from astropy.io import ascii
import pandas as pd
import numpy as np


if __name__ == '__main__':
    parser = argparse.ArgumentParser(
        prog=os.path.basename(sys.argv[0]),
        description='add entries to target generation table')
    parser.add_argument('-p', '--plan', dest='plan',
                        type=str, help='name of RS plan',
                        required=True)
    parser.add_argument('-f', '--first_release', dest='first_release',
                        help='planned first data release for this plan',
                        type=str, required=True)

    args = parser.parse_args()
    plan = args.plan
    first_release = args.first_release

    targetdb.database.connect_from_parameters(user='sdss',
                                              host='operations.sdss.utah.edu',
                                              port=5432)

    # grab the config file for this plan
    config = RobostrategyParams(plan=plan)
    label = config.cfg['Cartons']['version']
    label_str = 'v' + label

    # check if exists in db
    gen = targetdb.TargetingGeneration.select()\
                                      .where(targetdb.TargetingGeneration.label == label_str)
    if gen.exists():
        flag = 'Entry for targeting_generation already exists in targetdb'
        warnings.warn(flag, MugatuWarning)
        generation_pk = gen[0].pk
    else:
        targetgen = targetdb.TargetingGeneration.create(label=label_str,
                                                        first_release=first_release)
        targetgen.save()
        generation_pk = targetgen.pk

    # grab the cartons for this version if new targeting generation
    if ~gen.exists():
        cartons = ascii.read(os.getenv('RSCONFIG_DIR') + f'/etc/cartons-{label}.txt',
                             format='fixed_width',
                             delimiter='|')

        # add columns for targetdb.targeting_generation_to_carton
        cartons['rs_active'] = True
        cartons['rs_active'][cartons['active'] == 'n'] = False

        cartons['targeting_generation_pk'] = generation_pk

        cartons['stage'].name = 'rs_stage'

        cartons['carton_pk'] = 0
        for i in range(len(cartons)):
            cart = targetdb.Carton.select()\
                                  .join(targetdb.Version)\
                                  .where(targetdb.Carton.carton == cartons['carton'][i],
                                         targetdb.Version.plan == cartons['plan'][i])
            cartons['carton_pk'][i] = cart[0].pk

        rows = []
        cols = ['targeting_generation_pk', 'carton_pk', 'rs_stage', 'rs_active']
        for cart in cartons:
            row_dict = {}
            for c in cols:
                row_dict[c] =  cart[c]
            rows.append(row_dict)
        targetdb.TargetingGenerationToCarton.insert_many(rows).execute()

    # add generation_to_version entry
    version_pk = targetdb.Version.get(plan=plan).pk

    ver = targetdb.TargetingGenerationToVersion.select()\
                                               .where(targetdb.TargetingGenerationToVersion.targeting_generation_pk == generation_pk,
                                                      targetdb.TargetingGenerationToVersion.version_pk == version_pk)
    if ver.exists():
        flag = 'Entry for targeting_generation_to_version already exists in targetdb'
        warnings.warn(flag, MugatuWarning)
    else:
        targetgen_ver = targetdb.TargetingGenerationToVersion.create(targeting_generation_pk=generation_pk,
                                                                     version_pk=version_pk)
        targetgen_ver.save()
