from astropy_healpix import HEALPix
import astropy.units as u
from coordio.utils import offset_definition, Moffat2dInterp
import pickle
from multiprocessing import Pool
from functools import partial
from tqdm import tqdm
import numpy as np
import os.path

try:
    from sdssdb.peewee.sdss5db import database
    database.set_profile('operations')
    _database = True
except:
    _database = False

from sdssdb.peewee.sdss5db import targetdb
from mugatu.designmode import allDesignModes
from mugatu.designmode import build_brigh_neigh_query

fmagloss = Moffat2dInterp()


@database.connection_context()
def all_fields(plan):
    """
    return all of the fields
    """

    fields = targetdb.Field.select(targetdb.Field.racen,
                                   targetdb.Field.deccen,
                                   targetdb.Observatory.label)\
                           .join(targetdb.Version)\
                           .switch(targetdb.Field)\
                           .join(targetdb.Observatory)\
                           .where(targetdb.Version.plan == 'eta-6').dicts()
    return fields


@database.connection_context()
def return_desmodes():
    dms = allDesignModes()
    return dms


@database.connection_context()
def bn_field(field, mag_lim, instrument, mag_limits,
             lunation, hp):
    """
    get the bright neighbor exlusion regions
    """
    bright_idx = []
    db_query = build_brigh_neigh_query('designmode',
                                       instrument,
                                       mag_lim,
                                       field['racen'],
                                       field['deccen'],
                                       field['label'].strip().upper())
    if len(db_query) > 0:
        if isinstance(db_query, tuple):
            ras, decs, mags, catalogids, pmras, pmdecs = db_query
        else:
            ras, decs, mags, catalogids, pmras, pmdecs = map(list, zip(*list(db_query.tuples())))
        r_exclude, _ = offset_definition(mags,
                                         mag_limits,
                                         lunation=lunation,
                                         waveName=instrument.title(),
                                         fmagloss=fmagloss,
                                         obsSite=field['label'].strip().upper())
        for ra, dec, r in zip(ras, decs, r_exclude):
            bright_idx += list(hp.cone_search_lonlat(ra * u.deg, dec * u.deg, radius=r * u.arcsecond))
    return bright_idx


if __name__ == '__main__':
    # get all of the fields
    fields = all_fields('eta-6')

    # grab all designmode info
    desmodes = return_desmodes()

    ncores = 16

    desmode_check = []
    all_limits_boss = []
    all_limits_apogee = []

    # get unique designmodes
    for dm in desmodes.keys():
        if 'eng' not in dm:
            mag_limits_boss = desmodes[dm].bright_limit_targets['BOSS'][:, 0]
            mag_limits_apogee = desmodes[dm].bright_limit_targets['APOGEE'][:, 0]

            desmode_check.append(dm)
            all_limits_boss.append(mag_limits_boss)
            all_limits_apogee.append(mag_limits_apogee)

    finished_desmodes = []

    # do for all designmodes
    for dm in desmode_check:
        if dm not in finished_desmodes:
            hp = HEALPix(nside=2 ** 18, order='ring', frame='icrs')
            bn_maps_boss = []
            bn_maps_apogee = []

            mag_limits_boss = desmodes[dm].bright_limit_targets['BOSS'][:, 0]
            mag_limits_apogee = desmodes[dm].bright_limit_targets['APOGEE'][:, 0]

            if 'bright' in dm:
                mag_lim_boss = desmodes[dm].bright_limit_targets['BOSS'][5][0]
                lunation = 'bright'
            else:
                mag_lim_boss = desmodes[dm].bright_limit_targets['BOSS'][1][0]
                lunation = 'dark'

            mag_lim_apogee = desmodes[dm].bright_limit_targets['APOGEE'][8][0]

            dms_todo = [dm]
            for i in range(len(desmode_check)):
                if desmode_check[i] != dm:
                    if np.all(mag_limits_boss == all_limits_boss[i]):
                        dms_todo.append(desmode_check[i])

            # boss
            print('starting %s BOSS!' % dm)
            if not os.path.isfile('%s_boss_bn_healpix.pkl' % dm):
                with Pool(processes=ncores) as pool:
                    res = tqdm(pool.imap(partial(bn_field, mag_lim=mag_lim_boss,
                                                 instrument='BOSS', mag_limits=mag_limits_boss,
                                                 lunation=lunation, hp=hp), fields), total=len(fields))
                    res = [r for r in res]
                    for r in res:
                        bn_maps_boss += r
                # save the results
                bn_maps_boss = np.unique(bn_maps_boss)
                for dmt in dms_todo:
                    with open('%s_boss_bn_healpix.pkl' % dmt, 'wb') as f:
                        pickle.dump(bn_maps_boss, f)

            # apogee
            print('starting %s APOGEE!' % dm)
            if not os.path.isfile('%s_apogee_bn_healpix.pkl' % dm):
                with Pool(processes=ncores) as pool:
                    res = tqdm(pool.imap(partial(bn_field, mag_lim=mag_lim_apogee,
                                                 instrument='APOGEE', mag_limits=mag_limits_apogee,
                                                 lunation=lunation, hp=hp), fields), total=len(fields))
                    res = [r for r in res]
                    for r in res:
                        bn_maps_apogee += r
                # save the results
                bn_maps_apogee = np.unique(bn_maps_apogee)
                for dmt in dms_todo:
                    with open('%s_apogee_bn_healpix.pkl' % dmt, 'wb') as f:
                        pickle.dump(bn_maps_apogee, f)
            finished_desmodes += dms_todo
