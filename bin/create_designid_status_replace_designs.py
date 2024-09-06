import sys
import argparse
import os
import numpy as np
import glob
import datetime
from tqdm import trange, tqdm
from multiprocessing import Pool
from functools import partial

from astropy.io import fits
from astropy.table import Table

import sdss_access.path
import robostrategy.obstime as obstime
import coordio.time

from sdssdb.peewee.sdss5db import targetdb

from mugatu.fpsdesign import FPSDesign
from mugatu.designmode import find_designid_status
from mugatu.designs_to_targetdb import assignment_hash


def create_des_object(exp, design_file, obsTime):
    des = FPSDesign(design_pk=-1,
                    obsTime=obsTime,
                    design_file=design_file,
                    manual_design=True,
                    exp=exp)
    des.build_design_manual()
    return des.design


def get_designid_status(file, field_id, des_objs=None):
    """
    get the designid_status for a manual design
    """
    def designid_status(design_file, obsTime, exp, fexp, field_id, des_objs=None):
        if des_objs is None:
            des = FPSDesign(design_pk=-1,
                            obsTime=obsTime,
                            design_file=design_file,
                            manual_design=True,
                            exp=exp)
            des.build_design_manual()
            desob = des.design
        else:
            desob = des_objs[fexp]

        assign_hash = assignment_hash(desob['catalogID'][desob['robotID'] != -1],
                                      desob['holeID'][desob['robotID'] != -1])
        designid_status = find_designid_status(field_id, fexp, assign_hash=assign_hash)
        return designid_status

    with fits.open(file) as hdu:
        # get info on field
        n_exp = hdu[0].header['NEXP']
        racen = hdu[0].header['RACEN']
        ot = obstime.ObsTime(observatory=hdu[0].header['obs'].strip())
        obsTime = coordio.time.Time(ot.nominal(lst=racen)).jd

    # set up the opsdb connection
    from sdssdb.peewee.sdss5db import database
    database.set_profile('operations')
    from sdssdb.peewee.sdss5db import opsdb
    os.environ["OBSERVATORY"] = hdu[0].header['obs'].strip().upper()
    opsdb.database.connect()

    # get the designid_statuses
    designid = np.zeros(n_exp, dtype='>i4')
    status = np.zeros(n_exp, dtype='S20')
    if n_exp == 1:
        exp = 0
        designid[exp], status[exp] = designid_status(file, obsTime, exp, exp, field_id)
    else:
        for exp in trange(n_exp):
            designid[exp], status[exp] = designid_status(file, obsTime, exp + 1, exp,
                                                         field_id, des_objs=des_objs)
    return designid, status


if __name__ == '__main__':
    parser = argparse.ArgumentParser(
        prog=os.path.basename(sys.argv[0]),
        description='In a batch, validate a set of designs')
    parser.add_argument('-l', '--loc', dest='loc',
                        type=str, help='local or utah',
                        choices=['local', 'utah'], required=True)
    parser.add_argument('-p', '--plan', dest='plan',
                        type=str, help='name of RS plan',
                        required=True)
    parser.add_argument('-f', '--fieldids', dest='fieldids', nargs='+',
                        help='field_ids to replace)',
                        type=int, required=True)
    parser.add_argument('-n', '--Ncores', dest='Ncores',
                        type=int, help='number of cores to use. If Ncores=1, then not run in parallal.',
                        default=1, nargs='?')

    args = parser.parse_args()
    loc = args.loc
    plan = args.plan
    fieldids = args.fieldids
    Ncores = args.Ncores

    if loc == 'local':
        targetdb.database.connect_from_parameters(user='sdss',
                                                  host='localhost',
                                                  port=7502)
    else:
        targetdb.database.connect_from_parameters(user='sdss',
                                                  host='operations.sdss.utah.edu',
                                                  port=5432)
    # get the files with designs
    replace_path = ('/uufs/chpc.utah.edu/common/home/sdss50/sdsswork/'
                    'target/robostrategy_replacement/{plan}/'.format(
                        plan=plan))

    # grab the files
    files = []
    for fid in fieldids:
        files += [file for file in glob.glob(replace_path +
                                             '{plan}_{fid}*.fits'.format(
                                                 plan=plan,
                                                 fid=fid))]
    for f in files:
        if 'validation' in f or 'status' in f:
            files.remove(f)

    for file in files:
        with fits.open(file) as hdu:
            # get the metadata
            n_exp = hdu[0].header['NEXP']
            racen = hdu[0].header['RACEN']
            ot = obstime.ObsTime(observatory=hdu[0].header['obs'].strip())
            obsTime = coordio.time.Time(ot.nominal(lst=racen)).jd
            # make a STATUS HDU
            dtype = np.dtype([('fieldid', '>i4'),
                              ('designid', '>i4'),
                              ('status', 'S20')])
            status = np.zeros(n_exp, dtype=dtype)

            # find if the field already exists and get designids if it does
            same_field = targetdb.Field.select().where(hdu[0].header['RACEN'] == targetdb.Field.racen,
                                                       hdu[0].header['DECCEN'] == targetdb.Field.deccen,
                                                       hdu[0].header['PA'] == targetdb.Field.position_angle,
                                                       targetdb.Field.field_id >= 100000)
            if len(same_field) > 0:
                field_id = same_field[0].field_id
                # get the design objects if running in parallel
                if Ncores > 1 and n_exp > 1:
                    with Pool(processes=Ncores) as pool:
                        des_objs = tqdm(pool.imap(partial(create_des_object, design_file=file,
                                                          obsTime=obsTime),
                                                          range(1, n_exp + 1)),
                                                  total=n_exp)
                        des_objs = [r for r in des_objs]
                else:
                    des_objs = None
                designid, status = get_designid_status(file, field_id, des_objs=des_objs)
            else:
                field_id = -1
                designid = np.zeros(n_exp, dtype='>i4') - 1
                status = np.zeros(n_exp, dtype='S20')
                status[:] = 'not started'

            # add the new HDU
            status['fieldid'][:] = field_id
            status['designid'] = designid
            status['status'] = status

            hdu_status = fits.BinTableHDU(status, name='STATUS')
            hdu.append(hdu_status)
            hdu.write(file[:-5] + '_designid_status.fits')
