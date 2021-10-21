import os
import sys
import argparse
from astropy.io import fits

import sdss_access.path

from mugatu.fpsdesign import FPSDesign
from mugatu.designmode import DesignModeCheck
import robostrategy.obstime as obstime
import coordio.time

from sdssdb.peewee.sdss5db import targetdb
from sdssdb.peewee.sdss5db import catalogdb

targetdb.database.connect_from_parameters(user='sdss_user',
                                          host='localhost',
                                          port=7502)


def designmode_manual(file, exp, obsTime, observatory,
                      designmode_label):
    """
    test designmode check of manual design
    """
    des = FPSDesign(design_pk=-1,
                    obsTime=obsTime,
                    observatory=observatory.upper(),
                    design_file=file,
                    manual_design=True,
                    exp=exp,
                    idtype='carton_to_target')

    des.build_design_manual()
    des.validate_design(safety=False)

    mode = DesignModeCheck(FPSDesign=des,
                           desmode_label=designmode_label)

    mode.design_mode_check_all()


if __name__ == "__main__":
    sdss_path = sdss_access.path.Path(release='sdss5')

    parser = argparse.ArgumentParser(
        prog=os.path.basename(sys.argv[0]),
        description='Read Robostratgey design file to mugatu')

    parser.add_argument('-p', '--plan', dest='plan',
                        type=str, help='name of plan', required=True)
    parser.add_argument('-o', '--observatory', dest='observatory',
                        type=str, help='apo or lco',
                        choices=['apo', 'lco'], required=True)
    parser.add_argument('-f', '--fieldid', dest='fieldid',
                        type=int, help='fieldid',
                        required=True)
    parser.add_argument('-e', '--exp', dest='exp',
                        type=int, help='exposure of design',
                        required=True)
    parser.add_argument('-l', '--loc', dest='loc',
                        type=str, help='local or utah',
                        choices=['local', 'utah'], required=True)

    args = parser.parse_args()
    plan = args.plan
    observatory = args.observatory
    fieldid = args.fieldid
    exp = args.exp
    loc = args.loc

    if loc == 'local':
        file = 'rsFieldAssignments-%s-%s-%d.fits' % (plan, observatory, fieldid)
        targetdb.database.connect_from_parameters(user='sdss_user',
                                                  host='localhost',
                                                  port=7502)
        catalogdb.database.connect_from_parameters(user='sdss_user',
                                                  host='localhost',
                                                  port=7502)
    else:
        file = sdss_path.full('rsFieldAssignments',
                              plan=plan,
                              observatory=observatory,
                              fieldid=fieldid)
        if 'sas' in file:
            file = file[:32] + 'sdss50' + file[44:]
        targetdb.database.connect_from_parameters(user='sdss_user',
                                                  host='operations.sdss.utah.edu',
                                                  port=5432)
        catalogdb.database.connect_from_parameters(user='sdss_user',
                                                  host='operations.sdss.utah.edu',
                                                  port=5432)

    head = fits.open(file)[0].header
    racen = head['RACEN']
    designmode_label = head['DESMODE'].split(' ')[exp - 1]
    ot = obstime.ObsTime(observatory=observatory.lower())
    obsTime = coordio.time.Time(ot.nominal(lst=racen)).jd

    designmode_manual(file=file, exp=exp, obsTime=obsTime,
                     observatory=observatory,
                     designmode_label=designmode_label)
