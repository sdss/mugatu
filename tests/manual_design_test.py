# @Author: Ilija Medan
# @Date: February 10, 2021
# @Filename: manual_design_test.py
# @License: BSD 3-Clause
# @Copyright: Ilija Medan

import os
import sys
import argparse
from astropy.io import fits

import sdss_access.path

from mugatu.fpsdesign import FPSDesign
import robostrategy.obstime as obstime
import coordio.time


def test_warn_manual(file, exp, obsTime, observatory):
    """
    test a manual design with inputs that generate
    some warnings
    """
    des = FPSDesign(design_pk=-1,
                    obsTime=obsTime,
                    observatory=observatory.upper(),
                    design_file=file,
                    manual_design=True,
                    exp=exp)

    des.build_design_manual()
    des.validate_design()


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
    else:
        file = sdss_path.full('rsFieldAssignments',
                              plan=plan,
                              observatory=observatory,
                              fieldid=fieldid)
        if 'sas' in file:
            file = file[:32] + 'sdss50' + file[44:]

    head = fits.open(file)[0].header
    racen = head['RACEN']
    ot = obstime.ObsTime(observatory=observatory.lower())
    obsTime = coordio.time.Time(ot.nominal(lst=racen)).jd

    test_warn_manual(file=file, exp=exp, obsTime=obsTime,
                     observatory=observatory)
