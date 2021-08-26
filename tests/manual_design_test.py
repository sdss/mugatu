# @Author: Ilija Medan
# @Date: February 10, 2021
# @Filename: manual_design_test.py
# @License: BSD 3-Clause
# @Copyright: Ilija Medan

import os

import numpy as np

from mugatu.fpsdesign import FPSDesign
from coordio.utils import radec2wokxy, wokxy2radec
from astropy.io import fits
from astropy.time import Time
from coordio.exceptions import CoordinateError


def test_warn_manual():
    """
    test a manual design with inputs that generate
    some warnings
    """
    des = FPSDesign(design_pk=-1,
                    obsTime=2459410.7920949073,
                    design_file='rsFieldAssignments-test-designmode-2-apo-1.fits',
                    manual_design=True,
                    exp=1)

    des.build_design_manual()
    des.validate_design()
    print(len(des.targets_collided))


if __name__ == "__main__":
    test_warn_manual()
