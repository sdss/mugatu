# @Author: Ilija Medan
# @Date: February 10, 2021
# @Filename: manual_design_test.py
# @License: BSD 3-Clause
# @Copyright: Ilija Medan

import os

from mugatu.fpsdesign import FPSDesign


def test_warn_manual():
    """
    test a manual design with inputs that generate
    some warnings
    """
    des = FPSDesign(3, 2459519.51, racen=297.33348791823, deccen=-11.2681218565069,
                    position_angle=24.188576, observatory='APO',
                    design_file='tests/manual_design_test.fits',
                    manual_design=True)

    des.build_design_manual()
    des.validate_design()


if __name__ == "__main__":
    test_warn_manual()
