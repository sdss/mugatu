# @Author: Ilija Medan
# @Date: March 3, 2021
# @Filename: comm_design_test.py
# @License: BSD 3-Clause
# @Copyright: Ilija Medan

from sdssdb.peewee.sdss5db import targetdb
import os
import roboscheduler.cadence as cadence
from mugatu.comm_designs import all_sky_design_RS


targetdb.database.connect_from_parameters(user='sdss',
                                          host='localhost',
                                          port=7500)

os.environ['KAIJU_DIR'] = '/Users/imedan/Desktop/Graduate_Coursework/SDSS_V/Design_Package/kaiju-int64'


def test_comm():
    """
    test a all sky commisioning design
    """
    # build cadence list from file
    des = all_sky_design_RS(racen=20, deccen=20,
                           position_angle=24.188576,
                           observatory='APO',
                           obsTime=2459145.5,
                           n_sky_apogee=10,
                           n_sky_boss=490)


if __name__ == "__main__":
    cadence.CadenceList().fromfits(filename='tests/rsCadences-test-newfield-apo.fits')
    # cadence.CadenceList().fromdb()
    test_comm()
