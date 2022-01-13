import sys
import argparse
import os
from mugatu.designs_to_targetdb import make_desigmmode_results_targetdb
from astropy.io import fits
import numpy as np
from sdssdb.peewee.sdss5db import targetdb


if __name__ == '__main__':
    parser = argparse.ArgumentParser(
        prog=os.path.basename(sys.argv[0]),
        description='In a batch, validate a set of designs')
    parser.add_argument('-l', '--loc', dest='loc',
                        type=str, help='local or utah',
                        choices=['local', 'utah'], required=True)
    parser.add_argument('-d', '--dir', dest='dir',
                        type=str, help='directory with design validation_files',
                        required=True)

    args = parser.parse_args()
    loc = args.loc
    directory = args.dir

    if loc == 'local':
        targetdb.database.connect_from_parameters(user='sdss',
                                                  host='localhost',
                                                  port=7502)
    else:
        targetdb.database.connect_from_parameters(user='sdss',
                                                  host='operations.sdss.utah.edu',
                                                  port=5432)

    valid_data = fits.open(directory +
                           'design_validation_results.fits')[1].data
    valid_design_ids = fits.open(directory +
                                 'design_ids_for_design_files.fits')[1].data
    for valid in valid_data:
        ind = np.where((valid_design_ids['file_name'] == valid['file_name']) &
                       (valid_design_ids['exp'] == valid['exp']))[0][0]
        make_desigmmode_results_targetdb(
            design_id=valid_design_ids[ind]['design_id'],
            design_pass=True,
            design_valid_file_row=valid)
