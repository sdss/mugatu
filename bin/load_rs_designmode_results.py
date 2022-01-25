import sys
import argparse
import os
from mugatu.designs_to_targetdb import make_desigmmode_results_targetdb
from astropy.io import fits
from sdssdb.peewee.sdss5db import targetdb


if __name__ == '__main__':
    parser = argparse.ArgumentParser(
        prog=os.path.basename(sys.argv[0]),
        description='In a batch, validate a set of designs')
    parser.add_argument('-l', '--loc', dest='loc',
                        type=str, help='local or utah',
                        choices=['local', 'utah'], required=True)
    parser.add_argument('-p', '--plan', dest='plan',
                        type=str, help='name of plan', required=True)
    parser.add_argument('-o', '--observatory', dest='observatory',
                        type=str, help='apo or lco',
                        choices=['apo', 'lco'], required=True)

    args = parser.parse_args()
    loc = args.loc
    plan = args.plan
    observatory = args.observatory

    if loc == 'local':
        targetdb.database.connect_from_parameters(user='sdss',
                                                  host='localhost',
                                                  port=7502)
    else:
        targetdb.database.connect_from_parameters(user='sdss',
                                                  host='operations.sdss.utah.edu',
                                                  port=5432)

    valid_data = fits.open(('/uufs/chpc.utah.edu/common/home/sdss50/sdsswork/'
                            'sandbox/mugatu/rs_plan_validations/{plan}/'
                            'rs_{plan}_{obs}_design_validation_results.fits'.format(
                                plan=plan,
                                obs=observatory)))[1].data
    for valid in valid_data:
        fieldid = int(valid['file_name'].split('-')[-1])
        exp = int(valid['exp'])
        design = (targetdb.Design.select(targetdb.Design.design_id)
                                 .join(targetdb.Field)
                                 .where((targetdb.Field.field_id == fieldid) &
                                        (targetdb.Design.exposure == exp)))
        try:
            design_id = design[0].design_id
        except IndexError:
            print('No Design for field %d, exp %d' % (fieldid, exp))
            break
        make_desigmmode_results_targetdb(
            design_id=design_id,
            design_pass=True,
            design_valid_file_row=valid)
