import numpy as np
import matplotlib.pylab as plt
from os import listdir, mkdir
from os.path import isfile, join, isdir
from matplotlib.projections.geo import GeoAxes
import healpy as hp
import healpy
import jinja2
import sys
import argparse
import os
import fitsio
from astropy.io import fits
from mugatu.designmode import allDesignModes
import mugatu
from sdssdb.peewee.sdss5db import targetdb
from mugatu.exceptions import MugatuWarning, MugatuError
import warnings
import datetime


plt.rcParams.update({'font.size': 18})
plt.rcParams['savefig.facecolor'] = 'white'


class ThetaFormatterShiftPi(GeoAxes.ThetaFormatter):
    """Shifts labelling by pi
    Shifts labelling from -180,180 to 0-360"""
    def __call__(self, x, pos=None):
        if x != 0:
            x *= 1
        if x < 0:
            x += 2*np.pi
        return GeoAxes.ThetaFormatter.__call__(self, x, pos)


def get_healpix_N_dmode(valid_file, designmode, Hpix_num_obs,
                  m, met, nside, i):
    """
    get N per pix in healpix map per dmode
    """
    N_obs = np.zeros(healpy.nside2npix(nside)) - 1

    for j in np.unique(Hpix_num_obs):
        ev_obs = eval("(Hpix_num_obs == j) & (valid_file.designmode == m) & (valid_file[valid_file.columns.names[i]] >= 0.)")
        vals = valid_file[valid_file.columns.names[i]][ev_obs]
        if len(vals) > 0:
            N_obs[j]  = np.mean(vals)
        else:
            N_obs[j] = np.nan
    return N_obs


def get_healpix_N(valid_file, nside):
    """
    get healpix N
    """
    Hpix_num_obs = np.zeros(len(valid_file), dtype=int)

    for file in np.unique(valid_file.file_name):
        ra = valid_file['racen'][valid_file.file_name == file][0]
        dec = valid_file['deccen'][valid_file.file_name == file][0]
        hpix = healpy.pixelfunc.ang2pix(nside, ra,
                                        dec, lonlat=True)
        Hpix_num_obs[valid_file.file_name == file] = hpix
    return Hpix_num_obs


def plots_healpix(valid_apo, valid_lco, designmode):
    """
    Make sky plots with healpix
    """
    if valid_apo is not None:
        Hpix_num_apo = get_healpix_N(valid_apo, 24)
        column_names = valid_apo.columns.names
        valid_file = valid_apo
    if valid_lco is not None:
        Hpix_num_lco = get_healpix_N(valid_lco, 32)
        column_names = valid_lco.columns.names
        valid_file = valid_lco

    for i in range(len(column_names)):
        if 'value' in column_names[i]:
            pf = valid_file.columns.names[i - 1]
            mets = []
            for dmode in np.unique(valid_file['designmode']):
                evd = eval("designmode['label'] == dmode")
                dmode_val = designmode[valid_file.columns.names[i-1]][evd][0]
                if isinstance(dmode_val, np.ndarray):
                    mets.append(dmode_val[-1])
                else:
                    mets.append(dmode_val)
            for m, met in zip(np.unique(valid_file['designmode']), mets):
                if valid_apo is not None:
                    N_apo = get_healpix_N_dmode(valid_apo, designmode, Hpix_num_apo,
                                                m, met, 24, i)
                    N_apo[N_apo < 0] = np.nan
                if valid_lco is not None:
                    N_lco = get_healpix_N_dmode(valid_lco, designmode, Hpix_num_lco,
                                                m, met, 32, i)
                    N_lco[N_lco < 0] = np.nan

                xsize = int(2000)
                ysize = int(xsize/2)
                    
                theta = np.linspace(np.pi, 0, ysize)
                phi   = np.linspace(np.pi, -np.pi, xsize)
                longitude = np.radians(np.linspace(-180, 180, xsize))
                latitude = np.radians(np.linspace(-90, 90, ysize))

                # project the map to a rectangular matrix xsize x ysize
                PHI, THETA = np.meshgrid(phi, theta)
                grid_pix_apo = hp.ang2pix(24, THETA, PHI)
                grid_pix_lco = hp.ang2pix(32, THETA, PHI)

                vmin = 0
                vmax = 2 * met
                if vmax < vmin:
                    vmax = 1

                fig = plt.figure(figsize=(16, 8))
                # matplotlib is doing the mollveide projection
                ax = fig.add_subplot(111, projection='mollweide')


                # rasterized makes the map bitmap while the labels remain vectorial
                # flip longitude to the astro convention
                if valid_apo is not None:
                    image = plt.pcolormesh(longitude[::-1], latitude, N_apo[grid_pix_apo],
                                           vmin=vmin, vmax=vmax, rasterized=True, cmap='coolwarm')
                if valid_lco is not None:
                    image = plt.pcolormesh(longitude[::-1], latitude, N_lco[grid_pix_lco],
                                           vmin=vmin, vmax=vmax, rasterized=True, cmap='coolwarm')

                # graticule
                #ax.set_longitude_grid(30)
                ax.xaxis.set_major_formatter(ThetaFormatterShiftPi(30))
                ax.set_title('%s: %s' % (pf, m))

                plt.colorbar(image,
                             orientation='vertical', label=r'Mean %s for Field' % pf,
                             ax=ax)

                ax.set_xlabel('R.A.')
                ax.set_ylabel('Decl.')
                ax.grid()
                plt.savefig(path + '/sky_plots/%s_%s.png' % (pf, m), bbox_inches='tight',dpi=150)
                plt.figure().clear()
                plt.close()
                plt.cla()
                plt.clf()
                plt.close(fig)


def plot_hist(ax, values, cumulative, density, bins, title,
              label, xlabel, ylabel, dmode_val):
    ax.hist(values, cumulative=cumulative, density=density,
            bins=bins, alpha=0.3, label=label)
    ax.set_title(title)
    ax.set_xlabel(xlabel)
    ax.set_ylabel(ylabel)
    ax.grid()
    ax.legend()
    if isinstance(dmode_val, np.ndarray):
        ax.axvline(dmode_val[-1], linestyle='--', c='r',
                   label='DesignMode Value')
    else:
        ax.axvline(dmode_val, linestyle='--', c='r',
                   label='DesignMode Value')


def get_hist_values(valid_file, i, dmode, designmode):
    """
    get the histogram values
    """
    ev = eval("(valid_file['designmode'] == dmode) & (valid_file[valid_file.columns.names[i]] >= 0)")
    evd = eval("designmode['label'] == dmode")
    dmode_val = designmode[valid_file.columns.names[i-1]][evd][0]
    xval = valid_file[valid_file.columns.names[i]][ev]
    try:
        maxx = np.max(xval) * 1.1
        minn = np.min(xval) * 0.9
    except ValueError:
        maxx = 1
        minn = 0
    if minn < 0:
        minn = 0
    if maxx < 0:
        maxx = 1
    return xval, minn, maxx, dmode_val


def create_summary_dist_plots(valid_apo, valid_lco, designmode):
    """
    Create the summary distribution plots
    """
    if valid_apo is not None:
        column_names = valid_apo.columns.names
        valid_file = valid_apo
    if valid_lco is not None:
        column_names = valid_lco.columns.names
        valid_file = valid_lco
    for i in range(len(column_names)):
        if 'value' in column_names[i]:
            for dmode in np.unique(valid_file['designmode']):
                f, (ax1, ax2, ax3, ax4) = plt.subplots(1, 4, figsize=(40,10))

                if valid_apo is not None:
                    xval_apo, minn_apo, maxx_apo, dmode_val = get_hist_values(valid_apo, i, dmode, designmode)

                if valid_lco is not None:
                    xval_lco, minn_lco, maxx_lco, dmode_val = get_hist_values(valid_lco, i, dmode, designmode)

                if valid_apo is not None and valid_lco is not None:
                    if minn_lco < minn_apo:
                        minn = minn_lco
                    else:
                        minn = minn_apo
                    if maxx_lco > maxx_apo:
                        maxx = maxx_lco
                    else:
                        maxx = maxx_apo
                elif valid_apo is not None:
                    minn = minn_apo
                    maxx = maxx_apo
                else:
                    minn = minn_lco
                    maxx = maxx_lco

                if valid_apo is not None:
                    plot_hist(ax1, xval_apo,
                              True, True,
                              np.linspace(minn,
                                          maxx,
                                          50),
                              '%s: %s' % (dmode, 'APO'),
                              'APO', valid_apo.columns.names[i - 1],
                              'Cumulative Fraction', dmode_val)
                    plot_hist(ax2, xval_apo,
                              False, False,
                              np.linspace(minn,
                                          maxx,
                                          50),
                              '%s: %s' % (dmode, 'APO'),
                              'APO', valid_apo.columns.names[i - 1],
                              'N', dmode_val)
                if valid_lco is not None:
                    plot_hist(ax3, xval_lco,
                              True, True,
                              np.linspace(minn,
                                          maxx,
                                          50),
                              '%s: %s' % (dmode, 'LCO'),
                              'LCO', valid_lco.columns.names[i - 1],
                              'Cumulative Fraction', dmode_val)
                    plot_hist(ax4, xval_lco,
                              False, False,
                              np.linspace(minn,
                                          maxx,
                                          50),
                              '%s: %s' % (dmode, 'LCO'),
                              'LCO', valid_lco.columns.names[i - 1],
                              'N', dmode_val)

                plt.savefig(path + '/dist_plots/%s_%s.png' % (valid_file.columns.names[i-1],
                                                      dmode),
                            bbox_inches='tight', dpi=150)
                plt.figure().clear()
                plt.close()
                plt.cla()
                plt.clf()
                plt.close(f)
        if 'pass' in column_names[i]:
            for dmode in np.unique(valid_file['designmode']):
                f, (ax1, ax2, ax3, ax4) = plt.subplots(1, 4, figsize=(40,10))

                if valid_apo is not None:
                    ev = eval("valid_apo['designmode'] == dmode")
                    xval_apo = valid_apo[valid_apo.columns.names[i+1]][ev] - valid_apo[valid_apo.columns.names[i]][ev]
                    try:
                        maxx_apo = np.max(xval_apo)
                    except ValueError:
                        maxx_apo = 0
                if valid_lco is not None:
                    ev = eval("valid_lco['designmode'] == dmode")
                    xval_lco = valid_lco[valid_lco.columns.names[i+1]][ev] - valid_lco[valid_lco.columns.names[i]][ev]
                    try:
                        maxx_lco = np.max(xval_lco)
                    except ValueError:
                        maxx_lco = 0
                if valid_apo is not None and valid_lco is not None:
                    if maxx_lco > maxx_apo:
                        maxx = maxx_lco
                    else:
                        maxx = maxx_apo
                elif valid_apo is not None:
                    maxx = maxx_apo
                else:
                    maxx = maxx_lco
                if valid_apo is not None:
                    plot_hist(ax1, xval_apo,
                              True, True,
                              np.arange(0, maxx + 2, 1),
                              '%s: %s' % (dmode, 'APO'),
                              'APO', valid_apo.columns.names[i - 1] + ' (N Fibers in Design Failed)',
                              'Cumulative Fraction', np.nan)
                    plot_hist(ax2, xval_apo,
                              False, False,
                              np.arange(0, maxx + 2, 1),
                              '%s: %s' % (dmode, 'APO'),
                              'APO', valid_apo.columns.names[i - 1] + ' (N Fibers in Design Failed)',
                              'N', np.nan)
                if valid_lco is not None:
                    plot_hist(ax3, xval_lco,
                              True, True,
                              np.arange(0, maxx + 2, 1),
                              '%s: %s' % (dmode, 'LCO'),
                              'LCO', valid_lco.columns.names[i - 1] + ' (N Fibers in Design Failed)',
                              'Cumulative Fraction', np.nan)
                    plot_hist(ax4, xval_lco,
                              False, False,
                              np.arange(0, maxx + 2, 1),
                              '%s: %s' % (dmode, 'LCO'),
                              'LCO', valid_lco.columns.names[i - 1] + ' (N Fibers in Design Failed)',
                              'N', np.nan)
                plt.savefig(path + '/dist_plots/%s_%s.png' % (valid_file.columns.names[i-1],
                                                      dmode),
                            bbox_inches='tight', dpi=150)
                plt.figure().clear()
                plt.close()
                plt.cla()
                plt.clf()
                plt.close(f)


def write_html_jinja(valid_apo, valid_lco, designmode,
                     rs_run, mugatu_v, kaiju_v, coordio_v,
                     fps_calib_v, path):
    """
    Write HTML file using jinja
    """
    env = jinja2.Environment(
        loader=jinja2.PackageLoader("mugatu"),
        autoescape=jinja2.select_autoescape()
    )
    if not isdir('%s/indv_valid' % path):
        mkdir('%s/indv_valid' % path)
    if not isdir('%s/dist_plots' % path):
        mkdir('%s/dist_plots' % path)
    if not isdir('%s/sky_plots' % path):
        mkdir('%s/sky_plots' % path)
    plots_healpix(valid_apo, valid_lco, designmode)
    create_summary_dist_plots(valid_apo, valid_lco, designmode)
    mypath = path + '/dist_plots/'
    hist_files = [f for f in listdir(mypath) if isfile(join(mypath, f))]
    hist_files.sort()
    hist_files = np.array(hist_files)

    mypath_sky = path + '/sky_plots/'
    sky_files = [f for f in listdir(mypath_sky) if isfile(join(mypath_sky, f))]
    sky_files.sort()
    sky_files = np.array(sky_files)
    if valid_apo is not None:
        valid_file = valid_apo
    else:
        valid_file = valid_lco
    for col, ty in zip(valid_file.columns.names, valid_file.columns.formats):
        if ty == 'L' and 'apogee' not in col:
            html_dict = {}
            if 'boss' in col:
                valid_check = col[5:]

                boss_col = 'boss_' + valid_check
                ap_col = 'apogee_' + valid_check
                cols = [boss_col, boss_col, ap_col, ap_col]
                obs = ['APO', 'LCO', 'APO', 'LCO']

                for h in designmode.columns.names:
                    if valid_check in h:
                        for m in np.unique(valid_file['designmode']):
                            val = designmode[h][designmode.label == m][0]
                            if 'boss' in h:
                                html_dict[m + '_boss'] = '%s' % val
                            else:
                                html_dict[m + '_apogee'] = '%s' % val
            else:
                valid_check = col

                cols = [valid_check, valid_check]
                obs = ['APO', 'LCO']
            html_dict['path'] = path
            html_dict['valid_check'] = valid_check

            for o, col in zip(obs, cols):
                for m in np.unique(valid_file['designmode']):
                    val = 'NA'
                    if o == 'APO' and valid_apo is not None:
                        try:
                            value = (100 *
                                     len(valid_apo[col][(valid_apo[col]) & (valid_apo['designmode'] == m)]) /
                                     len(valid_apo[col][(valid_apo['designmode'] == m)]))
                            val = '%.2f' % value + '%'
                        except ZeroDivisionError:
                            val = 'NA'
                    elif o == 'LCO' and valid_lco is not None:
                        try:
                            value = (100 *
                                     len(valid_lco[col][(valid_lco[col]) & (valid_lco['designmode'] == m)]) /
                                     len(valid_lco[col][(valid_lco['designmode'] == m)]))
                            val = '%.2f' % value + '%'
                        except ZeroDivisionError:
                            val = 'NA'
                    else:
                        pass

                    if 'boss' in col:
                        html_dict[o + '_boss_' + m] = val
                    elif 'apogee' in col:
                        html_dict[o + '_apogee_' + m] = val
                    else:
                        html_dict[o + '_' + m] = val
            if len(cols) == 2:
                if any(valid_check in s for s in hist_files):
                    template = env.get_template('valid_check_w_dists.html')
                else:
                    template = env.get_template('valid_check.html')
            else:
                if any(valid_check in s for s in hist_files) and any(valid_check in s for s in sky_files):
                    template = env.get_template('dmode_check_w_dists_sky.html')
                else:
                    template = env.get_template('dmode_check_w_dists.html')
            page = template.render(html_dict)

            fp = open('%s/indv_valid/%s.html' % (path, valid_check), 'w')
            fp.write(page)
            fp.close()
    html_dict = {}
    html_dict['rs_run'] = rs_run
    html_dict['mugatu_v'] = mugatu_v
    html_dict['kaiju_v'] = kaiju_v
    html_dict['coordio_v'] = coordio_v
    html_dict['fps_calib_v'] = fps_calib_v
    template = env.get_template('main_validation_page.html')
    page = template.render(html_dict)
    fp = open('%s/%s_validation.html' % (path, rs_run), 'w')
    fp.write(page)
    fp.close()


if __name__ == '__main__':

    # grabbed the parser from robostratgey code to keep consistent
    # to initiate code use

    parser = argparse.ArgumentParser(
        prog=os.path.basename(sys.argv[0]),
        description='Write validaiton summary HTML file')

    parser.add_argument('-p', '--plan', dest='plan',
                        type=str, help='name of plan (for type=rs or type=rs_replace)',
                        required=False)
    parser.add_argument('-d', '--dir', dest='dir',
                        type=str, help='directory with design files (for type=dir)',
                        required=False)
    parser.add_argument('-da', '--date', dest='date',
                        type=str, help='date of validation in YYYY_MM_DD format (for type=dir)',
                        required=False)
    parser.add_argument('-t', '--type', dest='type',
                        type=str, help='Validating files in directory (dir) or robostrategy (rs)',
                        choices=['dir', 'rs', 'rs_catchup'], required=False,
                        default='rs')
    parser.add_argument('-v', '--ver_catch', dest='ver_catch',
                        type=str, help='version of catchup (for type=rs_catchup)', required=False)
    parser.add_argument('-i','--ignore_prev', help='True if want to ignore previously observed designs in validation',
                        type=bool, required=False, default=True)
    args = parser.parse_args()
    plan = args.plan
    directory = args.dir
    date = args.date
    vtype = args.type
    ver_catch = args.ver_catch
    ignore_prev = args.ignore_prev

    targetdb.database.connect_from_parameters(user='sdss_user',
                                              host='operations.sdss.utah.edu',
                                              port=5432)

    MUGATU_DATA = os.getenv('MUGATU_DATA')
    if vtype == 'dir':
        path = (directory + '/design_validation_{year}_{month}_{day}'
                  .format(year=date[:4],
                          month=date[5:7],
                          day=date[8:]))
    elif vtype == 'rs':
        path = MUGATU_DATA + '/rs_plan_validations/%s' % plan
    elif vtype == 'rs_catchup':
        path = MUGATU_DATA + '/rs_plan_validations/%s' % plan
    else:
        raise MugatuError(message='Improper Validation Type')

    # create designmode file for run
    designModeDict = allDesignModes()
    dmarr = None
    for i, d in enumerate(designModeDict):
        arr = designModeDict[d].toarray()
        if(dmarr is None):
            dmarr = np.zeros(len(designModeDict), dtype=arr.dtype)
        dmarr[i] = arr
    fitsio.write(path + '/designmodes_validation.fits', dmarr)
    # get validaiton results
    try:
        if vtype == 'dir':
            file = 'design_validation_results.fits'
        elif vtype == 'rs':
            file = 'rs_{plan}_apo_design_validation_results.fits'.format(plan=plan)
        elif vtype == 'rs_catchup':
            file = 'rs_Catchup{ver_catch}_{plan}_apo_design_validation_results.fits'.format(ver_catch=ver_catch,
                                                                                            plan=plan)
        else:
            raise MugatuError(message='Improper Validation Type')
        valid_apo = fits.open(path + '/' + file)[1].data
        header = fits.open(path + '/' + file)[0].header
        kaiju_v = header['kaiju_version']
        coordio_v = header['coordio_version']
        fps_calib_v = header['fps_calibrations_version']
        mugatu_version = header['mugatu_version']
        # ignore previously observed designs with designid_status != -1
        if ignore_prev:
            valid_apo = valid_apo[valid_apo['designid_status'] == -1]
    except FileNotFoundError:
        valid_apo = None
        flag = 'No validation for APO!'
        warnings.warn(flag, MugatuWarning)
    try:
        if vtype == 'dir':
            raise FileNotFoundError
        elif vtype == 'rs':
            file = 'rs_{plan}_lco_design_validation_results.fits'.format(plan=plan)
        elif vtype == 'rs_catchup':
            file = 'rs_Catchup{ver_catch}_{plan}_lco_design_validation_results.fits'.format(ver_catch=ver_catch,
                                                                                            plan=plan)
        else:
            raise MugatuError(message='Improper Validation Type')
        valid_lco = fits.open(path + '/' + file)[1].data
        header = fits.open(path + '/' + file)[0].header
        kaiju_v = header['kaiju_version']
        coordio_v = header['coordio_version']
        fps_calib_v = header['fps_calibrations_version']
        mugatu_version = header['mugatu_version']
        if ignore_prev:
            valid_lco = valid_lco[valid_lco['designid_status'] == -1]
    except FileNotFoundError:
        valid_lco = None
        flag = 'No validation for LCO!'
        warnings.warn(flag, MugatuWarning)
    if valid_apo is None and valid_lco is None:
        message = 'No validation files for this run'
        raise MugatuError(message=message)
    designmode = fits.open(path + '/designmodes_validation.fits')[1].data

    write_html_jinja(valid_apo, valid_lco, designmode,
                     plan, mugatu_version, kaiju_v, coordio_v, fps_calib_v,
                     path)
