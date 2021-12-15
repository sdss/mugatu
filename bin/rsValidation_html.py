import numpy as np
import matplotlib.pylab as plt
from os import listdir, mkdir
from os.path import isfile, join, isdir
from matplotlib.projections.geo import GeoAxes
import healpy as hp
import healpy
import jinja2

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


def plots_healpix(valid_apo, valid_lco, designmode):
    """
    Make sky plots with healpix
    """
    Hpix_num_apo = np.zeros(len(valid_apo), dtype=int)
    Hpix_num_lco = np.zeros(len(valid_lco), dtype=int)

    for file in np.unique(valid_apo.file_name):
        ra = valid_apo['racen'][valid_apo.file_name == file][0]
        dec = valid_apo['deccen'][valid_apo.file_name == file][0]
        hpix = healpy.pixelfunc.ang2pix(24, ra,
                                        dec, lonlat=True)
        Hpix_num_apo[valid_apo.file_name == file] = hpix
        
    for file in np.unique(valid_lco.file_name):
        ra = valid_lco['racen'][valid_lco.file_name == file][0]
        dec = valid_lco['deccen'][valid_lco.file_name == file][0]
        hpix = healpy.pixelfunc.ang2pix(32, ra,
                                        dec, lonlat=True)
        Hpix_num_lco[valid_lco.file_name == file] = hpix

    for i in range(len(valid_apo.columns.names)):
        if 'value' in valid_apo.columns.names[i]:
            pf = valid_apo.columns.names[i - 1]
            mets = []
            for dmode in np.unique(valid_apo['designmode']):
                evd = eval("designmode['label'] == dmode")
                dmode_val = designmode[valid_apo.columns.names[i-1]][evd][0]
                if isinstance(dmode_val, np.ndarray):
                    mets.append(dmode_val[-1])
                else:
                    mets.append(dmode_val)
            for m, met in zip(np.unique(valid_apo['designmode']), mets):
                N_apo = np.zeros(healpy.nside2npix(24)) - 1
                N_lco = np.zeros(healpy.nside2npix(32)) - 1

                for j in np.unique(Hpix_num_apo):
                    ev_apo =eval("(Hpix_num_apo == j) & (valid_apo.designmode == m) & (valid_apo[valid_apo.columns.names[i]] >= 0.)")
                    vals = valid_apo[valid_apo.columns.names[i]][ev_apo]
                    if len(vals) > 0:
                        N_apo[j]  = np.mean(vals)
                    else:
                        N_apo[j] = np.nan

                for j in np.unique(Hpix_num_lco):
                    ev_lco =eval("(Hpix_num_lco == j) & (valid_lco.designmode == m) & (valid_lco[valid_apo.columns.names[i]] >= 0.)")
                    vals = valid_lco[valid_apo.columns.names[i]][ev_lco]
                    if len(vals) > 0:
                        N_lco[j]  = np.mean(vals)
                    else:
                        N_lco[j] = np.nan

                N_apo[N_apo < 0] = np.nan
                N_lco[N_lco < 0] = np.nan

                xsize = int(2000)
                ysize = int(xsize/2)
                    
                theta = np.linspace(np.pi, 0, ysize)
                phi   = np.linspace(np.pi, -np.pi, xsize)
                longitude = np.radians(np.linspace(-180, 180, xsize))
                latitude = np.radians(np.linspace(-90, 90, ysize))

                # project the map to a rectangular matrix xsize x ysize
                PHI, THETA = np.meshgrid(phi, theta)
                grid_pix = hp.ang2pix(24, THETA, PHI)

                grid_map = N_apo[grid_pix]

                grid_pix2 = hp.ang2pix(32, THETA, PHI)

                grid_map2 = N_lco[grid_pix2]

                vmin = 0
                vmax = 2 * met
                if vmax < vmin:
                    vmax = 1

                fig = plt.figure(figsize=(16, 8))
                # matplotlib is doing the mollveide projection
                ax = fig.add_subplot(111, projection='mollweide')


                # rasterized makes the map bitmap while the labels remain vectorial
                # flip longitude to the astro convention
                image = plt.pcolormesh(longitude[::-1], latitude, grid_map,
                                       vmin=vmin, vmax=vmax, rasterized=True, cmap='coolwarm')
                image = plt.pcolormesh(longitude[::-1], latitude, grid_map2,
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
                plt.savefig('sky_plots/%s_%s.png' % (pf, m), bbox_inches='tight',dpi=150)
                plt.figure().clear()
                plt.close()
                plt.cla()
                plt.clf()


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


def create_summary_dist_plots(valid_apo, valid_lco, designmode):
    """
    Create the summary distribution plots
    """
    for i in range(len(valid_apo.columns.names)):
        if 'value' in valid_apo.columns.names[i]:
            for dmode in np.unique(valid_apo['designmode']):
                f, (ax1, ax2, ax3, ax4) = plt.subplots(1, 4, figsize=(40,10))

                ev = eval("(valid_apo['designmode'] == dmode) & (valid_apo[valid_apo.columns.names[i]] >= 0)")
                evd = eval("designmode['label'] == dmode")
                dmode_val = designmode[valid_apo.columns.names[i-1]][evd][0]
                xval_apo = valid_apo[valid_apo.columns.names[i]][ev]
                try:
                    maxx = np.max(xval_apo) * 1.1
                    minn = np.min(xval_apo) * 0.9
                except ValueError:
                    maxx = 1
                    minn = 0
                if minn < 0:
                    minn = 0
                if maxx < 0:
                    maxx = 1
                ev = eval("(valid_lco['designmode'] == dmode) & (valid_lco[valid_lco.columns.names[i]] >= 0)")
                xval_lco = valid_lco[valid_lco.columns.names[i]][ev]
                try:
                    maxx_lco = np.max(xval_lco) * 1.1
                    minn_lco = np.min(xval_lco) * 0.9
                except ValueError:
                    maxx_lco = 1
                    minn_lco = 0
                if minn_lco < 0:
                    minn_lco = 0
                if maxx_lco < 0:
                    maxx_lco = 1
                if minn_lco < minn:
                    minn = minn_lco
                if maxx_lco > maxx:
                    maxx = maxx_lco
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
                plot_hist(ax3, xval_lco,
                          True, True,
                          np.linspace(minn,
                                      maxx,
                                      50),
                          '%s: %s' % (dmode, 'LCO'),
                          'LCO', valid_apo.columns.names[i - 1],
                          'Cumulative Fraction', dmode_val)
                plot_hist(ax4, xval_lco,
                          False, False,
                          np.linspace(minn,
                                      maxx,
                                      50),
                          '%s: %s' % (dmode, 'LCO'),
                          'LCO', valid_apo.columns.names[i - 1],
                          'N', dmode_val)
                plt.savefig('dist_plots/%s_%s.png' % (valid_apo.columns.names[i-1],
                                                      dmode),
                            bbox_inches='tight', dpi=150)
                plt.figure().clear()
                plt.close()
                plt.cla()
                plt.clf()
        if 'pass' in valid_apo.columns.names[i]:
            for dmode in np.unique(valid_apo['designmode']):
                f, (ax1, ax2, ax3, ax4) = plt.subplots(1, 4, figsize=(40,10))

                ev = eval("valid_apo['designmode'] == dmode")
                xval_apo = valid_apo[valid_apo.columns.names[i+1]][ev] - valid_apo[valid_apo.columns.names[i]][ev]
                try:
                    maxx = np.max(xval_apo)
                except ValueError:
                    maxx = 0
                ev = eval("valid_lco['designmode'] == dmode")
                xval_lco = valid_lco[valid_lco.columns.names[i+1]][ev] - valid_lco[valid_lco.columns.names[i]][ev]
                try:
                    maxx_lco = np.max(xval_lco)
                except ValueError:
                    maxx_lco = 0
                if maxx_lco > maxx:
                    maxx = maxx_lco
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
                plot_hist(ax3, xval_lco,
                          True, True,
                          np.arange(0, maxx + 2, 1),
                          '%s: %s' % (dmode, 'LCO'),
                          'LCO', valid_apo.columns.names[i - 1] + ' (N Fibers in Design Failed)',
                          'Cumulative Fraction', np.nan)
                plot_hist(ax4, xval_lco,
                          False, False,
                          np.arange(0, maxx + 2, 1),
                          '%s: %s' % (dmode, 'LCO'),
                          'LCO', valid_apo.columns.names[i - 1] + ' (N Fibers in Design Failed)',
                          'N', np.nan)
                plt.savefig('dist_plots/%s_%s.png' % (valid_apo.columns.names[i-1],
                                                      dmode),
                            bbox_inches='tight', dpi=150)
                plt.figure().clear()
                plt.close()
                plt.cla()
                plt.clf()


def write_html_jinja(valid_apo, valid_lco, designmode,
                     rs_run, mugatu_v, kaiju_v, coordio_v,
                     path):
    """
    Write HTML file using jinja
    """
    env = jinja2.Environment(
        loader=jinja2.PackageLoader("mugatu"),
        autoescape=jinja2.select_autoescape()
    )
    if not isdir('%s/indv_valid' % path):
        mkdir('%s/indv_valid' % path)
    # plots_healpix(valid_apo, valid_lco, designmode)
    # create_summary_dist_plots(valid_apo, valid_lco, designmode)
    mypath = path + '/dist_plots/'
    hist_files = [f for f in listdir(mypath) if isfile(join(mypath, f))]
    hist_files.sort()
    hist_files = np.array(hist_files)

    mypath_sky = path + '/sky_plots/'
    sky_files = [f for f in listdir(mypath_sky) if isfile(join(mypath_sky, f))]
    sky_files.sort()
    sky_files = np.array(sky_files)
    for col, ty in zip(valid_apo.columns.names, valid_apo.columns.formats):
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
                        for m in np.unique(valid_apo['designmode']):
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
                for m in np.unique(valid_apo['designmode']):
                    if o == 'APO':
                        try:
                            value = (100 *
                                     len(valid_apo[col][(valid_apo[col]) & (valid_apo['designmode'] == m)]) /
                                     len(valid_apo[col][(valid_apo['designmode'] == m)]))
                            val = '%.2f' % value + '%'
                        except ZeroDivisionError:
                            val = 'NA'
                    else:
                        try:
                            value = (100 *
                                     len(valid_lco[col][(valid_lco[col]) & (valid_lco['designmode'] == m)]) /
                                     len(valid_lco[col][(valid_lco['designmode'] == m)]))
                            val = '%.2f' % value + '%'
                        except ZeroDivisionError:
                            val = 'NA'
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
    template = env.get_template('main_validation_page.html')
    page = template.render(html_dict)
    fp = open('%s/%s_validation.html' % (path, rs_run), 'w')
    fp.write(page)
    fp.close()
