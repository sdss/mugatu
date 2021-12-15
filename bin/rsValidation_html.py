import numpy as np
from astropy.io import fits
import matplotlib.pylab as plt
import matplotlib
from matplotlib.patches import Ellipse
from os import listdir
from os.path import isfile, join
from matplotlib.projections.geo import GeoAxes
import healpy as hp
import healpy

plt.rcParams.update({'font.size': 18})
plt.rcParams['savefig.facecolor'] = 'white'


def valid_check_html_table(valid_apo, valid_lco,
                           valid_check, boss_apogee=True):
    """
    Create a summary HTML table for a validation check.
    """
    if boss_apogee:
        boss_col = 'boss_' + valid_check
        ap_col = 'apogee_' + valid_check
        cols = [boss_col, boss_col, ap_col, ap_col]
        obs = ['APO', 'LCO', 'APO', 'LCO']
    else:
        cols = [valid_check, valid_check]
        obs = ['APO', 'LCO']

    header = [''] + list(np.unique(valid_apo['designmode']))
    header = tuple(header)

    rows = []
    for o, col in zip(obs, cols):
        row = [o + ' ' + col]
        for m in np.unique(valid_apo['designmode']):
            if o == 'APO':
                try:
                    value = (100 *
                             len(valid_apo[col][(valid_apo[col]) & (valid_apo['designmode'] == m)]) /
                             len(valid_apo[col][(valid_apo['designmode'] == m)]))
                    row.append('%.2f' % value + '%')
                except ZeroDivisionError:
                    row.append('NA')
            else:
                try:
                    value = (100 *
                             len(valid_lco[col][(valid_lco[col]) & (valid_lco['designmode'] == m)]) /
                             len(valid_lco[col][(valid_lco['designmode'] == m)]))
                    row.append('%.2f' % value + '%')
                except ZeroDivisionError:
                    row.append('NA')
        rows.append(row)
    html_table = """<table style="width:90%">
  <tr>"""
    for h in header:
        html_table += '<th>%s</th>\n' % h
    html_table += '</tr>\n'
    for row in rows:
        html_table += '<tr>\n'
        for r in row:
            html_table += '<td>%s</td>\n' % r
        html_table += '</tr>\n'
    html_table += '</table>\n'
    return html_table


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


def write_designmode_table(designmode,
                           dmode=None):
    """
    Write an HTML table for the designmodes for this
    run.
    """
    if dmode is None:
        html_table = """<table style="width:90%">
      <tr>"""
        for h in designmode.columns.names:
            html_table += '<th>%s</th>\n' % h
        html_table += '</tr>\n'
        for i in range(len(designmode)):
            html_table += '<tr>\n'
            for h in designmode.columns.names:
                html_table += '<td>%s</td>\n' % designmode[h][i]
            html_table += '</tr>\n'
        html_table += '</table>\n'
    else:
        header = [''] + list(designmode.label)
        rows = []
        for h in designmode.columns.names:
            if dmode in h:
                row = [h]
                row += list(designmode[h])
                rows.append(row)
        html_table = """<table style="width:90%">
      <tr>"""
        for h in header:
            html_table += '<th>%s</th>\n' % h
        html_table += '</tr>\n'
        for row in rows:
            html_table += '<tr>\n'
            for r in row:
                html_table += '<td>%s</td>\n' % r
            html_table += '</tr>\n'
        html_table += '</table>\n'
    return html_table


def write_html_file(valid_apo, valid_lco, designmode,
                    rs_run, mugatu_v, kaiju_v, coordio_v):
    """
    Create the summary distribution plots
    """
    dmode_descriptions = {}
    dmode_descriptions['decolide'] = ('A check on if the FPS grid was able '
                                      'to be decolided if any collisions existed.')
    dmode_descriptions['bright_safety'] = ('A check if any fibers in the grid were '
                                           'too close to a bright neighbor such '
                                           'that the total flux down the fiber would '
                                           'be greater than the bright_limit_targets '
                                           'for the designmode. This check only '
                                           'considers a subset of very bright stars '
                                           'with G < 13 (for BOSS) and H < 7 (for '
                                           'APOGEE) in the FOV.')
    dmode_descriptions['all_targets_assigned'] = ('A check if all of the requested target '
                                                  'assignments able to be made given '
                                                  'the physical constraints on the '
                                                  'robots.')
    dmode_descriptions['no_collisions'] = ('A check if any assignments were removed in order '
                                           'to allow for the decolision of the grid.')
    dmode_descriptions['n_skies_min'] = ('A check if the design has the minimum '
                                         'number of skies for the designmode.')
    dmode_descriptions['min_skies_fovmetric'] = ('A check if the skies in the design are '
                                                 'well distributed in the FOV of the FPS.')
    dmode_descriptions['n_stds_min'] = ('A check if the design has the minimum '
                                        'number of standards for the designmode.')
    dmode_descriptions['min_stds_fovmetric'] = ('A check if the standards in the design are '
                                                'well distributed in the FOV of the FPS.')
    dmode_descriptions['stds_mags'] = ('A check if the standards in the design are '
                                       'within the magntiude limits of a designmode.')
    dmode_descriptions['bright_limit_targets'] = ('A check if the science targets in the design are '
                                                  'within the magntiude limits of a designmode.')
    dmode_descriptions['sky_neighbors_targets'] = ('A check if any fibers in the grid were '
                                                   'too close to a bright neighbor such '
                                                   'that the total flux down the fiber would '
                                                   'be greater than the bright_limit_targets '
                                                   'for the designmode. This check '
                                                   'considers all stars down to the minimum '
                                                   'magnitude limit in the FOV.')
    with open('index.html', 'w') as f:
        f.write("""<!DOCTYPE html>
<html>
<style>
table, th, td {
  border:1px solid black;
}
</style>
<body>""")

        # plots_healpix(valid_apo, valid_lco, designmode)
        # create_summary_dist_plots(valid_apo, valid_lco, designmode)

        f.write('<h1>%s Validation Results</h1>' % rs_run)
        intro = ('Below are links to summaries of each design validation '
                 'criteria used for the %s run of robostrategy. These '
                 'designs were validated using v%s of mugatu, v%s '
                 'of kaiju and v%s of coordio.' % (rs_run,
                                                   mugatu_v,
                                                   kaiju_v,
                                                   coordio_v))
        f.write('<p>%s</p>' % intro)

        f.write('<h2>Design Validations</h2>')

        mypath = 'dist_plots/'
        hist_files = [f for f in listdir(mypath) if isfile(join(mypath, f))]
        hist_files.sort()
        hist_files = np.array(hist_files)

        mypath_sky = 'sky_plots/'
        sky_files = [f for f in listdir(mypath_sky) if isfile(join(mypath_sky, f))]
        sky_files.sort()
        sky_files = np.array(sky_files)
        for col, ty in zip(valid_apo.columns.names, valid_apo.columns.formats):
            if ty == 'L' and 'boss' in col:
                valid_check = col[5:]
                f.write('<p><a href="%s.html">%s</a>: %s</p>\n' % (valid_check,
                                                                   valid_check,
                                                                   dmode_descriptions[valid_check]))
                with open('%s.html' % valid_check, 'w') as fs:
                    fs.write("""<!DOCTYPE html>
<html>
<style>
table, th, td {
  border:1px solid black;
}
</style>
<body>""")
                    fs.write('<h1>%s</h1>' % valid_check)
                    html_tab = write_designmode_table(designmode,
                                                      dmode=valid_check)
                    fs.write('<h2>DesignModes for Run</h2>')
                    fs.write(html_tab)
                    fs.write('<h2>Percent Designs Pass</h2>')
                    html_tab = valid_check_html_table(valid_apo, valid_lco,
                                                      valid_check)
                    fs.write(html_tab)
                    for h in hist_files:
                        if valid_check in h:
                            x = str(h)
                            x = x.replace('.png', '')
                            x = x.replace('dist_plots/', '')
                            x = x.split('_')
                            y = ''
                            for i in range(len(x) - 2):
                                y += x[i]
                                if i < len(x) - 2 - 1:
                                    y += '_'
                            y += ' (%s)' % (x[-2] + '_' + x[-1])
                            fs.write('<h3>%s</h3>\n' % y)
                            if h in sky_files:
                                fs.write('<img src="%s" ' % (mypath_sky + h) + 'width="60%">\n')
                            fs.write('<img src="%s" ' % (mypath + h) + 'width="80%">\n')
                    fs.write("""</body>
</html>""")
            elif ty == 'L' and 'apogee' not in col:
                valid_check = col
                f.write('<p><a href="%s.html">%s</a>: %s</p>\n' % (valid_check,
                                                                   valid_check,
                                                                   dmode_descriptions[valid_check]))
                with open('%s.html' % valid_check, 'w') as fs:
                    fs.write("""<!DOCTYPE html>
<html>
<style>
table, th, td {
  border:1px solid black;
}
</style>
<body>""")
                    fs.write('<h1>%s</h1>' % valid_check)
                    fs.write('<h2>Percent Designs Pass</h2>')
                    html_tab = valid_check_html_table(valid_apo, valid_lco,
                                                      valid_check, boss_apogee=False)
                    fs.write(html_tab)
                    for h in hist_files:
                        if valid_check in h:
                            x = str(h)
                            x = x.replace('.png', '')
                            x = x.replace('dist_plots/', '')
                            x = x.split('_')
                            y = ''
                            for i in range(len(x) - 2):
                                y += x[i]
                                if i < len(x) - 2 - 1:
                                    y += '_'
                            y += '(%s)' % (x[-2] + '_' + x[-1])
                            fs.write('<h3>%s</h3>\n' % y)
                            if h in sky_files:
                                fs.write('<img src="%s" ' % (mypath_sky + h) + 'width="60%">\n')
                            fs.write('<img src="%s" ' % (mypath + h) + 'width="80%">\n')
                    fs.write("""</body>
</html>""")


        f.write("""</body>
</html>""")
