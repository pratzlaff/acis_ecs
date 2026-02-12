import argparse
from datetime import datetime,timedelta
from glob import glob
from matplotlib import pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages
import numpy as np
import os
import sys

ecsid=os.environ['ECSID']
srcdir=os.path.dirname(__file__)
datadir=os.popen(srcdir+'/datadir').read()

def plt_prefs():
    plt.rc('ytick', direction='in', color='gray', right=1)
    plt.rc('xtick', direction='in', color='gray', top=1)
    plt.rc('legend', fontsize='small', frameon=1, loc='center left')
    plt.rc('figure', labelsize='small', figsize=[9,6.5])
    plt.rc('axes', labelsize='small')
    plt.rc('lines', markerfacecolor='none')    ## always non-filled markers

# given year and month, return epoch number
def epoch(year, month):
    # first epoch started in 2000-02
    nmonths=(int(year)-2000)*12 + month-2
    return int(nmonths/3)+1

# given epoch number, return year and month of beginning
def epoch_start(epoch):
    nmonths = (epoch-1)*3
    year = 2000+int(nmonths/12)
    month = 2+nmonths%12
    return year, month

xtra=''; xtra_label='CALDB'
#xtra='_v772025'; xtra_label='v772025'
#xtra='_v7232025'; xtra_label='v7232025'
#xtra='_v7242025'; xtra_label='v7242025'
#xtra='_v7252025'; xtra_label='v7252025'

def mkplot(args, data, ccd):

    f,ax= plt.subplots(4,1,sharex=1,sharey=1)
    f.canvas.manager.set_window_title(ccd.upper())
    my_lef=0.1; my_rt=0.95; my_bot=0.08; my_top=0.93; my_hs=0; my_ws=0.2
    f.subplots_adjust(left=my_lef, right=my_rt, bottom=my_bot, top=my_top, hspace=my_hs, wspace=my_ws)
    c=['cyan','red','orange','lime']; fmt=['+','o','s','d']

    ylim = {
        'Al':(1.40, 1.53),
        'Ti':(4.41, 4.61),
        'Ml':(5.75, 5.95),
    }.get(args.line)

    td15= timedelta(days=5)
    td05= timedelta(days=2.5)
    xoffset = ( -td15, -td05, td05, td15 )

    labels = ('1:256y', '256:512y', '512:768y', '768:1024y')
    for j in range(data['ediff'].shape[0]):
        label = None
        for k in range(data['ediff'].shape[1]):
            mask = ~ np.isnan(data['ediff'][j,k])
            if j==0:
                label = labels[k]
            x = data['years'][mask] - xoffset[k]
            y = data['ediff'][j,k][mask]
            yerr = (data['lo'][j,k][mask], data['hi'][j,k][mask])

            ax[j].errorbar(x, y, yerr=yerr, fmt=fmt[k], ecolor=c[k], mec=c[k], mfc='none', ms=3, label=label)
            if ((ccd=='i3' and j==3 and k==3) or (ccd=='s3' and j==0 and k==2)):
                ax[j].errorbar(x, y, yerr=yerr, fmt='d', ecolor='black', mec='black', mfc='none', ms=6)

    ax[1].set_xlabel('ECS epoch midpoint')

    global nom
    for i in range(0,4):
        lbl1='+- 0.5%' if i==0 else ''
        ax[i].axhline(0,ls='--',lw=0.5,c='gray',label=lbl1)
        ax[i].axhline(-1e3*0.005*nom,ls='--',lw=0.5,c='gray')
        ax[i].axhline(1e3*0.005*nom,ls='--',lw=0.5,c='gray') 

    vfact=0.015 if args.line=='Mn' else 0.03
    vmin= nom - vfact*nom; vmax= nom + vfact*nom
    vmin = {
        'Al':-40,
        'Ti':-80,
        'Mn':-80,
    }.get(args.line)
    vmax = {
        'Al':40,
        'Ti':80,
        'Mn':80,
    }.get(args.line)

    vlines = ['2021-11-01', '2022-09-15', '2023-03-01', '2023-12-15', '2024-10-15']

    for i in range(4):
        for j, date in enumerate(vlines):
            label = None
            if j==0 and i==1:
                label='tgain begins'
            ax[i].vlines(datetime.strptime(date,'%Y-%m-%d'),vmin,vmax,color='blue',label=label,lw=0.5,ls='--')

    global temps
    f.suptitle(f'ECS {args.line}-Ka LineE -{temps.max()}.19:-{temps.min()-1}.19C -- {ccd.upper()} -- 256x256y \t\t{ecsid}'.expandtabs())
    f.supylabel('LineE shift [eV]')
    for i in range(4):
        ax[i].set_ylabel(f'node{i}')

    #ax[0].set_xlim(datetime.strptime('2020-10-01',dfmt), datetime.strptime('2025-05-01',dfmt))

    ylim = {
        'Al':[-50,50],
        'Ti':[-100,100],
        'Mn':[-100,100],
    }.get(args.line)

    #ax[0].set_ylim(ylim)
    z=ax[0].get_yticks()   ##; print(z)
    newz=[z[1],z[2],z[3]]; ax[0].set_yticks(newz)
    ax[0].legend()
    ax[1].legend(loc='lower left', frameon=False)

    plt.tight_layout()
    return f

def year_mid_epoch(e):
    year, month = epoch_start(e)
    date_str = f'{year:04d}-{month+1:02d}-{15}'
    return datetime.strptime(date_str, '%Y-%m-%d')

def plt_spec(args):
    global nom
    nom = {
        'Al':1.487,
        'Ti':4.511,
        'Mn':5.895,
    }.get(args.line)

    plt_prefs()

    global temps
    temps = np.array([ int(t) for t in args.temps.split(',')])

    fpt=args.temps.replace(',', '-')
    if args.pdf:
        pdf = PdfPages(args.pdf)

    for ccd in args.det:

        xi = [1,257,513,769]
        yi = [1,257,513,769]
        d = { 'e':np.empty((len(xi), len(yi), len(args.epochs))) }
        d['lo'] = d['e'].copy()
        d['hi'] = d['e'].copy()

        years = np.array([year_mid_epoch(e) for e in args.epochs])

        for i, e in enumerate(args.epochs):
            inf=f'{datadir}/e{e:03d}/fits/{ecsid}/fits/fpt_{fpt}_256x256y/{ccd}_ecs.txt'
            try:
                (xl,yl,al,al_lo,al_hi,ti,ti_lo,ti_hi,mn,mn_lo,mn_hi,stat)= np.loadtxt(inf, skiprows=2, unpack=1, usecols=[0,2,4,5,6,9,10,11,19,20,21,-1])

                li, li_lo, li_hi = {
                    'Al':(al, al_lo, al_hi),
                    'Ti':(ti, ti_lo, ti_hi),
                    'Mn':(mn, mn_lo, mn_hi),
                }.get(args.line)

                if (len(xl) != 16):
                    raise RuntimeError(f'len(xl)=f{len(xl)} in f{inf}')

                for j in range(len(yi)):
                    n = [int(np.where((xl==xi[k]) & (yl==yi[j]))[0][0]) for k in range(len(xi))]
                    for k in range(len(xi)):
                        d['e'][k,j,i] = li[n[k]]
                        d['lo'][k,j,i] = li_lo[n[k]]
                        d['hi'][k,j,i] = li_hi[n[k]]
            except:
                raise

        d['lo'] = np.abs(d['lo'])
        filt = d['lo'] >= 0.08
        d['ediff'] = 1e3 * (d['e'] - nom)
        d['lo'] *= 1e3
        d['hi'] *= 1e3

        for key in d:
            d[key][filt] = np.nan
        d['years'] = years

        fig = mkplot(args, d, ccd)
        if args.pdf:
            pdf.savefig(fig)

    if args.pdf:
        pdf.close()
    else:
        plt.show()

def main():
    parser = argparse.ArgumentParser(
        description='Plot fit results'
    )
    parser.add_argument('-p', '--pdf', help='Output PDF file.')
    parser.add_argument('temps', help='e.g., 120,119,118')
    parser.add_argument('line', choices=['Al','Ti','Mn'])
    parser.add_argument('epochs', type=int, nargs='+')
    parser.add_argument('--det', choices=['i0','i1','i2','i3','s0','s1','s2','s3','s4','s5'], nargs='+')
    args = parser.parse_args()
    if args.det is None:
        args.det = ['i0','i1','i2','i3','s0','s1','s2','s3','s4','s5']
    plt_spec(args)

if __name__ == '__main__':
    main()
