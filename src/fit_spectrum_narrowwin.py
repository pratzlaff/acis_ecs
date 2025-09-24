import argparse
import os
import sys
# if len(sys.argv) != 7: sys.exit('\nUsage: python *.py epoch ccd fp_temp[see script opts] yes/no[TG] binx biny\n\n')
# (sname,e,ccd,fptv,tg,binx,biny)= sys.argv
# e=f'{int(e):03d}'

if 'ECSID' not in os.environ:
    sys.stderr.write("no ECSID environment variable, exiting\n")
    sys.exit(1)
ecsid=os.environ['ECSID']

###### Editable ######################################

xtra=''  ## default
#xtra='_v772025'
#xtra='_v7232025'
#xtra='_v7242025'  ##1st attempt fix s0/4/5
#xtra='_v7252025'   ##2nd attempt fix s0/4/5


##### test run, only stdout output
test=0
if test: testccd='s0'

if test==1:
    clob='yes'   ## doesn't matter for test runs
    recpars=0    ## output .txt switch
    do2fit=1     ## 2nd fit for more fit accuracy
    xstart=1
    ystart=1
    testerr=1   ## turn conf on for test run
    verbose=9   ## 9=don't suppress stdout, 0= suppress all fit stdout
    allerr=0    ## always =0
else:
    clob='no'    ## checks if [ccd]_ecs.txt exists
    recpars=1    ## output .txt switch
    do2fit=1     ## 2nd fit for more fit accuracy
    xstart=1; ystart=1  ## chipX/Y loop start -- don't change!
    testerr=0    ## always keep this =0
    verbose=0   ## suppress all fit stdout
    allerr=1     ## 1== run errors on everything
    
noerr=0   ## 1= turn off conf completely
toperr=0  ## 1= turn on err only for 256x256y and yl=769, overrides all other settings

### use this for _very_ badly shifted spec
#xtrascl=0
#xtrascl=3  ## 7= 70-210eV, 5= 50-150eV, 3= 30-90eV, 1= 10-30eV

xcnt_thresh_min=2500  ## min counts to fit bkg widths, and run err on some bkg lines
xcnt_thresh_max=100000  ## max counts to fit bkg widths, and run err on some bkg lines, too high == ECS overwhelms bkg
cnt_thresh=200    ## min counts in the .pi spectrum, otherwise skip    
thresh=5e-5  ## norm theshold to run conf vs covar
numcores=16

cnt_thresh=100    ## min counts in the .pi spectrum, otherwise skip    

## don't edit below here ##
from astropy.time import Time as timetime
from astropy.io.fits import open as astro_open
from glob import glob
from matplotlib import pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages
import numpy as np
import shutil
from tempfile import TemporaryDirectory as mktempdir

plt.rc('ytick', direction='in', color='grey', right=1)
plt.rc('xtick', direction='in', color='grey', top=1)
plt.rc('legend', fontsize='small', frameon=1, loc='upper left')
plt.rc('figure', labelsize='small', figsize=[11,6.5])
plt.rc('axes', labelsize='small')
plt.rc('lines', markerfacecolor='none')    ## always non-filled markers
from matplotlib.ticker import ScalarFormatter

# if not test:
#     from matplotlib import use as pltuse
#     pltuse('Agg')

from sherpa.astro import ui
import logging
logger= logging.getLogger('sherpa')
def lerror(): logger.setLevel(logging.ERROR)
def lwarn(): logger.setLevel(logging.WARN)
def linfo(): logger.setLevel(logging.INFO)

pdir_tmp= mktempdir(dir='.')
pdir= pdir_tmp.name
os.environ['PFILES']=f'{pdir};/usr/local/ciao/contrib/param:/usr/local/ciao/param:'

if False:
    enum= int(e)
    fptv=int(fptv)

    if fptv==107: fpt_list=[107,106]
    if fptv==109: fpt_list=[109,108]
    if fptv==111: fpt_list=[111,110]
    if fptv==113: fpt_list=[113,112]
    if fptv==115: fpt_list=[115,114]
    if fptv==117: fpt_list=[117,116]
    if fptv==119: fpt_list=[119,118]  ## -119.19:-117.19
    if fptv==120: fpt_list=[120]      ## -120.19:-119.19
    if fptv==120118: fpt_list=[120,119,118]  ## -120.19:-117.19

dgdir= '/usr/local/ciao/CALDB/data/chandra/acis/det_gain/'
dg_yesCTI= dgdir+'acisD2000-01-29gain_ctiN0008.fits'
dg_noCTI= dgdir+'acisD2000-01-29gainN0005.fits'

pdf = None

#################################################
def do_fit(args):
    global pdf

    ccd = args.ccd
    xbin = args.binx
    ybin = args.biny
    # FIXME
    fptv = 120118
    enum = args.epoch
    e = f'{args.epoch:03d}'

    bin = int(f'{xbin}{ybin}')

    bscls = {
        32:{ 32:1/64, 64:1/32, 128:1/16, 256:1/8 },
        64:{ 256:1/4 },
        128:{ 256:1/2 },
        256:{ 256:1 }
    }
    bscl = bscls[xbin][ybin]

    ccd_list=['i0','i1','i2','i3','s0','s1','s2','s3','s4','s5']

    basedir = '/data/legs/rpete/data/ECS'
    tstr = args.temps.replace(',', '-')
    sfpt = tstr
    spec_dir = f'{basedir}/e{e}/fits/{ecsid}/spec/fpt_{tstr}_{xbin}x{ybin}y'
    fit_dir= f'{basedir}/e{e}/fits/{ecsid}/fits/fpt_{tstr}_{xbin}x{ybin}y'

    if not os.path.exists(spec_dir):
        sys.stderr.write(f'does not exist: {spec_dir}\n')
        sys.exit(1)

    ############ plot block
    def plt_me():
        global pdf

        lwarn(); ui.notice(); ui.group_snr(1,plt_grp)
        fplot = ui.get_fit_plot()
        dmx = fplot.dataplot.x
        daty = fplot.dataplot.y; datye = fplot.dataplot.yerr
        prat = ui.get_ratio_plot()
        ## ungrouped model & component plot
        ui.ungroup()

        pmdl= ui.get_model_component_plot(bkg_mdl+src_mdl)
        mdlx=(pmdl.xlo+pmdl.xhi)/2.; mdly= pmdl.y
        pbkg= ui.get_model_component_plot(bkg_mdl)
        xbkg=(pbkg.xlo+pbkg.xhi)/2.; ybkg= pbkg.y
        yal= ui.get_model_component_plot(alka+alkb).y
        ysi= ui.get_model_component_plot(sika+sikb).y
        ytika= ui.get_model_component_plot(tika1+tika2).y
        ytikb= ui.get_model_component_plot(tikb).y        
        pmnka= ui.get_model_component_plot(mnka1+mnka2); xcomp= (pmnka.xlo+pmnka.xhi)/2.; ymnka= pmnka.y
        ymnkb= ui.get_model_component_plot(mnkb).y
        yaum= ui.get_model_component_plot(auma+aumb).y
        ynik= ui.get_model_component_plot(nika+nikb).y
        yaul= ui.get_model_component_plot(aula+aulb+aula_fs).y                
        
        yaxl= 0.25*min(mdly[np.where( (3<mdlx) & (mdlx<9) & (mdly>0) )])
        yaxh=5*max(mdly)
        xax=[1,14]; yax=[yaxl,yaxh]

        f,ax= plt.subplots(2,1,sharex=1)

        if args.pdf and pdf is None:
            plt_dir= os.path.dirname(args.pdf)
            if plt_dir and not os.path.exists(plt_dir):
                os.makedirs(plt_dir)
            pdf = PdfPages(args.pdf)
        else:
            f.canvas.manager.set_window_title(ccd.upper()+f': x= {sxl}-{sxh}, y= {syl}-{syh}')

        my_lef=0.1; my_rt=0.95; my_bot=0.08; my_top=0.93; my_hs=0; my_ws=0.2
        f.subplots_adjust(left=my_lef, right=my_rt, bottom=my_bot, top=my_top, hspace=my_hs, wspace=my_ws)
        ## fit plot
        lw1=1; lw2=0.5
        ax[0].errorbar(dmx, daty, datye, fmt='+', ecolor='k', mec='k', ms=2, mew=0.5, elinewidth=0.5, zorder=0)
        ax[0].plot(mdlx, mdly, '-', c='orange', lw=lw1,zorder=100)
        ax[0].plot(xbkg, ybkg, '--', c='grey', lw=lw2)
        ax[0].plot(xcomp, yal, '--', c='blue', lw=lw2)
        ax[0].plot(xcomp, ysi, '--', c='grey', lw=lw2)
        ax[0].plot(xcomp, ytika, '--', c='lime', lw=lw2)
        ax[0].plot(xcomp, ytikb, '--', c='lime', lw=lw2)        
        ax[0].plot(xcomp, ymnka, '--', c='red', lw=lw2)
        ax[0].plot(xcomp, ymnkb, '--', c='red', lw=lw2)
        ax[0].plot(xcomp, yaum, '--', c='cyan', lw=lw2)
        ax[0].plot(xcomp, ynik, '--', c='green', lw=lw2)
        ax[0].plot(xcomp, yaul, '--', c='magenta', lw=lw2)                        

        ## plot narrow-window fit bounds
        ax[0].hlines(0.9*yaxh,mn_iglo,mn_ighi,color='blue',lw=2)
        ax[0].hlines(0.9*yaxh,mnb_iglo,mnb_ighi,color='blue',lw=2)
        ax[0].hlines(0.9*yaxh,ti_iglo,ti_ighi,color='blue',lw=2)
        ax[0].hlines(0.9*yaxh,tib_iglo,tib_ighi,color='blue',lw=2)
        ax[0].hlines(0.9*yaxh,al_iglo,al_ighi,color='blue',lw=2)
        ax[0].hlines(0.9*yaxh,si_iglo,si_ighi,color='blue',lw=2)

        
        ## ratio plot
        ax[1].step(prat.x, prat.y, '-', c='orange', lw=lw1,zorder=100, where='mid')
        ax[1].axhline(1, c='black', lw=0.5, ls='--')

        ## titles/axis/etc
        ax[0].set_yscale('log'); ax[0].set_xscale('log')
        ax[0].set_ylim(yax); ax[0].set_xlim(xax)
        ax[0].set_xticks([1,2,3,4,5,6,8,10,12,14])
        ax[0].get_xaxis().set_major_formatter(ScalarFormatter())
        ax[1].set_yscale('linear'); ax[1].set_ylim([0.5,1.5])

        ##
        f.suptitle(f'{ccd.upper()}\t{e}\t{dyear:.3f}+'.expandtabs())
        ax[0].text(0.02,0.94,f'x= {sxl}-{sxh} \ny= {syl}-{syh}', transform=ax[0].transAxes, va='top', ha='left')
        ax[0].text(0.98,0.94,f'-{sfpt} $^\circ$C \n{expo:.1f} ksec \nplot grouping SNR= {plt_grp}', transform=ax[0].transAxes, va='top', ha='right')
        ax[0].set(ylabel='cnt/sec/keV')
        ax[1].set(xlabel='keV', ylabel='data/model')
        ##
        ax[0].text(xax[1]*0.01,yaxl*1.2,'$\Delta$ eV', va='bottom', ha='left',fontsize=10)

        ##
        line1= 1.1; line2= 1.7
        ax[0].text(alka.LineE.val,yaxl*line2,'{:.0f} +{:.0f}/{:.0f}'.
                   format(1e3*(alka.LineE.val-alka_nom),1e3*alka_ehi,1e3*alka_elo), va='bottom', ha='center',fontsize=8)
        ax[0].text(sika.LineE.val,yaxl*line1,'{:.0f} +{:.0f}/{:.0f}'.
                   format(1e3*(sika.LineE.val-sika_nom),1e3*sika_ehi,1e3*sika_elo), va='bottom', ha='center',fontsize=8, c='grey')
        ax[0].text(sikb.LineE.val,yaxl*line2,'{:.0f} +{:.0f}/{:.0f}'.
                   format(1e3*(sikb.LineE.val-sikb_nom),1e3*sikb_ehi,1e3*sikb_elo), va='bottom', ha='center',fontsize=8, c='grey') 
        ##
        ax[0].text(tika1.LineE.val,yaxl*line2,'{:.0f} +{:.0f}/{:.0f}'.
                   format(1e3*(tika1.LineE.val-tika1_nom),1e3*tika1_ehi,1e3*tika1_elo), va='bottom', ha='center',fontsize=8)
        ax[0].text(tikb.LineE.val,yaxl*line1,'{:.0f} +{:.0f}/{:.0f}'.
                   format(1e3*(tikb.LineE.val-tikb_nom),1e3*tikb_ehi,1e3*tikb_elo), va='bottom', ha='center',fontsize=8, c='grey')
        ##
        ax[0].text(mnka1.LineE.val,yaxl*line2,'{:.0f} +{:.0f}/{:.0f}'.
                   format(1e3*(mnka1.LineE.val-mnka1_nom),1e3*mnka1_ehi,1e3*mnka1_elo), va='bottom', ha='center',fontsize=8)
        ax[0].text(mnkb.LineE.val,yaxl*line1,'{:.0f} +{:.0f}/{:.0f}'.
                   format(1e3*(mnkb.LineE.val-mnkb_nom),1e3*mnkb_ehi,1e3*mnkb_elo), va='bottom', ha='center',fontsize=8, c='grey')
        ##
        ax[0].text(auma.LineE.val,yaxl*line1,'{:.0f} +{:.0f}/{:.0f}'.
                   format(1e3*(auma.LineE.val-auma_nom),1e3*auma_ehi,1e3*auma_elo), va='bottom', ha='center',fontsize=8, c='grey')
        ax[0].text(aumb.LineE.val,yaxl*line2,'{:.0f} +{:.0f}/{:.0f}'.
                   format(1e3*(aumb.LineE.val-aumb_nom),1e3*aumb_ehi,1e3*aumb_elo), va='bottom', ha='center',fontsize=8, c='grey')
        ##
        ax[0].text(nika.LineE.val,yaxl*line2,'{:.0f} +{:.0f}/{:.0f}'.
                   format(1e3*(nika.LineE.val-nika_nom),1e3*nika_ehi,1e3*nika_elo), va='bottom', ha='center',fontsize=8, c='grey')
        ax[0].text(nikb.LineE.val,yaxl*line1,'{:.0f} +{:.0f}/{:.0f}'.
                   format(1e3*(nikb.LineE.val-nikb_nom),1e3*nikb_ehi,1e3*nikb_elo), va='bottom', ha='center',fontsize=8, c='grey')
        ##
        ax[0].text(aula.LineE.val,yaxl*line2,'{:.0f} +{:.0f}/{:.0f}'.
                   format(1e3*(aula.LineE.val-aula_nom),1e3*aula_ehi,1e3*aula_elo), va='bottom', ha='center',fontsize=8, c='grey')
        ax[0].text(aulb.LineE.val,yaxl*line1,'{:.0f} +{:.0f}/{:.0f}'.
                   format(1e3*(aulb.LineE.val-aulb_nom),1e3*aulb_ehi,1e3*aulb_elo), va='bottom', ha='center',fontsize=8, c='grey')
            
            
        if args.pdf:
            pdf.savefig(f)

    ##########################################
    linfo()
    ## ccd_id
    ccd_id=int(ccd_list.index(ccd))
    if ccd=='s1_noCTI': ccd_id=5
    if ccd=='s3_noCTI': ccd_id=7
    
    ## open det_gain info
    dg= dg_yesCTI
    if ccd=='s1_noCTI' or ccd=='s3_noCTI': dg=dg_noCTI
    hdulist = astro_open(dg)
    dg_ccd= hdulist[1].data.field('ccd_id')
    dg_xl= hdulist[1].data.field('chipx_min')
    dg_yl= hdulist[1].data.field('chipy_min')
    pha= hdulist[1].data.field('pha')
    nrg= hdulist[1].data.field('energy')

    opentxt=0   ## reset for each ccd run
        
    ####################################################################
    for xl in range(xstart,1025,xbin):
        linfo()
        xh= xl+xbin-1; sxl=str(xl).zfill(4); sxh=str(xh).zfill(4)

        sicol=0
        for yl in range(ystart,1025,ybin):
            linfo()
            yh= yl+ybin-1; syl=str(yl).zfill(4); syh=str(yh).zfill(4)

            pif= f'{spec_dir}/{ccd}_{sxl}-{sxh}x_{syl}-{syh}y.pi'
            tot_cnts=0

            if os.path.exists(pif):

                ui.clean()
                ui.set_method('neldermead'); ui.set_conf_opt('fast',True)
                ui.set_conf_opt('numcores',numcores); ui.set_conf_opt('max_rstat',1e6)
                ##ui.set_stat('cstat')  ## cstat too aggressive to include low count wings

                ui.load_pha(pif)

                ## spec info
                expo= 1e-3*ui.get_data().exposure
                date= ui.get_data().header['DATE-OBS']
                ##tt= time.Time(date); tt.format='decimalyear'; dyear= tt.value
                tt= timetime(date); tt.format='decimalyear'; dyear= tt.value                

                rspx= xbin; rspy= ybin; ccdrsp=ccd
                if ccd=='s1_noCTI': ccdrsp='s1'
                if ccd=='s3_noCTI': ccdrsp='s3'                
                if bin==64256 or bin==128256: rspx= 32

                tmin = np.array([int(t) for t in tstr.split('-')]).max()
                rmf_fpt = f'{tmin-1}-{tmin}' if tmin==120 else f'{tmin-2}-{tmin}'
                rmf_dir= f'{basedir}/acis_response/rmf/pi_{rspx}x{rspy}y_{rmf_fpt}'
                rmf,= glob(f'{rmf_dir}/{ccd}_{sxl}-*x_{syl}-*y.wrmf')

                arf_date= '{}-09-01'.format(round(dyear))
                if dyear>2024: arf_date='2024-09-01'
                arf_dir= '/data/hal9000/acis_response/arf/{}/{}x{}y/HRMA1'.format(arf_date,rspx,rspy)
                arf_dir= f'{basedir}/acis_response/arf/{arf_date}/{rspx}x{rspy}y/HRMA1'
                arf,= glob('{}/{}_{}-*x_{}-*y.warf'.format(arf_dir, ccdrsp, sxl, syl))

                ui.load_rmf(rmf); ui.load_arf(arf)

                lwarn()
                ui.ignore('8.0:')
                tot_cnts=sum(ui.get_counts(filter=1))  ## get_counts MUST be w/out snr grouping!
                print(f'\n\n{e}\t{sfpt}C\t{ccd}\t{xl}xl\t{yl}yl\t{bin}xy\ttot_cnts= {tot_cnts}')
                ui.notice()
                lwarn() if verbose<2 else linfo()

            ## turn on/off fit blocks, debugging/low-counts/BI
            fitbkginitlines=1; fitinit=1; fitmn=1; fitti=1; fital=1; fitsi=1; fitaum=1; fitnika=1; fitnikb=1; fitaula=1; fitaufs=1; fitaulb=1
            if tot_cnts<800 or ccd=='s1' or ccd=='s1_noCTI' or ccd=='s3' or ccd=='s3_noCTI':
                fitaum=0; fitnikb=0; fitaufs=0

            #fitbkginitlines=0; fitinit=0; fitmn=0; fitti=0; fital=0; fitsi=0; fitaum=0; fitnika=0; fitnikb=0; fitaula=0; fitaufs=0; fitaulb=0
            


            if tot_cnts>cnt_thresh:

                fit_grp=2; plt_grp=2
                if tot_cnts>1000: plt_grp=3
                if tot_cnts>4000: plt_grp=4

                
                ## run conf errors switches
                doerr=0  ## reset each run
                if test==0:
                    if (ccd=='i3') or (ccd=='s3'): doerr=1
                    if noerr==1: doerr=0
                if testerr==1: doerr=1
                if bin==256256 and yl==769 and toperr==1: doerr=1
                if allerr: doerr=1




                
                ##################
                ## apprx scaling for initial line norms
                l= np.log(2)/2.737
                nscl= np.exp(-l*(2.737*enum/10.948))

                ## kicks in if xtrascl != 0
                #xshift= -(yl*2e-5*xtrascl + 0.01*xtrascl)

                
                ##################
                # initital Mn-Ka fit:
                #   to find initial gain shift -- to set search windows
                #   to get mna fwhm -- for narrow fit windows
                lwarn()
                pl= ui.xspowerlaw.pl
                mnag = ui.xsgaussian.mnag; mnag_nom= 5.8988
                ui.set_par(mnag.LineE, mnag_nom, min= mnag_nom-0.25, max= mnag_nom+0.1)
                ui.set_par(mnag.sigma, 0.001, min= 0.001, max= 0.1)                
                mnag.norm= bscl*2.5*nscl
                mnbg = ui.xsgaussian.mnbg; mnbg_nom= 6.490
                mnbg.LineE= mnag.LineE + mnbg_nom-mnag_nom; mnbg.sigma= 0.001; mnbg.norm= mnag.norm*0.14
                rsp= ui.get_response(1); tst_mdl= rsp(mnag+mnbg+pl)
                ui.set_full_model(1, tst_mdl); ui.freeze(tst_mdl)
                ui.thaw(mnag,pl)
                ui.group_snr(1,fit_grp); ui.ignore(':5.5,7:')
                lwarn() if verbose<2 else linfo()
                ui.fit(1); lwarn(); ui.notice(); ui.ungroup()
                ig_shift= mnag.LineE.val-mnag_nom
                ## FWHM of Mna
                tmp= ui.get_model_component_plot(mnag)
                mdlx=tmp.xlo; mdly= tmp.y; idx,= np.where(mdly>0.5*max(mdly))
                mna_fwhm= mdlx[idx[-1]] - mdlx[idx[0]]
                ## debug plot
                show_mn_test= 0
                if show_mn_test:
                    ui.notice(); ui.ungroup(); ui.group_snr(1,plt_grp)
                    fplot = ui.get_fit_plot(); dmx= fplot.dataplot.x; daty= fplot.dataplot.y; datye= fplot.dataplot.yerr
                    pmdl= ui.get_model_component_plot(tst_mdl); mdlx=(pmdl.xlo+pmdl.xhi)/2.; mdly= pmdl.y                    
                    f,ax= plt.subplots()
                    ax.errorbar(dmx, daty, datye, fmt='+', ecolor='k', mec='k', ms=2, mew=0.5, elinewidth=0.5)
                    ax.plot(mdlx, mdly, '-', c='orange', lw=1)
                    ax.set_yscale('log'); ax.set_xscale('log')
                    ax.set_xlim([4,9]); plt.show()



                ##################                
                ## inst bkg model
                if dyear<=2005: byear=2000
                if 2005<dyear<2009: byear=2005
                if dyear>=2009: byear=2009
                bccd='i3'
                if ccd=='s3' or ccd=='s3_noCTI': bccd='s3'
                if ccd=='s1' or ccd=='s1_noCTI': bccd='s1'
                ## bkg chipy
                if yl<=256: byl=1
                if 256<yl<=512: byl=257
                if 512<yl<=768: byl=513
                if yl>768: byl=769
                ##
                ui.copy_data(1,10) ## avoids issues
                bgf,= glob(f'{basedir}/acis_bkg/s04_final_models_nolines/{bccd}_{byear}_{byl}-*y.dat')
                ui.load_table_model('bkg_arr',bgf)
                unit_arf = ui.get_arf(10)
                unit_arf.specresp = np.ones_like(unit_arf.specresp)
                rsp_bkg = ui.get_response(10)

                ## BKG lines
                sika = ui.xslorentz.sika; sika_nom= 1.740
                sika.LineE= sika_nom; sika.Width= 0.001
                sika.norm= 9e-3
                ## Si-Kb
                sikb = ui.xslorentz.sikb; sikb_nom= 1.836
                sikb.LineE= sikb_nom; sikb.Width= 0.001
                sikb.norm= 3e-3
                ##
                ui.xslorentz.auma; auma_nom= 2.120
                ui.set_par(auma.LineE, auma_nom, min= auma_nom-0.15, max= auma_nom+0.07)
                ui.set_par(auma.width, 0.01, min=0.005, max=0.1)
                auma.norm= 5e-4
                ui.xslorentz.aumb; aumb_nom= 2.205
                aumb.LineE= auma.LineE+ aumb_nom-auma_nom
                ui.set_par(aumb.width, 0.01, min=0.005, max=0.1)
                aumb.norm= 2e-4
                ##
                ui.xslorentz.nika; nika_nom= 7.478
                val= nika_nom+ig_shift; ui.set_par(nika.LineE, val, min= val-0.15, max= val+0.1)
                ui.set_par(nika.width, 0.01, min=0.005, max=0.1)
                nika.norm= 8e-4
                ui.xslorentz.nikb; nikb_nom= 8.265
                val= nikb_nom+ig_shift; ui.set_par(nikb.LineE, val, min= val-0.15, max= val+0.1)
                ui.set_par(nikb.width, 0.01, min=0.005, max=0.1)
                nikb.norm= 4e-4
                ##
                ui.xslorentz.aula; aula_nom= 9.713
                val= aula_nom+1.2*ig_shift; ui.set_par(aula.LineE, val, min= val-0.15, max= val+0.1)
                ui.set_par(aula.width, 0.01, min=0.005, max=0.1)
                aula.norm= 2e-3
                ##
                ui.xslorentz.aula_fs
                if yl<=256: aula_fs_shif= 0.65
                if 256<yl<=512: aula_fs_shif= 0.66
                if 512<yl<=768: aula_fs_shif= 0.7
                if yl>768: aula_fs_shif= 0.75
                aula_fs.LineE= aula_nom+ aula_fs_shif
                ui.set_par(aula_fs.width, 0.05, min=0.005, max=0.1)
                aula_fs.norm= aula.norm*0.3
                ##
                ui.xslorentz.aulb; aulb_nom= 11.442
                val= aulb_nom+1.3*ig_shift; ui.set_par(aulb.LineE, val, min= val-0.2, max= val+0.1)
                ui.set_par(aulb.width, 0.05, min=0.005, max=0.1)
                aulb.norm= 2e-3


                bkg_mdl=rsp_bkg(bkg_arr+ sika+ sikb+ auma+ aumb+ nika +nikb +aula +aula_fs +aulb)
                ## fit bkg scaling
                ui.set_full_model(10, bkg_mdl)

                ## Fit BKGND - above 8keV
                iglo=8.0; ighi=12.5; ign=f':{iglo},{ighi}:'
                lwarn(); ui.ungroup(); ui.ignore(ign); ui.group_snr(fit_grp)
                lwarn(); ui.freeze(bkg_mdl);
                if fitbkginitlines: ui.thaw(bkg_arr)
                if fitnikb: ui.thaw(nikb)
                if fitaula: ui.thaw(aula)
                if fitaufs: ui.thaw(aula_fs)
                if fitaulb: ui.thaw(aulb)
                ui.thaw(bkg_arr)
                if tot_cnts<xcnt_thresh_min or tot_cnts>xcnt_thresh_max:
                    ui.freeze(auma.width, aumb.width, aula.width, aula_fs.width, aulb.width, nika.width, nikb.width)
                lwarn() if verbose<2 else linfo()
                ui.fit(10); lwarn(); ui.notice(); ui.ungroup(); ui.freeze(bkg_mdl)

                mnmaxwid=0.08
                timaxwid= min([0.030, mnmaxwid])
                ##################                
                # Al-Ka
                alka = ui.xslorentz.alka; alka_nom= 1.4865
                val= alka_nom+0.5*ig_shift; ui.set_par(alka.LineE, val, min= val-0.1, max= val+0.1)
                ui.set_par(alka.Width, 0.005, min= 0.001, max= 0.02)
                val= bscl*1.19*nscl; ui.set_par(alka.norm, val, min= 0.1*val, max= 10*val)
                # Al-Kb
                alkb = ui.xslorentz.alkb
                alkb.LineE= alka.LineE+0.07095; alkb.Width= alka.Width; alkb.norm=  alka.norm*1e-2
                ## Si-Ka
                sika.LineE= alka.LineE+sika_nom-alka_nom; sika.Width= 0.001
                val= bscl*8e-3*nscl; ui.set_par(sika.norm, val, min= 1e-20, max= 0.5*alka.norm.val)
                ## Si-Kb
                sikb.LineE= alka.LineE+sikb_nom-alka_nom; sikb.Width= 0.001
                val= bscl*5e-3*nscl; ui.set_par(sikb.norm, val, min= 1e-20, max= 0.5*alka.norm.val)
                
                ##################
                # Ti-Ka1
                tika1 = ui.xslorentz.tika1; tika1_nom= 4.511
                val= tika1_nom+0.9*ig_shift; ui.set_par(tika1.LineE, val, min= val-0.075, max= val+0.075)
                ui.set_par(tika1.Width, 0.005, min= 0.001, max= timaxwid)
                val= bscl*0.71*nscl; ui.set_par(tika1.norm, val= val, min= 1e-20, max= 1e20)
                # Ti-Ka2
                tika2 = ui.xslorentz.tika2
                tika2.LineE= tika1.LineE-0.00598; tika2.Width= tika1.Width; tika2.norm=  tika1.norm*0.5
                # Ti-Kb
                tikb = ui.xslorentz.tikb; tikb_nom= 4.932
                tikb.LineE= tika1.LineE+tikb_nom-tika1_nom; tikb.Width= tika1.Width; tikb.norm= tika1.norm*0.20

                ##################
                # Mn-Ka1
                mnka1 = ui.xslorentz.mnka1; mnka1_nom= 5.8988
                val= mnka1_nom+ig_shift; ui.set_par(mnka1.LineE, val, min= val-0.075, max= val+0.075)
                val=0.005; ui.set_par(mnka1.Width, val, min= 0.001, max= mnmaxwid)
                val= bscl*2.5*nscl; ui.set_par(mnka1.norm, val= val, min= 1e-20, max= 1e20)
                # Mn-Ka2
                mnka2 = ui.xslorentz.mnka2; mnka2_nom= 5.8877
                mnka2.LineE= mnka1.LineE+mnka2_nom-mnka1_nom; mnka2.Width= mnka1.Width; mnka2.norm= mnka1.norm*0.5
                # Mn-Kb
                mnkb = ui.xslorentz.mnkb; mnkb_nom= 6.490
                mnkb.LineE= mnka1.LineE+mnkb_nom-mnka1_nom; mnkb.Width= mnka1.Width; mnkb.norm= mnka1.norm*0.23
                
                rsp= ui.get_response(1)
                src_mdl= rsp(alka+alkb + tika1+tika2+tikb + mnka1+mnka2+mnkb)
                ui.set_full_model(1, bkg_mdl + src_mdl); ui.freeze(src_mdl,bkg_mdl)

                

                #################                
                if fitinit:
                    ## initial rough fit
                    iglo=1.1; ighi=8.5; ign=':{},{}:'.format(iglo,ighi)
                    lwarn(); ui.notice(); ui.ungroup(); ui.ignore(ign); ui.group_snr(fit_grp)
                    ui.thaw(src_mdl); ui.thaw(sika, sikb, auma, aumb, nika); ui.freeze(sika.width, sikb.width, auma.width, aumb.width, nika.width)
                    lwarn() if verbose<2 else linfo()
                    ui.fit(1)
                    lwarn(); ui.notice(); ui.ungroup(); ui.freeze(src_mdl,bkg_mdl)

                #################
                if fitmn:
                    ##----------
                    ## untie Mn-Kb from Ka
                    val= mnka1.LineE.val+mnkb_nom-mnka1_nom
                    ui.set_par(mnkb.LineE, val, min= val-0.02, max= val+0.02)
                    val= mnka1.width.val; ui.set_par(mnkb.Width, val, min= 0.001, max= mnmaxwid)                    
                    mnkb.norm= mnkb.norm.val

                    ## fit untied Mna/b
                    mn_iglo= mnka1.LineE.val-0.1; mn_ighi=mnkb.LineE.val+0.1; ign=':{},{}:'.format(mn_iglo,mn_ighi)
                    lwarn(); ui.notice(); ui.ungroup(); ui.ignore(ign); ui.group_snr(fit_grp)
                    ui.freeze(src_mdl,bkg_mdl); ui.thaw(mnka1, mnka2, mnkb)
                    lwarn() if verbose<2 else linfo()
                    ui.fit(1); lwarn()
                    
                    ## FWHM of Mna
                    tmp= ui.get_model_component_plot(mnka1+mnka2)
                    mdlx=tmp.xlo; mdly= tmp.y; idx,= np.where(mdly>0.5*max(mdly))
                    mna_fwhm= mdlx[idx[-1]] - mdlx[idx[0]]

                    ## Mna-only tight window fit +- fw 1/4 m ## shape of rsp wonky for large noTG tiles... -- only for tot_cnts>4000
                    fwhm_fact=2
                    if tot_cnts>4000: fwhm_fact= 0.5
                    mn_iglo= mnka1.LineE.val-mna_fwhm*fwhm_fact; mn_ighi=mnka1.LineE.val+mna_fwhm*fwhm_fact; ign=':{},{}:'.format(mn_iglo,mn_ighi)
                    lwarn(); ui.notice(); ui.ungroup(); ui.ignore(ign); ui.group_snr(fit_grp)
                    ui.freeze(src_mdl,bkg_mdl); ui.thaw(mnka1, mnka2)
                    lwarn() if verbose<2 else linfo()
                    ui.fit(1); lwarn()
                    ## Mna errors
                    if doerr==1:
                        try:
                            ui.set_method('levmar')
                            iglo=mnka1.LineE.val-0.8; ighi=mnka1.LineE.val+1.2; ign=':{},{}:'.format(iglo,ighi)
                            lwarn(); ui.notice(); ui.ungroup(); ui.ignore(ign); ui.group_snr(fit_grp)
                            rstat= ui.get_fit_results().rstat; errl=None; errh=None
                            lerror() if verbose<2 else linfo()
                            if rstat < 5:
                                if test: print('\nLineE CONF...')
                                ui.conf(mnka1.LineE); tmp=ui.get_conf_results(); errl=tmp.parmins[0]; errh=tmp.parmaxes[0]
                            if rstat >= 5 or errl is None:
                                if test: print('\nLineE COVAR...')
                                ui.covar(mnka1.LineE); errl=ui.get_covar_results().parmaxes[0]
                                if errl is not None:
                                    errl=-errl; errh=-errl
                            if errl is None: errl=-0.099
                            if errh is None: errh=-errl
                            if errl > -0.001: errl=-0.001
                            if errh < 0.001: errh=0.001
                        except:
                            ui.set_method('neldermead'); errl=-0.099; errh=0.099
                    else:
                        errl=-0.099; errh=0.099
                    if errh==0.099: doerr=0  ## turn off err if Mn is bad
                    mnka1_elo= errl; mnka1_ehi=errh                    

                    ## Mnb-only fit, using fwhm wide window
                    val= mnka1.LineE.val+mnkb_nom-mnka1_nom
                    ui.set_par(mnkb.LineE, val, min= val-0.02, max= val+0.02)
                    val= mnka1.width.val; ui.set_par(mnkb.Width, val, min= 0.001, max= mnmaxwid)
                    val= mnka1.norm.val*0.21; ui.set_par(mnkb.norm, val, min=0.15*mnka1.norm.val, max=0.25*mnka1.norm.val)

                    mnb_iglo= mnkb.LineE.val-mna_fwhm*fwhm_fact*2; mnb_ighi=mnkb.LineE.val+mna_fwhm*fwhm_fact*2; ign=':{},{}:'.format(mnb_iglo,mnb_ighi)
                    lwarn(); ui.notice(); ui.ungroup(); ui.ignore(ign); ui.group_snr(fit_grp)
                    ui.freeze(src_mdl,bkg_mdl); ui.thaw(mnkb); ui.freeze(mnkb.width)
                    lwarn() if verbose<2 else linfo()
                    ui.fit(1); lwarn()
                    ## Mnb errors
                    if doerr==1:
                        try:
                            ui.set_method('levmar')
                            iglo=mnkb.LineE.val-1.2; ighi=mnkb.LineE.val+0.5; ign=':{},{}:'.format(iglo,ighi)
                            lwarn(); ui.notice(); ui.ungroup(); ui.ignore(ign); ui.group_snr(fit_grp)
                            rstat= ui.get_fit_results().rstat; errl=None; errh=None
                            lerror() if verbose<2 else linfo()
                            if rstat < 5:
                                if test: print('\nLineE CONF...')
                                ui.conf(mnkb.LineE); tmp=ui.get_conf_results(); errl=tmp.parmins[0]; errh=tmp.parmaxes[0]
                            if rstat >= 5 or errl is None:
                                if test: print('\nLineE COVAR...')
                                ui.covar(mnkb.LineE); errl=ui.get_covar_results().parmaxes[0]
                                if errl is not None:
                                    errl=-errl; errh=-errl
                            if errl is None: errl=-0.099
                            if errh is None: errh=-errl
                            if errl > -0.001: errl=-0.001
                            if errh < 0.001: errh=0.001
                        except:
                            ui.set_method('neldermead'); errl=-0.099; errh=0.099
                    else:
                        errl=-0.099; errh=0.099
                    mnkb_elo= errl; mnkb_ehi=errh
                    lwarn(); ui.notice(); ui.ungroup(); ui.freeze(src_mdl,bkg_mdl)
                else:
                    mnka1_elo= -0.099; mnka1_ehi=0.099; mnkb_elo= -0.099; mnkb_ehi=0.099                    


                

                #################
                if fitti:

                    ## untie Kb from Ka
                    val= tika1.LineE.val+tikb_nom-tika1_nom
                    ui.set_par(tikb.LineE, val, min= val-0.02, max= val+0.02)
                    val= tika1.Width.val; ui.set_par(tikb.Width, val, min= 0.001, max= timaxwid)
                    tikb.norm= tikb.norm.val

                    
                    ## fit untied Ka/b
                    ti_iglo= tika1.LineE.val-0.1; ti_ighi=tikb.LineE.val+0.1; ign=':{},{}:'.format(ti_iglo,ti_ighi)
                    lwarn(); ui.notice(); ui.ungroup(); ui.ignore(ign); ui.group_snr(fit_grp)
                    ui.freeze(src_mdl,bkg_mdl); ui.thaw(tika1, tika2, tikb)
                    lwarn() if verbose<2 else linfo()
                    ui.fit(1); lwarn()
                    
                    
                    ## Tia-only tight window fwhm
                    ti_iglo= tika1.LineE.val-mna_fwhm*fwhm_fact*2; ti_ighi=tika1.LineE.val+mna_fwhm*fwhm_fact*2; ign=':{},{}:'.format(ti_iglo,ti_ighi)
                    lwarn(); ui.notice(); ui.ungroup(); ui.ignore(ign); ui.group_snr(fit_grp)
                    ui.freeze(src_mdl,bkg_mdl); ui.thaw(tika1, tika2); ui.freeze(tika1.width)
                    lwarn() if verbose<2 else linfo()
                    ui.fit(1); lwarn()
                    ## Tia errors
                    if doerr==1:
                        try:
                            ui.set_method('levmar')
                            iglo=tika1.LineE.val-0.4; ighi=tika1.LineE.val+0.4; ign=':{},{}:'.format(iglo,ighi)
                            lwarn(); ui.notice(); ui.ungroup(); ui.ignore(ign); ui.group_snr(fit_grp)
                            rstat= ui.get_fit_results().rstat; errl=None; errh=None
                            lerror() if verbose<2 else linfo()
                            if rstat < 5:
                                if test: print('\nLineE CONF...')
                                ui.conf(tika1.LineE); tmp=ui.get_conf_results(); errl=tmp.parmins[0]; errh=tmp.parmaxes[0]
                            if rstat >= 5 or errl is None:
                                if test: print('\nLineE COVAR...')
                                ui.covar(tika1.LineE); errl=ui.get_covar_results().parmaxes[0]
                                if errl is not None:
                                    errl=-errl; errh=-errl
                            if errl is None: errl=-0.099
                            if errh is None: errh=-errl
                            if errl > -0.001: errl=-0.001
                            if errh < 0.001: errh=0.001
                        except:
                            ui.set_method('neldermead'); errl=-0.099; errh=0.099
                    else:
                        errl=-0.099; errh=0.099
                    tika1_elo= errl; tika1_ehi=errh                    
                                    
                    ## Tib-only fit, using fw1/4m wide window
                    val= tika1.LineE.val+tikb_nom-tika1_nom
                    ui.set_par(tikb.LineE, val, min= val-0.02, max= val+0.02)
                    tikb.width= tika1.width.val
                    val= 0.201*tika1.norm.val; ui.set_par(tikb.norm, val, min=0.15*tika1.norm.val, max=0.30*tika1.norm.val)

                    tib_iglo= tikb.LineE.val-mna_fwhm*fwhm_fact*2; tib_ighi=tikb.LineE.val+mna_fwhm*fwhm_fact*2; ign=':{},{}:'.format(tib_iglo,tib_ighi)
                    lwarn(); ui.notice(); ui.ungroup(); ui.ignore(ign); ui.group_snr(fit_grp)
                    ui.freeze(src_mdl,bkg_mdl); ui.thaw(tikb); ui.freeze(tikb.width)
                    lwarn() if verbose<2 else linfo()
                    ui.fit(1); lwarn()
                    ## Tib errors
                    if doerr==1:
                        try:
                            ui.set_method('levmar')
                            iglo=tikb.LineE.val-0.4; ighi=tikb.LineE.val+0.4; ign=':{},{}:'.format(iglo,ighi)
                            lwarn(); ui.notice(); ui.ungroup(); ui.ignore(ign); ui.group_snr(fit_grp)
                            rstat= ui.get_fit_results().rstat; errl=None; errh=None
                            lerror() if verbose<2 else linfo()
                            if rstat < 5:
                                if test: print('\nLineE CONF...')
                                ui.conf(tikb.LineE); tmp=ui.get_conf_results(); errl=tmp.parmins[0]; errh=tmp.parmaxes[0]
                            if rstat >= 5 or errl is None:
                                if test: print('\nLineE COVAR...')
                                ui.covar(tikb.LineE); errl=ui.get_covar_results().parmaxes[0]
                                if errl is not None:
                                    errl=-errl; errh=-errl
                            if errl is None: errl=-0.099
                            if errh is None: errh=-errl
                            if errl > -0.001: errl=-0.001
                            if errh < 0.001: errh=0.001
                        except:
                            ui.set_method('neldermead'); errl=-0.099; errh=0.099
                    else:
                        errl=-0.099; errh=0.099
                    tikb_elo= errl; tikb_ehi=errh
                    lwarn(); ui.notice(); ui.ungroup(); ui.freeze(src_mdl,bkg_mdl)
                else:
                    tika1_elo= -0.099; tika1_ehi=0.099; tikb_elo= -0.099; tikb_ehi=0.099


                
                #################
                if fital:
                    ##----------
                    ## untie Si from Al
                    sika.LineE= alka.LineE.val +sika_nom-alka_nom; sika.Width= alka.width.val; sika.norm= sika.norm.val
                    sikb.LineE= alka.LineE.val +sikb_nom-alka_nom; sikb.Width= alka.width.val; sikb.norm= sikb.norm.val                    
                   
                    ## Al-only tight window fwhm
                    al_iglo= alka.LineE.val-mna_fwhm; al_ighi=alka.LineE.val+mna_fwhm; ign=':{},{}:'.format(al_iglo,al_ighi)                    

                    lwarn(); ui.notice(); ui.ungroup(); ui.ignore(ign); ui.group_snr(fit_grp)
                    ui.freeze(src_mdl,bkg_mdl); ui.thaw(alka, alkb); ui.freeze(alka.width)
                    lwarn() if verbose<2 else linfo()
                    ui.fit(1); lwarn()
                    ## Al errors
                    if doerr==1:
                        try:
                            ui.set_method('levmar')
                            iglo=alka.LineE.val-0.4; ighi=alka.LineE.val+0.4; ign=':{},{}:'.format(iglo,ighi)
                            lwarn(); ui.notice(); ui.ungroup(); ui.ignore(ign); ui.group_snr(fit_grp)
                            rstat= ui.get_fit_results().rstat; errl=None; errh=None
                            lerror() if verbose<2 else linfo()
                            if rstat < 5:
                                if test: print('\nLineE CONF...')
                                ui.conf(alka.LineE); tmp=ui.get_conf_results(); errl=tmp.parmins[0]; errh=tmp.parmaxes[0]
                            if rstat >= 5 or errl is None:
                                if test: print('\nLineE COVAR...')
                                ui.covar(alka.LineE); errl=ui.get_covar_results().parmaxes[0]
                                if errl is not None:
                                    errl=-errl; errh=-errl
                            if errl is None: errl=-0.099
                            if errh is None: errh=-errl
                            if errl > -0.001: errl=-0.001
                            if errh < 0.001: errh=0.001
                        except:
                            ui.set_method('neldermead'); errl=-0.099; errh=0.099
                    else:
                        errl=-0.099; errh=0.099
                    alka_elo= errl; alka_ehi=errh
                else:
                    alka_elo= -0.099; alka_ehi=0.099

                    
                ## SiK fits                    
                if fitsi:
                    val= alka.LineE.val +sika_nom-alka_nom; ui.set_par(sika.LineE, val, min= val-0.04, max= val+0.04)
                    if sika.norm.val>1e-5*alka.norm.val: sika.norm.val= 2e-5*alka.norm.val
                    val= sika.norm.val; ui.set_par(sika.norm, val, min=1e-5*alka.norm.val, max=0.5*alka.norm.val)

                    val= alka.LineE.val +sikb_nom-alka_nom; ui.set_par(sikb.LineE, val, min= val-0.04, max= val+0.04)
                    if sikb.norm.val>1e-5*alka.norm.val: sikb.norm.val= 2e-5*alka.norm.val
                    val= sikb.norm.val; ui.set_par(sikb.norm, val, min=1e-5*alka.norm.val, max=0.5*alka.norm.val)
                    
                    si_iglo= sika.LineE.val-0.2; si_ighi=sika.LineE.val+0.2; ign=':{},{}:'.format(si_iglo,si_ighi)
                    lwarn(); ui.notice(); ui.ungroup(); ui.ignore(ign); ui.group_snr(fit_grp)
                    ui.freeze(src_mdl,bkg_mdl); ui.thaw(sika); ui.freeze(sika.width)
                    lwarn() if verbose<2 else linfo()
                    ui.fit(1); lwarn()
                    ## SiKa errors
                    if doerr==1 and tot_cnts>xcnt_thresh_min and tot_cnts<xcnt_thresh_max:
                        ui.set_method('levmar')
                        iglo=sika.LineE.val-0.4; ighi=sika.LineE.val+0.4; ign=':{},{}:'.format(iglo,ighi)
                        lwarn(); ui.notice(); ui.ungroup(); ui.ignore(ign); ui.group_snr(fit_grp)
                        rstat= ui.get_fit_results().rstat; errl=None; errh=None
                        lerror() if verbose<2 else linfo()
                        if rstat < 2.5:
                            if test: print('\nLineE CONF...')
                            ui.conf(sika.LineE); tmp=ui.get_conf_results(); errl=tmp.parmins[0]; errh=tmp.parmaxes[0]
                        else:
                            if test: print('\nLineE COVAR...')
                            ui.covar(sika.LineE); errl=ui.get_covar_results().parmaxes[0]
                            if errl is not None:
                                errl=-errl; errh=-errl
                        if errl is None: errl=-0.099
                        if errh is None: errh=-errl
                        if errl > -0.001: errl=-0.001
                        if errh < 0.001: errh=0.001
                        ui.set_method('neldermead')
                    else:
                        errl=-0.099; errh=0.099
                    sika_elo= errl; sika_ehi=errh

                    si_iglo= sikb.LineE.val-0.2; si_ighi=sikb.LineE.val+0.2; ign=':{},{}:'.format(si_iglo,si_ighi)
                    lwarn(); ui.notice(); ui.ungroup(); ui.ignore(ign); ui.group_snr(fit_grp)
                    ui.freeze(src_mdl,bkg_mdl); ui.thaw(sikb); ui.freeze(sikb.width)
                    lwarn() if verbose<2 else linfo()
                    ui.fit(1); lwarn()
                    ## SiKb errors
                    if doerr==1 and tot_cnts>xcnt_thresh_min and tot_cnts<xcnt_thresh_max:
                        ui.set_method('levmar')
                        iglo=sikb.LineE.val-0.4; ighi=sikb.LineE.val+0.4; ign=':{},{}:'.format(iglo,ighi)
                        lwarn(); ui.notice(); ui.ungroup(); ui.ignore(ign); ui.group_snr(fit_grp)
                        rstat= ui.get_fit_results().rstat; errl=None; errh=None
                        lerror() if verbose<2 else linfo()
                        if rstat < 2.5:
                            if test: print('\nLineE CONF...')
                            ui.conf(sikb.LineE); tmp=ui.get_conf_results(); errl=tmp.parmins[0]; errh=tmp.parmaxes[0]
                        else:
                            if test: print('\nLineE COVAR...')
                            ui.covar(sikb.LineE); errl=ui.get_covar_results().parmaxes[0]
                            if errl is not None:
                                errl=-errl; errh=-errl
                        if errl is None: errl=-0.099
                        if errh is None: errh=-errl
                        if errl > -0.001: errl=-0.001
                        if errh < 0.001: errh=0.001
                        ui.set_method('neldermead')
                    else:
                        errl=-0.099; errh=0.099
                    sikb_elo= errl; sikb_ehi=errh
                    lwarn(); ui.notice(); ui.ungroup(); ui.freeze(src_mdl,bkg_mdl)
                else:
                    sika_elo= -0.099; sika_ehi=0.099; sikb_elo= -0.099; sikb_ehi=0.099
                    


                #################
                if fitaum and tot_cnts>xcnt_thresh_min and tot_cnts<xcnt_thresh_max:
                    ##----------
                    ##untie Au-Mb
                    aumb.LineE= aumb.LineE.val
                    
                    ## Au-Ma
                    aum_iglo= 1.9; aum_ighi= 2.7; ign=':{},{}:'.format(aum_iglo,aum_ighi)
                    lwarn(); ui.notice(); ui.ungroup(); ui.ignore(ign); ui.group_snr(fit_grp)
                    ui.freeze(src_mdl,bkg_mdl); ui.thaw(auma); ui.freeze(auma.width)
                    lwarn() if verbose<2 else linfo()
                    ui.fit(1); lwarn()
                    ## Ni errors
                    if doerr==1 and alka_ehi!=0.099:
                        ui.set_method('levmar')
                        lerror() if verbose<2 else linfo()
                        rstat= ui.get_fit_results().rstat; errl=None; errh=None
                        if rstat < 2.5:
                            if test: print('\nLineE CONF...')
                            ui.conf(auma.LineE); tmp=ui.get_conf_results(); errl=tmp.parmins[0]; errh=tmp.parmaxes[0]
                        else:
                            if test: print('\nLineE COVAR...')
                            ui.covar(auma.LineE); errl=ui.get_covar_results().parmaxes[0]
                            if errl is not None:
                                errl=-errl; errh=-errl
                        if errl is None: errl=-0.099
                        if errh is None: errh=-errl
                        if errl > -0.001: errl=-0.001
                        if errh < 0.001: errh=0.001
                        ui.set_method('neldermead')
                    else:
                        errl=-0.099; errh=0.099
                    auma_elo= errl; auma_ehi=errh                  

                    ## Au-Mb
                    val= aumb.LineE.val; ui.set_par(aumb.LineE, val, min=val-0.04, max= val+0.1)
                    aum_iglo= 1.9; aum_ighi= 2.7; ign=':{},{}:'.format(aum_iglo,aum_ighi)
                    lwarn(); ui.notice(); ui.ungroup(); ui.ignore(ign); ui.group_snr(fit_grp)
                    ui.freeze(src_mdl,bkg_mdl); ui.thaw(aumb); ui.freeze(aumb.width)
                    lwarn() if verbose<2 else linfo()
                    ui.fit(1); lwarn()
                    ## Au-Mb errors
                    if doerr==1 and alka_ehi!=0.099:
                        ui.set_method('levmar')
                        lerror() if verbose<2 else linfo()
                        rstat= ui.get_fit_results().rstat; errl=None; errh=None
                        if rstat < 2.5:
                            if test: print('\nLineE CONF...')
                            ui.conf(aumb.LineE); tmp=ui.get_conf_results(); errl=tmp.parmins[0]; errh=tmp.parmaxes[0]
                        else:
                            if test: print('\nLineE COVAR...')
                            ui.covar(aumb.LineE); errl=ui.get_covar_results().parmaxes[0]
                            if errl is not None:
                                errl=-errl; errh=-errl
                        if errl is None: errl=-0.099
                        if errh is None: errh=-errl
                        if errl > -0.001: errl=-0.001
                        if errh < 0.001: errh=0.001
                        ui.set_method('neldermead')
                    else:
                        errl=-0.099; errh=0.099
                    aumb_elo= errl; aumb_ehi=errh
                    lwarn(); ui.notice(); ui.ungroup(); ui.freeze(src_mdl,bkg_mdl)
                else:
                    auma_elo= -0.099; auma_ehi=0.099; aumb_elo= -0.099; aumb_ehi=0.099                    


                    
                #################
                if fitnika:
                    ##----------
                    ## Ni-Ka +- fwhm
                    ni_iglo= nika.LineE.val-0.2; ni_ighi= nika.LineE.val+0.2; ign=':{},{}:'.format(ni_iglo,ni_ighi)
                    lwarn(); ui.notice(); ui.ungroup(); ui.ignore(ign); ui.group_snr(fit_grp)
                    ui.freeze(src_mdl,bkg_mdl); ui.thaw(nika); ui.freeze(nika.width)
                    lwarn() if verbose<2 else linfo()
                    ui.fit(1); lwarn()
                    ## NiKa errors
                    if doerr==1 and alka_ehi!=0.099 and tot_cnts>xcnt_thresh_min and tot_cnts<xcnt_thresh_max:
                        ui.set_method('levmar')
                        lerror() if verbose<2 else linfo()
                        rstat= ui.get_fit_results().rstat; errl=None; errh=None
                        if rstat < 2.5:
                            if test: print('\nLineE CONF...')
                            ui.conf(nika.LineE); tmp=ui.get_conf_results(); errl=tmp.parmins[0]; errh=tmp.parmaxes[0]
                        else:
                            if test: print('\nLineE COVAR...')
                            ui.covar(nika.LineE); errl=ui.get_covar_results().parmaxes[0]
                            if errl is not None:
                                errl=-errl; errh=-errl
                        if errl is None: errl=-0.099
                        if errh is None: errh=-errl
                        if errl > -0.001: errl=-0.001
                        if errh < 0.001: errh=0.001
                        ui.set_method('neldermead')
                    else:
                        errl=-0.099; errh=0.099
                    nika_elo= errl; nika_ehi=errh
                else:
                    nika_elo= -0.099; nika_ehi=0.099

                ## Ni-Kb +- fwhm                    
                if fitnikb:
                    ni_iglo= nikb.LineE.val-0.2; ni_ighi= nikb.LineE.val+0.2; ign=':{},{}:'.format(ni_iglo,ni_ighi)
                    lwarn(); ui.notice(); ui.ungroup(); ui.ignore(ign); ui.group_snr(fit_grp)
                    ui.freeze(src_mdl,bkg_mdl); ui.thaw(nikb); ui.freeze(nikb.width)
                    lwarn() if verbose<2 else linfo()
                    ui.fit(1); lwarn()
                    ## NiKb errors
                    if doerr==1 and alka_ehi!=0.099 and tot_cnts>xcnt_thresh_min and tot_cnts<xcnt_thresh_max:
                        ui.set_method('levmar')
                        lerror() if verbose<2 else linfo()
                        rstat= ui.get_fit_results().rstat; errl=None; errh=None
                        if rstat < 2.5:
                            if test: print('\nLineE CONF...')
                            ui.conf(nikb.LineE); tmp=ui.get_conf_results(); errl=tmp.parmins[0]; errh=tmp.parmaxes[0]
                        else:
                            if test: print('\nLineE COVAR...')
                            ui.covar(nikb.LineE); errl=ui.get_covar_results().parmaxes[0]
                            if errl is not None:
                                errl=-errl; errh=-errl
                        if errl is None: errl=-0.099
                        if errh is None: errh=-errl
                        if errl > -0.001: errl=-0.001
                        if errh < 0.001: errh=0.001
                        ui.set_method('neldermead')
                    else:
                        errl=-0.099; errh=0.099
                    nikb_elo= errl; nikb_ehi=errh
                    lwarn(); ui.notice(); ui.ungroup(); ui.freeze(src_mdl,bkg_mdl)
                else:
                    nikb_elo= -0.099; nikb_ehi=0.099      



                #################                    
                if fitaula:
                    ## Au-La +- fwhm
                    aula_iglo= aula.LineE.val-0.2; aula_ighi= aula.LineE.val+0.2; ign=':{},{}:'.format(aula_iglo,aula_ighi)
                    lwarn(); ui.notice(); ui.ungroup(); ui.ignore(ign); ui.group_snr(fit_grp)
                    ui.freeze(src_mdl,bkg_mdl); ui.thaw(aula.LineE)
                    lwarn() if verbose<2 else linfo()
                    ui.fit(1); lwarn()
                     ## Au-La errors
                    if doerr==1:
                        ui.set_method('levmar')
                        lerror() if verbose<2 else linfo()
                        rstat= ui.get_fit_results().rstat; errl=None; errh=None
                        if rstat < 2.5:
                            if test: print('\nLineE CONF...')
                            ui.conf(aula.LineE); tmp=ui.get_conf_results(); errl=tmp.parmins[0]; errh=tmp.parmaxes[0]
                        if rstat >= 2.5 or errl is None:
                            if test: print('\nLineE COVAR...')
                            ui.covar(aula.LineE); errl=ui.get_covar_results().parmaxes[0]
                            if errl is not None:
                                errl=-errl; errh=-errl
                        if errl is None: errl=-0.099
                        if errh is None: errh=-errl
                        if errl > -0.001: errl=-0.001
                        if errh < 0.001: errh=0.001
                        ui.set_method('neldermead')
                    else:
                        errl=-0.099; errh=0.099
                    aula_elo= errl; aula_ehi=errh                  
                else:
                    aula_elo= -0.099; aula_ehi=0.099

                    
                if fitaufs:
                    ## Au-La-FrameStore +- fwhm
                    val= aula_nom+ aula_fs_shif; ui.set_par(aula_fs.LineE, val, min= val-0.2, max= val+0.2)
                    val= aula_fs.norm.val; ui.set_par(aula_fs.norm, val, min= aula.norm.val*0.01, max= aula.norm.val)

                    aula_fs_iglo= aula_fs.LineE.val-0.4; aula_fs_ighi= aula_fs.LineE.val+0.4; ign=':{},{}:'.format(aula_fs_iglo,aula_fs_ighi)
                    lwarn(); ui.notice(); ui.ungroup(); ui.ignore(ign); ui.group_snr(fit_grp)
                    ui.freeze(src_mdl,bkg_mdl); ui.thaw(aula_fs); ui.freeze(aula_fs.width)
                    lwarn() if verbose<2 else linfo()
                    ui.fit(1); lwarn()
                     ## Au-La fs errors
                    if doerr==1 and tot_cnts>xcnt_thresh_min and tot_cnts<xcnt_thresh_max:
                        ui.set_method('levmar')
                        lerror() if verbose<2 else linfo()
                        rstat= ui.get_fit_results().rstat; errl=None; errh=None
                        if rstat < 2.5:
                            if test: print('\nLineE CONF...')
                            ui.conf(aula_fs.LineE); tmp=ui.get_conf_results(); errl=tmp.parmins[0]; errh=tmp.parmaxes[0]
                        else:
                            if test: print('\nLineE COVAR...')
                            ui.covar(aula_fs.LineE); errl=ui.get_covar_results().parmaxes[0]
                            if errl is not None:
                                errl=-errl; errh=-errl
                        if errl is None: errl=-0.099
                        if errh is None: errh=-errl
                        if errl > -0.001: errl=-0.001
                        if errh < 0.001: errh=0.001
                        ui.set_method('neldermead')
                    else:
                        errl=-0.099; errh=0.099
                    aula_fs_elo= errl; aula_fs_ehi=errh
                else:
                    aula_fs_elo= -0.099; aula_fs_ehi=0.099


                ## Au-Lb+Lb +- fwhm                    
                if fitaulb:    
                    aulb_iglo= aulb.LineE.val-0.3; aulb_ighi= aulb.LineE.val+0.3; ign=':{},{}:'.format(aulb_iglo,aulb_ighi)
                    lwarn(); ui.notice(); ui.ungroup(); ui.ignore(ign); ui.group_snr(fit_grp)
                    ui.freeze(src_mdl,bkg_mdl); ui.thaw(aulb); ui.freeze(aulb.width)
                    lwarn() if verbose<2 else linfo()
                    ui.fit(1); lwarn()
                    ## Au-Lb errors
                    if doerr==1 and tot_cnts>xcnt_thresh_min and tot_cnts<xcnt_thresh_max:
                        ui.set_method('levmar')
                        lerror() if verbose<2 else linfo()
                        rstat= ui.get_fit_results().rstat; errl=None; errh=None
                        if rstat < 2.5:
                            if test: print('\nLineE CONF...')
                            ui.conf(aulb.LineE); tmp=ui.get_conf_results(); errl=tmp.parmins[0]; errh=tmp.parmaxes[0]
                        else:
                            if test: print('\nLineE COVAR...')
                            ui.covar(aulb.LineE); errl=ui.get_covar_results().parmaxes[0]
                            if errl is not None:
                                errl=-errl; errh=-errl
                        if errl is None: errl=-0.099
                        if errh is None: errh=-errl
                        if errl > -0.001: errl=-0.001
                        if errh < 0.001: errh=0.001
                        ui.set_method('neldermead')
                    else:
                        errl=-0.099; errh=0.099
                    aulb_elo= errl; aulb_ehi=errh
                    lwarn(); ui.notice(); ui.ungroup(); ui.freeze(src_mdl,bkg_mdl)
                else:
                    aulb_elo= -0.099; aulb_ehi=0.099      


                    
                #################
                plt_me()
                
            ########
            ## Record fitpars
            if recpars==1:

                if opentxt==0:
                    ## open output fitpars txt
                    oufecs= f'{fit_dir}/{ccd}_ecs.txt'; fh_ecs= open(oufecs,'w')
                    oufbkg= f'{fit_dir}/{ccd}_bkg.txt'; fh_bkg= open(oufbkg,'w')
                    oufpha= f'{fit_dir}/{ccd}_pha.txt'; fh_pha= open(oufpha,'w')
                    
                    fmt_lines= '\t{:.3f}\t{:.3f}\t{:.3f}\t{:.3f}\t{:.1e}'
                    fmt_ecs='{}\t{}\t{}\t{}' + 5*fmt_lines + '\t{:f}\n'
                    fmt_bkg= '{}\t{}\t{}\t{}' + '\t{:.3f}' + 8*fmt_lines + '\n'
                    fmt_pha= '{}\t{}\t{}\t{}' + 12*'\t{:.1f}' + '\t{:.0f}\n'
                    
                    if not tot_cnts>cnt_thresh:
                        dyear=9999.9; expo=9999
                        
                    fh_ecs.write('#### epoch e{}\tstart= {:.3f}\tksec= {:.1f}\n'.format(e,dyear,expo))
                    fh_bkg.write('#### epoch e{}\tstart= {:.3f}\tksec= {:.1f}\n'.format(e,dyear,expo))
                    fh_pha.write('#### epoch e{}\tstart= {:.3f}\tksec= {:.1f}\n'.format(e,dyear,expo))                    
                    
                    fh_ecs.write('xl\txh\tyl\tyh\tAl_E\tAl_l\tAl_h\tAl_W\tAl_N\tTia_E\tTia_l\tTia_h\tTia_W\tTia_N\tTib_E\tTib_l\tTib_h\tTib_W\tTib_N\tMna_E\tMna_l\tMna_h\tMna_W\tMna_N\tMnb_E\tMnb_l\tMnb_h\tMnb_W\tMnb_N\tstat\n')
                    fh_ecs.flush()
                    fh_bkg.write('xl\txh\tyl\tyh\tamp\tSiKa_E\tSiKa_l\tSiKa_h\tSiKa_W\tSiKa_N\tAuMa_E\tAuMa_l\tAuMa_h\tAuMa_W\tAuMa_N\tAuMb_E\tAuMb_l\tAuMb_h\tAuMb_W\tAuMb_N\tNiKa_E\tNiKa_l\tNiKa_h\tNiKa_W\tNiKa_N\tNiKb_E\tNiKb_l\tNiKb_h\tNiKb_W\tNiKb_N\tAuLa_E\tAuLa_l\tAuLa_h\tAuLa_W\tAuLa_N\tAuFS_E\tAuFS_l\tAuFS_h\tAuFS_W\tAuFS_N\tAuLb_E\tAuLb_l\tAuLb_h\tAuLb_W\tAuLb_N\n')
                    fh_bkg.flush()
                    fh_pha.write('xl\txh\tyl\tyh\tAlKa\tSiKa\tAuMa\tAuMb\tTiKa1\tTiKb\tMnKa1\tMnKb\tNiKa\tAuLa\tAuFS\tAuLb\ttot_cnts\n')
                    
                    fh_pha.flush()                    
                    
                    opentxt+=1

                if tot_cnts>cnt_thresh:
                    ## get det_gain PHA conversion
                    idx_list,= np.where( (dg_ccd==ccd_id) & (dg_xl>=xl) & (dg_xl<xh) & (dg_yl>=yl) & (dg_yl<yh) )
                    dg_alka=[]; dg_tika1=[]; dg_tikb=[]; dg_mnka1=[]; dg_mnkb=[]
                    dg_sika=[]; dg_auma=[]; dg_aumb=[]; dg_nika=[]; dg_aula=[]; dg_aula_fs=[]; dg_aulb=[]
                    for i in range(0,idx_list.size):
                        idx=idx_list[i]
                        nrg0= nrg[idx,:]; pha0= pha[idx,:]
                        nz= np.where(nrg0>0)
                        nrg0= nrg0[nz]; pha0= pha0[nz]
                        ##
                        dg_alka.append(1e3*alka_nom/np.interp(1e3*alka_nom, nrg0, pha0))
                        dg_tika1.append(1e3*tika1_nom/np.interp(1e3*tika1_nom, nrg0, pha0))
                        dg_tikb.append(1e3*tikb_nom/np.interp(1e3*tikb_nom, nrg0, pha0))                        
                        dg_mnka1.append(1e3*mnka1_nom/np.interp(1e3*mnka1_nom, nrg0, pha0))
                        dg_mnkb.append(1e3*mnkb_nom/np.interp(1e3*mnkb_nom, nrg0, pha0))
                        dg_sika.append(1e3*sika_nom/np.interp(1e3*sika_nom, nrg0, pha0))
                        dg_auma.append(1e3*auma_nom/np.interp(1e3*auma_nom, nrg0, pha0))
                        dg_aumb.append(1e3*aumb_nom/np.interp(1e3*aumb_nom, nrg0, pha0))
                        dg_nika.append(1e3*nika_nom/np.interp(1e3*nika_nom, nrg0, pha0))
                        dg_aula.append(1e3*aula_nom/np.interp(1e3*aula_nom, nrg0, pha0))
                        dg_aula_fs.append(1e3*(aula_nom+aula_fs_shif)/np.interp(1e3*(aula_nom+aula_fs_shif), nrg0, pha0))
                        dg_aulb.append(1e3*aulb_nom/np.interp(1e3*aulb_nom, nrg0, pha0))                        

                        
                    dg_alka= np.mean(dg_alka); dg_tika1= np.mean(dg_tika1); dg_tikb= np.mean(dg_tikb); dg_mnka1= np.mean(dg_mnka1); dg_mnkb= np.mean(dg_mnkb)
                    dg_sika= np.mean(dg_sika); dg_auma= np.mean(dg_auma); dg_aumb= np.mean(dg_aumb); dg_nika= np.mean(dg_nika)
                    dg_aula= np.mean(dg_aula); dg_aula_fs= np.mean(dg_aula_fs); dg_aulb= np.mean(dg_aulb)

                    
                    fh_ecs.write(fmt_ecs.format(xl,xh,yl,yh,
                                                alka.LineE.val,alka_elo,alka_ehi,alka.width.val,alka.norm.val,
                                                tika1.LineE.val,tika1_elo,tika1_ehi,tika1.width.val,tika1.norm.val,
                                                tikb.LineE.val,tikb_elo,tikb_ehi,tikb.width.val,tikb.norm.val,
                                                mnka1.LineE.val,mnka1_elo,mnka1_ehi,mnka1.width.val,mnka1.norm.val,
                                                mnkb.LineE.val,mnkb_elo,mnkb_ehi,mnkb.width.val,mnkb.norm.val, ui.calc_stat()))
                    fh_bkg.write(fmt_bkg.format(xl,xh,yl,yh,bkg_arr.ampl.val,
                                                sika.LineE.val,sika_elo,sika_ehi,sika.width.val,sika.norm.val,
                                                auma.LineE.val,auma_elo,auma_ehi,auma.width.val,auma.norm.val,
                                                aumb.LineE.val,aumb_elo,aumb_ehi,aumb.width.val,aumb.norm.val,
                                                nika.LineE.val,nika_elo,nika_ehi,nika.width.val,nika.norm.val,
                                                nikb.LineE.val,nikb_elo,nikb_ehi,nikb.width.val,nikb.norm.val,                                                
                                                aula.LineE.val,aula_elo,aula_ehi,aula.width.val,aula.norm.val,
                                                aula_fs.LineE.val,aula_fs_elo,aula_fs_ehi,aula_fs.width.val,aula_fs.norm.val,
                                                aulb.LineE.val,aulb_elo,aulb_ehi,aulb.width.val,aulb.norm.val))
                    fh_pha.write(fmt_pha.format(xl,xh,yl,yh,
                                                1e3*alka.LineE.val/dg_alka, 1e3*sika.LineE.val/dg_sika, 1e3*auma.LineE.val/dg_auma, 1e3*aumb.LineE.val/dg_aumb, 
                                                1e3*tika1.LineE.val/dg_tika1, 1e3*tikb.LineE.val/dg_tikb, 1e3*mnka1.LineE.val/dg_mnka1, 1e3*mnkb.LineE.val/dg_mnkb,
                                                1e3*nika.LineE.val/dg_nika, 1e3*aula.LineE.val/dg_aula, 1e3*aula_fs.LineE.val/dg_aula_fs,
                                                1e3*aulb.LineE.val/dg_aulb, tot_cnts))
                    fh_ecs.flush(); fh_bkg.flush(); fh_pha.flush()

                else:
                    fh_ecs.write(fmt_ecs.format(xl,xh,yl,yh,
                                            9.999,-0.099,0.099,0.999,9.9e9,
                                            9.999,-0.099,0.099,0.999,9.9e9,
                                            9.999,-0.099,0.099,0.999,9.9e9,
                                            9.999,-0.099,0.099,0.999,9.9e9,
                                            9.999,-0.099,0.099,0.999,9.9e9, np.nan))
                    fh_bkg.write(fmt_bkg.format(xl,xh,yl,yh,999,
                                            9.999,-0.099,0.099,0.999,9.9e9,
                                            9.999,-0.099,0.099,0.999,9.9e9,
                                            9.999,-0.099,0.099,0.999,9.9e9,                                            
                                            9.999,-0.099,0.099,0.999,9.9e9,
                                            9.999,-0.099,0.099,0.999,9.9e9,
                                            9.999,-0.099,0.099,0.999,9.9e9,
                                            9.999,-0.099,0.099,0.999,9.9e9,
                                            9.999,-0.099,0.099,0.999,9.9e9))
                    fh_pha.write(fmt_pha.format(xl,xh,yl,yh,
                                            9999, 9999, 9999, 9999, 9999,
                                            9999, 9999, 9999, 9999, 9999, 9999, 9999, tot_cnts))
                    fh_ecs.flush(); fh_bkg.flush(); fh_pha.flush()


    
    if args.pdf and pdf is not None:
        pdf.close()
    else:
        plt.show()

############### MAIN RUN functions ################


def main():

    parser = argparse.ArgumentParser(
        description='Fit an ECS epoch.'
    )
    parser.add_argument('-p', '--pdf', help='Output PDF file.')
    parser.add_argument('temps', help='e.g., 120,119,118')
    parser.add_argument('binx', type=int)
    parser.add_argument('biny', type=int)
    #parser.add_argument('--foo', default=True, action=argparse.BooleanOptionalAction, help='--no-foo to change.')
    parser.add_argument('epoch', type=int, help='Epoch number.')
    parser.add_argument('ccd',
                        choices=[f'i{i}' for i in range(4)] +[f's{i}' for i in range(6)]
                        )
    args = parser.parse_args()

    do_fit(args)

if __name__ == '__main__':
    main()
