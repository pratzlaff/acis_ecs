import argparse
from datetime import datetime,timedelta
from glob import glob
import numpy as np
import os

from matplotlib import pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages

def plt_prefs():
    plt.rc('ytick', direction='in', color='gray', right=1)
    plt.rc('xtick', direction='in', color='gray', top=1)
    plt.rc('legend', fontsize='small', frameon=1, loc='upper left')
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

def plt_spec(args):

    plt_prefs()

    temps = np.array([ int(t) for t in args.temps.split(',')])
    e_list = args.epochs
    dfmt= '%Y-%m-%d'
    fpt=args.temps.replace(',', '-')
    if args.pdf:
        pdf = PdfPages(args.pdf)

    ## plot
    def plt_stuff():

        f,ax= plt.subplots(4,1,sharex=1,sharey=1)
        f.canvas.manager.set_window_title(ccd.upper())
        my_lef=0.1; my_rt=0.95; my_bot=0.08; my_top=0.93; my_hs=0; my_ws=0.2
        f.subplots_adjust(left=my_lef, right=my_rt, bottom=my_bot, top=my_top, hspace=my_hs, wspace=my_ws)
        c=['cyan','red','orange','lime']; fmt=['+','o','s','d']

        if mtl=='Al': ylim=[1.40,1.53]
        if mtl=='Mn': ylim=[5.75,5.95]

        td15= timedelta(days=5)
        td05= timedelta(days=2.5)
    
        i=0
        lbl1=''; lbl2=''; lbl3=''; lbl4=''
        if i==0:
            lbl1='1:256y'; lbl2='256:512y'; lbl3='512:768y'; lbl4='768:1024y'
        ax[i].errorbar(d_good[f0y1]-td15,li0y1,yerr=[li0y1_lo,li0y1_hi], fmt=fmt[0], ecolor=c[0], mec=c[0], mfc='none', ms=3, label=lbl1)
        ax[i].errorbar(d_good[f0y2]-td05,li0y2,yerr=[li0y2_lo,li0y2_hi], fmt=fmt[1], ecolor=c[1], mec=c[1], mfc='none', ms=3, label=lbl2)
        ax[i].errorbar(d_good[f0y3]+td05,li0y3,yerr=[li0y3_lo,li0y3_hi], fmt=fmt[2], ecolor=c[2], mec=c[2], mfc='none', ms=3, label=lbl3)
        ax[i].errorbar(d_good[f0y4]+td15,li0y4,yerr=[li0y4_lo,li0y4_hi], fmt=fmt[3], ecolor=c[3], mec=c[3], mfc='none', ms=3, label=lbl4)
        i=1
        ax[i].errorbar(d_good[f1y1]-td15,li1y1,yerr=[li1y1_lo,li1y1_hi], fmt=fmt[0], ecolor=c[0], mec=c[0], mfc='none', ms=3)
        ax[i].errorbar(d_good[f1y2]-td05,li1y2,yerr=[li1y2_lo,li1y2_hi], fmt=fmt[1], ecolor=c[1], mec=c[1], mfc='none', ms=3)
        ax[i].errorbar(d_good[f1y3]+td05,li1y3,yerr=[li1y3_lo,li1y3_hi], fmt=fmt[2], ecolor=c[2], mec=c[2], mfc='none', ms=3)
        ax[i].errorbar(d_good[f1y4]+td15,li1y4,yerr=[li1y4_lo,li1y4_hi], fmt=fmt[3], ecolor=c[3], mec=c[3], mfc='none', ms=3)
        i=2
        ax[i].errorbar(d_good[f2y1]-td15,li2y1,yerr=[li2y1_lo,li2y1_hi], fmt=fmt[0], ecolor=c[0], mec=c[0], mfc='none', ms=3)
        ax[i].errorbar(d_good[f2y2]-td05,li2y2,yerr=[li2y2_lo,li2y2_hi], fmt=fmt[1], ecolor=c[1], mec=c[1], mfc='none', ms=3)
        ax[i].errorbar(d_good[f2y3]+td05,li2y3,yerr=[li2y3_lo,li2y3_hi], fmt=fmt[2], ecolor=c[2], mec=c[2], mfc='none', ms=3)
        ax[i].errorbar(d_good[f2y4]+td15,li2y4,yerr=[li2y4_lo,li2y4_hi], fmt=fmt[3], ecolor=c[3], mec=c[3], mfc='none', ms=3)
        i=3
        ax[i].errorbar(d_good[f3y1]-td15,li3y1,yerr=[li3y1_lo,li3y1_hi], fmt=fmt[0], ecolor=c[0], mec=c[0], mfc='none', ms=3)
        ax[i].errorbar(d_good[f3y2]-td05,li3y2,yerr=[li3y2_lo,li3y2_hi], fmt=fmt[1], ecolor=c[1], mec=c[1], mfc='none', ms=3)
        ax[i].errorbar(d_good[f3y3]+td05,li3y3,yerr=[li3y3_lo,li3y3_hi], fmt=fmt[2], ecolor=c[2], mec=c[2], mfc='none', ms=3)
        ax[i].errorbar(d_good[f3y4]+td15,li3y4,yerr=[li3y4_lo,li3y4_hi], fmt=fmt[3], ecolor=c[3], mec=c[3], mfc='none', ms=3)

        for i in range(0,4):
            lbl1='+- 0.5%' if i==0 else ''
            ax[i].axhline(0,ls='--',lw=0.5,c='gray',label=lbl1)
            ax[i].axhline(-1e3*0.005*nom,ls='--',lw=0.5,c='gray')
            ax[i].axhline(1e3*0.005*nom,ls='--',lw=0.5,c='gray') 

        vfact=0.015 if mtl=='Mn' else 0.03
        vmin= nom - vfact*nom; vmax= nom + vfact*nom
        if mtl=='Al':
            vmin= -40; vmax=40
        if mtl=='Mn':
            vmin= -80; vmax=80

        vlines = {
            ''          : ['2021-11-01', '2022-08-02'],
            '_v772025'  : ['2021-11-01', '2022-09-15', '2024-01-15'],
            '_v7232025' : ['2021-11-01', '2022-09-15', '2023-03-01', '2023-12-15', '2024-10-15'],
        }
        for key in '_v7242025', '_v7252025':
            vlines[key] = vlines['_v7232025']

        for i in range(4):
            label=None
            for date in vlines[xtra]:
                if i==1:
                    label=date+' gain'
                ax[i].vlines(datetime.strptime(date,dfmt),vmin,vmax,color='blue',label=label,lw=0.5,ls='--')
            
        f.suptitle(f'ECS {mtl}-Ka LineE -{temps.max()}:-{temps.min()-1}C -- {ccd.upper()} -- 256x256y \t\t{xtra_label}'.expandtabs())
        for i in range(4):
            ax[i].set_ylabel(f'node{i}')

        #ax[0].set_xlim(datetime.strptime('2020-10-01',dfmt), datetime.strptime('2025-05-01',dfmt))

        if mtl=='Al': ylim=[-50,50]
        if mtl=='Mn': ylim=[-100,100]    

        #ax[0].set_ylim(ylim)
        z=ax[0].get_yticks()   ##; print(z)
        newz=[z[1],z[2],z[3]]; ax[0].set_yticks(newz)
        ax[0].legend()
        ax[1].legend()    

        plt.tight_layout()
        if args.pdf:
            pdf.savefig(f)

    ## end of plt_stuff()

    ccd_list=['i0','i1','i2','i3','s0','s1','s2','s3','s4','s5']
    for ccd in ccd_list:

        for mtl in [args.line]: #'Al']: #,'Mn']:
            al_nom= 1.487; mn_nom=5.895
            if mtl=='Al': nom= al_nom
            if mtl=='Mn': nom= mn_nom

            ##
            al0=[]; al_lo0=[]; al_hi0=[]
            al1=[]; al_lo1=[]; al_hi1=[]
            al2=[]; al_lo2=[]; al_hi2=[]
            al3=[]; al_lo3=[]; al_hi3=[]
            ##
            li0y1=[]; li0y1_lo=[]; li0y1_hi=[]; li1y1=[]; li1y1_lo=[]; li1y1_hi=[]
            li2y1=[]; li2y1_lo=[]; li2y1_hi=[]; li3y1=[]; li3y1_lo=[]; li3y1_hi=[]
            li0y2=[]; li0y2_lo=[]; li0y2_hi=[]; li1y2=[]; li1y2_lo=[]; li1y2_hi=[]
            li2y2=[]; li2y2_lo=[]; li2y2_hi=[]; li3y2=[]; li3y2_lo=[]; li3y2_hi=[]
            li0y3=[]; li0y3_lo=[]; li0y3_hi=[]; li1y3=[]; li1y3_lo=[]; li1y3_hi=[]
            li2y3=[]; li2y3_lo=[]; li2y3_hi=[]; li3y3=[]; li3y3_lo=[]; li3y3_hi=[]
            li0y4=[]; li0y4_lo=[]; li0y4_hi=[]; li1y4=[]; li1y4_lo=[]; li1y4_hi=[]
            li2y4=[]; li2y4_lo=[]; li2y4_hi=[]; li3y4=[]; li3y4_lo=[]; li3y4_hi=[]
            ##
            e_good=[]
            d_good=[]
        

            for i,e in enumerate(e_list):
                year, month = epoch_start(e)
                date_str = f'{year:04d}-{month+1:02d}-{15}'
                e_date = datetime.strptime(date_str, dfmt)
                ## grab Al/Mn
                inf= f'/data/legs/rpete/data/ECS/fits/ciao4.17.0_caldb4.12.2/e{e:03d}/fpt_{fpt}_256x256y_yesTG/{ccd}_ecs.txt'
                inf= f'/data/legs/rpete/data/ECS/e{e:03d}/fits/ciao4.17.0_caldb4.12.2/fits/fpt_{fpt}_256x256y/{ccd}_ecs.txt'
                if os.path.exists(inf):
                    (xl_tmp,yl_tmp,al,al_lo,al_hi,mn,mn_lo,mn_hi,stat)= np.loadtxt(inf, skiprows=2, unpack=1, usecols=[0,2,4,5,6,19,20,21,-1])

                    if mtl=='Al':
                        li= al; li_lo= al_lo; li_hi= al_hi
                    if mtl=='Mn':
                        li= mn; li_lo= mn_lo; li_hi= mn_hi

                    if len(xl_tmp)==16:
                        xl=xl_tmp; yl=yl_tmp
                        e_good= np.append(e_good,e)
                        d_good= np.append(d_good,e_date)

                        for yi in [1,257,513,769]:
                            n0,= np.where( (xl==1)   & (yl==yi) )
                            n1,= np.where( (xl==257) & (yl==yi) )
                            n2,= np.where( (xl==513) & (yl==yi) )
                            n3,= np.where( (xl==769) & (yl==yi) )
                            ##
                            if yi==1:
                                li0y1= np.append(li0y1, li[n0]); li0y1_lo= np.append(li0y1_lo, li_lo[n0]); li0y1_hi= np.append(li0y1_hi, li_hi[n0])
                                li1y1= np.append(li1y1, li[n1]); li1y1_lo= np.append(li1y1_lo, li_lo[n1]); li1y1_hi= np.append(li1y1_hi, li_hi[n1])
                                li2y1= np.append(li2y1, li[n2]); li2y1_lo= np.append(li2y1_lo, li_lo[n2]); li2y1_hi= np.append(li2y1_hi, li_hi[n2])
                                li3y1= np.append(li3y1, li[n3]); li3y1_lo= np.append(li3y1_lo, li_lo[n3]); li3y1_hi= np.append(li3y1_hi, li_hi[n3])        
                            if yi==257:
                                li0y2= np.append(li0y2, li[n0]); li0y2_lo= np.append(li0y2_lo, li_lo[n0]); li0y2_hi= np.append(li0y2_hi, li_hi[n0])
                                li1y2= np.append(li1y2, li[n1]); li1y2_lo= np.append(li1y2_lo, li_lo[n1]); li1y2_hi= np.append(li1y2_hi, li_hi[n1])
                                li2y2= np.append(li2y2, li[n2]); li2y2_lo= np.append(li2y2_lo, li_lo[n2]); li2y2_hi= np.append(li2y2_hi, li_hi[n2])
                                li3y2= np.append(li3y2, li[n3]); li3y2_lo= np.append(li3y2_lo, li_lo[n3]); li3y2_hi= np.append(li3y2_hi, li_hi[n3])       
                            if yi==513:
                                li0y3= np.append(li0y3, li[n0]); li0y3_lo= np.append(li0y3_lo, li_lo[n0]); li0y3_hi= np.append(li0y3_hi, li_hi[n0])
                                li1y3= np.append(li1y3, li[n1]); li1y3_lo= np.append(li1y3_lo, li_lo[n1]); li1y3_hi= np.append(li1y3_hi, li_hi[n1])
                                li2y3= np.append(li2y3, li[n2]); li2y3_lo= np.append(li2y3_lo, li_lo[n2]); li2y3_hi= np.append(li2y3_hi, li_hi[n2])
                                li3y3= np.append(li3y3, li[n3]); li3y3_lo= np.append(li3y3_lo, li_lo[n3]); li3y3_hi= np.append(li3y3_hi, li_hi[n3])   
                            if yi==769:
                                li0y4= np.append(li0y4, li[n0]); li0y4_lo= np.append(li0y4_lo, li_lo[n0]); li0y4_hi= np.append(li0y4_hi, li_hi[n0])
                                li1y4= np.append(li1y4, li[n1]); li1y4_lo= np.append(li1y4_lo, li_lo[n1]); li1y4_hi= np.append(li1y4_hi, li_hi[n1])
                                li2y4= np.append(li2y4, li[n2]); li2y4_lo= np.append(li2y4_lo, li_lo[n2]); li2y4_hi= np.append(li2y4_hi, li_hi[n2])
                                li3y4= np.append(li3y4, li[n3]); li3y4_lo= np.append(li3y4_lo, li_lo[n3]); li3y4_hi= np.append(li3y4_hi, li_hi[n3])   

            ## abs for lower err
            li0y1_lo= np.abs(li0y1_lo); li1y1_lo= np.abs(li1y1_lo); li2y1_lo= np.abs(li2y1_lo); li3y1_lo= np.abs(li3y1_lo)
            li0y2_lo= np.abs(li0y2_lo); li1y2_lo= np.abs(li1y2_lo); li2y2_lo= np.abs(li2y2_lo); li3y2_lo= np.abs(li3y2_lo)
            li0y3_lo= np.abs(li0y3_lo); li1y3_lo= np.abs(li1y3_lo); li2y3_lo= np.abs(li2y3_lo); li3y3_lo= np.abs(li3y3_lo)
            li0y4_lo= np.abs(li0y4_lo); li1y4_lo= np.abs(li1y4_lo); li2y4_lo= np.abs(li2y4_lo); li3y4_lo= np.abs(li3y4_lo)


            ## filt for err < 80eV
            f0y1= np.where(li0y1_lo < 0.08); f1y1= np.where(li1y1_lo < 0.08); f2y1= np.where(li2y1_lo < 0.08); f3y1= np.where(li3y1_lo < 0.08)
            f0y2= np.where(li0y2_lo < 0.08); f1y2= np.where(li1y2_lo < 0.08); f2y2= np.where(li2y2_lo < 0.08); f3y2= np.where(li3y2_lo < 0.08)
            f0y3= np.where(li0y3_lo < 0.08); f1y3= np.where(li1y3_lo < 0.08); f2y3= np.where(li2y3_lo < 0.08); f3y3= np.where(li3y3_lo < 0.08)
            f0y4= np.where(li0y4_lo < 0.08); f1y4= np.where(li1y4_lo < 0.08); f2y4= np.where(li2y4_lo < 0.08); f3y4= np.where(li3y4_lo < 0.08)
            ##
            li0y1= 1e3*(li0y1[f0y1]-nom); li0y2= 1e3*(li0y2[f0y2]-nom); li0y3= 1e3*(li0y3[f0y3]-nom); li0y4= 1e3*(li0y4[f0y4]-nom)
            li1y1= 1e3*(li1y1[f1y1]-nom); li1y2= 1e3*(li1y2[f1y2]-nom); li1y3= 1e3*(li1y3[f1y3]-nom); li1y4= 1e3*(li1y4[f1y4]-nom)
            li2y1= 1e3*(li2y1[f2y1]-nom); li2y2= 1e3*(li2y2[f2y2]-nom); li2y3= 1e3*(li2y3[f2y3]-nom); li2y4= 1e3*(li2y4[f2y4]-nom)
            li3y1= 1e3*(li3y1[f3y1]-nom); li3y2= 1e3*(li3y2[f3y2]-nom); li3y3= 1e3*(li3y3[f3y3]-nom); li3y4= 1e3*(li3y4[f3y4]-nom)    
            li0y1_lo= 1e3*li0y1_lo[f0y1]; li0y2_lo= 1e3*li0y2_lo[f0y2]; li0y3_lo= 1e3*li0y3_lo[f0y3]; li0y4_lo= 1e3*li0y4_lo[f0y4]
            li1y1_lo= 1e3*li1y1_lo[f1y1]; li1y2_lo= 1e3*li1y2_lo[f1y2]; li1y3_lo= 1e3*li1y3_lo[f1y3]; li1y4_lo= 1e3*li1y4_lo[f1y4]
            li2y1_lo= 1e3*li2y1_lo[f2y1]; li2y2_lo= 1e3*li2y2_lo[f2y2]; li2y3_lo= 1e3*li2y3_lo[f2y3]; li2y4_lo= 1e3*li2y4_lo[f2y4]
            li3y1_lo= 1e3*li3y1_lo[f3y1]; li3y2_lo= 1e3*li3y2_lo[f3y2]; li3y3_lo= 1e3*li3y3_lo[f3y3]; li3y4_lo= 1e3*li3y4_lo[f3y4]    
            li0y1_hi= 1e3*li0y1_hi[f0y1]; li0y2_hi= 1e3*li0y2_hi[f0y2]; li0y3_hi= 1e3*li0y3_hi[f0y3]; li0y4_hi= 1e3*li0y4_hi[f0y4]
            li1y1_hi= 1e3*li1y1_hi[f1y1]; li1y2_hi= 1e3*li1y2_hi[f1y2]; li1y3_hi= 1e3*li1y3_hi[f1y3]; li1y4_hi= 1e3*li1y4_hi[f1y4]
            li2y1_hi= 1e3*li2y1_hi[f2y1]; li2y2_hi= 1e3*li2y2_hi[f2y2]; li2y3_hi= 1e3*li2y3_hi[f2y3]; li2y4_hi= 1e3*li2y4_hi[f2y4]
            li3y1_hi= 1e3*li3y1_hi[f3y1]; li3y2_hi= 1e3*li3y2_hi[f3y2]; li3y3_hi= 1e3*li3y3_hi[f3y3]; li3y4_hi= 1e3*li3y4_hi[f3y4]    

            plt_stuff()

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
    parser.add_argument('line', choices=['Al', 'Mn'])
    parser.add_argument('epochs', type=int, nargs='+')
    args = parser.parse_args()

    plt_spec(args)

if __name__ == '__main__':
    main()


