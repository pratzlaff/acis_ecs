
import numpy as np
import os
from tempfile import TemporaryDirectory as mktempdir
from glob import glob

from sherpa.astro import ui,xspec
from matplotlib import pyplot as plt

plt.rc('ytick', direction='in', color='gray', right=1)
plt.rc('xtick', direction='in', color='gray', top=1)
plt.rc('legend', fontsize='small', frameon=1, loc='upper left')
plt.rc('figure', labelsize='small', figsize=[9,6.5])
plt.rc('axes', labelsize='small')
plt.rc('lines', markerfacecolor='none')    ## always non-filled markers






xtra=''; xtra_label='CALDB'
#xtra='_v772025'; xtra_label='v772025'
#xtra='_v7232025'; xtra_label='v7232025'
#xtra='_v7242025'; xtra_label='v7242025'
#xtra='_v7252025'; xtra_label='v7252025'





fpt='120-119-118'
estart=87
e_list= np.arange(estart,101+1)


## dates weighted to exptime (eye-balled based on obsids per epoch) / mid-epoch
##x9= [          83        84           85           86           87           88           89           90           91           92           93           94           95           96           97           98           99          100          101 ]
##date_list= ['2020-11'... 
date_list_84_101= [     '2021-01-01','2021-04-01','2021-06-01','2021-10-01','2021-12-25','2022-03-15','2022-06-15','2022-09-15','2022-12-15','2023-03-01','2023-06-15','2023-09-15','2024-01-25','2024-02-05','2024-07-01','2024-10-11','2024-12-15','2025-03-15']
date_list_87_101= [                                            '2021-10-01','2021-12-25','2022-03-15','2022-06-15','2022-09-15','2022-12-15','2023-03-01','2023-06-15','2023-09-15','2024-01-25','2024-02-05','2024-07-01','2024-10-11','2024-12-15','2025-03-15']


if estart==87: date_list= date_list_87_101



from datetime import datetime,timedelta
dfmt= '%Y-%m-%d'; tmp=[]
for d in date_list: tmp.append(datetime.strptime(d, dfmt))

date_list= tmp




############################## plot
def plt_stuff():
    f,ax= plt.subplots(4,1,sharex=1,sharey=1)
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


    if xtra=='_v7232025' or xtra=='_v7242025' or xtra=='_v7252025':
        for i in range(0,4):
            lbl1=''; lbl2=''; lbl3=''; lbl4=''; lbl5=''
            if i==1:
                lbl1='11/1/2021 gain'; lbl2='9/15/2022 gain'; lbl3='3/1/2023 gain'; lbl4='12/15/2023 gain'; lbl5='10/15/2024 gain'
            ax[i].vlines(datetime.strptime('2021-11-01',dfmt),vmin,vmax,color='blue',label=lbl1,lw=0.5,ls='--')
            ax[i].vlines(datetime.strptime('2022-09-15',dfmt),vmin,vmax,color='blue',label=lbl2,lw=0.5,ls='--')
            ax[i].vlines(datetime.strptime('2023-03-01',dfmt),vmin,vmax,color='blue',label=lbl3,lw=0.5,ls='--')
            ax[i].vlines(datetime.strptime('2023-12-15',dfmt),vmin,vmax,color='blue',label=lbl4,lw=0.5,ls='--')                
            ax[i].vlines(datetime.strptime('2024-10-15',dfmt),vmin,vmax,color='blue',label=lbl5,lw=0.5,ls='--')        

    if xtra=='_v772025':
        for i in range(0,4):
            lbl1=''; lbl2=''; lbl3=''; lbl4=''; lbl5=''
            if i==1:
                lbl1='11/1/2021 gain'; lbl2='9/15/2022 gain'; lbl5='1/15/2024 gain'
            ax[i].vlines(datetime.strptime('2021-11-01',dfmt),vmin,vmax,color='blue',label=lbl1,lw=0.5,ls='--')
            ax[i].vlines(datetime.strptime('2022-09-15',dfmt),vmin,vmax,color='blue',label=lbl2,lw=0.5,ls='--')
            ax[i].vlines(datetime.strptime('2024-1-15',dfmt),vmin,vmax,color='blue',label=lbl5,lw=0.5,ls='--')        

    if xtra=='':
        for i in range(0,4):
            lbl1=''; lbl2=''; lbl3=''; lbl4=''; lbl5=''
            if i==1:
                lbl1='11/1/2021 gain'; lbl2='8/2/2022 gain'
            ax[i].vlines(datetime.strptime('2021-11-01',dfmt),vmin,vmax,color='blue',label=lbl1,lw=0.5,ls='--')
            ax[i].vlines(datetime.strptime('2022-08-02',dfmt),vmin,vmax,color='blue',label=lbl2,lw=0.5,ls='--')



            
    f.suptitle('ECS {}-Ka LineE -120:-117C -- {} -- 256x256y \t\t{}'.expandtabs().format(mtl,ccd.upper(),xtra_label))
    ax[0].set_ylabel('node0'); ax[1].set_ylabel('node1'); ax[2].set_ylabel('node2'); ax[3].set_ylabel('node3')



    
    ax[0].set_xlim(datetime.strptime('2020-10-01',dfmt), datetime.strptime('2025-05-01',dfmt))

    if mtl=='Al': ylim=[-50,50]
    if mtl=='Mn': ylim=[-100,100]    

    ax[0].set_ylim(ylim)
    z=ax[0].get_yticks()   ##; print(z)
    newz=[z[1],z[2],z[3]]; ax[0].set_yticks(newz)
    ax[0].legend()
    ax[1].legend()    

    f.savefig('{}/tst_{}_{}_{}.pdf'.format(tdir,ccd,xtra_label,mtl))
    plt.close()
    #plt.show()

    




######################################################

tdir_tmp= mktempdir(dir='.')  ## tmpdir created in ./  otherwise /tmp/
tdir= tdir_tmp.name
    
ccd_list=['i0','i1','i2','i3','s0','s1','s1_noCTI','s2','s3','s3_noCTI','s4','s5']
for ccd in ccd_list:

    for mtl in ['Al','Mn']:
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
            e_date= date_list[i]
            ## grab Al/Mn
            inf= '/data/legs/rpete/data/ECS/fits/e{}/fpt_{}_256x256y_yesTG{}/{}_ecs.txt'.format(e,fpt,xtra,ccd)
            if os.path.exists(inf):
                (xl_tmp,yl_tmp,al,al_lo,al_hi,mn,mn_lo,mn_hi)= np.loadtxt(inf, skiprows=2, unpack=1, usecols=[0,2,4,5,6,19,20,21])

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


    #os.system('pdfjam --landscape {}/*pdf -o s999_plts/{}_120-117_e84-101_{}.pdf'.format(tdir,ccd,xtra_label))

    

for mtl in ['Al','Mn']:
    os.system('pdfjam --landscape {}/tst_*_{}.pdf -o s999_plts/{}_120-117_e{}-101_{}.pdf'.format(tdir,mtl,mtl,estart,xtra_label))    




