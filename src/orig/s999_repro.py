
#
#
# 2025/1/15 creation
#
#
#
#

## ARGS
import shutil
import sys
if len(sys.argv) != 2: sys.exit('\nUsage: python s*.py epochnum\n\n')
(sname,e)= sys.argv
e=f'{int(e):03d}'

clobber= 1   # 0= prompt if repro directory exists
cleanup= 0   # 0= keep/delete intermediate files
zip= 0       # 1= gzip all download files when finished
notgnocti=1  # 0= turn off noTG and noCTI output

dock=0       # careful turning this off!! on= confirms non-standard options (test/dir/etc)


test=0; xtra=''  ## default caldb

##test=1; testdir='haltg_772025'; xtra=''   ###### custom reprodir output name
#test=1; testdir='haltg_7232025'; xtra=''   ###### custom reprodir output name
#test=1; testdir='haltg_7242025'; xtra=''   ###### custom reprodir output name
#test=1; testdir='haltg_7252025'; xtra=''   ###### custom reprodir output name




## t_gain / cti file defaults
yescti_file='CALDB'; nocti_file='CALDB'
tg_file_yescti='CALDB'; tg_file_nocti='CALDB'



### Edit these lines if needed ###
#xtra='_tgv1'    ## change to add suffix to output repro dir name

## e99 based tg test
#tg_file_yescti='/data/hal9000/2025-02_pyTGain/s04_tgain/acisD2024-10-11t_gain_biN0000.fits'
#tg_file_nocti='/data/hal9000/2025-02_pyTGain/s04_tgain/acisD2024-10-11t_gainN0000.fits'

## e91 based tg test
#tg_file_yescti='/data/hal9000/2025-02_pyTGain/s04_tgain/acisD2022-09-01t_gain_biN0001.fits'
#tg_file_nocti='/data/hal9000/2025-02_pyTGain/s04_tgain/acisD2022-09-01t_gainN0001.fits'


#yescti_file='/data/chandra_caldb/ciao/data/chandra/acis/cti/acisD2020-01-01ctiN0009.fits'
#tg_file_yescti='/data/chandra_caldb/ciao/data/chandra/acis/t_gain/acisD2013-11-01t_gainN0008.fits'
#tg_file_nocti='/data/chandra_caldb/ciao/data/chandra/acis/t_gain/acisD2013-11-01t_gainN0008B.fits'


########## be careful turning this off! #############
if dock:
    ## confirm changes from defaults before running
    ck=''
    if test:
        ck= input('\n\n*** Running "test" version, output: {} \n*** proceed? [y/n]\t'.format(testdir))
    if ck=='n': sys.exit('\naborting\n')

    if notgnocti==0:
        ck= input('\n\n*** noTG and noCTI output turned off \n*** proceed? [y/n]\t')
    if ck=='n': sys.exit('\naborting\n')


    if not test:
        if yescti_file !='CALDB':
            ck= input('\n\n*** yescti_file not CALDB! *** proceed? [y/n]\t')
        if ck=='n': sys.exit('\naborting\n')
        if nocti_file !='CALDB':
            ck= input('\n\n*** nocti_file not CALDB! *** proceed? [y/n]\t')
        if ck=='n': sys.exit('\naborting\n')
        if tg_file_yescti !='CALDB':
            ck= input('\n\n*** tg_file_yescti not CALDB! *** proceed? [y/n]\t')
        if ck=='n': sys.exit('\naborting\n')
        if tg_file_nocti !='CALDB':
            ck= input('\n\n*** tg_file_nocti not CALDB! *** proceed? [y/n]\t')
        if ck=='n': sys.exit('\naborting\n')
        if xtra != '':
            ck= input('\n\n*** repro dir contains additional suffix ---{}---! *** proceed? [y/n]\t'.format(xtra))
        if ck=='n': sys.exit('\naborting\n')






from glob import glob
import os
import sys
import shutil
import numpy as np
from subprocess import *
from ciao_contrib.runtool import dmlist,dmstat,dmkeypar


############# function #### run / log / print info 
def logtool(par_list, silent, log):   ##always prints stderr to LOG
    p=run(par_list, stdout=PIPE, stderr=PIPE, universal_newlines=True)
    print(par_list)
    if p.returncode== 1: sys.exit('\n\tERROR: {}\n'.format(p.stderr))
    if p.returncode== 0 and silent== 0:
        if p.stdout: print('\n\t'+p.stdout)
        if p.stderr: print('\n\t'+p.stderr)
    if log:
        if p.stdout: logger.info('\n\t'+p.stdout)
    if p.stderr: logger.info('\n\t'+p.stderr)
    rtstdout= p.stdout; return rtstdout
#############



## create output repro dir / exit if clobber= 0
from subprocess import getoutput
ciaover= getoutput(['ciaover']).splitlines()
ciaov= ciaover[0].split()[1]
calv= ciaover[2].split()[2]
repro_dir= '/data/hal9000/ecs/data/e{}/ciao{}_caldb{}{}'.format(e,ciaov,calv,xtra)
repro_dir= f'/data/legs/rpete/data/ECS/e{e}/ciao{ciaov}_caldb{calv}'
if test: repro_dir= '/data/hal9000/ecs/data/e{}/{}'.format(e,testdir)

if os.path.exists(repro_dir):
    if clobber: shutil.rmtree(repro_dir)
    if not clobber: sys.exit('\n\n{} exists!\n\n'.format(repro_dir))
os.makedirs(repro_dir)
cwd=os.getcwd()
os.chdir(repro_dir)

import logging
logger = logging.getLogger()
fh_log= logging.FileHandler('repro.log','w'); logger.addHandler(fh_log)
logger.setLevel(logging.INFO)

from tempfile import TemporaryDirectory as mktempdir
tdir_tmp= mktempdir(dir='.'); tdir= tdir_tmp.name
os.environ['PFILES']='{};/usr/local/ciao/contrib/param:/usr/local/ciao/param:'.format(tdir)
lookuptab= '/data/legs/rpete/flight/acis_ecs/data/dmmerge_lookupTab.txt'


## read obsid list
inf= '../support/obsids.dat'
obs_list= np.loadtxt(inf, unpack=1, dtype='str')
nobs= obs_list.size
if nobs==1: obs_list=[obs_list]  ## prevents crash for single obs lists

ccd_list=['i0','i1','i2','i3','s0','s1','s2','s3','s4','s5']


## start / announce info
print('\n**************************************')
print(f'ECS DATA PROCESSING - Epoch e{e}\n')


## setup output EXPOSURE record
hdr1= ['OBSID','evt1','evt2']
for t in range(120,101-1,-1): hdr1.append('{}C'.format(t))
hdr2= ['#####']
for i in range(0,22): hdr2.append('##')
fmt_hdr='{:10s}'+22*'{:9s}'+'\n'
fmt= '{:5d}'+22*'{:9.1f}'+'\n'
for ccd in ccd_list:
    ouf= f'e{e}_{ccd}_exposures.rdb'
    fh= open(ouf,'w'); fh.write(fmt_hdr.format(*hdr1)); fh.write(fmt_hdr.format(*hdr2))
    fh.close()


## unzip all download files
tmp= glob('../download/*gz')
if tmp:
    print('uncompressing download files...\n')
    os.system('gunzip ../download/*gz')

for obs in obs_list:
    obsnum,= np.where(obs_list==obs)[0]
    msg=f'e{e} ObsID {obs}\t{obsnum+1}/{nobs}'; print(msg); logger.info(msg)
    
    ## input file names
    mtl,= glob(f'../download/*{obs}*mtl1.fits*')
    evt1,= glob(f'../download/*{obs}*evt1.fits*')
    msk,= glob(f'../download/*{obs}*msk1.fits*')
    stat,= glob(f'../download/*{obs}*stat1.fits*')
    flt,= glob(f'../download/*{obs}*flt1.fits*')
    pbk= f'../download/{obs}_pbk.fits'; bias= f'../download/{obs}_bias.lis'

    ## output file names    
    fpt_gti=[]; evt2_notg=[]; evt2_notgnocti=[]; evt2_yestg=[]; evt2_yestgnocti=[]
    for t in range(101, 120):
        fpt_gti.append(f'obs{obs}_fpt_gti_{t}.fits')
        evt2_notg.append(f'obs{obs}_evt2_noTG_{t}.fits')
        evt2_notgnocti.append(f'obs{obs}_evt2_noTGnoCTI_{t}.fits')
        evt2_yestg.append(f'obs{obs}_evt2_yesTG_{t}.fits')
        evt2_yestgnocti.append(f'obs{obs}_evt2_yesTGnoCTI_{t}.fits')
    reset= f'{obs}.reset'
    dstrk= f'{obs}.dstrk'; dstrk_noCTI= f'{obs}_noCTI.dstrk'
    aglow= f'{obs}.aglow'
    abb_bpix= f'{obs}.abb_bpix'
    bpix= f'{obs}.bpix'
    obspar= f'obs{obs}.par'
    
    ## file combos: 'noTG', 'yesTG', 'noTGnoCTI', 'yesTGnoCTI'
    tg_list= ['no','yes','no','yes']
    cti_list= ['','','no','no']
    # tg_list= ['yes','yes']
    # cti_list= ['','no']
    if notgnocti==0:
        tg_list= ['yes']; cti_list=['']

    ## always start with punlearn ardlib
    logtool(['punlearn', 'ardlib'], silent=0, log=1)
    
    ## reset evt1
    shutil.copyfile(evt1, reset)    
    logtool(['acis_clear_status_bits', reset], silent=0, log=1)

    ## destreak 
    logtool(['punlearn', 'destreak'], silent=0, log=1)
    logtool(['destreak', reset, dstrk, 'mask=NONE', 'filter=no'], silent=0, log=1)
    ## new badpix
    logtool(['dmmakepar', dstrk, obspar], silent=0, log=1)
    ##msg='\n***Ignore any \'CALDB search returned 10 files\' warning***'; print(msg)
    #
    logtool(['punlearn', 'acis_build_badpix'], silent=0, log=1)
    logtool(['acis_build_badpix', 'obsfile='+obspar, 'pbkfile='+pbk, 'biasfile=@'+bias, 'outfile='+abb_bpix, 'bitflag=00000000000000120021100020022222', 'calibfile=CALDB','mode=h'], silent=0, log=1)
    #
    ##msg='\n*** Ignore any \'Could not get region\' warnings ***'; print(msg)
    logtool(['punlearn', 'acis_find_afterglow'], silent=0, log=1)
    logtool(['acis_find_afterglow', 'infile='+dstrk, 'outfile='+aglow, 'badpixfile='+abb_bpix, 'maskfile='+msk, 'statfile='+stat], silent=0, log=1)
    #
    logtool(['punlearn', 'acis_build_badpix'], silent=0, log=1)
    logtool(['acis_build_badpix', 'obsfile='+obspar, 'pbkfile='+pbk, 'biasfile=NONE', 'outfile='+bpix, 'calibfile='+aglow, 'procbias=no', 'mode=h', 'clobber=yes'], silent=0, log=1)
    # separate s1/s3 for noCTI repro
    logtool(['dmcopy', dstrk+'[ccd_id=5,7][opt mem=999]', dstrk_noCTI], silent=0, log=1)
    

    ## file combos: 'noTG', 'yesTG', 'noTGnoCTI', 'yesTGnoCTI'
    ## a_p_e, (yes/no TG correction)
    for i in range(0,len(tg_list)):
        tg= tg_list[i]; cti= cti_list[i]
        tgcti=f'{tg}TG{cti}CTI' if cti=='no' else f'{tg}TG'
        appcti= 'yes' if cti=='' else 'no'

        if cti=='':
            tg_file= tg_file_yescti; cti_file= yescti_file; ape_inf= dstrk
        if cti=='no':
            tg_file= tg_file_nocti; cti_file= nocti_file; ape_inf= dstrk_noCTI


        ## a_p_e
        logtool(['punlearn', 'acis_process_events'], silent=0, log=1)
        ape= f'{obs}_{tgcti}.ape'
        logtool(['acis_process_events', 'infile='+ape_inf, 'outfile='+ape, 'badpixfile='+bpix, 'acaofffile=NONE', 'mtlfile='+mtl, 'apply_cti='+appcti, 'apply_tgain='+tg, 'ctifile='+cti_file, 'tgainfile='+tg_file, 'check_vf_pha=no', 'stop=tdet', 'pix_adj=NONE', 'eventdef={d:time,l:expno,s:ccd_id,s:node_id,s:chip,s:tdet,d:phas,l:pha,l:pha_ro,f:energy,l:pi,s:fltgrade,s:grade,x:status}'], silent=0, log=1)


        ## new evt2: filter grade/status/ pipeline GTI
        evt= f'{obs}_{tgcti}.evt2'
        logtool(['punlearn', 'dmcopy'], silent=0, log=1)
        logtool(['dmcopy', ape+'[EVENTS][grade=0,2,3,4,6,status=0]', evt], silent=0, log=1)
        logtool(['dmcopy', evt+f'[@{flt}]', evt, 'clobber=yes'], silent=0, log=1)  ## update file in place


    ## filter into FP_TEMP bins
    ## create gti files
    kzero=273.15
    for t in range(120,100,-1):
        fpt_lo= -t-0.19+kzero; fpt_hi= fpt_lo+1
        gti=f'{obs}_{t}.gti'
        logtool(['punlearn', 'dmgti'], silent=0, log=1)
        logtool(['dmgti', mtl, f'userlimit=( (FP_TEMP>={fpt_lo:.2f})&&(FP_TEMP<{fpt_hi:.2f}))', 'outfile='+gti], silent=0, log=1)

    ## check for exptime, filt on FP_TEMP, record EXPOS
    for ccd in ccd_list:
        ccd_id= str(ccd_list.index(ccd))
        expos_arr=[int(obs)]
        dmkeypar(f'{dstrk}[ccd_id={ccd_id}][cols time]', 'EXPOSURE')
        if float(dmkeypar.value) > 10:
            expos_arr.append(float(dmkeypar.value))
            dmkeypar(f'{obs}_yesTG.evt2[ccd_id={ccd_id}][cols time]', 'EXPOSURE')
            if float(dmkeypar.value) > 10:
                expos_arr.append(float(dmkeypar.value))
                for t in range(120,100,-1):
                    gti=f'{obs}_{t}.gti'
                    try:
                        dmstat(gti+'[GTI]')
                        dmkeypar(f'{obs}_yesTG.evt2[@{gti}][ccd_id={ccd_id}][cols time]', 'EXPOSURE')
                        if float(dmkeypar.value) > 10:
                            expos_arr.append(float(dmkeypar.value))                        
                            rangehi=2
                            if (ccd=='s1' or ccd=='s3'): rangehi=4
                            for i in range(0,rangehi):
                                tg= tg_list[i]; cti= cti_list[i]
                                tgcti=f'{tg}TG{cti}CTI' if cti=='no' else f'{tg}TG'
                                evt= f'{obs}_{tgcti}.evt2[@{gti}][ccd_id={ccd_id}]'
                                evt_fpt= f'{obs}_{tgcti}_{t}_{ccd}.evt2'
                                logtool(['dmcopy', evt, evt_fpt], silent=0, log=1)
                        else: expos_arr.append(0.)
                    except: expos_arr.append(0.)
                ouf= f'e{e}_{ccd}_exposures.rdb'
                fh= open(ouf,'a'); fh.write(fmt.format(*expos_arr)); fh.close()                        


                
    ## intermediate cleanup
    if cleanup:
        rm_list= [reset, dstrk, obspar, abb_bpix, aglow, bpix]
        for tmp in rm_list: os.remove(tmp)
        tmp_list= glob('*gti')+glob('*ape')
        for tmp in tmp_list: os.remove(tmp)
        if notgnocti==0: 
            for f in ['yesTG']: os.remove(f'{obs}_{f}.evt2')
        else:
            for f in ['yesTG','noTG','noTGnoCTI','yesTGnoCTI']: os.remove(f'{obs}_{f}.evt2')                



## merge
msg=f'e{e} merging...'; print(msg); logger.info(msg)
for ccd in ccd_list:
    for t in range(120,100,-1):
        rangehi=2
        if (ccd=='s1' or ccd=='s3'): rangehi=4
        if notgnocti==0: rangehi=1
        for i in range(0,rangehi):
            tg= tg_list[i]; cti= cti_list[i]
            tgcti=f'{tg}TG{cti}CTI' if cti=='no' else f'{tg}TG'
            ouf= f'e{e}_{ccd}_{tgcti}_{t}.evt2'
            tmp= glob(f'*_{tgcti}_{t}_{ccd}.evt2')
            if len(tmp)>=1:
                if tgcti=='yesTG':  ## informative, only needed for one tgcti list
                    msg=f'\te{e}_{ccd}_*_{t}.evt2'; print(msg); logger.info(msg)                
                mergef= f'merge_{ccd}_{tgcti}_{t}.lis'
                mergefh= open(mergef,'w')
                for inf in tmp: mergefh.write(inf+'\n')
                mergefh.close()
                logtool(['punlearn', 'dmmerge'], silent=0, log=1)
                logtool(['dmmerge', f'@{mergef}[events][subspace -expno]', 'outfile='+ouf, 'lookupTab='+lookuptab, 'mode=h'], silent=1, log=0)
            if len(tmp)==1: shutil.copy2(tmp[0], ouf)


## move merge lists to separate dir for cleanup            
os.makedirs('merge_lis')
os.system('mv merge_*.lis merge_lis')

## final cleanup
if cleanup:
    for obs in obs_list:
        tmp_list= glob(f'{obs}*')
        for tmp in tmp_list: os.remove(tmp)
if zip:
    print('\ncompressing download files...')
    os.chdir('../download')
    z_list= glob('*evt1.fits') + glob('*bpix1.fits') + glob('*flt1.fits') + glob('*msk1.fits') + glob('*mtl1.fits') + glob('*stat1.fits') + glob('*bias0.fits') + glob('*pbk0.fits')
    for z in z_list: os.system(f'gzip {z}')
    os.chdir(cwd)  ## change back so tmpdir gets deleted
    
msg= f'\ne{e}\tDONE!\nrepro directory: {repro_dir}\n'
print(msg); logger.info(msg)
