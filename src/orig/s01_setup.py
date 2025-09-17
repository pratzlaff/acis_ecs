import sys
import numpy as np

if len(sys.argv) != 4: sys.exit('\nUsage: python *.py fp_temp[see script opts] yes/no[TG] [bin]\n\n')
(sname,fptv,tg,bin)= sys.argv
fptv=int(fptv); bin= int(bin)


import os
from subprocess import getoutput
from ciao_contrib.runtool import dmmerge,dmextract,dmstat
from glob import glob
from tempfile import TemporaryDirectory as mktempdir

cwd=os.getcwd()



############## EDITABLE ###################################
## mk_merge
## OUTPUT evt dir naming ex. ==> evt2/e1e2e3/fpt_120-119-118/
## EVT2 file naming: [ccd].evt2
##
## mk_spec
## OUTPUT evt dir naming ex. ==> spec/e1e2e3/fpt_120-119-118/32x256y/
## EVT2 file naming: [ccd].evt2


rdir='LATEST'; xtra='' ## DEFAULT
#rdir='haltg_772025'; xtra='_v772025'  ## xtra= suffix for merge/spec dirs
#rdir='haltg_7232025'; xtra='_v7232025'  ## xtra= suffix for merge/spec dirs
#rdir='haltg_7242025'; xtra='_v7242025'  ## xtra= suffix for merge/spec dirs
#rdir='haltg_7252025'; xtra='_v7252025'  ## xtra= suffix for merge/spec dirs



## individual epochs e100-e1, combine fp_temp as needed
#main_e_list=[]
#for i in np.arange(100,17-1,-1): main_e_list.append('e{}'.format(i))
#main_e_list= main_e_list + ['e16c', 'e16b', 'e16a']
#for i in np.arange(15,1-1,-1): main_e_list.append('e{}'.format(i))

main_e_list=[]
#for i in np.arange(91,101+1,1): main_e_list.append('e{}'.format(i))
for i in np.arange(87,101+1,1): main_e_list.append('e{}'.format(i))
main_e_list=['e083', 'e094', 'e099']
main_e_list=['e094', 'e099']
print(main_e_list)



## individual epochs, manual lists
#main_e_list=['e97','e98']
#main_e_list=['e101']


## combine epoch combos, manual lists
#main_e_list=[['e25','e26'],['e27','e28']]





###### DO NOT EDIT #####
if fptv==107: fpt_list=[107,106]
if fptv==109: fpt_list=[109,108]
if fptv==111: fpt_list=[111,110]
if fptv==113: fpt_list=[113,112]
if fptv==115: fpt_list=[115,114]
if fptv==117: fpt_list=[117,116]
if fptv==119: fpt_list=[119,118]  ## -119.19:-117.19
if fptv==120: fpt_list=[120]      ## -120.19:-119.19
if fptv==120118: fpt_list=[120,119,118]  ## -120.19:-117.19


if bin==256256:
    xbin=256; ybin=256; sbin='{}x{}y'.format(xbin,ybin)
if bin==32256:
    xbin=32; ybin=256; sbin='{}x{}y'.format(xbin,ybin)
if bin==32128:
    xbin=32; ybin=128; sbin='{}x{}y'.format(xbin,ybin)
if bin==3264: 
    xbin=32; ybin=64; sbin='{}x{}y'.format(xbin,ybin)


print(main_e_list); print('\n')
## fpt
sfpt='{}'.format(fpt_list[0])
if len(fpt_list)>1:
    for fpt in fpt_list[1:]: sfpt= sfpt+'-{}'.format(fpt)

basedir= '/data/hal9000/ECS_fits'
basedir= '/data/legs/rpete/data/ECS'
ccd_list=['i0','i1','i2','i3','s0','s1','s2','s3','s4','s5','s1_noCTI','s3_noCTI']
#ccd_list=['i0','i1','i2','i3']   ## okay to use subset, ccd_id not used below







#################################################
def do_merge(e_list):
    if type(e_list)==str: e_list= [e_list]
       
    ## ename -- for directory/naming
    ename=''
    for e in e_list: ename+=e
    
    os.chdir(cwd)    
    pdir_tmp= mktempdir(dir='.', ignore_cleanup_errors=True)
    pdir= pdir_tmp.name
    os.environ["PFILES"]="{};/usr/local/ciao/contrib/param:/usr/local/ciao/param:".format(pdir)
    ltab= '/data/legs/rpete/flight/acis_ecs/data/dmmerge_lookupTab.txt'
    merge_dir= '{}/evt2/{}/fpt_{}{}/'.format(basedir,ename,sfpt,xtra)
    if not os.path.exists(merge_dir): os.makedirs(merge_dir)
    os.chdir(merge_dir)
    print('\nOUTPUT dir: '+merge_dir)

    ciaover= getoutput(['ciaover']).splitlines()
    ciaov= ciaover[0].split()[1]
    calv= ciaover[2].split()[2]

    ## create merge list for each ccd + gunzip if needed
    for ccd in ccd_list:
        evt_ouf= '{}_{}TG.evt2'.format(ccd,tg)

        if not os.path.exists(evt_ouf):
            merge_file= 'merge_{}_{}TG.lis'.format(ccd,tg); fh= open(merge_file, 'w')
               
            for e in e_list:  ## if combining multi epochs
                edir='/data/hal9000/ecs/data/{}/{}'.format(e,rdir)
                edir=f'{basedir}/{e}/ciao{ciaov}_caldb{calv}'

                for fpt in fpt_list:  ## if combining multi fpt bins
                    if ccd=='s1_noCTI' or ccd=='s3_noCTI':
                        evt='{}/*{}_{}TGnoCTI_{}.evt2'.format(edir,ccd[0:2],tg,fpt)
                    else:
                        evt='{}/*{}_{}TG_{}.evt2'.format(edir,ccd,tg,fpt)                            
                    if glob(evt+'.gz'): os.system('gunzip '+gz)
                    if glob(evt):
                        tmp,= glob(evt)
                        fh.write(tmp+'\n')

            fh.close()
            f= open(merge_file,'r'); lines= f.readlines()
            if len(lines)==1:
                os.system('ln -s {} {}'.format(lines[0].rstrip(),evt_ouf))
            if len(lines)>1:
                dmmerge('@{}[cols -phas,-node_id][subspace -expno]'.format(merge_file), evt_ouf, lookupTab=ltab)
            #os.remove(merge_file)

        else: print('evt2 already exists '+merge_dir+evt_ouf)
    os.chdir(cwd)



    
    
#################################################
def do_spec(e_list):
    if type(e_list)==str: e_list= [e_list]

    ## ename -- for directory/naming
    ename=''
    for e in e_list: ename+=e

    os.chdir(cwd)    
    pdir_tmp= mktempdir(dir='.')
    pdir= pdir_tmp.name
    os.environ["PFILES"]="{};/usr/local/ciao/contrib/param:/usr/local/ciao/param:".format(pdir)

    merge_dir= '{}/evt2/{}/fpt_{}{}'.format(basedir,ename,sfpt,xtra)    
    spec_dir= '{}/spec/{}/fpt_{}_{}_{}TG{}'.format(basedir,ename,sfpt,sbin,tg,xtra)
    if not os.path.exists(spec_dir): os.makedirs(spec_dir)
    os.chdir(spec_dir)
    print('\nOUTPUT dir: '+spec_dir)

    ## dmextract
    for ccd in ccd_list:
        evt= merge_dir+'/{}_{}TG.evt2'.format(ccd,tg)
        tst_glob= glob(ccd+'*.pi')
        msg=''
        if not os.path.exists(evt): msg= '\nEVT2 NOT FOUND {}\n'.format(evt)
        if glob(ccd+'*.pi'): msg= '\nOUTPUT SPEC EXISTS FOR {} {} {} {}TG\n'.format(ccd,e,fpt,tg)
        if not msg=='':
            print(msg)
        else:
            print('extracting {}...'.format(ccd))
            for x in range(1,1025,xbin):
                xl=str(x).zfill(4); xh=str(x+xbin-1).zfill(4)

                for y in range(1,1025,ybin):
                    yl=str(y).zfill(4); yh=str(y+ybin-1).zfill(4)

                    ouf='{}_{}-{}x_{}-{}y.pi'.format(ccd,xl,xh,yl,yh)
                    if not os.path.exists(ouf):
                        print(spec_dir+'\t'+ouf,end='\r',flush=1)

                        ## min 200 events in region
                        inf='{}[chipx={}:{},chipy={}:{}][cols time]'.format(evt,xl,xh,yl,yh)
                        tmp= dmstat(inf)
                        tmp= int(dmstat.out_good)
                        if tmp > 300:
                            inf='{}[chipx={}:{},chipy={}:{}][bin pi=1:1024:1]'.format(evt,xl,xh,yl,yh)
                            dmextract(inf,ouf)
                        else: print('low count region, counts= {}, skipping...\n'.format(tmp))
            print('')
    os.chdir(cwd)

            


    
########################### MAIN ###############################            

for e_list in main_e_list: do_merge(e_list)
for e_list in main_e_list: do_spec(e_list)
sys.exit(0)




import multiprocessing
pool = multiprocessing.Pool(1) #specify #cpu's
dothis=1
if dothis==1:
    for e_list in main_e_list:
        pool.apply_async(do_spec, args=(e_list,))
    pool.close(); pool.join()
