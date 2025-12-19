ACIS ECS Fitting
========

Downloading of new epochs currently requires that Nick has listed the
appropriate obsids in `/data/hal9000/ecs/data/e${epoch}/support/obsids.dat`,
and then, e.g.,
```
# get data for several epochs
epochs='94 95 96'
src/download_epochs $epochs
```

Downloaded files end up in `$datadir/e{epoch:03d}`, where `$datadir` is set in
`src/process_epoch_functions.sh`

To reprocess obsids in each epoch with current CalDB,
```
src/process_epochs $epochs
```
`src/process_epochs` looks at environment variables `TGAIN`, `TGAINFILE`, `CTI`, `CTIFILE` and applies the appropriate `apply_{tgain,cti}` and `{tgain,cti}file` parameters to `acis_process_events`.

To merge the various temperatures and CCDs from each ObsID in the epochs,
```
ciao
export ECSID=$(src/ciaostr)
src/merge_epochs $epochs
```

The output files end up in `$datadir/e{epoch:03d}/merge/$ECSID`

After this step, the `repro` files for each ObsID are unneeded, and occupy a
large amount of disk space, so they can usually be deleted.

Further, compressing the merge `evt2` files can reduce their size by
70% or so.

To extract spectra for a specific temperature range and chip binning, e.g.,
```
tstrs='120,119,118 117,116 115,114 113,112 111,110 109,108 107,106'
binx=256
biny=256
time for tstr in $tstrs; do
    src/extract_spectra_parallel $tstr $binx $biny $epochs
done
```
and those spectra end up in `$datadir/e{epoch:03d}/fits/$ECSID/spec/` subdirectories.

To fit those spectra,
```
time for tstr in $tstrs; do
    src/fit_epochs $tstr $binx $biny $epochs
done
```

Finally, to plot fitted line energy deviations (currently 256x256y),
```
ciao
python3 src/plt_spec.py $tstr Al|Mn $epochs -p out.pdf
```
