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
`src/process_epoch_functions.sh'

To reprocess obsids in each epoch with current CALDB,
```
src/process_epochs $epochs
```
`src/process_epochs` looks at environment variables `TGAIN`, `TGAINFILE`, `CTI`, `CTIFILE` and applies
options `apply_(tgain|cti)` and `(tgain|cti)file` appropriately to `acis_process_events`.

To merge the various temperatures and CCDs from each ObsID in the
epochs
```
ciao
export ECSID=$(src/ciaostr)
src/merge_epochs $epochs
```

The output files in up in `$datadir/e{epoch:03d}/merge/$ECSID`

After this the `repro` files for each ObsID are unneeded, and require a
good deal of disk space, so they can usually be deleted.

Further, comressing the merge `evt2` files can reduce their footprint
by 70% or so.

To extract spectra for a specific temperature range and chip binning, e.g.,
```
tstr=120,119,118
binx=256
biny=256
src/extract_spectra_parallel $tstr $binx $biny $epochs
```
and those spectra end up in `$datadir/e{epoch:03d}/fits/$ECSID/spec/` subdirectories.

To fit those spectra,
```
src/fit_epochs $tstr $binx $biny $epochs
```

Finally, to plot fitted line energy deviations (currently 256x256y),
```
ciao
python3 src/plt_spec.py $tstr Al|Mn $epochs -p out.pdf
```
