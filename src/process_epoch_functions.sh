#! /bin/bash

SCRIPTDIR=$( cd -- "$( dirname -- "${BASH_SOURCE[0]}" )" &> /dev/null && pwd )

datadir=$("$SCRIPTDIR"/datadir)
lookuptab="$SCRIPTDIR"/../data/dmmerge_lookupTab.txt

# given year and month, return the epoch number
epoch() {
    [ $# -eq 2 ] || {
	\echo "Usage: $0 year month" 1>&2
	return 1
    }
    local year="$1"
    local month=$(( 10#$2 ))

    # first epoch started in 2000-02
    nmonths=$(( ($year-2000)*12 + $month-2 ))
    echo $(( $nmonths/3 + 1))
}

# given epoch number, return year and month of beginning
epoch_start() {
    [ $# -eq 1 ] || {
	\echo "Usage: $0 epoch" 1>&2
	return 1
    }
    local epoch=$(( 10#$1 ))
    nmonths=$(( ($epoch-1)*3 ))
    year=$(( 2000+$nmonths/12 ))
    month=$(( 2+$nmonths%12 ))
    echo $year $month
}

# comparison of two floating points, tidies things a bit later
gt() {
    (( $(echo "$1>$2" | bc) )) && return;
    return 1
}

# pbk files are links, sometimes those link targets must be unzipped
unlink_pbk() {
    [ $# -eq 1 ] || {
	\echo "Usage: $0 pbk" 1>&2
	return 1
    }

    local dir=$(dirname "$pbk")
    cd "$dir"
    local link=$(readlink $(basename "$pbk"))

    [ -z "$link" ] && { cd -; return; }

    [ ! -f "$link" -a -f "$link".gz ] && gunzip "$link".gz
    cd -

    [ -f "$dir/$link" ] || {
	echo "Failed to gunzip '$dir/$link.gz'" 1>&2
	return 1
    }
}

shopt -s extglob
process_obsid() {
    [ $# -eq 1 ] || {
	\echo "Usage: $0 indir" 1>&2
	return 1
    }
    local indir="$1"
    local outdir="$indir/repro"
    mkdir -p "$outdir"

    local evt1=$(echo "$indir"/secondary/*evt1.fits?(.gz))
    local mtl=$(\ls "$indir"/secondary/*mtl1.fits*)
    local msk=$(\ls "$indir"/secondary/*msk1.fits*)
    local stat=$(\ls "$indir"/secondary/*stat1.fits*)
    local flt=$(\ls "$indir"/secondary/*flt1.fits*)
    local pbk=$(\ls "$indir"/secondary/*pbk0.fits*)
    local bias_stack=$(\ls -w 0 -m $indir/$obsid/secondary/*bias* | sed 's/ //g')

    local obsid=$(dmkeypar "$evt1" obs_id ec+)

    local reset="$outdir/${obsid}.reset"
    local dstrk="$outdir/${obsid}.dstrk"
    local aglow="$outdir/${obsid}.aglow"
    local abb_bpix="$outdir/${obsid}.abb_bpix"
    local bpix="$outdir/${obsid}.bpix"
    local obspar="$outdir/obs${obsid}.par"

    if [[ "$evt1" =~ \.gz ]]; then
	gzip -dc "$evt1" > "$reset"
    else
	cp -a "$evt1" "$reset"
    fi

    punlearn ardlib

    acis_clear_status_bits "$reset"

    punlearn destreak
    destreak "$reset" "$dstrk" mask=NONE filter=no cl+
    dmmakepar "$dstrk" "$obspar" cl+

    punlearn acis_build_badpix
    acis_build_badpix \
	obsfile="$obspar" \
	pbkfile="$pbk" \
	biasfile="$bias_stack" \
	outfile="$abb_bpix" \
	bitflag=00000000000000120021100020022222 \
	calibfile=CALDB \
	mode=h \
	cl+

    punlearn acis_find_afterglow
    acis_find_afterglow \
	infile="$dstrk" \
	outfile="$aglow" \
	badpixfile="$abb_bpix" \
	maskfile="$msk" \
	statfile="$stat" \
	cl+

    punlearn acis_build_badpix
    acis_build_badpix \
	obsfile="$obspar" \
	pbkfile="$pbk" \
	biasfile=NONE \
	outfile="$bpix" \
	calibfile="$aglow" \
	procbias=no \
	mode=h \
	cl+

    tgain_opts=
    cti_opts=
    [ -n "$TGAIN" ] && tgain_opts .= " apply_tgain=$TGAIN"
    [ -n "$TGAINFILE" ] && tgain_opts .= " tgainfile=$TGAINFILE"
    [ -n "$CTI" ] && cti_opts .= " apply_cti=$CTI"
    [ -n "$CTIFILE" ] && cti_opts .= " ctifile=$CTIFILE"

    ape_infile="$dstrk"
    ape="$outdir/${obsid}.ape"

    #logtool(['acis_process_events', 'infile='+ape_inf, 'outfile='+ape, 'badpixfile='+bpix, 'acaofffile=NONE', 'mtlfile='+mtl, 'apply_cti='+appcti, 'apply_tgain='+tg, 'ctifile='+cti_file, 'tgainfile='+tg_file, 'check_vf_pha=no', 'stop=tdet', 'pix_adj=NONE', 'eventdef={d:time,l:expno,s:ccd_id,s:node_id,s:chip,s:tdet,d:phas,l:pha,l:pha_ro,f:energy,l:pi,s:fltgrade,s:grade,x:status}'], silent=0, log=1)

    punlearn acis_process_events
    acis_process_events \
	infile="$ape_infile" \
	outfile="$ape" \
	acaofffile=NONE \
	$ape_opt \
	$tgain_opts \
	$cti_opts \
	badpixfile="$bpix" \
	mtlfile="$mtl" \
	eventdef='{d:time,l:expno,s:ccd_id,s:node_id,s:chip,s:tdet,d:phas,l:pha,l:pha_ro,f:energy,l:pi,s:fltgrade,s:grade,x:status}' \
	stop=tdet \
	pix_adj=NONE \
	cl+

    local evt2="$outdir/${obsid}.evt2"
    punlearn dmcopy
    dmcopy "$ape"'[events][grade=0,2,3,4,6,status=0]' "$evt2" cl+
    dmcopy "$evt2[@$flt]" "$evt2" cl+

    kzero=273.15
    for t in $(seq 101 120); do
	fpt_lo=$(echo "-($t+0.19)+$kzero" | bc)
	fpt_hi=$(echo "${fpt_lo}+1" | bc)
	gti="$outdir/${obsid}_${t}.gti"
	# FIXME: cl+ gives
	#        ERROR: both parameters userlimit and lkupfile do not
	#        have a value specified.  At least one must have a value.
	rm -f "$gti"
	punlearn dmgti
	dmgti \
	    infile="$mtl" \
	    "userlimit=((fp_temp>=${fpt_lo})&&(fp_temp<${fpt_hi}))" \
	    outfile="$gti"
    done

    ccd_id=0
    for ccd in i{0,1,2,3} s{0,1,2,3,4,5}; do
	exp=$(dmkeypar "$dstrk[ccd_id=${ccd_id}][cols time]" exposure ec+)
	if gt $exp 10; then
	    exp=$(dmkeypar "$evt2[ccd_id=${ccd_id}][cols time]" exposure ec+)
	    if gt $exp 10; then
		for t in $(seq 101 120); do
		    gti="$outdir/${obsid}_${t}.gti"
		    dmlist "$gti" blocks | grep -q '3: GTI' || continue
		    exp=$(dmkeypar "$evt2[@${gti}][ccd_id=${ccd_id}][cols time]" exposure ec+)
		    if gt $exp 10; then
			evt_fpt="$outdir/${obsid}_${t}_${ccd}.evt2"
			dmcopy \
			    "$evt2[@${gti}][ccd_id=${ccd_id}]" \
			    "$evt_fpt" \
			    cl+
		    fi
		done
	    fi
	fi

	(( ++ccd_id ))
    done
}

merge_epoch() {
    [ $# -eq 2 ] || {
	\echo "Usage: $0 outdir epoch" 1>&2
	return 1
    }
    local outdir="$1"
    local e=$(printf %03d $((10#"$2")))

    mkdir -p "$outdir/merge_lis"
    cd "$outdir"

    for ccd in i{0,1,2,3} s{0,1,2,3,4,5}; do
	for t in $(seq 101 120); do
	    outf="e${e}_${ccd}_${t}.evt2"
	    globstr="$datadir/e${e}/[0-9][0-9][0-9][0-9][0-9]/repro/[0-9][0-9][0-9][0-9][0-9]_${t}_${ccd}.evt2"
	    inevt2=$(\ls $globstr 2>/dev/null || :)
	    [ -z "$inevt2" ] || {
		nfiles=$(wc -l <<<"$inevt2")
		mergef="merge_lis/merge_${ccd}_${t}.lis"
		rm -f "$mergef"
		for f in $inevt2; do echo "$f" >> "$mergef"; done
		[ $nfiles -eq 1 ] && {
		    cp -a "$inevt2" "$outf"
		} || {
		    punlearn dmmerge
		    dmmerge \
			"@${mergef}[events][subspace -expno]" \
			outfile="$outf" \
			lookupTab="$lookuptab" \
			mode=h \
			cl+
		}
	    }
	done
    done

    #mkdir -p merge_lis
    #mv merge_*.lis merge_lis

    cd -
}

cmp_files() {
    [ $# -eq 2 ] || {
	echo "Usage: $0 dir1 dir2" 1>&2
	return 1
    }
    local dir1="$1"
    local dir2="$2"
    for f1 in "$dir1"/*.evt2; do
	bf1=$(basename "$f1")
	f2="$dir2/$bf1"
	[ -f "$f2" ] || { echo "'$bf1' not found in '$dir2'" 1>&2; }
    done
}

cmp_datasum() {
    [ $# -eq 3 ] || {
	echo "Usage: $0 dir1 dir2 ext" 1>&2
	return 1
    }
    local dir1="$1"
    local dir2="$2"
    local ext="$3"
    for f1 in "$dir1"/*.$ext; do
	f2="$dir2"/$(basename "$f1")
	[ -f "$f2" ] || { echo "'$f1' not found in '$dir2'" 1>&2; continue; }
	ds1=$(dmkeypar "$f1" datasum ec+)
	ds2=$(dmkeypar "$f2" datasum ec+)
	[ $ds1 -eq $ds2 ] || echo $f1
    done
}

cmp_fits() {
    # fits dir is /data/legs/rpete/data/ECS/e${epoch}/fits/$ECSID/fits
    [ -z "$ECSID" ] && {
	\echo "\$ECSID is not set, exiting" 1>&2
	return 1
    }

    [ $# -eq 1 ] || {
	echo "Usage: $0 epoch" 1>&2
	return 1
    }
    local epoch=$(printf %03d "$1")
    for f in "$datadir/e$epoch/fits/$ECSID/fits/fpt_120-119-118_256x256y"/[is][0-5]_{ecs,bkg,pha}.txt; do
	bname=$(basename "$f")
	\diff -u \
	     $f \
	     "$datadir/fits/ciao4.17.0_caldb4.12.2/e$epoch/fpt_120-119-118_256x256y_yesTG/$bname"
	\diff -u \
	     "$datadir/fits/ciao4.17.0_caldb4.12.2/e$epoch/fpt_120-119-118_256x256y_yesTG/$bname" \
	     $f
    done
}

link_new_old() {
    base=/data/legs/rpete/data/ECS
    ciao=ciao4.17.0_caldb4.12.2
    epochs="$@"
    for epoch in $epochs; do
	epoch=e$(printf %03d $((10#$epoch)))
	for d1 in evt2 spec fits; do
	    for d2 in $base/$d1/$ciao/$epoch/fpt_120-119-118*; do
		bname=$(basename $d2)
		echo ln -fsT $base/$epoch/fits/$ciao/$d1/$bname $base/$d1/$ciao/$epoch/$bname
	    done
	done
    done
}

cmp_merges() {
    local epoch=e$(printf %03d $((10#$1)))

    cd /data/legs/rpete/data/ECS/$epoch
    for f1 in merge/test2/*.evt2; do
	b=$(basename $f1)
	f2=merge/test2.bak/$(perl -ple 's/_(\d{3})/_yesTG_$1/' <<< $b)
	ds1=$(dmkeypar $f1 datasum ec+)
	ds2=$(dmkeypar $f2 datasum ec+)
        [ $ds1 -eq $ds2 ] || echo $b;
    done
    cd -
}

cmp_repros() {
    local epoch=e$(printf %03d $((10#$1)))

    cd /data/legs/rpete/data/ECS/$epoch
    for o in [0-9]*; do
        for f1 in $o/repro/*.evt2; do
            b=$(basename $f1);
            f2=$o/repro.bak/$b;
            f2=$(sed s/\\/"$o"/\\/"$o"_yesTG/ <<<$f2);
            ds1=$(dmkeypar $f1 datasum ec+);
            ds2=$(dmkeypar $f2 datasum ec+);
            [ $ds1 -eq $ds2 ] || echo $b;
        done;
    done
    cd -
}

