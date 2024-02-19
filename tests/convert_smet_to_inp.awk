#!/usr/bin/awk -f 
BEGIN {
	header=0;
	data=0;
}
{
	if(header && /fields/) {
		for(i=3; i<=NF; i++) {
			if($i=="T_2M") {t2m=i-2};
			if($i=="PS") {ps=i-2};
			if($i=="QV") {qv=i-2};
			if($i=="U_10M") {u=i-2};
			if($i=="V_10M") {v=i-2};
			if($i=="ASWDIR_S") {iswr_dir=i-2};
			if($i=="ASWDIFD_S") {iswr_dif=i-2};
			if($i=="ATHD_S") {ilwr=i-2};
			if($i=="TOT_PREC") {psum=i-2}
		}
	}
	if(data) {
		# HACK: We need to only print every second line because of a timestep doubling in MeteoIO output
		d++
		if(d%2==1) print $t2m "," $ps "," $qv "," sqrt($u*$u+$v*$v) "," $iswr_dir "," $iswr_dif "," $ilwr "," $psum "," 273.15;
	}
	if(/SMET 1.1 ASCII/) {header=1; data=0};
	if(/\[DATA\]/) {header=0; data=1};
}
