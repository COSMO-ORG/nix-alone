#!/usr/bin/awk -f 
BEGIN {
	header=0;
	data=0;
	sep=" ";
}
{
	if(header && /fields/) {
		for(i=3; i<=NF; i++) {
			if($i=="timestamp") {timestamp=i-2};
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
		print "time T PS QV U V ASWDIR_S ASWDIFD_S ATHD_S TOT_PREC RAIN_GSP SNOW_GSP GRAU_GSP T_SO"
	}
	if(data) {
		print $timestamp sep $t2m sep $ps sep $qv sep $u sep $v sep $iswr_dir sep $iswr_dif sep $ilwr sep $psum sep "0" sep "0" sep "0" sep 273.15;
	}
	if(/SMET 1.1 ASCII/) {header=1; data=0};
	if(/\[DATA\]/) {header=0; data=1};
}
