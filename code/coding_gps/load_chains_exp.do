
 * This do file: Loads and codes trips data 
 * This version: 25 July, 2017
 * Authors: Gabriel Kreindler

clear all
pause off
set more off

if "`1'" == ""{
	local chain_th "60"	
}
else{
	local chain_th "`1'"		
}

* Path based on global set in profile.do
	local path "$congestion_pricing_path"
	cd "`path'"

**************
*** Load trips
	import delimited using "data/coded_gps/sample_exp_chains_`chain_th' no PII.csv", clear

	* checks
	mdesc
	isid uidp date chain
	sort uidp date chain
	// assert chain==0 if center_orig_lat == . 
	// by uidp date: assert chain==_N if center_dest_lat == . 

	*** date
 	gen date_ = date(date,"YMD")
 	drop date
 	rename date_ date
 	format %td date
 	gen dow = dow(date) // 0 = Sunday ... 6 = Saturday

 *** Special drops
	* DROP all data before 04-03 (inclusive) for L0314124709
	drop if uidp == "L0314124709" & date <= mdy(4,3,2017)

	* DROP all data before 04-06 (inclusive) for L0221092234
	drop if uidp == "L0221092234" & date <= mdy(4,3,2017)

	* DROP last day (Sept 11) due to lack of quality measures for that day
	assert date < mdy(9,12,2017)
	drop if date == mdy(9,11,2017)

* CODING

*** HWW2 definitions 
	tab oh, m
	count if oh + ow > 1
	assert r(N) <= 1

	// gen trip_dir_hw = oh + dw == 1
	// gen trip_dir_wh = ow + dh == 1

	gen trip_hw = (oh == 1 & dw == 1) | (ow == 1 & dh == 1)
	gen trip_hw2 = (oh == 1 & dw2 == 1) | (ow2 == 1 & dh == 1)
	gen trip_ww2 = (ow == 1  | ow2 == 1) & (dw == 1  | dw2 == 1)  // covers ALL combinations of w-w, w-w2, w2-w and w2-w2

	* home - work or home - work2
	gen trip_hwp = trip_hw == 1 | trip_hw2 == 1

	* home - work or home - work2 or any work*-work*
	gen trip_hwpwp = trip_hw == 1 | trip_hw2 == 1 | trip_ww2 == 1

*** generate start and end times in fractional hours
	gen t_start = t_start_bg_ / 3600 
	gen t_stop = t_stop_bg_ / 3600 

*** Define in-bangalore (opposite of outstation)
	assert !missing(inb_frac)
	gen sample_inb_70p = inb_frac >= 0.7
	gen sample_inb_60p = inb_frac >= 0.6

	tab sample_inb_70p, m
	tab sample_inb_60p, m

*** Departure time
	gen sample_dt722 = inrange(t_start, 7, 22)

*** Departure time precision
	// sum t_start_prec t_start_prec_max t_stop_prec t_stop_prec_max ,d
	assert t_start_prec_max >= t_start_prec
	assert t_stop_prec_max >= t_stop_prec
	assert !missing(t_start_prec_max)
	assert !missing(t_stop_prec_max)
	gen sample_dep_precise = t_start_prec_max <= 15
	gen sample_arr_precise = t_stop_prec_max <= 15

	tab sample_dep* sample_arr*, m

*** Very short/ very long trips
	sum plb, d
	sum dur_trips, d
	assert plb!=. & dur_trips!=.

	gen sample_dur_long = dur_trips > 180
	gen sample_dur_short = dur_trips < 5
	
	gen sample_dist_long = plb > 35000 // 99th percentile is actually 35223
	gen sample_dist_short = plb < 500

*** Jump distance and duration
	sum trips_dur_jump if trips_dur_jump > 0, d
	sum trips_plb_jump if trips_plb_jump > 0, d

	gen  jump_frac_dur = trips_dur_jump / dur_trips
	egen jump_frac_dur_cat = cut(jump_frac_dur), at(0,0.3,5)
	bys jump_frac_dur_cat: sum jump_frac_dur

	gen  jump_frac_plb = trips_plb_jump / plb
	egen jump_frac_plb_cat = cut(jump_frac_plb), at(0,0.3,5)
	bys jump_frac_plb_cat: sum jump_frac_plb

	// count if trips_dur_jump > 0
	// count if trips_plb_jump > 0
	// count if trips_dur_jump > 0 & trips_plb_jump == 0  // 3811 examples!
	// count if trips_dur_jump ==0 & trips_plb_jump > 0

	cap drop jump_dur_cat
	egen jump_dur_cat = cut(trips_dur_jump), at(0,1,30,2000) // cut at less than half
	assert jump_dur_cat!=.
	bys jump_dur_cat: sum trips_dur_jump

	cap drop jump_plb_cat
	egen jump_plb_cat = cut(trips_plb_jump), at(0,1,2000,100000) // cut at less than half
	assert jump_plb_cat!=.
	bys jump_plb_cat: sum trips_plb_jump	

	// tab jump_dur_cat jump_plb_cat, m

	gen sample_jumpy_1 = trips_dur_jump >= 30 | trips_plb_jump >= 2000

	* more conservative - removes all trips with jumpps > 30% of duration OR > 30% of distance (plb)
	gen sample_jumpy_2 = trips_dur_jump >= 30 | trips_plb_jump >= 2000 | abs(jump_frac_dur_cat - .3) <0.01 | abs(jump_frac_plb_cat - .3) <0.01

	tab sample_jumpy_*, m

*** Swiggly factor - the diameter of the trip is no more than 30% of the path length. 
	gen diameter_ratio = chain_diameter / plb
	gen sample_swiggly = diameter_ratio < 0.3 & !missing(diameter_ratio)

*** Loopy trips factor - very small origin-dest distance relative to total path, and the entire trip is less than 2km long
	gen loopy_ratio = center_line_dist / plb
	gen sample_loopy = loopy_ratio < 0.3 & !missing(loopy_ratio) & plb < 2000

	tab sample_loopy, m

************************
*** SAMPLE SELECTION ***
************************
	// label define sample_drop ///
	// 		01 "01 weekend" ///
	// 		02 "02 outstation" 	 ///
	// 		03 "03 dt night" ///
	// 		04 "04 dep imprec" ///
	// 		05 "05 arr imprec" ///
	// 		06 "06 dur_long" ///
	// 		07 "07 dist_long" ///
	// 		08 "08 dur_short" ///
	// 		09 "09 dist_short" ///
	// 		10 "10 swiggly" 

	gen 	sample_drop = ""
	replace sample_drop = "01 weekend" 		if  inlist(dow,0,6)
	replace sample_drop = "02 outstation" 	if  sample_inb_70p 	== 0 & sample_drop == ""
	replace sample_drop = "03 dt night" 	if  sample_dt722 	== 0 & sample_drop == ""
	
	replace sample_drop = "04 dep imprec" 	if  sample_dep_precise == 0 & sample_drop == ""
	replace sample_drop = "05 arr imprec" 	if  sample_arr_precise == 0 & sample_drop == ""

	replace sample_drop = "06 dur_long" 	if  sample_dur_long == 1 	& sample_drop == ""
	replace sample_drop = "07 dist_long" 	if  sample_dist_long == 1 	& sample_drop == ""

	replace sample_drop = "08 dur_short" 	if  sample_dur_short == 1 	& sample_drop == ""
	replace sample_drop = "09 dist_short" 	if  sample_dist_short == 1 	& sample_drop == ""

	replace sample_drop = "10 jump"			if sample_jumpy_2 == 1 		& sample_drop == ""

	replace sample_drop = "11 swiggly" 		if 	sample_swiggly == 1 	& sample_drop == ""

	replace sample_drop = "12 loopy short" 	if 	sample_loopy == 1 		& sample_drop == ""

	replace sample_drop = "13 KEEP" if sample_drop == ""

	encode sample_drop, gen(sample_drop_encoded)
	la li sample_drop_encoded // check

*** Stats
	tab sample_drop_encoded
	tab sample_drop_encoded if sample_drop_encoded > 3
	// 75% trips keep
	// 12% departure time imprecise
	// 3% arrival time imprecise
	// 7.5% short duration (<5min)

	gen sample_trip_ok = sample_drop == "13 KEEP"

* save raw chains
	compress
	save "data/coded_gps_dta/chains_exp_`chain_th'.dta", replace



