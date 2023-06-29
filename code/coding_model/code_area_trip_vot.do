
 * This do file: coding for demand estimation
 * This version: Nov 23, 2021
 * Authors: Gabriel Kreindler

clear all
pause off
set more off
set matsize 10000

* Path based on global set in profile.do
	local path "$congestion_pricing_path"
	cd "`path'"

*** LOAD 
	use "data/coded_gps_dta/coded_trips_15.dta", clear

****************
**** SAMPLE ****
****************

*** Drop dates after the study ends (=9) and during the study but not in the experiment (=99)
	assert inlist(study_cycle,0,1,2,3,4,5,9,99)
	drop if inlist(study_cycle,9,99) // drop post-experiment
	// assert sample_study_proper == inlist(study_cycle,0,1,2,3,4,5)

*** SAMPLE keep only AREA respondents
	drop if in_area_treatment == 0

	*** Merge area candidates 
	merge m:1 uidp using "data\coded_gps_dta\area_candidates.dta"
	assert _m==3
	drop _m

	* sample where both >=1 short and >=1 long candidates exist
	gen sample_short_long_subtr = ac_valid_short==1 & ac_valid_long==1

*** SAMPLE drop 11 respondents with faulty AREA locations
	// gunique uidp if sample_area_drop == 1 
	drop if sample_area_drop == 1 

	* Drop DAYS without trips
	assert trip == 1 - date_wo_trips
	keep if trip == 1

*** Drop April 26, which was a holiday and wrongly included in participants' timelines
	assert date != mdy(4,26,2017)

*
	isid uidp date chain
	sort uidp date chain

****************
*** SAMPLING ***
****************

* trip quality sample 
	local trip_qual_sample "ok" //"all" 
	if "`trip_qual_sample'" == "ok"{
		keep if sample_trip_ok==1
	}	

* Sample of weeks during the experiment
	local sample_weeks "_wall" // "_w1"  "_w14" "_w4"
	local ifsample_weeks "if sample`sample_weeks'==1 "

* day data quality sample (GOOD+BAD)
	assert sample_analysis == qual_gb
	local quality_sample "& sample_analysis==1 "

* Regular commuters only	
	keep if inlist(regular_commuter,2,3)

* Restrict study cycle, and day quality
	keep `ifsample_weeks' `quality_sample'

*** Trip direction coding
		* ONLY h-w and w-h chains
	gen trip2_hw = (oh == 1 & dw == 1) 
	gen trip2_wh = (ow == 1 & dh == 1)
	assert trip2_hw + trip2_wh == trip_hw


*** SAMPLE - only H - W and back trips
	keep if trip_hw == 1
	assert trip2_hw + trip2_wh == 1

	gegen route_FE = group(uidp trip2_hw)


*** Interactions for persistence
	assert area_treat == ((study_cycle == 1) & a_early == 1) | (((study_cycle == 4) & a_late == 1))

	gen area_treat_w1 = area_treat * (study_cycle == 1)
	gen area_treat_w4 = area_treat * (study_cycle == 4)
	assert area_treat_w1 + area_treat_w4 == area_treat

	gen early_w2 = a_early * (study_cycle == 2)
	gen early_w3 = a_early * (study_cycle == 3)
	gen early_w4 = a_early * (study_cycle == 4)
	
	gen area_treat_early_w23 = a_early * inlist(study_cycle,2,3)


***
	assert inlist(c_area_charge,0,100)

	gen is_long_route = (c_area_charge == 0)
	gen is_long_route_pre = (c_area_charge == 0) * (study_cycle == 0)

	bys uid: gegen has_long_route_pre	= max(is_long_route_pre)
	bys uid: gen u1=_n==1


*** Treatment interactions with detour length
	sum adetourtime if u1==1, d
	assert !missing(adetourtime)
	gen high_dttime = adetourtime > r(p50)   // 5.8 minutes
	gen area_treat_high_dttime = area_treat * high_dttime
	gen early_w2_high_dttime = early_w2 * high_dttime
	gen early_w3_high_dttime = early_w3 * high_dttime
	gen early_w4_high_dttime = early_w4 * high_dttime


	gen area_treat_adetourtime = area_treat * (adetourtime - r(mean)) / r(sd)



*** 
	gen area_treat_long_pre = area_treat * has_long_route_pre
	gen early_w2_long_pre = early_w2 * has_long_route_pre
	gen early_w3_long_pre = early_w3 * has_long_route_pre
	gen early_w4_long_pre = early_w4 * has_long_route_pre



*** Export for Julia GMM analysis
	gen acharge_weekmean = (ad1 + ad2 + ad3 + ad4 + ad5)/5

	keep uidp date route_FE study_cycle is_long_route a_early a_late ///
				has_long_route_pre high_dttime adetourtime acharge_weekmean

	label drop study_cycle
	format %tdCCYY-MM-DD date

	by uidp: gegen max_cycle=max(study_cycle)
	tab max_cycle
	gunique uidp if max_cycle == 0

	drop if max_cycle == 0

	gunique uidp

	export delimited using "data/coded_model/area_trip_vot.csv", replace

