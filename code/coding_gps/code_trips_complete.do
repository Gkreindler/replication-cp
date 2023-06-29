
 * This do file: Codes chain trips at daily level (v0 - no HW information, no charges)
 * This version: 30 July, 2017
 * Authors: Gabriel Kreindler

clear all
pause on
set more off

if "`1'" == ""{
	local chain_th "15"	
}
else{
	local chain_th "`1'"		
}

* Path based on global set in profile.do
	local path "$congestion_pricing_path"
	cd "`path'"


******************
*** LOAD TRIPS ***
******************

	use "data/coded_gps_dta/chains_complete_`chain_th'.dta", clear

************************
******* SAMPLE *********
************************
	* drop chains from outstation (based on trip data), weekends and nights (10pm-5am)
	// keep if sample_drop_encoded >= 3
	keep if !inlist(sample_drop_encoded,2,3)

************************

* merge with data quality to check number of days (good+bad) with missing trips (or with outstation)
	merge m:1 uidp date using "data/coded_gps_dta/quality_complete.dta" , ///
	keepusing(qual_* data_quality outstation_qual jump_dist tot_gap_hr)
	count if _m==1
	assert r(N) <= 6
	drop if _m==1
	
	* if _m==2, either no data on that day (70%), or good/bad data but no trips
	gen date_wo_trips = _m==2 
	drop _m

	* NOTE: we have DROPPED outstation chains, so the outstation_qual = 1 (N=320) that are merged are discrepancies between chains and qual/
	tab date_wo_trips data_quality, m
	tab date_wo_trips data_quality if sample_trip_ok==1, m
	gunique uidp date if date_wo_trips == 0 & qual_gb == 1 & outstation_qual == 0
	* --> 
	count 			 if date_wo_trips == 1 & qual_gb == 1 & outstation_qual == 0
	* --> 
	* ratio = 597 / (597 + 21154) = .02744701  <- includes night trips
	* ratio = 832 / (832 + 20919) = .03825111  < - excludes night trips (10pm-7am)
	* small number of - stay at home days

	* no days with trips and no_data, EXCEPT April 26
	count if qual_no_data == 1 & chain !=.
	assert r(N) == 0

*** CODING
	* fix chain ID on days without trips
	assert chain == . if date_wo_trips == 1
	replace chain = 1 if date_wo_trips == 1

	* INCLUDE days without trips in the sample of good trips
	assert  sample_trip_ok == . if date_wo_trips == 1
	replace sample_trip_ok = 1 if date_wo_trips == 1
	assert  sample_trip_ok != .

*** Clean Vars
	replace plb = plb / 1000

	foreach va of varlist plb dur_trips dur_bg_mm{
		rename `va' `va'_old
		winsor `va'_old, gen(`va') highonly p(0.01)
	}

*** Day-level coding
	isid uidp date chain

	* minor cleaning
	replace dow = dow(date) if dow==.
 	order dow, after(date)


	** Replace KM, duration, DT charges = 0 on days with >= halfday data but no trips (also for all locations)
	foreach va of varlist plb dur_bg_mm dur_trips{
		assert `va' == . if date_wo_trips == 1
		replace `va' = 0 if date_wo_trips == 1
		assert `va' != .
	}

	recode trip trips_dur_jump trips_plb_jump (.=0)

	* On days without trips, any_trip = 0. =1 otherwise
	gen any_trip = 1 - date_wo_trips
	assert inlist(any_trip,0,1)


*** Check missing values (see notes below)
	mdesc

	/* Notes:
	- generally missing 12,178 obs for days without trips
	- try to add d_oh and others on days with location (but no trip) data using location location
	- same for inb_frac
	- segs/sids?

	- INFO: all sample_ have missing for date_wo_trips==1. (except sample_trip_ok)
	- INFO: DT ramp has missing values
	- INFO: Area charges missing except intersect and charge, 

	- INFO: quality: jump_dist not defined for no_data
	 */

	 order uidp date chain trip date_wo_trips 

********************
** Save Trip Data **
********************
	isid date uidp chain
	save "data/coded_gps_dta/coded_trips_complete_`chain_th'.dta", replace
