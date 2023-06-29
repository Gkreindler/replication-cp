
 * This do file: Loads and codes quality data 

 * This version: 29 April, 2017
 * Authors: Gabriel Kreindler

clear all
pause off
set more off

* Path based on global set in profile.do

	local path "$congestion_pricing_path"
	cd "`path'"
	adopath ++ "code/ado/"


*** Load QUALITY data - uidp-date level
	// import delimited using "data/treatment/quality_all_new.csv", clear varnames(1)
	import delimited using "data/coded_gps/quality_all.csv", clear varnames(1)

	* date
	gen date_ = date(date,"YMD")
	drop date 
	rename date_ date
	format %td date
	gen dow = dow(date)

	* SAMPLE - weekdays only
	// keep if inrange(dow,1,5)

	count if inlist(data_quality, "good", "bad") // 37,489
	count if inlist(data_quality, "good", "bad") & outstation == 0 // 35,223

	rename outstation outstation_qual

*** SAMPLE
	keep uidp date dow data_quality gap_hr tot_gap_hr jump_dist outstation_qual   

*** CODING

	* at the individual level
	gen qual_good 		= data_quality == "good"
	gen qual_bad 		= data_quality == "bad"
	gen qual_gb 		= qual_good + qual_bad
	gen qual_halfday 	= data_quality == "halfday"
	gen qual_no_data 	= data_quality == "no data"
	assert qual_good + qual_bad + qual_halfday + qual_no_data == 1

	* missing gap means full missing (14 hours)
	recode tot_gap_hr (-1=14)
	recode jump_dist (-1=.)
	sum jump_dist, d 
	replace jump_dist = 15 if jump_dist > 15 & !missing(jump_dist)


	* define time relative to experiment start
	isid uidp date
	sort uidp date


*** Define no-data spells
	gen nodata_spell = 0
	by uidp (date):	replace nodata_spell = qual_no_data if _n == 1
	by uidp (date): replace nodata_spell = nodata_spell[_n - 1] * qual_no_data + qual_no_data if _n > 1


	* check that all days are present
	by uidp: egen date_max = max(date)
	by uidp: egen date_min = min(date)
	by uidp: egen n_dates = count(date)
	weekdays_bw date_min date_max, generate(n_weekdays_)
	// assert n_weekdays_ == n_dates
	drop date_max date_min n_dates n_weekdays_
	* done

*** Final coding
	by uidp: egen max_nodata_spell = max(nodata_spell)

*** SAVE
	compress
	save "data/coded_gps_dta/quality_complete.dta", replace
