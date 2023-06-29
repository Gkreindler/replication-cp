
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


*** STUDY DATES - uidp-date level, all dates (including outstation, special, trial)
	use "data/coded_gps_dta/study_dates_full.dta", clear
	keep uidp date outstation charge bal_current cycle_type study_cycle study_day time_in_exp
	isid uidp date
	rename charge report_charge
	rename bal_current report_bal_current
	rename outstation report_outstation
	tempfile study_dates
	save 	`study_dates'

*** Study sample - uidp level, experiment start and end dates
	* merge based on uidp (so it directly covers all observations)
	use "data/coded_gps_dta/study_dates_full.dta", clear
	keep uidp date_exp_m* 
	duplicates drop
	isid uidp
	tempfile study_sample
	save 	`study_sample'

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
	keep if inrange(dow,1,5)

	count if inlist(data_quality, "good", "bad") // 37,489
	count if inlist(data_quality, "good", "bad") & outstation == 0 // 35,223

	rename outstation outstation_qual

	keep uidp date dow data_quality gap_hr tot_gap_hr jump_dist outstation_qual

*** merge info on dates in the experiment period that are not in the experiment (trial, outstation, special days off)
	merge m:1 uidp date using  "data/coded_gps_dta/study_dates_full.dta", keepusing(date_type)
	assert _m!=2 if date < mdy(9,11,2017)
	drop _m

* after we stopped the study
	replace data_quality = "no data"		if date >= mdy(9,11,2017)
	replace gap_hr = -1						if date >= mdy(9,11,2017)
	replace tot_gap_hr = -1					if date >= mdy(9,11,2017)
	replace jump_dist = -1				if date >= mdy(9,11,2017)
	replace outstation_qual = 0				if date >= mdy(9,11,2017)
	replace dow = dow(date) 				if date >= mdy(9,11,2017)
	// mdesc

	* only UIDPs from the study sample
	merge m:1 uidp using `study_sample'
	assert _m!=2
	keep if _m==3
	drop _m

	** Special drops
	* DROP all data before 04-03 (inclusive) for L0314124709
	// replace data_quality = "no data" if uidp == "L0314124709" & date <= mdy(4,3,2017)
	drop if uidp == "L0314124709" & date <= mdy(4,3,2017)

	* DROP all data before 04-06 (inclusive) for L0221092234
	// replace data_quality = "no data" if uidp == "L0221092234" & date <= mdy(4,3,2017)
	drop if uidp == "L0221092234" & date <= mdy(4,3,2017)

*** SAMPLE
* Drop various data quality variables
	keep uidp date dow data_quality gap_hr tot_gap_hr jump_dist outstation_qual date_exp_min date_exp_max date_type

*** CODING
	* assume we don't have data on April 26
	replace data_quality = "no data" if date==mdy(4,26,2017)
	foreach var of varlist gap_hr jump_dist tot_gap_hr{

		* covers 4/26 as well as the two H changes above
		count if `var' != -1 & data_quality == "no data"
		assert r(N) == 257

		replace `var' = -1 if data_quality == "no data" 
	}

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

*** ADD STUDY DATE INFO (this includes the timeline during the study (experiment+trial+outstation+special dates))
	* merge study dates
	merge 1:1 uidp date using `study_dates'
	assert _m!=2 
	drop _m

	* define time relative to experiment start
	isid uidp date
	sort uidp date
	assert date_exp_min!=.
	gen calendar_time = date - date_exp_min + 1

	* study_cycle is : 
	*  = 0 for pre meeting
	*  = 1-5 during experiment
	*  = 9 after the experiment
	replace study_cycle = 0 if calendar_time <= 0

	replace study_cycle = 9 if date > date_exp_max

	* 3 exceptions (why? Aug 8 GK)
	replace study_cycle = 99 if inlist(date_type,"outstation", "special", "trial") | ///
			(date == mdy(4,26,2017) & inlist(uidp,"L0208195435", "L0307094926", "L0308170807"))

	label define study_cycle 0 "before exp" 1 "week 1 (exp)" ///
											2 "week 2 (exp)" ///
											3 "week 3 (exp)" ///
											4 "week 4 (exp)" ///
											5 "week 5 (exp)" 9 "after exp" 99 "trial, out, special"
	label values study_cycle study_cycle
	tab study_cycle, m
	assert study_cycle != .

*** Define "notyet" = 1 if no data at all earlier dates
	sort uidp date
	gen notyet = 0
	by uidp (date): replace notyet = qual_no_data * (date != mdy(4,26,2017)) if _n==1
	by uidp (date): replace notyet = notyet[_n-1]*qual_no_data if _n > 1
	tab notyet, m
	assert qual_no_data == 1 if notyet == 1


*** Define "dropped" = 1 if no data at all subsequent dates
	gen date_reversed = -date
	sort uidp date_reversed
	gen dropped = 0
	by uidp (date_reversed): replace dropped = qual_no_data if _n==1
	by uidp (date_reversed): replace dropped = dropped[_n-1]*qual_no_data if _n > 1
	drop date_reversed

	sort uidp date 

	tab dropped if inlist(study_cycle,1,2,3,4,5)

*** Define no-data spells
	gen nodata_spell = 0
	by uidp (date):	replace nodata_spell = qual_no_data if _n == 1
	by uidp (date): replace nodata_spell = nodata_spell[_n - 1] * qual_no_data + qual_no_data if _n > 1

*** Code other quality stats
	sum tot_gap_hr , d
	sum tot_gap_hr if inrange(study_cycle,1,5), d

	*** TESTING
	// tab date if uidp =="L0525103515_DA90"
	// list date study_cycle dropped if uidp =="L0525103515_DA90"
	// sdfdsf

*** SAMPLE
	* drop dates after the experiment if the respondents dropped
	// drop if study_cycle == 9 & dropped == 1
	drop if notyet == 1 

	* check that all days are present
	by uidp: egen date_max = max(date)
	by uidp: egen date_min = min(date)
	by uidp: egen n_dates = count(date)
	weekdays_bw date_min date_max, generate(n_weekdays_)
	assert n_weekdays_ == n_dates
	drop date_max date_min n_dates n_weekdays_
	* done

*** Final coding
	by uidp: egen max_nodata_spell = max(nodata_spell)
	// hist max_nodata_spell if date_exp_min==date, d
	by uidp: egen n_days_after = sum(study_cycle == 9) // number of days after
	by uidp: egen n_gb_after = sum((study_cycle == 9) * qual_gb) // number of days after

*** SAVE
	compress
	save "data/coded_gps_dta/quality_exp.dta", replace
