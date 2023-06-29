
 * This do file: Prepares for demand estimation
 * Authors: Gabriel Kreindler

** RESTRICTED TO AM **

 /* Plan:
at trip level (?)
- departure time  			OK
- intersect with area 		OK

at ID level
- treatment status  		OK
-- DT treat 				OK
-- Area treat 				OK
-- Area control 			OK
- peak value 				OK	
- area detour 				OK
- typical KM 				OK
- normal mu and sigma for h_arrival

at Google Map level
- travel time at each HD, one way

  */

clear all
pause on
set more off

* Path based on global set in profile.do
	local path "$congestion_pricing_path"
	cd "`path'"
	adopath ++ "code/ado/"

* Income data
	use "data/coded_cto/recruitment coded.dta", clear
	keep uidp_original log_income log_income_with_pred income_outlier
	rename uidp_original uidp
	tempfile income
	save 	`income'

************************
*** TRIPS micro data ***
************************
*** Visit H/W/W2
	use "data/coded_gps_dta/coded_trips_15.dta", clear

	drop if inlist(study_cycle,9,99) // drop post-experiment
	assert trip == 1 - date_wo_trips
	keep if trip == 1

* checks
	assert sample_analysis == qual_gb
	assert date != mdy(4,26,2017)
	assert qual_good + qual_bad + qual_halfday + qual_no_data == 1
	isid uidp date chain
	sort uidp date chain

	*** SAMPLE drop 11 respondents with faulty AREA locations
	* tag as not part of the AREA treatment
		replace in_area_treatment = 0 if sample_area_drop == 1 
		replace area_treat = 0        if sample_area_drop == 1 

****************
**** CODING ****
****************
	*** add daily rate (120 vs 80 and 240 vs 160)
	preserve
		keep if area_treat == 1
		keep uidp study_cycle ad?
		duplicates drop
		isid uidp
		reshape long ad, i(uidp study_cycle) j(study_day)

		assert inlist(ad,80,120,160,240)

		bys uidp: egen min_charge = min(ad)
		bys uidp: egen max_charge = max(ad)

		assert (min_charge == 80 & max_charge==120) | (min_charge==160 & max_charge==240)

		gen high_rate_day = ad == max_charge
		rename ad rate_today

		sort uidp study_day
		by uidp (study_day): gen hrd_l1 = high_rate_day[_n-1]
		by uidp (study_day): gen hrd_l2 = high_rate_day[_n-2]
		recode hrd_l1 hrd_l2 (.=0)

		keeporder uidp study_cycle study_day rate_today 

		tempfile daily_charges
		save 	`daily_charges'
	restore

	merge m:1 uidp study_cycle study_day using `daily_charges'
	count if _m==2
	assert r(N)==101
	drop if _m==2
	* we have daily rate data for every day in treatment
	assert area_treat == (_m==3)
	drop _m

	* average area charge
	gen area_charge_mean = (ad1 + ad2 + ad3 + ad4 + ad5) / 5
	// recode area_charge_mean rate_today (.=0)

*** Add Google Maps data
	merge m:1 uidp using "data/google_maps/gmaps_hw_predicted.dta"
	assert _m!=2
	drop _m


************************
**** GENERAL VARS ****
************************
	
	gen peak_rate = 0
	replace peak_rate = 12 if dt_low  == 1
	replace peak_rate = 24 if dt_high == 1

	*** typical KM from home to work
	* define home to work trips
	gen hwp = (oh==1) & (dw==1 | dw2 == 1)
	cap drop temp
	bys uidp: gegen temp = mean(plb) if hwp == 1
	bys uidp: gegen mean_km = mean(temp)

	// bys uidp: gen o1 = _n==1
	reg mean_km dist_mode
	predict mean_km_predicted
	replace mean_km = mean_km_predicted if mean_km == .

********************
*** Final coding ***
********************

	keeporder uidp mean_km dist_mode dist_mean dt_wcha peak_rate pk_am pk_pm ///
			  adetourtime adetourdist in_area_treatment area_charge_mean dur_gm_*
	duplicates drop
	gisid uidp

***************************
*** Merge AREA outcomes ***
***************************
	merge 1:1 uidp using "data/coded_model/area_het.dta"
	assert _m!=2
	gen has_area_data = _m==3
	drop _m


***************************
*** Merge DT   outcomes ***
***************************
	merge 1:1 uidp using "data/coded_model/dt_het.dta"
	assert _m!=2
	gen has_dt_data = _m==3
	drop _m

	tab has_area_data has_dt_data, m
	// 	has_area_d |      has_dt_data
	//        ata |         0          1 |     Total
	// -----------+----------------------+----------
	//          0 |       193         99 |       292 
	//          1 |        49        156 |       205 
	// -----------+----------------------+----------
	//      Total |       242        255 |       497 


	drop if has_area_data == 0 & has_dt_data == 0
	recode sample_dt_pre sample_dt_pos (.=0)
	recode idfxbe idfx idfxnpre idfxnpos (.=0)

	* only missing norm_mean and norm_sd for missing DT data
	// assert (sample_dt_pre == 0 & norm_mean == .) | (sample_dt_pre == 1 & norm_mean != .)
	assert norm_mean != . if sample_dt_pre == 1

	* for convenience (to avoid errors in matlab, replace with mean = 0 and sd = 1)
	replace norm_mean = 0 if sample_dt_pre == 0
	replace norm_sd   = 1 if sample_dt_pre == 0

	*** Merge in income data
	merge 1:1 uidp using `income' // /* log_income income_outlier */ log_income_with_pred
	assert _m!=1
	drop if _m==2
	drop _m

	* mean-1 inverse income
	gen inverse_income_with_pred = 1 / exp(log_income_with_pred)
	sum inverse_income_with_pred
	replace inverse_income_with_pred = inverse_income_with_pred / r(mean)
	sum inverse_income_with_pred
	sum inverse_income_with_pred if peak_rate != 0
	sum inverse_income_with_pred if in_area_treatment == 1

*** SORTING ***
	sort in_area_treatment uidp, stable 

* encode uidp
	// sencode uidp, generate(uidpn)
	// label drop uidpn
	// gen n=_n
	// sort uidpn
	// assert n==_n
	// drop n
	// drop uidp

*************
*** Order ***
*************
	keeporder uidp mean_km dist_mode dist_mean ///
			  norm_mean norm_sd ///
			  dt_wcha peak_rate pk_am pk_pm ///
			  adetourtime adetourdist in_area_treatment area_charge_mean ///
			  nct23 nct ntr cha_ct23 cha_ct cha_tr ///
			  sample_dt_pre sample_dt_pos mean_dt_pre* mean_dt_pos* ///
			  idfxnpre idfxnpos idfxbe idfx  ///
			  dur_gm_* n_pre n_pos inverse_income_with_pred

	mdesc
	count if ntr > 0 & nct23 >0 & ~missing(ntr) &~missing(nct23) // v1 (186), v2 (154)

**************
*** Export ***
**************
	assert inlist(in_area_treatment,0,1)
	export delimited using "data/coded_model/coded_uidpn_a.csv" if in_area_treatment == 1, replace
	export delimited using "data/coded_model/coded_uidpn_na.csv" if in_area_treatment == 0, replace


