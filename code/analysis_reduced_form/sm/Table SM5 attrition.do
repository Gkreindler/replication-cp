
 * This do file: Attrition Results (Attrition = Data Quality)
 * This version: August 15, 2018
 * Authors: Gabriel Kreindler

 /* T2 - DT responses at day level (full sample)
   - Columns: all/AM/PM, DD/FE. (6 columns)
   - panel A. Attrition (days with useable data)
   - panel B. Number of Trips
   - panel C. Charges */

clear all
pause off
set more off
set matsize 10000

*** TABLE number
local tfolname "smtable5"
local tname "table"
local outputfolder "paper/tables/`tfolname'/"
cap mkdir "`outputfolder'"

* Path based on global set in profile.do
	local path "$congestion_pricing_path"
	cd "`path'"

*** LOAD 
	use "data/coded_gps_dta/coded_days_15.dta", clear

****************
**** SAMPLE ****
****************

*** Drop dates after the study ends (=9) and during the study but not in the experiment (=99)
	assert inlist(study_cycle,0,1,2,3,4,5,9,99)
	drop if inlist(study_cycle,9,99) // drop post-experiment
	// assert sample_study_proper == inlist(study_cycle,0,1,2,3,4,5)

	* sample if all good+bad days
	assert sample_analysis == qual_gb

*** Drop April 26, which was a holiday and wrongly included in participants' timelines
	assert date != mdy(4,26,2017)

*
	gisid uidp date
	sort uidp date

* checks, stats

	** sample sizes
	tab sample_analysis post, m col

	* checks on number of data points per participant 
	drop o1
	bys uidp: egen n_in_exp = sum(post)
	bys uidp: egen n_pre = sum((study_cycle==0)*sample_analysis)
	bys uidp:  gen o1=_n==1
	tab n_in_exp if o1==1
	tab n_pre if o1==1
	// hist n_pre if o1==1, d
	tab n_in_exp if o1==1 
	assert inlist(n_in_exp,19,20,24,25) if o1==1  // perfect (not exact because of 4/26)
	drop o1 n_in_exp n_pre

**************
*** CODING ***
**************
	assert qual_good + qual_bad + qual_halfday + qual_no_data == 1
	assert regular_commuter != .

	rename qual_gb qual_gb_

	// fvset base 8 strat_cell

	tempfile DT_data
	save 	`DT_data', replace

*******************
*** AREA CODING ***
*******************
* keep only AREA respondents, and drop 11 respondents with faulty AREA locations
	drop if in_area_treatment == 0

*** Merge area candidates 
	merge m:1 uidp using "data\coded_gps_dta\area_candidates.dta"
	assert _m==3
	drop _m

* demean chosen detour
	cap drop o1
	bys uidp: gen o1 = _n==1

	* sample where both >=1 short and >=1 long candidates exist
	gen sample_short_long_subtr = ac_valid_short==1 & ac_valid_long==1


*** SAMPLE CLEAN
	drop if sample_area_drop == 1 

	*** add daily rate (120 vs 80 and 240 vs 160)
	preserve
		keep if area_treat == 1
		keep uidp study_cycle ad?
		gduplicates drop
		gisid uidp
		greshape long ad, i(uidp study_cycle) j(study_day)

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

		keeporder uidp study_cycle study_day high_rate_day rate_today hrd_l?

		tempfile daily_charges
		save 	`daily_charges'
	restore

	merge m:1 uidp study_cycle study_day using `daily_charges'

	* some area days are not yet in the GPS data
	count if _m==2
	// br uidp date study_cycle study_day if _m==2
	assert r(N)==9
	drop if _m==2

	* we have daily rate data for every day in treatment
	assert area_treat == (_m==3)
	drop _m

	* 
	recode high_rate_day rate_today hrd_l? (.=0)

	* by day
	gen area_treat_1 = area_treat * (study_day == 1)
	gen area_treat_2 = area_treat * (study_day == 2)
	gen area_treat_3 = area_treat * (study_day == 3)
	gen area_treat_4 = area_treat * (study_day == 4)
	gen area_treat_5 = area_treat * (study_day == 5)

	* recode time in experiment
	assert study_cycle == 0 if time_in_exp == .
	recode time_in_exp (.=0)

	* separate dummies for weeks 1 and 4
	gen area_treat_w1 = area_treat * (study_cycle == 1)
	gen area_treat_w4 = area_treat * (study_cycle == 4)
	assert area_treat_w1 + area_treat_w4 == area_treat


	tempfile AREA_data
	save 	`AREA_data', replace

*********************************************
*** Panel A - DT_TREAT (all 4)
*********************************************

* Sample of weeks during the experiment
	local sample_weeks "_3w" // "_3w_fullc"
	local ifsample_weeks "if sample`sample_weeks'==1 "		

*** file names
	cap mkdir "`outputfolder'"
	cap erase "`outputfolder'`tname'_panel_A.tex"
	cap erase "`outputfolder'`tname'_panel_B.tex"
	cap erase "`outputfolder'`tname'.tex" // main table file

	* estimmate decimal precision
	local dec = 2
		

	di "************************** Running iteration `varp' `time_sample'"

	* load data, decide sample
	use `DT_data', clear
	keep `ifsample_weeks' 

	* sanity checks
	qui assert qual_gb != .
	qui count if qual_gb == 0
	qui assert r(N) != 0

	assert inlist(study_cycle,0,1,2,3,4,5)

**********************
*** individual FEs ***

	display "individual FE"
	areg qual_gb dt_high_post dt_low_post dt_info_post post ///
				i.study_cycle , ab(uidp) vce(cluster uidp)
				
* control mean
	qui sum qual_gb if dt_ctrl_post == 1 & e(sample) == 1 // control mean at endline
	qui estadd scalar control_mean = r(mean)
	qui estadd scalar control_sd = r(sd)

* number of participants in sample
	cap drop esample
	qui gen esample = e(sample)
	cap drop o1
	qui bys uidp esample: gen o1=_n==1
				
	qui count if o1==1 & esample==1
	qui estadd scalar n_total = r(N)

*** Store
	estimates store dt_fe

	estimates describe dt_fe

*********************************************
*** Panel B - AREA and subtreatments 
*********************************************

	local dec=2

* Sample of weeks during the experiment
	local sample_weeks "_wall" //"_w1" // "_w14" "_w4"
	local ifsample_weeks "if sample`sample_weeks'==1 "		

		* load data, decide sample
		use `AREA_data', clear
		keep `ifsample_weeks' 
		keep if study_cycle <= 1

		* sanity checks
		qui assert qual_gb != .
		qui count if qual_gb == 0
		qui assert r(N) != 0

		assert inlist(study_cycle,0,1,2,3,4)

		**********************
		*** individual FEs ***

			display "individual FE"
			areg qual_gb area_treat post, ab(uidp) vce(cluster uidp)
						
		* control mean
			sum qual_gb if area_treat == 0 & sample_w14==1 & e(sample) == 1 // control mean at endline
			estadd scalar control_mean = r(mean)
			estadd scalar control_sd = r(sd)

		* number of participants in sample
			cap drop esample
			gen esample = e(sample)
			cap drop o1
			bys uidp esample: gen o1=_n==1
						
			qui count if o1==1 & esample==1
			estadd scalar n_total = r(N)

		*** Store
			estimates store route_fe


	** Output
	esttab dt_fe route_fe using "`outputfolder'`tname'_panel.tex", replace ///
			 keep(dt_high_post dt_low_post dt_info_post area_treat post) ///
			order(dt_high_post dt_low_post dt_info_post area_treat post) ///
			b(%12.`dec'f) se(%12.`dec'f) ///
			coeflabels( ///
						dt_high_post "High Rate $\times$ Post" dt_low_post "Low Rate $\times$ Post" dt_info_post "Information $\times$ Post" ///
				    area_treat "Route Charges" /// 
						post "Post") ///
			starlevels(* .10 ** .05 *** .01) ///
			nomtitles nonumbers ///
			stats(N control_mean, labels("Observations" "Control Mean") fmt("%12.0fc" "%12.`dec'f")) booktabs ///
			substitute(_ \_ "<" "$<$" "\midrule" "\addlinespace\addlinespace") ///
			fragment 

*** Write main table tex file
	file open myfile using "`outputfolder'`tname'.tex", write replace
	file write myfile "\begin{tabular}{lcc}" _n "\toprule" _n  ///
					  " & (1) & (2) \\" _n ///
					  "Treatment & \textit{Departure Time} & \textit{Route} \\" _n ///
					  "Commuter FE & X & X \\" _n ///
				   	  "\ExpandableInput{tables/`tfolname'/`tname'_panel}" _n ///
					  "\bottomrule" _n "\end{tabular}" _n 

	file close myfile
