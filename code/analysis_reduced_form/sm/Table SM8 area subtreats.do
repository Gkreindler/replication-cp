
 * This do file: Area Subtreatments
 * This version: Jan 2021
 * Authors: Gabriel Kreindler

clear all
pause off
set more off
set matsize 10000

*** TABLE number
local tfolname "smtable8"
local tname "table"
local outputfolder "paper/tables/`tfolname'/"

* Path based on global set in profile.do
	local path "$congestion_pricing_path"
	cd "`path'"

*** Experiment
	use "data/coded_gps_dta/treat_assign.dta", clear
	keep uidp dttreat atreat dt_high dt_low dt_info dt_ctrl ///
		 strat_cell a_early a_late a_short a_low adetourdist adetourtime
	gen in_area_treat = atreat != "0 No treatment"
	assert _N==497

	keep if in_area_treat == 1

	// rename uidp uidp_original in_area_treat
	keep uidp in_area_treat a_short a_low adetourdist adetourtime
	tempfile treat_status
	save 	`treat_status'


*** LOAD 
	use "data/coded_gps_dta/coded_days_15.dta", clear

****************
**** SAMPLE ****
****************

*** Drop dates after the study ends (=9) and during the study but not in the experiment (=99)
	assert inlist(study_cycle,0,1,2,3,4,5,9,99)
	drop if inlist(study_cycle,9,99) // drop post-experiment
	// assert sample_study_proper == inlist(study_cycle,0,1,2,3,4,5)

	* keep only AREA respondents, and drop 11 respondents with faulty AREA locations
	drop if in_area_treatment == 0

*** Drop April 26, which was a holiday and wrongly included in participants' timelines
	assert date != mdy(4,26,2017)

	gisid uidp date
	sort uidp date

* checks, stats

	** sample sizes
	tab sample_analysis post, m col

**************
*** CODING ***
**************
	assert post == inlist(study_cycle,1,2,3,4,5,99)  	// check that the definitions are equivalent
	assert qual_good + qual_bad + qual_halfday + qual_no_data == 1
	assert regular_commuter != .
	* sample if all good+bad days
	assert sample_analysis == qual_gb

	*** Merge area candidates 
	merge m:1 uidp using "data\coded_gps_dta\area_candidates.dta"
	assert _m==3
	drop _m

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

		bys uidp: gegen min_charge = min(ad)
		bys uidp: gegen max_charge = max(ad)

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

	* some area days are not in the GPS data
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

	tempfile AREA_data
	save 	`AREA_data', replace


*********************************************
*** AREA subtreatments IV
*********************************************

* Sample of weeks during the experiment
	local sample_weeks "_wall" //"_w1" // "_w14" "_w4"
	local ifsample_weeks "if sample`sample_weeks'==1 "			

	keep `ifsample_weeks' & sample_analysis==1

	tempfile alldata
	save 	`alldata'

*** HIGH/LOW
	gegen uidp_ = group(uidp)
	gen post_high = post * a_high
	gegen study_cycle_high = group(study_cycle a_high)

	cap drop o1
	bys uidp: gen o1=_n==1


	* RF (full)
	// areg cha_ok_ area_treat_high area_treat_low i.study_cycle_high, ab(uidp) vce(cluster uidp)
	// 	sum cha_ok_ if area_treat == 0 & post == 1 & e(sample) == 1	

	areg cha_ok_ area_treat_high area_treat_low i.study_cycle_high if study_cycle <= 1, ab(uidp) vce(cluster uidp)
	sum cha_ok_ if area_treat == 0 & post == 1 & e(sample) == 1	

	estadd scalar control_mean = r(mean)

		sum o1 if e(sample)
		di r(sum)
	estadd scalar nuidp = r(sum)

		test area_treat_high == area_treat_low
	estadd scalar pval = r(p)

	estimates store hl_rf





**************
**************
*** SHORT/LONG
	use `alldata', clear 

	* sample where both >=1 short and >=1 long candidates exist
	keep if sample_short_long_subtr == 1

	gegen uidp_ = group(uidp)
	gen post_short = post * a_short
	gegen study_cycle_short = group(study_cycle a_short)

	cap drop o1
	bys uidp: gen o1=_n==1

*** reduced form
	areg cha_ok_ area_treat_short area_treat_long i.study_cycle_short if study_cycle <= 1, ab(uidp) vce(cluster uidp)

	sum cha_ok_ if area_treat == 0 & post == 1 & e(sample) == 1

	estadd scalar control_mean = r(mean)
		sum o1 if e(sample)
		di r(sum)
	estadd scalar nuidp = r(sum)

		test area_treat_short == area_treat_long
	estadd scalar pval = r(p)

	estimates store sl_rf

**** 
	cap mkdir "`outputfolder'"
	// cap erase "`outputfolder'`tname')sma.tex"

	** Output Panel A
	esttab hl_rf sl_rf using "`outputfolder'`tname'_small.tex", ///
			 replace keep(area_treat_high area_treat_low area_treat_short area_treat_long) ///
				b(%12.1f) se(%12.1f) ///
				nomtitles mgroups("Hypothetical Route Charges", pattern(1 0) ///
				prefix(\multicolumn{@span}{c}{) suffix(}) span erepeat(\cmidrule(lr){@span})) ///
				coeflabels(	area_treat_high "Treated $\times$ High Rate" ///
						 	area_treat_low  "Treated $\times$ Low Rate" ///
							area_treat_short "Treated $\times$ Short Detour" ///
						 	area_treat_long  "Treated $\times$ Long Detour") ///
				starlevels(* .10 ** .05 *** .01) ///
				stats(N nuidp control_mean pval, ///
					labels("Observations" "Commuters" "Control Mean" "P-val Equal Sub-treatment Effects") ///
					 fmt("%12.0fc" "%12.0fc" "%12.1f" "%12.2fc")) booktabs ///
				substitute(_ \_ "<" "$<$" "\midrule" "\addlinespace\addlinespace")
