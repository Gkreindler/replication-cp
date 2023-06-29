
 * This do file: main figure 
 * This version: Nov 23, 2021
 * Authors: Gabriel Kreindler

clear all
pause off
set more off
set matsize 10000

*** TABLE number
local tfolname "smtable7"
local tname "table"
local outputfolder "paper/tables/`tfolname'/"
	cap mkdir "`outputfolder'"
	cap erase "`outputfolder'`tname'_panel_A.tex"
	cap erase "`outputfolder'`tname'_panel_B.tex"
	cap erase "`outputfolder'`tname'.tex" // main table file

* Path based on global set in profile.do
	local path "$congestion_pricing_path"
	cd "`path'"


*** LOAD 
	use "data\coded_gps_dta\area_coded", clear


*** more coding
	forv i=4(-1)0{
		gen sc_`i'_early = (study_cycle == `i') * a_early
		gen sc_`i'_late  = (study_cycle == `i') * a_late
	}


************************
**** Trip-level detour rate
************************


*** collapse to mean detour route rate per uidp -- 
	// gcollapse (mean) is_long_route, by(uidp study_cycle sc_* a_early a_late has_long_route_pre)

	gunique uidp
	gunique uidp if has_long_route_pre == 1

*** panel A - Sample = all commuters

	eststo A1: reghdfe is_long_route sc_1_early ///
			if study_cycle <= 1, ab(uidp study_cycle) vce(cl uidp)
	sum is_long_route if sc_1_late == 1 
	estadd scalar control_mean = r(mean)


	eststo A2: reghdfe is_long_route sc_1_early sc_2_early ///
			if study_cycle <= 2, ab(uidp study_cycle) vce(cl uidp)
	sum is_long_route if sc_1_late == 1 
	estadd scalar control_mean = r(mean)


*** panel B. Sample = commuters who used the detour route at baseline
	keep if has_long_route_pre == 1

	gunique uidp

	eststo B1: reghdfe is_long_route sc_1_early ///
			if study_cycle <= 1, ab(uidp study_cycle) vce(cl uidp)
	sum is_long_route if sc_1_late == 1 
	estadd scalar control_mean = r(mean)

	eststo B2: reghdfe is_long_route sc_1_early sc_2_early ///
			if study_cycle <= 2, ab(uidp study_cycle) vce(cl uidp)
	sum is_long_route if sc_1_late == 1 
	estadd scalar control_mean = r(mean)







****************
**** DAILY SAMPLE ****
****************

*** baseline use
	use "data\coded_gps_dta\area_coded", clear

	keep uidp has_long_route_pre
	gduplicates drop
	gisid uidp

	tempfile baselineuse
	save 		`baselineuse'

*** LOAD 
	use "data/coded_gps_dta/coded_days_15.dta", clear

*** SAMPLE

*** Drop dates after the study ends (=9) and during the study but not in the experiment (=99)
	assert inlist(study_cycle,0,1,2,3,4,5,9,99)
	drop if inlist(study_cycle,9,99) // drop post-experiment
	// assert sample_study_proper == inlist(study_cycle,0,1,2,3,4,5)

	* keep only AREA respondents, and drop 11 respondents with faulty AREA locations
	drop if in_area_treatment == 0
	drop if sample_area_drop == 1 

*** Drop April 26, which was a holiday and wrongly included in participants' timelines
	assert date != mdy(4,26,2017)

	isid uidp date
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

	gen area_treat_w1 = area_treat * (study_cycle == 1)
	gen area_treat_w4 = area_treat * (study_cycle == 4)
	assert area_treat_w1 + area_treat_w4 == area_treat

	rename qual_gb qual_gb_

	// fvset base 8 strat_cell

	// tempfile alldata
	// save 	`alldata', replace

	forv i=4(-1)0{
		gen sc_`i'_early = (study_cycle == `i') * a_early
		gen sc_`i'_late  = (study_cycle == `i') * a_late
	}


*** Sample
	keep if sample_wall == 1


*** Add data on baseline detour route
	merge m:1 uidp using `baselineuse'
	assert _m!=2
	// keep if _m==3
	drop _m


************************
**** Trip-level detour rate
************************

*** panel A - Sample = all commuters

	eststo A3: reghdfe nt_ok_hwp_ sc_1_early ///
			if study_cycle <= 1, ab(uidp study_cycle) vce(cl uidp)
	sum nt_ok_hwp_ if sc_1_late == 1 
	estadd scalar control_mean = r(mean)


	eststo A4: reghdfe nt_ok_hwp_ sc_1_early sc_2_early ///
			if study_cycle <= 2, ab(uidp study_cycle) vce(cl uidp)
	sum nt_ok_hwp_ if sc_1_late == 1 
	estadd scalar control_mean = r(mean)


*** panel B. Sample = commuters who used the detour route at baseline
	keep if has_long_route_pre == 1

	eststo B3: reghdfe nt_ok_hwp_ sc_1_early ///
			if study_cycle <= 1, ab(uidp study_cycle) vce(cl uidp)
	sum nt_ok_hwp_ if sc_1_late == 1 
	estadd scalar control_mean = r(mean)

	eststo B4: reghdfe nt_ok_hwp_ sc_1_early sc_2_early ///
			if study_cycle <= 2, ab(uidp study_cycle) vce(cl uidp)
	sum nt_ok_hwp_ if sc_1_late == 1 
	estadd scalar control_mean = r(mean)

	

** Output Panel A
local dec=2

esttab A1 A2 A3 A4 using "`outputfolder'`tname'_panel_A.tex", ///
		append keep(sc_1_early sc_2_early ) ///
		b(%12.`dec'f) se(%12.`dec'f) ///
		coeflabels(	sc_1_early "Treatment: Early $ \times $ week 1" ///
					sc_2_early "Persistence: Early $ \times $ week 2") ///
		starlevels(* .10 ** .05 *** .01) ///
		nomtitles nonumbers ///
		stats(N control_mean, labels("Observations" "Control Mean (week 1)") ///
		 fmt("%12.0fc" "%12.`dec'f")) booktabs ///
		substitute(_ \_ "<" "$<$" "\midrule" "\addlinespace\addlinespace") ///
		fragment 

** Output Panel B
esttab B1 B2 B3 B4 using "`outputfolder'`tname'_panel_B.tex", ///
		append keep(sc_1_early sc_2_early ) ///
		b(%12.`dec'f) se(%12.`dec'f) ///
		coeflabels(	sc_1_early "Treatment: Early $ \times $ week 1" ///
					sc_2_early "Persistence: Early $ \times $ week 2") ///
		starlevels(* .10 ** .05 *** .01) ///
		nomtitles nonumbers ///
		stats(N control_mean, labels("Observations" "Control Mean (week 1)") ///
			fmt("%12.0fc" "%12.`dec'f")) booktabs ///
		substitute(_ \_ "<" "$<$" "\midrule" "\addlinespace\addlinespace") ///
		fragment 

*** Write main table tex file
file open myfile using "`outputfolder'`tname'.tex", write replace
file write myfile /* "\resizebox*{!}{0.8\textheight}{" _n */ "\begin{tabular}{lcccc}" _n "\toprule" _n  ///
				  " & (1) & (2) & (3) & (4) \\" _n ///
				  "Outcome & \multicolumn{2}{c}{\textit{Use Detour Route}} & \multicolumn{2}{c}{\textit{Number of Trips Today}}  \\" _n ///
				  "Commuter FE & X & X & X & X  \\" _n ///
				  "\addlinespace\addlinespace\multicolumn{5}{l}{\emph{Panel A. All Commuters}} \\" _n ///
			   	  "\ExpandableInput{tables/`tfolname'/`tname'_panel_A}" _n ///
			   	  "\addlinespace\addlinespace\multicolumn{5}{l}{\emph{Panel B. Commuters Who Used Detour at Baseline}} \\" _n ///
			   	  "\ExpandableInput{tables/`tfolname'/`tname'_panel_B}" _n ///
				  "\bottomrule" _n "\end{tabular}" _n /* "}" _n  */

file close myfile
