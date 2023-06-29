
 * This do file: Table 1 (departure time impact, daily)
 * This version: Feb 11 2020
 * Authors: Gabriel Kreindler

clear all
pause off
set more off
set matsize 10000

* What type of regression? 
	* uidp = individual FE
	* strate_cell_new = strata fixed effects
local regFE "uidp" // "strat_cell_new"

*** TABLE number
local tfolname "table1"
local outputfolder "paper/tables/`tfolname'/"

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
	isid uidp date
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

	tempfile alldata
	save 	`alldata', replace

*********************************************
*** DT_TREAT (all 4)
*********************************************

* Sample of weeks during the experiment
	local sample_weeks "_3w" // "_3w_fullc"
	local ifsample_weeks "if sample`sample_weeks'==1 "		

*** file names
	cap mkdir "`outputfolder'"
	cap erase "`outputfolder'table_daily_panel_A.tex"
	cap erase "`outputfolder'table_daily_panel_B.tex"
	cap erase "`outputfolder'table_daily.tex" // main table file

*** panel A
forv i=1/2{
	local varp: word `i' of "ra_ok"	"nt_ok"

	forv l=1/3{
		local time_sample: word `l' of "" "am" "pm"

		* outcome variable
		local yvar "`varp'_`time_sample'"
		local estp "`varp'_`time_sample'"

		di "************************** Running iteration `varp' `time_sample'"

		* load data, decide sample
		use `alldata', clear
		keep `ifsample_weeks' & sample_analysis==1

		* sanity checks
		qui assert `yvar' != .
		qui count if `yvar' == 0
		qui assert r(N) != 0

		assert inlist(study_cycle,0,1,2,3,4,5)

	**********************
	*** individual FEs ***

		display "individual FE"
		qui areg `yvar' dt_high_post dt_low_post dt_info_post post ///
					i.study_cycle dt_high dt_low dt_info, ab(`regFE') vce(cluster uidp)
					
	* control mean
		qui sum `yvar' if dt_ctrl_post == 1 & e(sample) == 1 // control mean at endline
		qui estadd scalar control_mean = r(mean)
		qui estadd scalar control_sd = r(sd)

	*** Store
		qui estimates store `estp'_fe

	} //l
} //i

	** Output
	esttab ra_* nt_* using "`outputfolder'table_daily_panel_A.tex", ///
			append keep(post dt_info_post dt_low_post dt_high_post) ///
			b(%12.2f) se(%12.2f) ///
			coeflabels(post "Post" dt_info_post "Information $\times$ Post" dt_low_post "Low Rate $\times$ Post" dt_high_post "High Rate $\times$ Post") ///
			starlevels(* .10 ** .05 *** .01) nomtitles nonumbers ///
			stats(N control_mean, labels("Observations" "Control Mean") fmt("%12.0fc" "%12.2f")) booktabs ///
			substitute(_ \_ "<" "$<$" "\midrule" "\addlinespace") ///
			fragment 

	estimates clear


*** panel B
forv i=1/2{
	local varp: word `i' of "ra_ok"	"nt_ok"

	forv l=1/3{
		local time_sample: word `l' of "" "c1am" "c1pm"

		* outcome variable
		local yvar "`varp'_`time_sample'"
		local estp "`varp'_`time_sample'"

		di "************************** Running iteration `varp' `time_sample'"

		* load data, decide sample
		use `alldata', clear
		keep `ifsample_weeks' & sample_analysis==1

		* sanity checks
		qui assert `yvar' != .
		qui count if `yvar' == 0
		qui assert r(N) != 0

		assert inlist(study_cycle,0,1,2,3,4,5)

	**********************
	*** individual FEs ***

		display "individual FE"
		qui areg `yvar' dt_wcha_post post ///
					i.study_cycle dt_wcha, ab(`regFE') vce(cluster uidp)
					
	* control mean
		qui sum `yvar' if dt_noch_post == 1 & e(sample) == 1 // control mean at endline
		qui estadd scalar control_mean = r(mean)
		qui estadd scalar control_sd = r(sd)

	*** Store
		qui estimates store `estp'_fe

	} //l
} //i

	** Output
	esttab ra_* nt_* using "`outputfolder'table_daily_panel_B.tex", ///
			append keep(post dt_wcha_post) ///
			b(%12.2f) se(%12.2f) ///
			coeflabels(post "Post" dt_wcha_post "Charges $\times$ Post" ) ///
			starlevels(* .10 ** .05 *** .01) nomtitles nonumbers ///
			stats(N control_mean, labels("Observations" "Control Mean") fmt("%12.0fc" "%12.2f")) booktabs ///
			substitute(_ \_ "<" "$<$" "\midrule" "\addlinespace") ///
			fragment 

	estimates clear


*** Write main table tex file
	file open myfile using "`outputfolder'table_daily.tex", write replace
	file write myfile "\begin{tabular}{lccc@{\hskip 0.25in}ccc}" _n "\toprule" _n  ///
					  " & (1) & (2) & (3) & (4) & (5) & (6) \\" _n ///
					  "Outcome & \multicolumn{3}{c}{\textit{Total Hypothetical Rates Today}} & \multicolumn{3}{c}{\textit{Number of Trips Today}}  \\" _n ///
					  "Time of Day & AM \& PM & AM & PM   &   AM \& PM & AM & PM \\" _n ///
					  "Commuter FE & X & X & X & X & X & X \\" _n ///
					  "\addlinespace\addlinespace\multicolumn{7}{l}{\emph{Panel A. All Departure Time Sub-Treatments}} \\" _n ///
				   	  "\ExpandableInput{tables/`tfolname'/table_daily_panel_A}" _n ///
				   	  "\addlinespace\addlinespace\multicolumn{7}{l}{\emph{Panel B. Any Departure Time Charge vs. Control or Information}} \\" _n ///
				   	  "\ExpandableInput{tables/`tfolname'/table_daily_panel_B}" _n ///
					  "\bottomrule" _n "\end{tabular}" _n ///

	file close myfile
