
 * This do file: Figure 2
 * This version: 30 July, 2017
 * Authors: Gabriel Kreindler

clear all
pause on
set more off
version 15
set seed 23432


local chain_th "15"	

* Path based on global set in profile.do
	local path "$congestion_pricing_path"
	cd "`path'"
	adopath ++ "code/ado/"

*** FIGURE number
	local fnum "figure2"
	local outputfolder "paper/figures/`fnum'/"
	cap mkdir "`outputfolder'"


*** loop over four types of figures: AM/PM and home-work/all-trips
scalar hwsample 	 = "all" // "all" or "hw"
scalar period_of_day = "am"  // "am" or "pm" 


********************
** Load Trip Data **
********************
	use "data/coded_gps_dta/coded_trips_`chain_th'.dta", clear

**************************
** DAILY OUTCOME CODING **
**************************
	isid uidp date chain
	sort uidp date chain
	
	** SUM OVER TRIPS	
	gen one = 1
	gen sample_trip_all = 1
	assert inlist(sample_trip_ok,0,1)
	assert am_sample * pm_sample * am_pk_sample * pm_pk_sample * am_pk1_sample * pm_pk1_sample != .

	merge m:1 uidp using "data/coded_gps_dta/treat_assign.dta"
	assert _m==3
	drop _m

****************
**** CODING  ***
****************
	*** Drop dates after the study ends (=9) and during the study but not in the experiment (=99)
	assert inlist(study_cycle,0,1,2,3,4,5,9,99)
	drop if inlist(study_cycle,9,99) // drop post-experiment
	* sample if all good+bad days
	assert sample_analysis == qual_gb
	*** Drop April 26, which was a holiday and wrongly included in participants' timelines
	assert date != mdy(4,26,2017)
	isid uidp date chain
	sort uidp date chain

	assert qual_good + qual_bad + qual_halfday + qual_no_data == 1
	assert regular_commuter != .
	assert inlist(study_cycle,0,1,2,3,4,5)
	gen dt_any = dt_high + dt_low
	gen dt_any_post = dt_any * post

	cap drop charge
	local scale = 0.2
	gen charge = max(0, min(t_rel + 1.5, 1, 1.5 - t_rel)) * `scale'

	tempfile alldata
	save 	`alldata', replace

*********************************************
*** Sample
*********************************************

* Sample of weeks during the experiment
	local sample_weeks "_3w" // "_3w_fullc"
	local ifsample_weeks "if sample`sample_weeks'==1 "		

* day data quality sample (GOOD+BAD)
	local quality_sample "& sample_analysis==1 "

*** SAMPLE
	keep `ifsample_weeks' `quality_sample'

*** Home work sample

	local hwsample = "`=hwsample'"
	local ymax = 0.20
	if hwsample == "hw"{
		keep if trip_hwp == 1
	}
	

*** morning or evening
	local period_of_day = "`=period_of_day'"
	if period_of_day == "am"{
		keep if am_sample == 1 // 1 or 0 (for am / pm)
	}
	else{
		assert period_of_day == "pm"
		keep if am_sample == 0 // 1 or 0 (for am / pm)	
	}
	
	
*********************************************
*** Analysis
*********************************************
	keep if inrange(t_rel,-3,3)
	sort t_rel

	*** Graph the four densities
	// local bw = 0.15
	// twoway 	(kdensity t_rel if post == 0 & dt_any == 0, bw(`bw') lcolor(blue) lpattern(dash)) ///
	// 		(kdensity t_rel if post == 0 & dt_any == 1, bw(`bw') lcolor(red)  lpattern(dash)) ///
	// 		(kdensity t_rel if post == 1 & dt_any == 0, bw(`bw') lcolor(blue) lpattern(solid)) ///
	// 		(kdensity t_rel if post == 1 & dt_any == 1, bw(`bw') lcolor(red)  lpattern(solid)) ///
	// 		(line charge t_rel if inrange(t_rel,-3,3), lcolor(blue) lpattern(dash) lwidth(medthick)) ///
	// 		, xline(-1.5 -0.5 0.5 1.5)

	// gen nbins = 72

	keep t_rel dt_any post uidp date

	tempfile mydata
	save 	`mydata'


*** Bootstrap the differences-in-differences in densities
scalar nboot = 1000
local nboot = nboot

forv iboot = 1/`nboot'{
	qui use `mydata', clear

	* bootstrap: sample with replacement at the individual (uidp) level
	qui bsample , cluster(uidp) strata(dt_any) idcluster(new_uidp)

	* Departure time grid
	qui gen tre = -3 + (_n-1)/72 * 6  if _n <= 73

	* the four densities
	local bw = 0.15
	qui kdensity t_rel if post == 0 & dt_any == 0, bw(`bw') gen(kden_post0_treat0) at(tre) nograph
	qui kdensity t_rel if post == 0 & dt_any == 1, bw(`bw') gen(kden_post0_treat1) at(tre) nograph
	qui kdensity t_rel if post == 1 & dt_any == 0, bw(`bw') gen(kden_post1_treat0) at(tre) nograph
	qui kdensity t_rel if post == 1 & dt_any == 1, bw(`bw') gen(kden_post1_treat1) at(tre) nograph

	* adjust by the number of trips per person-day
	qui gunique new_uidp date if post == 0 & dt_any == 0
	qui local n_post0_treat0 = r(N) / r(J)
	qui gunique new_uidp date if post == 0 & dt_any == 1
	qui local n_post0_treat1 = r(N) / r(J)
	qui gunique new_uidp date if post == 1 & dt_any == 0
	qui local n_post1_treat0 = r(N) / r(J)
	qui gunique new_uidp date if post == 1 & dt_any == 1
	qui local n_post1_treat1 = r(N) / r(J)

	// di `n_post0_treat0'
	// di `n_post0_treat1'
	// di `n_post1_treat0'
	// di `n_post1_treat1'

	* the differences-in-differences change in density
	qui gen impact = (kden_post1_treat1 * `n_post1_treat1' - kden_post1_treat0 * `n_post1_treat0') - ///
					 (kden_post0_treat1 * `n_post0_treat1' - kden_post0_treat0 * `n_post0_treat0')

	qui keep if tre != .
	qui keep tre impact

	gen iboot = `iboot'

	tempfile data`iboot'
	save 	`data`iboot''
}

*** Same analysis for main sample
	use `mydata', clear

	* Departure time grid
	gen tre = -3 + (_n-1)/72 * 6  if _n <= 73

	* the four densities
	local bw = 0.15
	kdensity t_rel if post == 0 & dt_any == 0, bw(`bw') gen(kden_post0_treat0) at(tre) nograph
	kdensity t_rel if post == 0 & dt_any == 1, bw(`bw') gen(kden_post0_treat1) at(tre) nograph
	kdensity t_rel if post == 1 & dt_any == 0, bw(`bw') gen(kden_post1_treat0) at(tre) nograph
	kdensity t_rel if post == 1 & dt_any == 1, bw(`bw') gen(kden_post1_treat1) at(tre) nograph

	* adjust by the number of trips per person-day
	qui gunique uidp date if post == 0 & dt_any == 0
	qui local n_post0_treat0 = r(N) / r(J)
	qui gunique uidp date if post == 0 & dt_any == 1
	qui local n_post0_treat1 = r(N) / r(J)
	qui gunique uidp date if post == 1 & dt_any == 0
	qui local n_post1_treat0 = r(N) / r(J)
	qui gunique uidp date if post == 1 & dt_any == 1
	qui local n_post1_treat1 = r(N) / r(J)

	* the differences-in-differences change in density
	qui gen impact = (kden_post1_treat1 * `n_post1_treat1' - kden_post1_treat0 * `n_post1_treat0') - ///
					 (kden_post0_treat1 * `n_post0_treat1' - kden_post0_treat0 * `n_post0_treat0')

	keep if tre != .
	keep tre impact
	rename impact real_impact

	gen iboot = -1

	tempfile data_main
	save 	`data_main'

*** Append together bootstrap runs to get 2.5 and 97.5% CI limits
	use `data1', clear
	forv iboot = 2/`nboot'{
		append using `data`iboot''
	}
	gcollapse (p2.5) impact_2_5=impact (p97.5) impact_97_5=impact, by(tre)

	append using `data_main'

*** Graph
	// local ymax = 0.30
	local ymax = 0.20

	cap drop charge
	local scale = 0.1
	gen charge = max(0, min(tre + 1.5, 1, 1.5 - tre)) * `scale'

	sort tre
	twoway 	(line charge tre if inrange(tre,-3,3), lcolor(blue) lpattern(dash)) ///
			(line impact_2_5 impact_97_5 tre, lcolor(gs8 gs8) lpattern(dash dash)) ///
			(line real_impact tre if iboot == -1, lcolor(red)) ///
			, legend(order( ///
				4 "Change in number of trips" ///
				2 "95% CI" ///
				1 "Departure Time Rate" ) cols(1) position(1) ring(0) ) ///
		 	graphregion(color(white)) ///
		 	xlabel(-2.5 -1.5 -0.5 0 0.5 1.5 2.5, grid) ///
			xtitle("Departure time (hours relative to peak)") ///
			ylabel(-`ymax'(0.10)`ymax', angle(0)) ///
			ytitle("Trips" "per day" "per 1 hr", orientation(horiz))

		 	// yscale(titlegap(*-30)) ///


*** Save 
	local outputfile "`outputfolder'figure_bw15_boot`=nboot'_`=hwsample'_`=period_of_day'"
	graph export "`outputfile'.eps", replace
	graph export "`outputfile'.png", replace
	graph export "`outputfile'.pdf", replace
