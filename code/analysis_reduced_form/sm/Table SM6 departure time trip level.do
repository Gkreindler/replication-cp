
 * This do file: SM version of Table 1, daily level
 * This version: August 8, 2017
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
local tfolname "smtable6"
local tname "table"
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

* sample if all good+bad days
	assert sample_analysis == qual_gb

*** Drop April 26, which was a holiday and wrongly included in participants' timelines
	assert date != mdy(4,26,2017)

*
	gisid uidp date
	sort uidp date

**************
*** CODING ***
**************
	assert qual_good + qual_bad + qual_halfday + qual_no_data == 1
	assert regular_commuter != .

	rename qual_gb qual_gb_

	tempfile alldata
	save 	`alldata', replace

*********************************************
*** DT_TREAT (all 4)
*********************************************

* Sample of weeks during the experiment
	local sample_weeks "_3w" // "_3w_fullc"
	local ifsample_weeks "if sample`sample_weeks'==1 "		

* day data quality sample (GOOD+BAD)
	local quality_sample "& sample_analysis==1 "

* trip quality sample 
	local trip_qual_sample "ok" //"all" 

*** file names
	cap mkdir "`outputfolder'"
	cap erase "`outputfolder'`tname'_panel_A.tex"
	cap erase "`outputfolder'`tname'_panel_B.tex"
	cap erase "`outputfolder'`tname'_panel_C.tex"
	cap erase "`outputfolder'`tname'.tex" // main table file


	local varp "ra_ok" // "ch_"

	* User and trip type (regular / irregular commuter, all trips vs HW)
	local maxj = 3
	forv j=1/`maxj'{
		local tnam: word `j' of "panel_A" "panel_B" "panel_C"
		local trip_sample: word `j' of 	"" ///
										"hwp_" ///
										""

		local user_sample: word `j' of 	"" /// 
										"keep if inlist(regular_commuter,2,3)" ///
										"keep if regular_commuter==4"

		* Time period: all, AM, AM-pre, AM-post, PM, PM-pre, PM-post
		local maxl = 7
		forv l=1/`maxl'{
			local time_sample: word `l' of "" "am" "pream" "postam" "pm" "prepm" "postpm"

			* outcome variable
			local yvar "`varp'_`trip_sample'`time_sample'"
			local estp "`tnam'_`trip_sample'_`time_sample'"

			di "************************** Running iteration `estp'"

		* load data, decide sample
			use `alldata', clear

			* date and data quality sample
			keep `ifsample_weeks' `quality_sample'

			* keep all/regular/regular/irregular commuters
			`user_sample'

			* sanity checks
			qui assert `yvar' != .	
			assert inlist(study_cycle,0,1,2,3,4)

		**********************
		*** individual FEs ***

			display "individual FE"
			qui areg `yvar' dt_wcha_post post ///
						i.study_cycle dt_wcha, ab(`regFE') vce(cluster uidp)
						
		* control mean
			qui sum `yvar' if dt_noch_post == 1 & e(sample) == 1 // control mean at endline
			qui estadd scalar control_mean = r(mean)
			qui estadd scalar control_sd = r(sd)

		* number of participants in sample
			cap drop esample
			qui gen esample = e(sample)
			cap drop o1
			qui bys uidp esample: gen o1=_n==1
						
			qui count if o1==1 & esample==1
			qui estadd scalar n_total = r(N)

			foreach group in "low" "high" "info" "ctrl"{
				qui count if dt_`group' == 1 & o1==1 & esample==1
				qui estadd scalar n_`group' = r(N)	
			}

		*** Store
			qui estimates store `estp'_fe

		} //l


	** Output
		local dec = 2
		esttab * using "`outputfolder'`tname'_`tnam'.tex", ///
				append keep(post dt_wcha_post) ///
				b(%12.`dec'f) se(%12.`dec'f) nomtitles nonumbers ///
				coeflabels(post "Post" dt_wcha_post "Charges $\times$ Post") ///
				starlevels(* .10 ** .05 *** .01) ///
				stats(N control_mean, labels("Observations" "Control Mean") ///
					fmt("%12.0fc" "%12.`dec'f")) booktabs ///
				substitute(_ \_ "<" "$<$" "\midrule" "\addlinespace") ///
				fragment 

	estimates clear

	} //j


*** Write main table tex file
	file open myfile using "`outputfolder'`tname'.tex", write replace
	file write myfile /* "\resizebox*{\textwidth}{!}{" _n */ "\begin{tabular}{lc@{\hskip 0.25in}ccc@{\hskip 0.3in}ccc}" _n "\toprule" _n  ///
					  " & (1) & (2) & (3) & (4) & (5) & (6) & (7)\\" _n ///
					  "Time of Day & AM \& PM & \multicolumn{3}{c}{AM}     & \multicolumn{3}{c}{PM}  \\" _n ///
					  "            &          & all & \thead{pre\\ peak} & \thead{post\\ peak} & all & \thead{pre\\ peak} & \thead{post\\ peak} \\" _n ///
					  "Commuter FE & X & X & X & X & X & X & X \\" _n ///
					  "\addlinespace\addlinespace\multicolumn{6}{l}{\emph{Panel A. Full Sample}} \\" _n ///
				   	  "\ExpandableInput{tables/`tfolname'/`tname'_panel_A}" _n ///
				   	  "\addlinespace\addlinespace\multicolumn{6}{l}{\emph{Panel B. Regular Commuters, Home-Work and Work-Home Trips}} \\" _n ///
				   	  "\ExpandableInput{tables/`tfolname'/`tname'_panel_B}" _n ///
				   	  "\addlinespace\addlinespace\multicolumn{6}{l}{\emph{Panel C. Variable Commuters, All Trips}} \\" _n ///
				   	  "\ExpandableInput{tables/`tfolname'/`tname'_panel_C}" _n ///
					  "\bottomrule" _n "\end{tabular}" _n  /* "}" _n  */

	file close myfile


*** Write main table tex file --- SHORT version
	file open myfile using "`outputfolder'`tname'_short.tex", write replace
	file write myfile /* "\resizebox*{\textwidth}{!}{" _n */ "\begin{tabular}{lc@{\hskip 0.25in}ccc@{\hskip 0.3in}ccc}" _n "\toprule" _n  ///
					  " & (1) & (2) & (3) & (4) & (5) & (6) & (7)\\" _n ///
					  "Time of Day & AM \& PM & \multicolumn{3}{c}{AM}     & \multicolumn{3}{c}{PM}  \\" _n ///
					  "            &          & all & \thead{pre\\ peak} & \thead{post\\ peak} & all & \thead{pre\\ peak} & \thead{post\\ peak} \\" _n ///
					  "Commuter FE & X & X & X & X & X & X & X \\" _n ///
				   	  "Sample:     & \multicolumn{7}{c}{ \textit{Regular Commuters, Home-Work and Work-Home Trips}} \\" _n ///
				   	  "\addlinespace\ExpandableInput{tables/`tfolname'/`tname'_panel_B}" _n ///
					  "\bottomrule" _n "\end{tabular}" _n  /* "}" _n  */

	file close myfile
