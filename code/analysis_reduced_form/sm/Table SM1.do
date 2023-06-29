
 * This do file: Descriptive stats of GPS trip data
 * This version: June 10, 2019
 * Authors: Gabriel Kreindler

clear all
pause off
set more off
set matsize 10000

*** Table number
	local tnum "smtable1"
	local outputfolder "paper/tables/`tnum'/"
	cap mkdir "`outputfolder'"

* Path based on global set in profile.do
	local path "$congestion_pricing_path"
	cd "`path'"
	adopath ++ "code/ado/"

*** Dummy for being at home/work
	use "data/coded_gps_dta/coded_days_15.dta", clear
	keep uidp date at_h at_w at_w_all
	isid uidp date
	tempfile at_hw
	save 	`at_hw'


*** LOAD full gps data
	use "data/coded_gps_dta/coded_trips_15.dta", clear

	merge m:1 uidp date using `at_hw'
	assert _m==3
	drop _m

****************
**** SAMPLE ****
****************

*** Drop dates after the study ends (=9) and during the study but not in the experiment (=99)
	assert inlist(study_cycle,0,1,2,3,4,5,9,99)
	drop if inlist(study_cycle,9,99) // drop post-experiment

	* Drop DAYS without trips
	assert trip == 1 - date_wo_trips
	keep if trip == 1

* sample if all good+bad days
	assert sample_analysis == qual_gb

*** Drop April 26, which was a holiday and wrongly included in participants' timelines
	assert date != mdy(4,26,2017)

* checks
	isid uidp date chain
	sort uidp date chain

**************
*** CODING ***
**************
	assert qual_good + qual_bad + qual_halfday + qual_no_data == 1
	assert regular_commuter != .

	fvset base 8 strat_cell

	tempfile alldata
	save 	`alldata', replace

*********************************************
*** Panel A
*********************************************

	use `alldata', clear
	// sum plb
	// di %12.0fc r(sum)

	keep if qual_gb == 1
	keep if sample_trip_ok == 1

	gen regular = inlist(regular_commuter,2,3)

**** CODING
	sort uidp date pm_a_sample t_start

	* check that implicitly sorted by t_start
	cap drop temp
	by uidp date: gen temp = _n
	sort uidp date t_start
	by uidp date: assert temp == _n
	drop temp


	by uidp: gen o1 = _n == 1
	by uidp date: gen od1 = _n == 1


	*** Average number of trips per day
	cap drop temp
	by uidp date: egen temp = sum(sample_trip_ok)
	by uidp: egen n_trips = mean(temp) if od1==1
	sum n_trips if o1==1
	sum n_trips if o1==1 & regular==1
	sum n_trips if o1==1 & regular==0

	* pre period only - identical
	// cap drop temp
	// by uidp date: egen temp = sum((post == 0))
	// by uidp: egen n_trips_pre = mean(temp) if od1==1 & post == 0
	// sum n_trips_pre if o1==1


** AT WORK 
	* Fraction of trips that are HW or HW2
	by uidp: egen trips_hwp_frac = mean(trip_hwp)

	* Fraction of trips that are WW or WW2 or W2W2
	by uidp: egen trips_wwp_frac = mean(trip_ww2)

	// sum trips_hwp_frac trips_wwp_frac if o1==1, d


	*** Average number of HW trips per day
	// cap drop temp
	// by uidp date: egen temp = sum(trip_hw == 1 | trip_hw2 == 1) if regular==1
	// by uidp: egen n_trips_hw = mean(temp)			   			if regular==1 & od1==1
	// sum n_trips_hw if o1==1 & regular==1

	*** Fraction of days when visit work
	by uidp: egen at_w_frac_days = mean(at_w_all) if od1==1
	sum at_w_frac_days if o1==1 & regular ==1, d

** Median trip duration for each commuter
	by uidp: egen dur_median = median(dur_trips)	
	by uidp: egen dur_median_hwp = median(dur_trips)	if trip_hwp == 1

** Median trip distance for each commuter
	by uidp: egen plb_median = median(plb)	


*** SDs
	by uidp date: gen trip1_am	   = _n==1 if am_a_sample == 1
	by uidp date: gen trip1_am_hwp = _n==1 if am_a_sample == 1 & oh == 1 & dw ==1 

	by uidp date: gen trip1_pm	   = _n==_N if pm_a_sample == 1
	by uidp date: gen trip1_pm_hwp = _n==_N if pm_a_sample == 1 & ow == 1 & dh ==1 

	tab trip1_am_hwp
	tab trip1_pm_hwp

* all trips
	by uidp: egen dtsd_am = sd(t_start) if am_a_sample == 1
	by uidp: egen dtsd_pm = sd(t_start) if pm_a_sample == 1

* first AM trip, last PM trip
	by uidp: egen dtsd_t1_am = sd(t_start) if trip1_am == 1
	by uidp: egen dtsd_t1_pm = sd(t_start) if trip1_pm == 1

* first AM trip H-> W, last PM trip W->H
	by uidp: egen dtsd_t1hw_am = sd(t_start) if trip1_am_hwp == 1
	by uidp: egen dtsd_t1wh_pm = sd(t_start) if trip1_pm_hwp == 1

	fill_by_group dtsd_* dur_median_hwp, fillby(uidp) replace

	// sum dur_median if o1==1, d
	// sum dur_median_hwp if o1==1, d

*******************
*** file names
	cap mkdir "`outputfolder'"
	cap erase "`outputfolder'panel_A.tex"
	cap erase "`outputfolder'panel_B.tex"
	cap erase "`outputfolder'panel_C.tex"
	cap erase "`outputfolder'table_1.tex" // main table file

*** Write to table 

*** PANEL A
	replace n_trips 	= . if o1==0
	replace dur_median 	= . if o1==0
	replace plb_median 	= . if o1==0
	qui estpost summarize sample_trip_ok n_trips dur_median plb_median, d
	qui estimates store ests_A

	esttab ests_A using "`outputfolder'panel_A.tex", append ///
			cells((p50(label("Median") fmt("%4.2f")) mean(label("Mean") fmt("%4.2f")) sd(par([ ]) label("Std. Dev.") fmt("%4.2f")) ///
			p10(label("10 Perc.") fmt("%4.2f")) p90(label("90 Perc.") fmt("%4.2f")) count(label("Obs.") fmt("%7.0fc")))) nomtitles nonumbers ///
			coeflabels(sample_trip_ok "Total Number of Trips" n_trips "Number of Trips per Day" ///
				dur_median "Median trip duration (minutes)" plb_median "Median trip length (Km.)") ///
			substitute(_ \_ "<" "$<$" "\midrule" "\addlinespace\addlinespace" "\hline" "\addlinespace\addlinespace") ///
			noobs fragment //prehead( "\begin{tabular}{lcccccc}" "\toprule") //  "\resizebox{!}{\textheight}{"


*** PANEL B
	replace trips_hwp_frac = . if regular == 0
	replace trips_wwp_frac = . if regular == 0
	replace at_w_frac_days = . if regular == 0
	qui estpost summarize regular trips_hwp_frac trips_wwp_frac at_w_frac_days if o1==1, d
	qui estimates store ests_B

	esttab ests_B using "`outputfolder'panel_B.tex", append ///
			cells((p50(label("Median") fmt("%4.2f")) mean(label("Mean") fmt("%4.2f")) sd(par([ ]) label("Std. Dev.") fmt("%4.2f")) ///
			p10(label("10 Perc.") fmt("%4.2f")) p90(label("90 Perc.") fmt("%4.2f")) count(label("Obs.") fmt("%4.0fc")))) nomtitles nonumbers  ///
			coeflabels(regular "Regular Commuter" trips_hwp_frac "Frac. trips Home-Work, Work-Home" trips_wwp_frac "Frac. of trips Work-Work" ///
				at_w_frac_days "Frac. of days present at Work") ///
			substitute(_ \_ "<" "$<$" "\midrule" "\addlinespace\addlinespace" "\hline" "\addlinespace\addlinespace") ///
			noobs fragment 
			

*** PANEL C
	replace dtsd_t1hw_am = . if regular == 0
	replace dtsd_t1wh_pm = . if regular == 0
	qui estpost summarize dtsd_t1_am dtsd_t1_pm dtsd_t1hw_am dtsd_t1wh_pm if o1==1, d
	qui estimates store ests_C

	esttab ests_C using "`outputfolder'panel_C.tex", append ///
			cells((p50(label("Median") fmt("%4.2f")) mean(label("Mean") fmt("%4.2f")) sd(par([ ]) label("Std. Dev.") fmt("%4.2f")) ///
			p10(label("10 Perc.") fmt("%4.2f")) p90(label("90 Perc.") fmt("%4.2f")) count(label("Obs.") fmt("%4.0fc")))) nomtitles nonumbers  ///
			coeflabels(	dtsd_t1_am 		"First Trip (AM) " ///
						dtsd_t1_pm 		"Last Trip (PM) " ///
						dtsd_t1hw_am 	"First Home to Work Trip (AM)" ///
						dtsd_t1wh_pm 	"Last Work to Home Trip (PM)") ///
			substitute(_ \_ "<" "$<$" "\midrule" "\addlinespace\addlinespace" "\hline" "\addlinespace\addlinespace") ///
			noobs fragment //postfoot("\bottomrule" "\end{tabular}" )  // "}"


*** Write main table tex file
	file open myfile using "`outputfolder'table_1.tex", write replace
	file write myfile "\begin{tabular}{lcccccc}" _n "\toprule" _n  ///
					  "\addlinespace\multicolumn{7}{l}{\emph{Panel A. Trip Characteristics}} \\" _n ///
				   	  "\ExpandableInput{tables/`tnum'/panel_A}" _n ///
				   	  "\addlinespace\multicolumn{7}{l}{\emph{Panel B. Commute Destination Variability}} \\" _n ///
				   	  "\ExpandableInput{tables/`tnum'/panel_B}" _n ///
				   	  "\addlinespace\multicolumn{7}{l}{\emph{Panel C. Departure Time Variability}} \\ " _n ///
				   	  "\addlinespace\multicolumn{7}{l}{\emph{(Standard Deviation of the Departure Time in hours)}} \\" _n ///
				   	  "\ExpandableInput{tables/`tnum'/panel_C}" _n ///
					  "\bottomrule" _n "\end{tabular}" _n ///

	file close myfile


