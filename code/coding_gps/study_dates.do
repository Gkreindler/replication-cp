
 * This do file: Dates in study, from daily logs
 * This version: 8 August 2017
 * Authors: Gabriel Kreindler

clear all
pause on
set more off

* Path based on global set in profile.do
	local path "$congestion_pricing_path"
	cd "`path'"


*** Load the type of cycles and days (panel) for each treatment cell
	import excel using "data/treatment/time_in_exp_wtrial.xlsx", clear firstrow
	tempfile time_in_exp_wtrial
	save 	`time_in_exp_wtrial'

********************************************************
*** Respondent status today: ongoing/notyet/ended
********************************************************
	import delimited using "data/treatment/respondent_stages noPII.csv", clear varnames(1)
	keep uidp today_stage treat_cell

	replace today_stage = "ended" if inlist(uidp, "L0425164649_732B", "L0425194311_634C")

	tempfile resp_stages
	save 	`resp_stages'
	
********************************************************
*** Load all daily logs. This file keeps track of communication with participants during the experiment
********************************************************
	import delimited using "data/treatment/full_daily_logs.csv", clear varnames(1)

*** merge respondent status and throw an error is there are any "not yet" respondents
	merge m:1 uidp using `resp_stages'
	assert _m==3
	drop _m

	gunique uidp

	assert today_stage == "ended"
	// gen study_ongoing = 0
	drop today_stage
	gunique uidp

********************************************************
*** SAMPLE 
********************************************************
	* Keep only reports (not reminders, welcome and summary messages)
	tab type, m
	assert  inlist(type, "report", "special", "reminder", "summary", "welcome")
	keep if inlist(type, "report", "special")

********************************************************
*** CODING
********************************************************

	* reference date (usualy the sending date (timestamp) is one day later for reports, one day earlier for welcome and other)
	gen date = date("2017-" + substr(reference_date,2,5), "YMD")
	format date %td
	gen dow = dow(date)
	assert inrange(dow,1,5)

	* checks
	gisid uidp date
	sort uidp date

	* outstation dummy
	replace outstation = proper(outstation)
	replace outstation = "1" if outstation == "True"
	replace outstation = "0" if outstation == "False"
	assert inlist(outstation, "1", "0")
	destring outstation, replace
	assert inlist(outstation, 0, 1)

	gunique uidp

/*** 
There are 5 special dates in the data set:
- 3,29,2017 with OUTSTATION==1 and type == "special"
- 4,26,2017 with OUTSTATION==0 and type == "report"  <-- this is a problem because we lost a day from everyone's cycle
- 6,26,2017 with OUTSTATION==1 and type == "report" 
- 8,15,2017 with OUTSTATION==1 and type == "report"
- 8,25,2017 with OUTSTATION==1 and type == "report"

***/
	* 3-29-2017
	* 6-26-2017 - replace with outstation == 0 unless padded by outstation
	* 8-15-2017 - replace with outstation == 0 unless padded by outstation
	* 8-25-2017 - replace with outstation == 0 unless padded by outstation

	assert date == mdy(3,29,2017) if type == "special"

forv i=1/4{
	local month: word `i' of 3  6  8  8
	local day:   word `i' of 29 26 15 25
	local type_kind: word `i' of "special" "report" "report" "report"
	assert type == "`type_kind'" 		if date==mdy(`month',`day',2017)
	assert 	outstation == 1 		if date==mdy(`month',`day',2017)
	replace outstation = 0 			if date==mdy(`month',`day',2017)
	by uidp (date): replace outstation = 1 ///
									if date==mdy(`month',`day',2017) & outstation[_n-1] == 1 & outstation[_n+1] == 1
	replace type = "special" 		if date==mdy(`month',`day',2017)
	replace cycle_type = "special" 	if date==mdy(`month',`day',2017)
}

	* 4-26-2017
	* we will deal with this later.
	
	tab date if cycle_type == "special"
	tab cycle_type type, m
	gen special = cycle_type == "special"

**** there are 6 types of cycles OR special
	assert 	inlist(cycle_type,"control_qual","control_info","dt","area", "trial_a", "trial_dt") | ///
			(type == "special" & inlist(date, mdy(3,29,2017), mdy(6,26,2017), mdy(8,15,2017), mdy(8,25,2017)) & cycle_type == "special")
	
*** TIME
	* Generate indicator progress in the experiment timeline (range 1-20 or 1-25)
	* This excludes trial periods, outstation days, and special days (4-26 is included for now)
	gen time_in_exp = (study_cycle - 1)*5 + study_day
	gen in_experiment = outstation == 0 & inlist(cycle_type, "area", "dt", "control_qual", "control_info")
	assert inlist(cycle_type, "trial_a", "trial_dt", "special") | outstation == 1 if in_experiment == 0

	* 
	bys uidp: egen ndays_in_exp    = sum(in_experiment)
	bys uidp: egen max_time_in_exp = max(in_experiment * time_in_exp)


	* Assert all cycles are exactly 20 or 25 days long
	assert inlist(max_time_in_exp,20,25)


*** First and last days in study (includes outstation) and in experiment
	by uidp (date): gen date_first = _n== 1
	by uidp (date): gen date_last  = _n==_N

	*
	by uidp: egen date_study_min = min(date)
	by uidp: egen date_study_max = max(date)
	gen n_calendar_days_in_study = date_study_max - date_study_min

	* first and last dates in the ***experiment***
	cap drop temp
	by uidp: egen temp = min(date) if in_experiment == 1
	by uidp: egen date_exp_min = mean(temp)
	cap drop temp
	by uidp: egen temp = max(date) if in_experiment == 1
	by uidp: egen date_exp_max = mean(temp)
	cap drop temp

	format %td date_study_* date_exp_* 

*** DATE checks
	* check that in_experiment data is later than the first day in study only if outstation == 1 or trial
	assert (date_study_min == date_exp_min) | (outstation == 1)  | inlist(cycle_type, "trial_a", "trial_dt", "special") if date_first==1

*** OUTSTATION spells: identify runs of consecutive days outstation
	gen outstation_spell = 0
	by uidp (date):	replace outstation_spell = outstation if _n == 1
	by uidp (date): replace outstation_spell = outstation_spell[_n - 1] * outstation + outstation if _n > 1
	by uidp: egen max_outstation_spell = max(outstation_spell)
	tab max_outstation_spell 

	* only two respondents with more than 11 consecutive outstation days 
	* L0419153703 - 13 days, ongoing (Aug 8)
	* L0411175102 - 21 days, ended
	tab uidp if max_outstation_spell > 11 & date_first == 1
	assert inlist(uidp, "L0411175102", "L0419153703") if max_outstation_spell > 11 & date_first == 1

*** Create Unique ID
	* adjust for trial periods
	gen 	study_cycle_unique = study_cycle
	replace study_cycle_unique = study_cycle - 0.5 if inlist(cycle_type,"trial_a", "trial_dt")

	gisid uidp study_cycle_unique special study_day outstation_spell

	** define alternative time_in_exp including trial periods
	merge m:1 treat_cell study_cycle_unique study_day using `time_in_exp_wtrial'
	assert _m==3
	drop _m

	assert cycle_type_check == cycle_type if outstation == 0 & special == 0
	drop cycle_type_check

	* check that time_in_exp_wtrial uniquely identifies observations
	preserve
		keep if outstation == 0 & special == 0
		gisid uidp time_in_exp_wtrial
		sort uidp date
		by uidp: egen max_time_in_exp_wtrial = max(time_in_exp_wtrial)
		by uidp:  gen ndays_in_exp_wtrial = _N
		// assert max_time_in_exp_wtrial == ndays_in_exp_wtrial
		* FAILS because initially 14 respondents had 2-day, 4-day, or no trial periods sometimes
		gunique uidp if max_time_in_exp_wtrial != ndays_in_exp_wtrial
	restore

*** Check that ALL days in the experiment between 1 and 20 or 25 are included
	assert ndays_in_exp == max_time_in_exp

* Histogram of reports per day
	// hist date , d freq

*** One date that was included in the timeline but it was a holiday - should be dropped
	* April 26 was a holiday
	gen includes_apr26 = date_exp_max >= mdy(4, 26, 2017) & mdy(4, 26, 2017) >= date_exp_min

*** Final coding
	gen date_type = ""
	replace date_type = "outstation" 	if outstation == 1 & special == 0
	replace date_type = "special" 		if special == 1
	replace date_type = "trial" 		if outstation == 0 & special == 0 & inlist(cycle_type, "trial_a", "trial_dt")
	replace date_type = "experiment" 	if outstation == 0 & special == 0 & in_experiment == 1
	assert date_type != ""

********************************************************

*** Final cleaning
	sort uidp date
	order uidp treat_cell includes_apr26 reference_date date_type date dow type special ///  
			cycle_type study_cycle study_cycle_unique study_day time_in_exp in_experiment ///
			outstation outstation_spell max_outstation_spell

********************************
*** All reports
********************************
	compress
	save "data/coded_gps_dta/study_dates_full.dta", replace

********************************
*** except outstation and special
********************************
	drop if inlist(date_type, "outstation", "special")
	assert type == "report"
	drop notes reference_date timestamp special date_first date_last

	sort uidp date
	by uidp (date): gen time_in_exp_wtrial_solid = _n
	by uidp: egen max_time_in_exp_wtrial = max(time_in_exp_wtrial)
	by uidp: egen max_time_in_exp_wtrial_solid = max(time_in_exp_wtrial_solid)
	by uidp:  gen ndays_in_exp_wtrial = _N
	assert max_time_in_exp_wtrial_solid == ndays_in_exp_wtrial
	assert inlist(max_time_in_exp_wtrial - ndays_in_exp_wtrial,0,1,2,3) // not super high discrepancies

* save
	compress
	save "data/coded_gps_dta/study_dates_wtrial.dta", replace	

********************************
*** save only experiment data (drop trial)
********************************
	keep if date_type == "experiment"
	assert type == "report"

	drop time_in_exp_wtrial_solid max_time_in_exp_wtrial* ndays_in_exp_wtrial

	* check that ALL days in the experiment between 1 and 20 or 25 are included
	by uidp: assert _N == max_time_in_exp[1]

	* all COMPLETED cycles are 20 unless INFO or CONTROL, in which case they are 25
	assert max_time_in_exp == 25 if  inlist(treat_cell, "CONTROL", "INFO")
	assert (max_time_in_exp == 20) if !inlist(treat_cell, "CONTROL", "INFO")			

	* confirm unique ID
	gisid uidp study_cycle study_day

	compress
	save "data/coded_gps_dta/study_dates.dta", replace
