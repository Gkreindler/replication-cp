
 * This do file:  Treatment assignment (uidp level)
 * This version: 11 June, 2017
 * Authors: Gabriel Kreindler

clear all
pause off
set more off

* Path based on global set in profile.do

	local path "$congestion_pricing_path"
	cd "`path'"

/* *** Regular / irregular commuters
	// Note: this data is not included in the replication package because it contains PII

	import excel using "data/raw_gps_homework/hw_locations_v3.xlsx", clear cellrange(A3) case(lower) firstrow
	keep uidp dest dest2 status_2
	isid uidp

	assert inlist(dest, "yes", "no", "dk")
	assert inlist(dest2, "yes", "", "dk")
	assert inlist(dest2, "yes", "") if dest == "yes"
	assert inlist(status_2, "PITDEST", "STUB") if dest == "dk"
	assert dest2 == "dk" if dest == "dk"

	gen 	reg_com = ""
	replace reg_com = "1 dk" 	if dest == "dk"
	replace reg_com = "2 reg"  	if dest == "yes" & dest2 == "" & reg_com == ""
	replace reg_com = "3 reg2" 	if dest == "yes" & dest2 == "yes" & reg_com == ""
	replace reg_com = "4 irr"  	if dest == "no" & reg_com == ""
	assert reg_com != ""

	encode reg_com, gen(regular_commuter)

	keep uidp regular_commuter

	save "data/raw_gps_homework/regular_commuter", replace
 */

*** Drops
	import excel using "data/raw_gps_homework/uidp_drop_list.xlsx", clear firstrow case(lower)
	keep if action == "DROP" & category == "AREA"
	count
	assert r(N) == 11
	isid uidp
	keep uidp
	tempfile area_drop_list
	save	`area_drop_list'

*** Study stage
	use "data/coded_gps_dta/study_dates.dta", clear
	keep uidp treat_cell
	duplicates drop
	tempfile uidp_sample
	save 	`uidp_sample'

*** Treatment assignment
	// import delimited using "data/treatment/TREATMENT_ROSTER.csv", clear varnames(1)
	use "data/treatment/treatment roster noPII.dta", clear
	merge m:1 idx strat_cell using "data/treatment/prerandomized_assignment.dta", keep(match master)
	assert _merge!=1
	isid uidp
	drop _merge
	
* Keep only those enlisted in the experiment
	assert inlist(meeting,"","done","dropped")
	keep if meeting == "done"

	merge 1:1 uidp using `uidp_sample'
	assert _m==3
	drop _m

	merge 1:1 uidp using `area_drop_list'
	assert _m!=2
	gen sample_area_drop = _m==3
	drop _m

	merge 1:1 uidp using "data/raw_gps_homework/regular_commuter"
	assert _m==3
	drop _m

* generate separate strata for first "set" (first 8 respondents in strata 1-4, and first 12 in strata 5-8)
	gen 	strat_cell_new = strat_cell
	replace strat_cell_new = strat_cell + 10 if set == 1
	tab strat_cell_new strat_cell, m

* coding
	* meet date and start date
	foreach va in dateadded datemeeting datestart{
		cap drop temp
		gen temp = "2017-" + substr(`va',2,5)
		drop `va'
		gen `va' = date(temp, "YMD")
		format `va' %td
	}

	save "data/coded_gps_dta/treat_assign.dta", replace
