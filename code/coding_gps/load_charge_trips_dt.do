
 * This do file: Loads and codes trips data 
 * This version: 25 July, 2017
 * Authors: Gabriel Kreindler

clear all
pause off
set more off

if "`1'" == ""{
	local chain_th "15"	
}
else{
	local chain_th "`1'"		
}

* Path based on global set in profile.do
	local path "$congestion_pricing_path"
	cd "`path'"


*** Load trips with charges
	import delimited using "data/coded_gps/charges_trips_dt_`chain_th'.csv", clear
	rename id_new chain

	* start time
	assert length(start) == 5
	gen t_start_hh = substr(start,1,2)
	gen t_start_mm = substr(start,4,2)
	destring t_start*, replace
	gen t_start = t_start_hh + t_start_mm / 60
	drop t_start_?? start

	foreach va of varlist charge distance duration ramp rate t_start{
		rename `va' c_dt_`va'
	}

	gen date_ = date(date,"YMD")
	format %td date_
	drop date
	rename date_ date

*** Special drops
	* DROP all data before 04-03 (inclusive) for L0314124709
	drop if uidp == "L0314124709" & date <= mdy(4,3,2017)

	* DROP all data before 04-06 (inclusive) for L0221092234
	drop if uidp == "L0221092234" & date <= mdy(4,3,2017)

	* DROP last day (Sept 11) due to lack of quality measures for that day
	assert date < mdy(9,12,2017)
	drop if date == mdy(9,11,2017)


*** SAVE
	isid uidp date chain

	compress
	save "data/coded_gps_dta/charges_trips_dt_`chain_th'.dta", replace


