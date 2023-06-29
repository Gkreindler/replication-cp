
 * This do file: Loads and codes DT charges at daily level
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


*** Load trips with charges (DT)
	import delimited using "data/coded_gps/charges_daily_all_`chain_th'.csv", clear
	
*** DT 
	keep if type == "dt"
	drop type *intersect*

*** uniform high and high
	keep if bal_target == "high" & h=="high"
	drop bal_target h

* 
	assert flag_out_am == outstation
	assert flag_out_pm == outstation
	drop flag_out*

* rename except // flag_no_data_am flag_no_data_pm flag_no_trip_am flag_no_trip_pm
	foreach va of varlist charge_am charge_maxday charge_nodata charge_notrip charge_pm charge_total charge_total_raw outstation{
		rename `va' cd_`va'
	}

	gen date_ = date(date,"YMD")
	format %td date_
	drop date
	rename date_ date
	gen dow = dow(date)

*** Special drops
	* DROP all data before 04-03 (inclusive) for L0314124709
	drop if uidp == "L0314124709" & date <= mdy(4,3,2017)

	* DROP all data before 04-06 (inclusive) for L0221092234
	drop if uidp == "L0221092234" & date <= mdy(4,3,2017)

	* DROP last day (Sept 11) due to lack of quality measures for that day
	assert date < mdy(9,12,2017)
	drop if date == mdy(9,11,2017)

*** Drop weekends
	drop if inlist(dow,0,6)

*** SAVE
	gisid uidp date
	order uidp date cd_charge_am cd_charge_pm cd_charge_total cd_charge_total_raw 

	compress
	save "data/coded_gps_dta/charges_days_dt_`chain_th'.dta", replace

