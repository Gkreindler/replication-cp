
 * This do file: Loads and codes AREA chages at the daily level
 * This version: 12 August, 2017
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

*** Load trips with charges (AREA)
	import delimited using "data/coded_gps/charges_daily_all_`chain_th'.csv", clear
	
*** AREA
	keep if type == "area_new"
	cap drop trips_am trips_pm
	drop bal_target charge_maxday charge_nodata charge_notrip charge_total_raw h /* time_intersect_am time_intersect_pm */ type
* 
	assert flag_out_am == outstation
	assert flag_out_pm == outstation
	drop flag_out*

* rename all
	foreach va of varlist charge_am charge_intersect_am charge_intersect_pm charge_pm charge_total ///
						  flag_no_data_am flag_no_data_pm flag_no_trip_am flag_no_trip_pm ///
						  intersect_am intersect_pm outstation{
		assert `va' != .
		rename `va' cad_`va' // Charge Area Daily (CAD)
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
	isid uidp date
	order uidp date cad_charge_am cad_charge_pm cad_charge_total cad_charge_intersect_am cad_charge_intersect_pm cad_intersect_am cad_intersect_pm cad_outstation

	compress
	save "data/coded_gps_dta/charges_days_area_`chain_th'.dta", replace


	
