
 * This do file: Implied beta/alpha ratios
 * Authors:Gabriel Kreindler
 
clear all
pause off
set more off 

* Path based on global set in profile.do
	local path "$congestion_pricing_path"
	cd "`path'"
		
***************	
** LOAD DATA **
***************	
	use "data/google_maps/gmaps_citywide_live_coded.dta", clear

	merge m:1 odid using "data\google_maps\gmaps_citywide_query_btm.dta"
	assert _m!=2
	drop if _m==1
	drop _m

************
** SAMPLE **
************

*** Sample
	* drop if weekend, drop if 3/29 or 5/1 (holidays)
	drop if dow==0 | dow==6
	drop if inlist(date, mdy(3,29,2017), mdy(5,1,2017), mdy(9,15,2017))
	
	// large outlier: Varasiddhi Vinayaka Vratha -- Aug 25th 2017
	drop if inlist(date, `=mdy(8,25,2017)')

	tab date
	
	gcollapse (mean) delay, by(odid hff_clean dow btm)

	assert btm==1
	gunique odid
	gunique hff_clean // 72 balanced groups of 140

	* Collapse at weekday x time level
	collapse (mean) delay*, by(hff_clean) 

	sort hff_clean
	gen delay_next = delay[_n+1]
	gen slope = (delay_next - delay) / 0.3333333

	sum slope, d // max = 0.89 min/km/hour
	local max_slope = r(max)

*** 
	use "data\google_maps\gmaps_hw_predicted.dta", clear
	sum dist_mean, d
	local dist_p50 = r(p50)
	local dist_p95 = r(p95)


*** Numbers cited in the text:
	
	* \beta_early / alpha for the median home-work distance (9.8 km)
		* max_slope is in minutes per kilometer per hour, so divide by 60 to get hour/km/hour
		* dist_p50 is in kilometers
	di (`max_slope' / 60) * `dist_p50' // .14546616

	* \beta_early / alpha for the p95 home-work distance (27.8 km)
	di (`max_slope' / 60) * `dist_p95' // .41145334

	di 552  / 609






