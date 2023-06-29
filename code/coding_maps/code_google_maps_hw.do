
 * This do file: Code Google Maps Data (Home-Work) - used for structural analysis
 * This version: 13 Aug 2018
 * Authors: Gabriel Kreindler
 
clear all
pause off
set more off

* Path based on global set in profile.do
	local path "$congestion_pricing_path"
	cd "`path'"
		
*******************************
*** Google Travel Time data ***
*******************************

*** query definitions
	import delimited using "data/google_maps/gmaps_hw_query.csv", clear
	rename id odid
	isid name
	gen uidp =  trim(substr(name,1,strpos(name," ")))
	
	gen 	dir = .
	replace dir = 0 if origin == "home" & destination == "work"
	replace dir = 1 if origin == "home" & destination == "work1"
	replace dir = 2 if origin == "home" & destination == "work2"

	replace dir = 10 if origin == "work"  & destination == "home"
	replace dir = 11 if origin == "work1" & destination == "home"
	replace dir = 12 if origin == "work2" & destination == "home"
	assert dir!=.

	label define direction 0 "home-work" 1 "home-work1" 2 "home-work2" 10 "work-home" 11 "work1-home" 12 "work2-home"
	label values dir direction

	isid uidp dir
	keep uidp dir odid

	save "data\google_maps\gmaps_hw_query.dta", replace

*** data
	import delimited using "data/google_maps/gmaps_hw_predicted.csv", clear
	assert v7 == " "
	drop v7

	merge m:1 odid using "data\google_maps\gmaps_hw_query.dta"
	assert _m==3
	drop _m

	isid odid dep_time_local
	sort odid dep_time_local

* clean
	drop dur_bg querytime_utc
	destring dur_traffic_bg, replace

* decode time
	gen double dt_local = clock(dep_time_local,"YMD hms")
	format dt_local %tcCCYY-NN-DD_HH:MM:SS

	* fractional hour
	gen double hf = hh(dt_local) + mm(dt_local) / 60 + ss(dt_local) / 3600
	* integer from 1 to 144 (multiple of 10 minutes)
	gen time = hf * 6 
	assert time == floor(time)

	mdesc dur_traffic_bg

*** Check Km - looks great
	// bys uidpn: gen o1=_n==1
	// replace distance_bg = distance_bg / 1000
	// reg mean_km distance_bg if o1== 1 & distance_bg < 30
	// scatter mean_km distance_bg if o1==1 & distance_bg < 30

*** more coding
	keep uidp dir time dur_traffic_bg distance_bg
	rename dur_traffic_bg dur_gm_
	rename distance_bg dist
	replace dist = dist / 1000

	sort uidp dir time
	by uidp dir: egen dist_mean = mean(dist)
	by uidp dir: egen dist_mode = mode(dist)
	count if dist_mode == .
	assert r(N) == 144
	replace dist_mode = dist_mean if dist_mode == .

*** interpolate missing values
	count if dur == .
	list uidp dir if dur == ., nol
	by uidp dir: ipolate dur time, gen(dur_gm_full)

	// br time dur_gm_ dur_full if uidp == "L0509084542_17D2" & dir == 0
	replace dur_gm_ = dur_gm_full
	drop dur_gm_full


						// gen delay = dur_gm_ / 60 / dist_mode
						// sum delay if inrange(time,`=7*6', `=22*6'), d
						// // kdensity delay
						// reg delay dist_mode
						// sdfsdf

*** Reshape into wide
	drop dist
	reshape wide dur_gm_, i(uidp dir) j(time)
	isid uidp dir

*** Save for dir = 0 or 1

* export hom - work or home - work1
	gunique uidp if inlist(dir, 0,1)
	keep if inlist(dir,0,1)
	drop dir
	isid uidp
	order uidp dist_mean dist_mode
	save "data\google_maps\gmaps_hw_predicted.dta", replace
// */
