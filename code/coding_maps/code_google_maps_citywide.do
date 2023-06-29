
 * This do file: Code Google Maps Data
 * This version: 13 Aug 2018
 * Authors: Tammy Tseng, Gabriel Kreindler
 
clear all
pause off
set more off

* Path based on global set in profile.do
	local path "$congestion_pricing_path"
	cd "`path'"
		
***************	
** LOAD DATA **
***************	
	
* Load live data
	import delimited using "data/google_maps/gmaps_citywide_live.csv", clear
	drop mode v8 dep_time_utc
	
* format time
	gen double dt_local = clock(dep_time_local,"YMD hms")
	format dt_local %tcCCYY-NN-DD_HH:MM:SS
	
* define clean departure time that merges easily across dates
	gen double dt_local_clean = round(dt_local, 20 * 60 * 1000)
	format dt_local_clean %tcCCYY-NN-DD_HH:MM:SS
	drop dep_time_local

* clean and generate nice clean vars - LOCAL
	gen date = dofc(dt_local_clean)
	gen day = day(date)
	gen month = month(date)
	
	gen hour = hh(dt_local_clean)
	gen min  = mm(dt_local_clean)
	gen hff_clean = hour + min/60

	* check we are aleady rounding to closest 20 minutes
	gen hff20 = round(hff_clean * 3) / 3
	assert hff20 == hff_clean
	drop hff20
	
	gen dow = dow(date)
	
	format date %td
		
* gen delay
	gen delay = (dur/60) / (dist/1000)
	label var delay "Delay: minutes per KM"

	bys odid: gegen dist_mod = mode(dist)
	gen delay_mode = (dur/60) / (dist_mod/1000)
	
* clean up labels
	order odid dist_bg dur dt_local month day hour min
	label var odid "Edge ID"
	label var dist_bg "Distance of best guess shortest path"
	label var dur "Duration in minutes"
	label var dur_bg "Duration in minutes of best guess shortest path"
	label var dt_local "Date and time of query (local)"
	label var date "Date of query (local)"
	label var hff_clean "Clean time of query (rounded to nearest 20 min)"
	label var dt_local_clean "Clean date and time of query (rounded to nearest 20 min)"
	label var dow "Day of week (0 = Sunday, 6 = Saturday)"

	* drop duplicates
	gduplicates drop odid date hff_clean, force
	gisid odid date hff_clean
	
* Save
	save "data/google_maps/gmaps_citywide_live_coded.dta", replace


***************************
*** BTM Area Edges only ***
***************************
	*** load origin and destination
	import delimited "data\google_maps\gmaps_citywide_query.csv", clear
	rename id odid
	keep origin destination odid
	replace origin = trim(origin)
	replace destination = trim(destination)
	gen l1 = substr(origin,1,1)
	assert l1=="p"
	drop l1
	gen l1 = substr(destination,1,1)
	assert l1=="p"
	drop l1
	gen l1 = substr(destination,length(destination),1)
	assert inlist(l1,"a", "b")
	drop l1
	gen l1 = substr(origin,length(origin),1)
	assert inlist(l1,"a", "b")
	drop l1

	gen orig = substr(origin,2,length(origin)-2)
	gen dest = substr(destination,2,length(destination)-2)

	destring orig dest, replace

	// gen btm = inlist(orig, 27, 28, 29, 30, 31, 32, 33, 34, 35, 36, 37, 38, 39, 40, 41, 42, 43, 44, 45, 46, 97, 101, 160, 161) | ///
	// 		  inlist(dest, 27, 28, 29, 30, 31, 32, 33, 34, 35, 36, 37, 38, 39, 40, 41, 42, 43, 44, 45, 46, 97, 101, 160, 161)

	// gen btm = inlist(orig, 27, 28, 29, 30, 31, 32, 33, 34, 35, 36, 37, 38, 39, 40, 41, 42, 43, 44, 45, 46, 97, 101, 160, 161) | ///
	// 		  inlist(orig, 1,2,3,4,12,13,14,15,20,21,26,27,28,29,96,97,106,107,108,109,   92,93,102,103,110,111) | ///
	// 		  inlist(dest, 27, 28, 29, 30, 31, 32, 33, 34, 35, 36, 37, 38, 39, 40, 41, 42, 43, 44, 45, 46, 97, 101, 160, 161) | ///
	// 		  inlist(dest, 1,2,3,4,12,13,14,15,20,21,26,27,28,29,96,97,106,107,108,109,   92,93,102,103,110,111)

	gen btm = inlist(orig, 27, 28, 29, 30, 31, 32, 33, 34, 35, 36, 37, 38, 39, 40, 41, 42, 43, 44, 45, 46, 97, 101, 160, 161) | ///
			  inlist(dest, 27, 28, 29, 30, 31, 32, 33, 34, 35, 36, 37, 38, 39, 40, 41, 42, 43, 44, 45, 46, 97, 101, 160, 161)

	replace btm = 1 if inlist(odid, 1,2,3,4,12,13,14,15,20,21,26,27,28,29,96,97,106,107,108,109,164,165,   92,93,102,103,110,111, 56, 57, 62, 63)

	keep if btm == 1
	keep odid btm
	
	save "data\google_maps\gmaps_citywide_query_btm.dta", replace
