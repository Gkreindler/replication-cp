
 * This do file: Codes the speeds (GPS + Google Maps) in prep for road technology estimation
 * This version: 14 Sept 2017
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


*** Load Google Maps travel delay data
	* code at t-h (date-hour) level
	use "data/google_maps/gmaps_citywide_live_coded.dta", clear
	merge m:1 odid using "data/google_maps/gmaps_citywide_query_btm.dta"
	assert _m!=2
	keep if _m==3
	drop _m


	gen hour3 = round(hff_clean * 3)
	assert inrange(hour3,0,71)

	cap rename distance_bg dist_bg

	*** date-hour level
	preserve
		gcollapse (mean) dur dist=dist_bg, by(date hour)
		gen delay_gm = (dur / 60) / (dist / 1000)
		save "data/coded_road_tech/google_maps_th_level", replace
	restore

	*** date-hour level
	preserve
		gcollapse (mean) dur dist=dist_bg, by(date hour3)
		gen delay_gm = (dur / 60) / (dist / 1000)
		save "data/coded_road_tech/google_maps_th3_level", replace
	restore

	*** date level
	preserve
		gcollapse (mean) dur dist=dist_bg, by(date)
		gen delay_gm = (dur / 60) / (dist / 1000)
		save "data/coded_road_tech/google_maps_t_level", replace
	restore

	*** hour level (weekdays only)
	preserve
		keep if inrange(dow,1,5)

		gcollapse (mean) dur dist=dist_bg, by(hour)
		gen delay_gm = (dur / 60) / (dist / 1000)
		save "data/coded_road_tech/google_maps_h_level", replace
		export delimited "data/coded_road_tech/google_maps_h_level.csv", replace
	restore

	preserve
		keep if inrange(dow,1,5)

		gcollapse (mean) dur dist=dist_bg, by(hour3)
		gen delay_gm = (dur / 60) / (dist / 1000)
		save "data/coded_road_tech/google_maps_h3_level", replace
		export delimited "data/coded_road_tech/google_maps_h3_level.csv", replace
	restore



*** Load GPS trips data
	use "data/coded_gps_dta/coded_trips_complete_`chain_th'.dta", clear

	*** SAMPLE 
	keep if date_wo_trips == 0 & sample_trip_ok == 1
	drop if plb < 2

	count
	gunique uidp
	gunique uidp date
	gunique date

	keep if n_struct_loc==0 & inrange(diameter_ratio,0.6,1.0)

	count
	gunique uidp
	gunique uidp date
	gunique date

*** Coding 
	gen hour = floor(t_start)
	gen hour3 = floor(t_start * 3)
	assert inrange(hour3, 0, 71)

*** Travel Delay Measure from GPS data
	gen delay =  dur_bg_mm / plb
	winsor delay, gen(delay_w) p(0.01)


*** number of UIDPS by date
preserve
	assert inlist(qual_gb,0,1)
	cap drop qual_gb
	gen qual_gb = 1

	gen week = week(date)

	gcollapse (mean) qual_gb, by(uidp week)
	gcollapse (sum) n_uidps=qual_gb, by(week)

	tempfile number_users_t
	save	`number_users_t'
restore

*** date-hour level
preserve
	gen volume = 1
	gcollapse (sum) volume (mean) delay_gps=delay_w, by(date hour)

	gen week = week(date)

	merge m:1 week using `number_users_t'
	assert _m==3
	drop _m

	drop week

	gen volume_pc = volume / n_uidps
	sum volume_pc
	replace volume_pc = volume_pc / r(mean)

	keep if inrange(date, mdy(3,1,2017), mdy(9,1,2017))  // ???????????????????????????????????

	save "data/coded_road_tech/volumes_th_level", replace
restore

*** date-hour level
preserve
	gen volume = 1
	gcollapse (sum) volume (mean) delay_gps=delay_w, by(date hour3)

	gen week = week(date)

	merge m:1 week using `number_users_t'
	assert _m==3
	drop _m

	drop week

	gen volume_pc = volume / n_uidps
	sum volume_pc
	replace volume_pc = volume_pc / r(mean)

	keep if inrange(date, mdy(3,1,2017), mdy(9,1,2017))  // ???????????????????????????????????

	save "data/coded_road_tech/volumes_th3_level", replace
restore

*** date level
preserve
	
	// bys date: gen d1 = _n==1
	// bys date: gen n  = _N
	// sum n if d1==1, d

	// keep if n >= r(p50)
	// sum n if d1==1, d
	// count

	count
	gunique uidp
	gunique uidp date
	gunique date

	gen volume = 1
	gcollapse (sum) volume (mean) delay_gps=delay_w, by(date)

	gen week = week(date)

	merge m:1 week using `number_users_t'
	assert _m==3
	drop _m

	drop week

	gen volume_pc = volume / n_uidps
	sum volume_pc
	replace volume_pc = volume_pc / r(mean)

	keep if inrange(date, mdy(3,1,2017), mdy(9,1,2017))

	save "data/coded_road_tech/volumes_t_level", replace
restore

*** hour level -- only weekdays
preserve
	keep if inrange(dow,1,5)

	gen volume = 1
	gcollapse (sum) volume ///
			(p10) delay_gps_p10=delay_w ///
			(p50) delay_gps_p50=delay_w ///
			(p90) delay_gps_p90=delay_w ///
			(mean) delay_gps=delay_w, by(hour)

	save "data/coded_road_tech/volumes_h_level", replace
restore

preserve
	keep if inrange(dow,1,5)

	gen volume = 1
	gcollapse (sum) volume ///
			(p10) delay_gps_p10=delay_w ///
			(p50) delay_gps_p50=delay_w ///
			(p90) delay_gps_p90=delay_w ///
			(mean) delay_gps=delay_w, by(hour3)

	save "data/coded_road_tech/volumes_h3_level", replace
restore


