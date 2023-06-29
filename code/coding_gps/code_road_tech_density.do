
 * This do file: Codes the quantity (GPS) and speeds (GPS + Google Maps) in prep for road technology estimation
 * This version: 14 Sept 2017
 * Authors: Gabriel Kreindler
 
clear all
pause off
set more off

* Path based on global set in profile.do
	local path "$congestion_pricing_path"
	cd "`path'"


*** TABLE number
// local tnum "T_RT"
// local outputfolder_table "$cp_git_path/paper/tables/`tnum'/"
// cap mkdir "`outputfolder_table'"


*** Load GPS trips data
	use "data/coded_gps_dta/coded_trips_complete_15.dta", clear

	*** SAMPLE 
	keep if date_wo_trips == 0 
	keep if sample_trip_ok == 1
	
	*** Keep all high-quality trips, no restrictions
	*** Dropping these leads to a more _concave_ density - instant speed curve
	// drop if plb < 2
	// keep if n_struct_loc==0 
	// keep if inrange(diameter_ratio,0.6,1.0)

*** Coding 
	gen hour = floor(t_start)
	gen hour3 = floor(t_start * 3)
	assert inrange(hour3, 0, 71)

*** Travel Delay Measure from GPS data
	gen delay =  dur_bg_mm / plb
	winsor delay, gen(delay_w) p(0.01)


*** Road DENSITY
	// use "data/coded_trips_dta/coded_trips_complete_15.dta", clear

	* SAMPLE
	keep if inrange(dow,1,5)

	scalar scalar_n_trips = _N

	gen m_start = round(t_start*60)
	gen m_stop = round(t_stop*60)

	gen trip_duration = t_stop - t_start
	sum trip_duration
	scalar scalar_trip_duration = r(mean)

	keep m_start m_stop

	forv i=0/1439{
		di "." _continue
		gen d`i' = m_start <= `i' & `i' <= m_stop
	}

	gcollapse (sum) d*, fast

	gen one=1
	greshape long d, i(one) j(mmd)

	line d mmd
	drop one
	drop if mmd == 0	
	rename d density

	gen n_trips = scalar_n_trips
	gen mean_trip_dur = scalar_trip_duration

	sum density
	sum n_trips
	di 1328.821 * 24 / 0.5 / 61893

	gen density_norm = density * 24 / 0.5 / n_trips
	sum density_norm

	save "data/coded_road_tech/density_mdd", replace

preserve
	gen hour = floor(mmd/60)
	assert inrange(hour, 0, 23)
	gcollapse (mean) density_norm, by(hour)
	save "data/coded_road_tech/density_h", replace	
restore 

preserve
	gen hour3 = floor(mmd / 60 * 3)
	assert inrange(hour3, 0, 71)
	gcollapse (mean) density_norm, by(hour3)
	save "data/coded_road_tech/density_h3", replace	
restore 
