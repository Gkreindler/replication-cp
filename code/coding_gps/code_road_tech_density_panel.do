
 * This do file: Codes the quantity (GPS) and speeds (GPS + Google Maps) in prep for road technology estimation
 * This version: 14 Sept 2017
 * Authors: Gabriel Kreindler
 
clear all
pause off
set more off

* Path based on global set in profile.do
	local path "$congestion_pricing_path"
	cd "`path'"



*** Load GPS trips data
	use "data/coded_gps_dta/coded_trips_complete_15.dta", clear

	*** SAMPLE 
	keep if date_wo_trips == 0 
	keep if sample_trip_ok == 1


*** Numbers cited in the text on sample size (section "Supply Estimation")
preserve

	// keep if inrange(dow, 1,5)
	keep if inrange(date, mdy(3,1,2017), mdy(9,1,2017)) 

	count
	gunique uidp
	gunique date
	gunique uidp date

restore
	
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


*** number of UIDPS by date
preserve
	assert inlist(qual_gb,0,1)
	cap drop qual_gb
	gen qual_gb = 1

	gen week = week(date)

	gcollapse (mean) qual_gb, by(uidp week)
	gcollapse (sum) n_uidps=qual_gb, by(week)

	sum n_uidps
	replace n_uidps = n_uidps/r(mean)

	* mean 1 by construction
	sum n_uidps

	tempfile number_users_t
	save	`number_users_t'
restore


*** Road DENSITY at the date-hour level

	* SAMPLE
	// keep if inrange(dow,1,5)

	scalar scalar_n_trips = _N

	gen m_start = round(t_start*60)
	gen m_stop = round(t_stop*60)

	gen trip_duration = t_stop - t_start
	sum trip_duration
	scalar scalar_trip_duration = r(mean)

	keep date m_start m_stop

*** generate density at the daily level
	forv i=0/1439{
		di "." _continue
		gen d_`i' = m_start <= `i' & `i' <= m_stop
	}

	gen n_trips=1
	gcollapse (sum) d_* n_trips, by(date) fast
	gisid date

	greshape long d_, i(date n_trips) j(mmd)
	rename d_ density
	gisid date mmd

*** merge information on number of active users at the weekly level
	gen week = week(date)

	merge m:1 week using `number_users_t'
	assert _m==3
	drop _m

	* should be close to 1
	sum n_uidps
	//     Variable |        Obs        Mean    Std. dev.       Min        Max
	// -------------+---------------------------------------------------------
	//      n_uidps |    312,480     1.02771    .5333463   .1363978   1.947012

	drop week

*** Normalization
	gen n_trips_norm = n_trips / n_uidps
	gen density_norm = density / n_uidps

	egen mean_n_trips_norm = mean(n_trips_norm)

	replace density_norm = density_norm * 24 / 0.5 / mean_n_trips_norm

	gen hour = floor(mmd/60)
	assert inrange(hour, 0, 23)

	gen hour3 = floor(mmd / 60 * 3)
	assert inrange(hour3, 0, 71)

	keep if inrange(date, mdy(3,1,2017), mdy(9,1,2017)) 

	save "data/coded_road_tech/density_tm", replace


*** Analysis files are other levels
	// use "data/coded_road_tech/density_tm", clear

preserve
	// gen hour = floor(mmd/60)
	// assert inrange(hour, 0, 23)
	gcollapse (mean) density_norm, by(date hour)
	save "data/coded_road_tech/density_th", replace	
restore 

preserve
	// gen hour3 = floor(mmd / 60 * 3)
	// assert inrange(hour3, 0, 71)
	gcollapse (mean) density_norm, by(date hour3)
	save "data/coded_road_tech/density_th3", replace	
restore 

preserve 
	gcollapse (mean) density_norm, by(date)
	save "data/coded_road_tech/density_t", replace	
restore
