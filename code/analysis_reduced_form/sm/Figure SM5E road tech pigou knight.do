
 * This do file: Figure SM.5.E robustness of road tech, two routes
 * This version: Sept 2021
 * Authors: Gabriel Kreindler
 
clear all
pause off
set more off

* Path based on global set in profile.do
	local path "$congestion_pricing_path"
	cd "`path'"

*** FIGURE number
	local fnum "smfigure5"
	local outputfolder "paper/figures/`fnum'/"
	cap mkdir "`outputfolder'"

*******************************
*** Google Travel Time data ***
*******************************

*** Load Google Maps data (collected in 2021)
	import delimited using "data\google_maps\all_unique_routes_by24h.csv", clear
	save "data\google_maps\all_unique_routes_by24h.dta", replace

	use "data\google_maps\all_unique_routes_by24h.dta", clear

	*** Convert from Boston time (query location) to Bangalore time
	rename hour hour_boston
	gen 	hour = hour_boston + 10.5
	replace hour = hour - 24 if hour >= 24
	// drop hour_boston

	* check looks sensible!
	// gcollapse (mean) duration, by(hour)
	// line duration hour

	tab result if hour_boston == 9  // 5% of routes do not match
	drop if result != "ok"

	sort route_index hour
	gisid route_index hour

*** check number of routes per commuter
	preserve
		gcollapse (nunique) route_index, by(uidp)
		tab route_index
		sum route_index
	restore


*** only keep non-dominated
	bys uidp hour: gegen duration_min = min(duration)
	gen fastest_route = duration == duration_min
	bys route_index: gegen ever_fastest_route = max(fastest_route)

	tab ever_fastest_route

	drop if ever_fastest_route == 0

*** (re-)check number of routes per commuter
	preserve
		gcollapse (nunique) route_index, by(uidp)
		tab route_index
		sum route_index
	restore

	bys uidp: gegen n_unique_routes = nunique(route_index)
	bys uidp: gen u1=_n==1
	tab n_unique_routes if u1==1
	sum n_unique_routes if u1==1, d

*** Only keep UIDPs with at least two (non-dominated) routes
	gunique uidp
	// keep if n_unique_routes > 1
	gunique uidp

*** Coding: duration in minutes
	replace duration = duration / 60

*** pick highest externality route
	tempfile fulldata
	save 	`fulldata'


	*** code slope
		* morning interval
		keep if inlist(hour, 6.5, 9.5)

		* second data point
		gen late = hour == 9.5

		gen delay = duration/(distance/1000)

		keep duration delay route_index uidp late 
		greshape wide duration delay, i(uidp route_index) j(late)

		gen slope = delay1 - delay0
		// gen slope = (duration1-duration0)
		// gen slope = (duration1-duration0)/duration0
		sum slope

	// sort uidp 
	by uidp: egen slope_max = max(slope)
	by uidp: egen slope_min = min(slope)
	by uidp: gen u1=_n==1

	gen slope_minmax_diff = slope_max - slope_min
	sum slope_minmax_diff if u1==1, d

	gen has_max_slope = abs(slope_max - slope) < 0.0001
	
	gisid uidp route_index

	keep route_index has_max_slope slope 

	tempfile has_max_slope
	save 	`has_max_slope'


*** Graph	
	use `fulldata', clear

	merge m:1 route_index using `has_max_slope'
	assert _m==3
	drop _m

	*** keep only best route
	sort uidp hour duration
	by uidp hour: keep if _n==1

	gcollapse (mean) has_max_slope slope, by(hour)

*** Intensive margin
	qui sum slope, d
	di (r(max)-r(min))/r(max)

	twoway 	(line slope hour, lcolor(blue) lwidth(0.5) ) ///
			(line has_max_slope hour, yaxis(2) lcolor(black) lpattern(dash)) ///
			, graphregion(color(white)) ///
			xtitle("Departure time") xlabel(0(2)24, grid) ///
			ytitle("Route" "slope" "(min/km" "/3hr)", color(blue) orientation(horizontal)) ///
			ytitle("Use" "route" "with" "highest" "slope" , orientation(horizontal) axis(2)) ///
			yscale(lcolor(blue)) ///
			ylabel(0.0(0.1)0.9, angle(0) labcolor(blue) tlcolor(blue)) ///
			ylabel(0(0.2)1.0, angle(0) axis(2)) ///
			legend(order(2 "Use highest slope route" ///
						 1 "Slope of optimal route") cols(1) position(5) ring(0)) ///
			scale(1.1)


	local outputfile "`outputfolder'panel_E_two_routes"
	graph export "`outputfile'.pdf", replace
	graph export "`outputfile'.png", replace

