
 * This do file: Figure SM.5.A robustness of road tech to using traffic volume instead of traffic delay
 * This version: 14 Sept 2017
 * Authors: GK
 
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


***************	
** ANALYSIS Main Figure and appendix figure (including daily variation too)
***************	


*** hour level data
	use "data/coded_road_tech/volumes_h_level", clear

	merge 1:1 hour using "data/coded_road_tech/google_maps_h_level"
	keep if _m==3
	drop _m

	tsset hour

	// gen hour = floor(hour/3)

	sum volume
	replace volume = volume / r(mean)

*** Line Fit
	reg delay_gm volume
	predict delay_gm_predicted

*** GRAPH
	// replace hour = . if hour3 != 3 * hour

	cap drop hour_str 
	tostring hour, gen(hour_str)
	replace hour_str = hour_str + ":00" if hour_str != "."
	replace hour_str = "" if !inlist(hour,0,7,8,9,11,15,19,21,23)
	cap drop mlabposition
	gen 	mlabposition = 11
	// replace mlabposition = 11 if inlist(hour, 7)
	// replace mlabposition = 11 if inlist(hour, 9)
	// replace mlabposition = 11 if inlist(hour, 19)
	// replace mlabposition = 11 if inlist(hour, 21)
	// replace mlabposition = 11 if inlist(hour, 23)

	// (line delay_live_south vol_smooth, lwidth(thin) lcolor(gs6)) ///
	sort volume
	twoway 	(line    delay_gm_predicted volume, /* lwidth(medium) */ lcolor(gs10) lpattern(dash)) ///
			(scatter delay_gm 			volume, mlabel(hour_str) mcolor(gs4) mlabcolor(gs4) msize(small) ///
												 mlabvposition(mlabposition)) ///
			, graphregion(color(white)) xtitle("Traffic Volume (normalized)") ///
			ytitle("Travel Delay (minutes / km)") ///
			legend(off) scale(1.2) 

	local outputfile "`outputfolder'panel_B_volume_h"
	graph save "`outputfile'", replace
	graph export "`outputfile'.eps", replace
	graph export "`outputfile'.png", replace width(1800)
	graph export "`outputfile'.pdf", replace
