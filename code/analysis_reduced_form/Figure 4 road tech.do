
 * This do file: 
 * This version: 14 Sept 2017
 * Authors: GK
 
clear all
pause off
set more off

* Path based on global set in profile.do
	local path "$congestion_pricing_path"
	cd "`path'"

*** FIGURE number
local fnum "figure4"
local outputfolder "paper/figures/`fnum'/"
cap mkdir "`outputfolder'"


***************	
** ANALYSIS Main Figure and appendix figure (including daily variation too)
***************	


*** hour level data
	use "data/coded_road_tech/density_h", clear

	merge 1:1 hour using "data/coded_road_tech/google_maps_h_level"
	keep if _m==3
	drop _m

	tsset hour

*** Line Fit
	newey2 delay_gm density_norm, lag(3)
	// predict delay_gm_predicted
	predict delay_gm_pred_lin_sd, stdp
	predict delay_gm_pred_lin   

	gen ci_low = delay_gm_pred_lin - 1.96 * delay_gm_pred_lin_sd
	gen ci_hig = delay_gm_pred_lin + 1.96 * delay_gm_pred_lin_sd


	// scatter delay_gm density_norm, mlabel(hour)

	newey2 delay_gm density_norm if inrange(hour, 8,20), lag(3)
	predict delay_gm_pred_lin_daytime_sd if inrange(hour, 8,20), stdp
	predict delay_gm_pred_lin_daytime    if inrange(hour, 8,20)

	gen ci_daytime_low = delay_gm_pred_lin_daytime - 1.96 * delay_gm_pred_lin_daytime_sd
	gen ci_daytime_hig = delay_gm_pred_lin_daytime + 1.96 * delay_gm_pred_lin_daytime_sd

*** GRAPH
	// replace hour = . if hour3 != 3 * hour

	cap drop hour_str 
	tostring hour, gen(hour_str)
	replace hour_str = hour_str + ":00" if hour_str != "."
	replace hour_str = "" if !inlist(hour,0,7,8,9,11,15,19,21,23)
	cap drop mlabposition
	gen 	mlabposition = 11
	replace mlabposition = 11 if inlist(hour, 7)
	replace mlabposition = 11 if inlist(hour, 8)
	replace mlabposition = 11 if inlist(hour, 9)
	replace mlabposition = 5 if inlist(hour, 15)
	replace mlabposition = 11 if inlist(hour, 19)
	replace mlabposition = 11 if inlist(hour, 21)
	replace mlabposition = 11 if inlist(hour, 23)


	// (line delay_live_south vol_smooth, lwidth(thin) lcolor(gs6)) ///
	sort density_norm
	twoway 	(line    delay_gm_pred_lin  density_norm, lcolor(gs8) lpattern(solid)) ///
			(scatter delay_gm 			density_norm, mlabel(hour_str) mcolor(gs4) mlabcolor(gs4) msize(small) ///
												 mlabvposition(mlabposition)) ///
			(line    ci_low density_norm, lcolor(gs10) lpattern(dash)) ///
			(line    ci_hig density_norm, lcolor(gs10) lpattern(dash)) ///
			(line    delay_gm_pred_lin_daytime  density_norm, lcolor(blue) lpattern(longdash)) ///
			(line    ci_daytime_low density_norm, lcolor(blue%50) lpattern(dash)) ///
			(line    ci_daytime_hig density_norm, lcolor(blue%50) lpattern(dash)) ///
			, graphregion(color(white)) xtitle("Traffic Density (normalized)") ytitle("Travel Delay (minutes / km)") ///
			legend(order(1 "Linear full sample" 5 "Linear 8am-8pm" ) cols(1) position(11) ring(0) ) ///
			scale(1.2) 

	local outputfile "`outputfolder'figure_`fnum'_h"
	graph save "`outputfile'", replace
	graph export "`outputfile'.eps", replace
	graph export "`outputfile'.png", replace width(1800)
	graph export "`outputfile'.pdf", replace

	gen hour_level=1

	tempfile hour3level
	save 	`hour3level'


*** GPS date level
	use "data/coded_road_tech/density_t", clear

	merge 1:1 date using "data/coded_road_tech/google_maps_t_level"
	keep if _m==3
	drop _m

	tsset date

	gen dow = dow(date)

	*** Column 4 OLS over dates
	reg delay_gm density_norm, r
	predict delay_gm_predicted

	gen hour_level = 0

	append using `hour3level'

	sort hour_level density_norm
	twoway 	(line    delay_gm_predicted density_norm if hour_level == 1, lcolor(gs10) lpattern(dash)) ///
			(scatter delay_gm 			density_norm if hour_level == 1, msymbol(square) mcolor(gs4) msize(small)) ///
			(line    delay_gm_predicted density_norm if hour_level == 0, lcolor(eltblue) lpattern(dash)) ///
			(scatter delay_gm 			density_norm if hour_level == 0, mcolor(eltblue) msize(small )) ///
			(scatter delay_gm 			density_norm if hour_level == 0 & dow == 0, mcolor(pink) msize(small )) ///
			, graphregion(color(white)) xtitle("Traffic Density (normalized)") ytitle("Travel Delay (minutes / km)") ///
			scale(1.2) xlabel(0(0.5)2) ylabel(2(0.5)4.5) ///
			legend(order(2 "Departure time (h)" 4 "Calendar Date (t)" 5 "Sundays") position(11) ring(0))

local outputfile "`outputfolder'figure_`fnum'_h_and_t"
graph save "`outputfile'", replace
graph export "`outputfile'.eps", replace
graph export "`outputfile'.png", replace width(1800)
graph export "`outputfile'.pdf", replace
