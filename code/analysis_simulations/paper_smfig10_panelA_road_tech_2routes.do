
* This do file: SM Figure 10, panel A, plot the linear road technology  with 15% and 30% steeper slope
 * This version: Oct 2020
* Author: Gabriel Kreindler
 
clear all
pause off
set more off

* Path based on global set in profile.do
	local path "$congestion_pricing_path"
	cd "`path'"

*** FIGURE number
local outputfolder "paper/figures/smfigure10/"
cap mkdir "`outputfolder'"

***************	
** Main relationship: density on instant delay
***************	


*** Hour level data
	use "data/coded_road_tech/volumes_h_level", clear

	merge 1:1 hour using "data/coded_road_tech/density_h"
	assert _m==3
	drop _m

	merge 1:1 hour using "data/coded_road_tech/google_maps_h_level"
	keep if _m==3
	drop _m

	tsset hour

	// sum volume
	// replace volume = volume / r(mean)

	*** do NOT normalize density -- already normalized
	sum density_norm, d 

*** Force higher power
	newey2 delay_gm density_norm, lag(3)
	local my_cons = _b[_cons]
	local my_slop = _b[density_norm]
	predict delay_gm_pred_lin_sd, stdp
	predict delay_gm_pred_lin   

	gen ci_low = delay_gm_pred_lin - 1.96 * delay_gm_pred_lin_sd
	gen ci_hig = delay_gm_pred_lin + 1.96 * delay_gm_pred_lin_sd

	gen delay_gm_s15 = `=`my_cons'*(1-0.15/2)' + `=`my_slop'*1.15' * density_norm
	gen delay_gm_s30 = `=`my_cons'*(1-0.30/2)' + `=`my_slop'*1.30' * density_norm

	sort density_norm
	twoway 	(scatter delay_gm  		density_norm, mlabel(hour) msize(small) mcolor(gs6)) ///
			(line delay_gm_pred_lin	density_norm, lcolor(gs8) lwidth(0.5) ) ///
			(line delay_gm_s15	density_norm, lcolor(blue) lwidth(0.5) lpattern(longdash)) ///
			(line delay_gm_s30	density_norm, lcolor(red) lwidth(0.5) ) ///
			(line ci_low density_norm, lcolor(gs8) lpattern(dash)) ///
			(line ci_hig density_norm, lcolor(gs8) lpattern(dash)) ///			
			, legend(order(2 "linear fit" 5 "linear 95% CI" ///
						   3 "15% higher slope" ///
						   4 "30% higher slope") ///
					 cols(1) position(11) ring(0) ) ///
			graphregion(color(white)) ///
			xtitle("Traffic Density (relative)") ///
			ytitle("Travel Delay (minutes / km)") scale(1.2)

	local outputfile "`outputfolder'panel_A_2routes_road_tech"
	graph export "`outputfile'.eps", replace
	graph export "`outputfile'.pdf", replace
	graph export "`outputfile'.png", replace
