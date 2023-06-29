
* This do file: SM Figure 8, panel C, plot the non-linear road technology with gamma=1.5 and gamma=0.5
 * This version: Oct 2020
* Author: Gabriel Kreindler
 
clear all
pause off
set more off

* Path based on global set in profile.do
	local path "$congestion_pricing_path"
	cd "`path'"

*** FIGURE number
local outputfolder "paper/figures/smfigure8/"
cap mkdir "`outputfolder'"


***************	
** Main relationship: density on instant delay
***************	


*** Column 1,2,5,6: hour level data
	use "data/coded_road_tech/volumes_h_level", clear

	merge 1:1 hour using "data/coded_road_tech/density_h"
	assert _m==3
	drop _m

	merge 1:1 hour using "data/coded_road_tech/google_maps_h_level"
	keep if _m==3
	drop _m

	tsset hour

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

	forv i=5/15{
		gen delay_gm_power_`i' = `my_cons' + `my_slop' * density_norm^(`i'/10)
	}

	sum density_norm, d
	local maximum_density = r(max)

	di "Power 0.5, slope at max = `=`my_slop' * 0.5 * (`maximum_density')^(-0.5)'"
	di "Power 1.0, slope at max = `=`my_slop''"
	di "Power 1.5, slope at max = `=`my_slop' * 1.5 * (`maximum_density')^(0.5)'"

	sort density_norm
	twoway 	(scatter delay_gm  		density_norm, mlabel(hour) msize(small) mcolor(gs6)) ///
			(line delay_gm_power_5	density_norm, lcolor(blue) lwidth(0.5) lpattern(longdash)) ///
			(line delay_gm_power_10	density_norm, lcolor(gs8)  lwidth(0.5)) ///
			(line delay_gm_power_15	density_norm, lcolor(red)  lwidth(0.5) ) ///
			(line ci_low density_norm, lcolor(gs8) lpattern(dash)) ///
			(line ci_hig density_norm, lcolor(gs8) lpattern(dash)) ///			
			, legend(order(3 "linear fit" 5 "linear 95% CI" ///
						   2 "(density){superscript:{bf:0.5}}" ///
						   4 "(density){superscript:{bf:1.5}}") ///
					 position(11) ring(0) cols(1) ) ///
			graphregion(color(white)) ///
			xtitle("Traffic Density (relative)") ///
			ytitle("Travel Delay (minutes / km)") scale(1.2)

	local outputfile "`outputfolder'panel_C_density_power"
	graph export "`outputfile'.eps", replace
	graph export "`outputfile'.pdf", replace
	graph export "`outputfile'.png", replace
