
 * This do file: Figure SM.5.A robustness of road tech to GPS data and percentiles
 * This version: 1 Dec 2021
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
** ANALYSIS DT LEVEL
***************	
	use "data/coded_road_tech/density_h3", clear

	merge 1:1 hour3 using "data/coded_road_tech/google_maps_h3_level"
	assert _m==3
	drop _m

	merge 1:1 hour3 using "data/coded_road_tech/volumes_h3_level", keepusing(delay_gps*)
	drop if _m==2
	drop _m

	tsset hour3
	gen hour = floor(hour3 / 3)


*** Analysis
		reg delay_gps density_norm
	predict delay_gps_predicted

foreach va of varlist delay_gps delay_gm delay_gps_p10 delay_gps_p90 {
	nl (`va' = {l0=2.0} + {l1=1000.0} * (density_norm ^ {gamma=1.0})), vce(hac nwest 9)
	predict `va'_nl
}
	
*** GRAPH
	cap drop hour_str 
	tostring hour, gen(hour_str)
	replace hour_str = "" if hour3/3 != hour
	replace hour_str = hour_str + ":00" if hour_str != ""
	replace hour_str = "." if !inlist(hour,0,7,8,9,11,15,19,21,23)
	cap drop mlabposition
	gen 	mlabposition = 5
	replace mlabposition = 11 if inlist(hour, 7)
	replace mlabposition = 6 if inlist(hour, 15,9)
	replace mlabposition = 11 if inlist(hour, 19)
	replace mlabposition = 11 if inlist(hour, 21)

	// medium

	sort density_norm
	twoway 	(line 	 delay_gm_nl 		density_norm, lcolor(blue%50)  lpattern(dash)) ///
			(line 	 delay_gps_p10_nl  	density_norm, lcolor(green%50) lpattern(dash)) ///
			(line 	 delay_gps_p90_nl  	density_norm, lcolor(red%50)   lpattern(dash)) ///
			(line 	 delay_gps_nl    	density_norm, lcolor(gs6%50)   lpattern(dash)) ///
			(scatter delay_gm			density_norm, msymbol(circle)   msize(medsmall) mcolor(blue%50)) /// 
			(scatter delay_gps_p10		density_norm, msymbol(diamond)  msize(medsmall) mcolor(green%50)) /// 
			(scatter delay_gps_p90		density_norm, msymbol(triangle) msize(medsmall) mcolor(red%50)) /// 
			(scatter delay_gps			density_norm, msymbol(square)   msize(medsmall) mcolor(gs4%50) mlabel(hour_str) mlabcolor(gs4) mlabvposition(mlabposition)) ///
			, graphregion(color(white)) xtitle("Traffic Density (normalized)") ytitle("Travel Delay (minutes / km)") ///
			scale(1.2) ylabel(2/6.5) ///
			 legend(order( ///
			 	7 "GPS 90th p'tile" ///
			 	8 "GPS average" ///
			 	5 "Google Maps average" ///
			 	6 "GPS 10th p'tile" ) cols(1) size(small) position(11) ring(0))


	local outputfile "`outputfolder'panel_A_gps_gm_p10_p90_time"
	graph save "`outputfile'", replace
	graph export "`outputfile'.eps", replace
	graph export "`outputfile'.png", replace width(1800)
	graph export "`outputfile'.pdf", replace

	
