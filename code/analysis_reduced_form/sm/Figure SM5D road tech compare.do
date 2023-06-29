
 * This do file: Figure SM.5.D compare estimation results with Akbar and Duranton (2017)
 * This version: 15 Jan 2018
 * Authors: Gabriel Kreindler
 
clear all
pause off
set more off
version 17

* Path based on global set in profile.do
	local path "$congestion_pricing_path"
	cd "`path'"

*** FIGURE number
	local fnum "smfigure5"
	local outputfolder "paper/figures/`fnum'/"
	cap mkdir "`outputfolder'"

*** Load Duranton
	import delimited using "data/raw_other/Akbar and Duranton (2017)/f4_pc_data.csv", clear 
	rename v1 log_volume
	rename v2 log_time
	sort log_volume

	reg log_time log_volume
	line log_time log_volume

	*** Normalize
	sum log_volume
	replace log_volume = log_volume - r(min)
	sum log_time
	replace log_time = log_time - r(min)

	rename log_time log_delay_pred

	gen study = "Akbar_Duranton"

	tempfile ad
	save 	`ad'


***************	
** ANALYSIS DT LEVEL
***************	
	use "data/coded_road_tech/volumes_h3_level", clear

	merge 1:1 hour3 using "data/coded_road_tech/google_maps_h3_level"
	assert _m==3
	drop _m

	sum volume
	replace volume = volume / r(mean)

	*** logs
	gen log_volume = log(volume)
	gen log_delay = log(delay_gm)

	gen hour = floor(hour3/3)

*** Normalize
	sum log_volume
	replace log_volume = log_volume - r(min)
	sum log_delay
	replace log_delay = log_delay - r(min)

	lpoly log_delay log_volume, ci generate(log_delay_pred) at(log_volume)

	gen study = "GK"

	append using `ad'

*** Graph
	sort log_volume
	twoway 	(line log_delay_pred	log_volume if study == "Akbar_Duranton", lcolor(red) lwidth(medthick)) ///
			(scatter log_delay 		log_volume if study == "GK", mcolor(gs8) msize(medium)) ///
			(line log_delay_pred 	log_volume if study == "GK", lcolor(blue) lpattern(dash) lwidth(medthick)) ///
			, legend(order( 3 "Bangalore" 1 "Bogota")  pos(6) cols(2)) ///
			xtitle("Log Traffic Volume (Normalized, min=0)") ytitle("Log Travel Delay (min=0)") ///
			graphregion(color(white)) scale(1.2) xlabel(0/6)

	local outputfile "`outputfolder'panel_D_roadtech_comparison"
	graph save "`outputfile'", replace
	graph export "`outputfile'.eps", replace
	graph export "`outputfile'.pdf", replace
	graph export "`outputfile'.png", replace width(1800)


*** Maximum slopes
	// gen logv = _n / 100 if _n < 500
	// ipolate log_delay log_volume if study == "GK" , generate() [epolate]
	// sort logv
	sort study log_volume
	by study: gen d_log_delay = (log_delay_pred - log_delay_pred[_n-1]) / (log_volume - log_volume[_n-1])

	twoway (line d_log_delay log_volume if study == "GK") (line log_delay_pred log_volume if study == "GK")

	lpoly d_log_delay log_volume if study == "GK", gen(dlog_gk) at(log_volume)
	lpoly d_log_delay log_volume if study == "Akbar_Duranton", gen(dlog_ad) at(log_volume)

	sum dlog_ad dlog_gk
