
 * This do file: Figure SM.5.C robustness of road tech. Histogram of recruitment times
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

*** Load data 
	use "data/coded_road_tech/volumes_h3_level.dta", clear

	// 10 minute bandwidth
	cap drop vol_smooth
	lpoly volume hour3, bw(1) nosc gen(vol_smooth) at(hour3) nograph

	gen hff = hour3/3
	gen mmd = hff * 60

	keep mmd hour3 hff vol_smooth
	tempfile vol_trips
	save 	`vol_trips'

*** Sample
	use "data/coded_gps_dta/coded_trips_complete_15.dta", clear
	keep uidp
	duplicates drop
	tempfile sample 
	save 	`sample'


*** histogram
	use "data/coded_cto/recruitment coded.dta", clear

	* interview  start time (hour)
	gen start_hff = starthour + startmin / 60
	drop uidp
	rename uidp_original uidp
	merge 1:1 uidp using `sample'
	count if _m==2
	assert r(N) <= 2
	drop if _m==2
	// assert _m!=2
	keep if _m==3
	drop _m

	keep uidp start_hff
	tempfile recruit_time
	save 	`recruit_time'

	*** Compare with volume
	keep start_hff

	cap set obs `=24*60'
	gen mmd = _n
	merge 1:1 mmd using `vol_trips'

	sort hff
	twoway 	(hist start_hff, s(6) w(0.5) fcolor(gs12) lcolor(gs12)) ///
	 		(line vol_smooth hff if inrange(hff,5,22), yaxis(2) lwidth(medthick) lcolor(blue)) ///
			, graphregion(color(white)) ///
			 legend(order(1 "Recruit time" 2 "Traffic Volume") pos(6) cols(2)) ///
			xtitle("Departure Time") xlabel(5(1)22, grid) ytitle("") ytitle("", axis(2)) ylabel(none) ylabel(none, axis(2)) scale(1.2)

	local outputfile "`outputfolder'panel_C_recruitment"
	graph export "`outputfile'.eps", replace
	graph export "`outputfile'.png", replace
	graph export "`outputfile'.pdf", replace

