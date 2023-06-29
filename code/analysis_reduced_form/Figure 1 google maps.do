
 * This do file: Code Google Maps Data
 * This version: 13 Aug 2018
 * Authors: Tammy Tseng, Gabriel Kreindler
 
clear all
pause off
set more off
version 15

* Path based on global set in profile.do
	local path "$congestion_pricing_path"
	cd "`path'"

*** FIGURE number
	local fnum "figure1"
	local outputfolder "paper/figures/`fnum'/"
	cap mkdir "`outputfolder'"

*** Check secular trend
	use "data/google_maps/gmaps_citywide_live_coded.dta", clear
	drop if dow==0 | dow==6
	drop if inlist(date, mdy(3,29,2017), mdy(5,1,2017), mdy(9,15,2017))
	drop if inlist(date, `=mdy(8,25,2017)')
	keep if inrange(hff_clean,5,22)

	gcollapse (mean) delay, by(date)

	sort date

*** Load data
	use "data/google_maps/gmaps_citywide_live_coded.dta", clear

	merge m:1 odid using "data/google_maps/gmaps_citywide_query_btm.dta"
	assert _m!=2
	drop if _m==1
	drop _m

*** Sample
	* drop if weekend, drop if 3/29 or 5/1 (holidays)
	drop if dow==0 | dow==6
	drop if inlist(date, mdy(3,29,2017), mdy(5,1,2017), mdy(9,15,2017))
	
	* also drop holiday (puja) on Varasiddhi Vinayaka Vratha -- Aug 25th 2017
	drop if inlist(date, `=mdy(8,25,2017)')

	tab date

	*** Check for outliers
	// preserve
	// 	gcollapse (mean) delay, by(date)
	// 	sort date
	// 	line delay date, xline(`=mdy(8,25,2017)')
	// 	dsfds
	// restore

	gcollapse (mean) delay, by(odid hff_clean dow btm)

**** Analysis
	* stats used in the text
	// On average between 7 am and 10 pm on weekdays and across all routes,
	// it takes 3.7 minutes to advance one kilometer.
	// 60 / 3.7 km / h
	// 60 / 3.7 / 1.60934 miles / h   =   10.08
	// 1 mile = 1.60934 km
	sum delay if inrange(hff_clean,7,22) // mean: 3.707016 
	di 60/r(mean) / 1.60934

*** only the 28 routes in South Bangalore
	assert btm==1
	gunique odid
	gunique hff_clean // 72 balanced groups of 140

	* Collapse at weekday x time level
	collapse (mean) delay*, by(hff_clean dow) 

	*** USED in the text:
	preserve 
		gcollapse (mean) delay*, by(hff_clean)
		sum delay if hff_clean == 7
		local d7 = r(mean)
		sum delay if hff_clean == 9
		local d9 = r(mean)

		di " Difference `=`d9'-`d7''"
		di " Ratio      `=(`d9'-`d7') / `d7''"

		// .                 di " Difference `=`d9'-`d7''"
		// Difference 1.169727325439453

		// .                 di " Ratio      `=(`d9'-`d7') / `d7''"
		//  Ratio      .4837814809267294

	restore 

	twoway 	///
			(line delay hff_clean if dow==0 & inrange(hff_clean,0,24), lcolor(gs2))  ///
			(line delay hff_clean if dow==1 & inrange(hff_clean,0,24), lcolor(gs5))  ///
			(line delay hff_clean if dow==2 & inrange(hff_clean,0,24), lcolor(gs8))  ///
			(line delay hff_clean if dow==3 & inrange(hff_clean,0,24), lcolor(gs11))  ///
			(line delay hff_clean if dow==4 & inrange(hff_clean,0,24), lcolor(gs14))  ///
			, xlabel(0(2)24, grid) ///
			legend(order(1 "Mon" 2 "Tue" 3 "Wed" 4 "Thu" 5 "Fr")  cols(5) region(lwidth(none)) pos(6) ) ///
			graphregion(color(white)) bgcolor(white) ///
			xtitle("Departure Time (hour)") ///
			ytitle("Delay" "(min/km)", orientation(horizontal)) //ylabel(0/5)

	graph save   "`outputfolder'avg_delay_south_bangalore.gph", replace
	graph export "`outputfolder'avg_delay_south_bangalore.pdf", replace
	graph export "`outputfolder'avg_delay_south_bangalore.eps", replace
	graph export "`outputfolder'avg_delay_south_bangalore.png", replace width(1800)
