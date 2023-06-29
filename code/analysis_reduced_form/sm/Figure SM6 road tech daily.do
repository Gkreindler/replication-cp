
 * This do file: Road tech at daily level
 * This version: 1 Dec 2021
 * Authors: Gabriel Kreindler
 
clear all
pause off
set more off
set seed 329743


* Path based on global set in profile.do
	local path "$congestion_pricing_path"
	cd "`path'"

*** Load data
	use "data/coded_road_tech/density_th", clear

	merge 1:1 date hour using "data/coded_road_tech/google_maps_th_level"
	keep if _m==3
	drop _m

*** Coding
	nl (delay_gm = {l0=2.0} + {l1=1000.0} * (density_norm ^ {gamma=1.0}))
	predict delay_mean_hour_predicted
	
*** Non-linear least squares
	gen delay_pred = .

	keep if inrange(date, 20920, 21009)

	levelsof date, local(all_possible_dates)

	foreach my_date in `all_possible_dates'{
		di %td `my_date'

		qui nl (delay_gm = {l0=2.0} + {l1=1.5} * (density_norm ^ {gamma=1.0})) if date == `my_date' // , vce(hac nwest 180)
		qui cap drop delay_temp 
		qui predict delay_temp
		qui replace delay_pred = delay_temp if date == `my_date'

	}

*** Prep graphs	
	local ffolname "smfigure6"
	local outputfolder "paper/figures/`ffolname'/"
	cap mkdir "`outputfolder'"
	cap mkdir "`outputfolder'bydate/"

*** All days, only NL fit
	gen dow = dow(date)

	format %4.0f date
	tab date if dow != 0
	tab date if dow == 0

	sort date density_norm
	local lcolorspec_w_day "blue%20"
	local lcolorspec_w_end "green%75"

	twoway 	(line delay_pred density_norm if date == 20920, lcolor(`lcolorspec_w_day')) /// 
			(line delay_pred density_norm if date == 20921, lcolor(`lcolorspec_w_day')) /// 
			(line delay_pred density_norm if date == 20922, lcolor(`lcolorspec_w_day')) /// 
			(line delay_pred density_norm if date == 20923, lcolor(`lcolorspec_w_day')) /// 
			(line delay_pred density_norm if date == 20924, lcolor(`lcolorspec_w_day')) /// 
			(line delay_pred density_norm if date == 20927, lcolor(`lcolorspec_w_day')) /// 
			(line delay_pred density_norm if date == 20928, lcolor(`lcolorspec_w_day')) /// 
			(line delay_pred density_norm if date == 20929, lcolor(`lcolorspec_w_day')) /// 
			(line delay_pred density_norm if date == 20930, lcolor(`lcolorspec_w_day')) /// 
			(line delay_pred density_norm if date == 20931, lcolor(`lcolorspec_w_day')) /// 
			(line delay_pred density_norm if date == 20933, lcolor(`lcolorspec_w_day')) /// 
			(line delay_pred density_norm if date == 20934, lcolor(`lcolorspec_w_day')) /// 
			(line delay_pred density_norm if date == 20935, lcolor(`lcolorspec_w_day')) /// 
			(line delay_pred density_norm if date == 20936, lcolor(`lcolorspec_w_day')) /// 
			(line delay_pred density_norm if date == 20937, lcolor(`lcolorspec_w_day')) /// 
			(line delay_pred density_norm if date == 20938, lcolor(`lcolorspec_w_day')) /// 
			(line delay_pred density_norm if date == 20941, lcolor(`lcolorspec_w_day')) /// 
			(line delay_pred density_norm if date == 20942, lcolor(`lcolorspec_w_day')) /// 
			(line delay_pred density_norm if date == 20943, lcolor(`lcolorspec_w_day')) /// 
			(line delay_pred density_norm if date == 20944, lcolor(`lcolorspec_w_day')) /// 
			(line delay_pred density_norm if date == 20945, lcolor(`lcolorspec_w_day')) /// 
			(line delay_pred density_norm if date == 20947, lcolor(`lcolorspec_w_day')) /// 
			(line delay_pred density_norm if date == 20948, lcolor(`lcolorspec_w_day')) /// 
			(line delay_pred density_norm if date == 20949, lcolor(`lcolorspec_w_day')) /// 
			(line delay_pred density_norm if date == 20950, lcolor(`lcolorspec_w_day')) /// 
			(line delay_pred density_norm if date == 20951, lcolor(`lcolorspec_w_day')) /// 
			(line delay_pred density_norm if date == 20952, lcolor(`lcolorspec_w_day')) /// 
			(line delay_pred density_norm if date == 20954, lcolor(`lcolorspec_w_day')) /// 
			(line delay_pred density_norm if date == 20955, lcolor(`lcolorspec_w_day')) /// 
			(line delay_pred density_norm if date == 20956, lcolor(`lcolorspec_w_day')) /// 
			(line delay_pred density_norm if date == 20957, lcolor(`lcolorspec_w_day')) /// 
			(line delay_pred density_norm if date == 20958, lcolor(`lcolorspec_w_day')) /// 
			(line delay_pred density_norm if date == 20959, lcolor(`lcolorspec_w_day')) /// 
			(line delay_pred density_norm if date == 20961, lcolor(`lcolorspec_w_day')) /// 
			(line delay_pred density_norm if date == 20962, lcolor(`lcolorspec_w_day')) /// 
			(line delay_pred density_norm if date == 20963, lcolor(`lcolorspec_w_day')) /// 
			(line delay_pred density_norm if date == 20964, lcolor(`lcolorspec_w_day')) /// 
			(line delay_pred density_norm if date == 20965, lcolor(`lcolorspec_w_day')) /// 
			(line delay_pred density_norm if date == 20966, lcolor(`lcolorspec_w_day')) /// 
			(line delay_pred density_norm if date == 20968, lcolor(`lcolorspec_w_day')) /// 
			(line delay_pred density_norm if date == 20969, lcolor(`lcolorspec_w_day')) /// 
			(line delay_pred density_norm if date == 20970, lcolor(`lcolorspec_w_day')) /// 
			(line delay_pred density_norm if date == 20971, lcolor(`lcolorspec_w_day')) /// 
			(line delay_pred density_norm if date == 20972, lcolor(`lcolorspec_w_day')) /// 
			(line delay_pred density_norm if date == 20973, lcolor(`lcolorspec_w_day')) /// 
			(line delay_pred density_norm if date == 20975, lcolor(`lcolorspec_w_day')) /// 
			(line delay_pred density_norm if date == 20976, lcolor(`lcolorspec_w_day')) /// 
			(line delay_pred density_norm if date == 20977, lcolor(`lcolorspec_w_day')) /// 
			(line delay_pred density_norm if date == 20978, lcolor(`lcolorspec_w_day')) /// 
			(line delay_pred density_norm if date == 20979, lcolor(`lcolorspec_w_day')) /// 
			(line delay_pred density_norm if date == 20980, lcolor(`lcolorspec_w_day')) /// 
			(line delay_pred density_norm if date == 20982, lcolor(`lcolorspec_w_day')) /// 
			(line delay_pred density_norm if date == 20983, lcolor(`lcolorspec_w_day')) /// 
			(line delay_pred density_norm if date == 20984, lcolor(`lcolorspec_w_day')) /// 
			(line delay_pred density_norm if date == 20985, lcolor(`lcolorspec_w_day')) /// 
			(line delay_pred density_norm if date == 20986, lcolor(`lcolorspec_w_day')) /// 
			(line delay_pred density_norm if date == 20987, lcolor(`lcolorspec_w_day')) /// 
			(line delay_pred density_norm if date == 20989, lcolor(`lcolorspec_w_day')) /// 
			(line delay_pred density_norm if date == 20990, lcolor(`lcolorspec_w_day')) /// 
			(line delay_pred density_norm if date == 20991, lcolor(`lcolorspec_w_day')) /// 
			(line delay_pred density_norm if date == 20992, lcolor(`lcolorspec_w_day')) /// 
			(line delay_pred density_norm if date == 20993, lcolor(`lcolorspec_w_day')) /// 
			(line delay_pred density_norm if date == 20994, lcolor(`lcolorspec_w_day')) /// 
			(line delay_pred density_norm if date == 20996, lcolor(`lcolorspec_w_day')) /// 
			(line delay_pred density_norm if date == 20997, lcolor(`lcolorspec_w_day')) /// 
			(line delay_pred density_norm if date == 20998, lcolor(`lcolorspec_w_day')) /// 
			(line delay_pred density_norm if date == 20999, lcolor(`lcolorspec_w_day')) /// 
			(line delay_pred density_norm if date == 21000, lcolor(`lcolorspec_w_day')) /// 
			(line delay_pred density_norm if date == 21001, lcolor(`lcolorspec_w_day')) /// 
			(line delay_pred density_norm if date == 21003, lcolor(`lcolorspec_w_day')) /// 
			(line delay_pred density_norm if date == 21004, lcolor(`lcolorspec_w_day')) /// 
			(line delay_pred density_norm if date == 21005, lcolor(`lcolorspec_w_day')) /// 
			(line delay_pred density_norm if date == 21006, lcolor(`lcolorspec_w_day')) /// 
			(line delay_pred density_norm if date == 21007, lcolor(`lcolorspec_w_day')) /// 
			(line delay_pred density_norm if date == 21008, lcolor(`lcolorspec_w_day')) /// 
			(line delay_pred density_norm if date == 21010, lcolor(`lcolorspec_w_day')) /// 
			(line delay_pred density_norm if date == 21011, lcolor(`lcolorspec_w_day')) /// 
			(line delay_pred density_norm if date == 21012, lcolor(`lcolorspec_w_day')) /// 
			(line delay_pred density_norm if date == 21013, lcolor(`lcolorspec_w_day')) /// 
			(line delay_pred density_norm if date == 21014, lcolor(`lcolorspec_w_day')) /// 
			(line delay_pred density_norm if date == 21015, lcolor(`lcolorspec_w_day')) /// 
			(line delay_pred density_norm if date == 21017, lcolor(`lcolorspec_w_day')) /// 
			(line delay_pred density_norm if date == 21018, lcolor(`lcolorspec_w_day')) /// 
			(line delay_pred density_norm if date == 21019, lcolor(`lcolorspec_w_day')) /// 
			(line delay_pred density_norm if date == 21020, lcolor(`lcolorspec_w_day')) /// 
			(line delay_pred density_norm if date == 21021, lcolor(`lcolorspec_w_day')) /// 
			(line delay_pred density_norm if date == 21022, lcolor(`lcolorspec_w_day')) /// 
			(line delay_pred density_norm if date == 21024, lcolor(`lcolorspec_w_day')) /// 
			(line delay_pred density_norm if date == 21025, lcolor(`lcolorspec_w_day')) /// 
			(line delay_pred density_norm if date == 21026, lcolor(`lcolorspec_w_day')) /// 
			(line delay_pred density_norm if date == 21027, lcolor(`lcolorspec_w_day')) /// 
			(line delay_pred density_norm if date == 21028, lcolor(`lcolorspec_w_day')) /// 
			(line delay_pred density_norm if date == 21029, lcolor(`lcolorspec_w_day')) /// 
			(line delay_pred density_norm if date == 21031, lcolor(`lcolorspec_w_day')) /// 
			(line delay_pred density_norm if date == 21032, lcolor(`lcolorspec_w_day')) /// 
			(line delay_pred density_norm if date == 21033, lcolor(`lcolorspec_w_day')) /// 
			(line delay_pred density_norm if date == 21034, lcolor(`lcolorspec_w_day')) /// 
			(line delay_pred density_norm if date == 21035, lcolor(`lcolorspec_w_day')) /// 
			(line delay_pred density_norm if date == 21040, lcolor(`lcolorspec_w_day')) /// 
			(line delay_pred density_norm if date == 20939, lcolor(`lcolorspec_w_end')) /// 
			(line delay_pred density_norm if date == 20946, lcolor(`lcolorspec_w_end')) /// 
			(line delay_pred density_norm if date == 20953, lcolor(`lcolorspec_w_end')) /// 
			(line delay_pred density_norm if date == 20960, lcolor(`lcolorspec_w_end')) /// 
			(line delay_pred density_norm if date == 20967, lcolor(`lcolorspec_w_end')) /// 
			(line delay_pred density_norm if date == 20974, lcolor(`lcolorspec_w_end')) /// 
			(line delay_pred density_norm if date == 20981, lcolor(`lcolorspec_w_end')) /// 
			(line delay_pred density_norm if date == 20988, lcolor(`lcolorspec_w_end')) /// 
			(line delay_pred density_norm if date == 20995, lcolor(`lcolorspec_w_end')) /// 
			(line delay_pred density_norm if date == 21009, lcolor(`lcolorspec_w_end')) /// 
			(line 	 delay_mean_hour_predicted 	density_norm, lcolor(red) lpattern(dash) sort ) ///
			, legend(order(1 "Mon-Sat" 108 "Sunday" 110 "all dates fit") cols(1) position(5) ring(0) ) graphregion(color(white)) ///
			xlabel(0(0.5)2.5, grid) xtitle("Traffic density") ylabel(2(1)4.5, nolabels grid) ytitle(" ") ///
			title(" Fit by calendar date ", position(12) ring(0) size(medium))

		graph export "`outputfolder'bydate\all_lines.png", replace
		graph save   "`outputfolder'bydate\all_lines.gph", replace

*** Randomly select date
	preserve
		keep date dow
		duplicates drop
		gen ur = runiform()
		// gen sunday = dow == 6
		sort dow ur, stable
		by dow: keep if _n == 1
		keep date
		gen to_include = 1

		tempfile to_include
		save 	`to_include'
	restore

	merge m:1 date using `to_include'
	assert _m!=2
	drop _m

	format %td date
	levelsof date if to_include == 1, local(random_list_of_dates)
	di `random_list_of_dates'
	// levelsof date, local(random_list_of_dates)
	// tab date if to_include == 1


*** Non-linear least squares
	sort density_norm

	cap drop delay_pred
	gen delay_pred = .
	
	local i = 1
	// local random_list_of_dates "`=mdy(5,11,2017)'" // testing
	foreach my_date in `random_list_of_dates'{
		local my_date_str: di %td `my_date'
		local my_date_tit: di %td_dd_Mon_CCYY,_DAYNAME `my_date'
		// format %td_dd_Mon_CCYY,_DAYNAME date

		qui nl (delay_gm = {l0=2.0} + {l1=1.5} * (density_norm ^ {gamma=1.0})) if date == `my_date' // , vce(hac nwest 180)
		qui cap drop delay_temp 
		qui predict delay_temp
		qui replace delay_pred = delay_temp if date == `my_date'

		local ytitle = " "
		local xtitle = " "
		local ylabel = "nolabels"
		local xlabel = "nolabels"
		local mylinecolor = "blue"
		local mydotcolor = "red"

		if (mod(`i',4) == 1){
			local ytitle = "Travel Delay (min/km)"
			local ylabel = ""
		}
		if (`i'>4){
			local xtitle = "Traffic density"
			local xlabel = ""
		}
		if (dow(`my_date') == 0){
			local mylinecolor = "green"
			local mydotcolor = "orange"
		}


		twoway 	(line 	 delay_mean_hour_predicted 	density_norm, lcolor(gs8) lpattern(dash) sort) ///
				(line 	 delay_pred 				density_norm if date == `my_date', lcolor(`mylinecolor')) ///
				(scatter delay_gm 					density_norm if date == `my_date', mcolor(`mydotcolor') msize(small) mlabel(hour) mlabcolor(red)) ///
				, legend(order(3 "hour" 2 "daily fit" 1 "all dates fit") cols(1) /* size(tiny) symysize(2pt) symxsize(2pt) */ position(5) ring(0) ) ///
				graphregion(color(white)) ///
				xtitle("`xtitle'") ytitle("`ytitle'") xlabel(0(0.5)2.5, grid `xlabel') ylabel(2(1)4.5, grid `ylabel') ///
				title(" `my_date_tit'", position(12) ring(0) size(medium))

		graph export "`outputfolder'bydate/scatter_`my_date_str'.png", replace
		graph save   "`outputfolder'bydate/scatter_`my_date_str'.gph", replace

		local i = `i' + 1
	}

*** Combined graph 
	
	cd "`outputfolder'bydate/"
	
	local list_of_graphs = ""
	foreach my_date in `random_list_of_dates'{
		local my_date_str: di %td `my_date'
		local list_of_graphs = "`list_of_graphs' scatter_`my_date_str'.gph"
	}
	di "`list_of_graphs'"

	graph combine `list_of_graphs' all_lines.gph ///
			, scheme(s1color) commonscheme ///
			col(4) ycommon xcommon imargin(zero) iscale(0.6)

	graph display, ysize(10) xsize(20) 

cd "../"
	local figure "combined_all"
	graph export "`figure'_1800.png", replace width(1800)
	graph export "`figure'_3600.png", replace width(3600)
	graph save "`figure'", replace
	graph export "`figure'.eps", replace
	graph export "`figure'.pdf", replace
