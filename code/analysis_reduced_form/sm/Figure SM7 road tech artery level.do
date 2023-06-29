
 * This do file: Road tech at the artery level
 * This version: January 8, 2022
 * Authors: Gabriel Kreindler

clear all
pause on
set more off
set matsize 10000

* Path based on global set in profile.do
	local path "$congestion_pricing_path"
	cd "`path'"

*** FIGURE number
	local ffolname "smfigure7"
	local outputfolder "paper/figures/`ffolname'/"
	cap mkdir "`outputfolder'"
	cap mkdir "`outputfolder'byodid/"

*** 
	import delimited using "data\raw_other\road tech artery level\road_segment_btm.csv", clear

	gisid id
	list

	local n_arteries = _N

	forv i=1/`n_arteries'{

		local e = gm_east[`i']
		local w = gm_west[`i']
		local n = gm_north[`i']
		local s = gm_south[`i']

		* route Google Maps ID and condition -- 1
		if (east_west[`i'] == 1){
			local gm_id_1_`i' = `e'
			local gm_id_2_`i' = `w'
			di `gm_id_1_`i''	
			
			*** both N-S and E-W
			if (north_sout[`i'] == 1){
				
				if( `e' == `n'){ 
					*** NE
					local dir_1_`i' = "NE"
					local dir_2_`i' = "SW"

					local cond_1_`i' = " dir_lng > 0 & dir_lat > 0 "
					local cond_2_`i' = " dir_lng < 0 & dir_lat < 0 "
				}
				else{
					*** SE
					local dir_1_`i' = "SE"
					local dir_2_`i' = "NW"

					local cond_1_`i' = " dir_lng > 0 & dir_lat < 0 "
					local cond_2_`i' = " dir_lng < 0 & dir_lat > 0 "
				}
			}
			else{
				*** Only E-W
				local dir_1_`i' = "E"
				local dir_2_`i' = "W"

				local cond_1_`i' = " dir_lng > 0 "
				local cond_2_`i' = " dir_lng < 0 "
			}
		}
		else{
			assert north_sout[`i'] == 1
			local gm_id_1_`i' = `n'
			local gm_id_2_`i' = `s'
			di `gm_id_1_`i''	

			*** Only N-S
			local dir_1_`i' = "N"
			local dir_2_`i' = "S"

			local cond_1_`i' " dir_lat > 0 "
			local cond_2_`i' " dir_lat < 0 "
		}

		* route notes 
		local note`i': display notes[`i']
		// di "`note`i''"
	}


*** LOAD GPS and Google by artery
forv i=1/`n_arteries'{
	local note `note`i''

	*** LOAD DATA SETS 
	// import delimited using "data/road_network/intermediate/all_trips_btm_pts.csv", clear
	local i1 = `i' - 1
	import delimited using "data/raw_other/road tech artery level/btm_volumes/all_trips_btm_pts_`i1'.csv", clear

	gen speed = dist / dur
	// sum speed, d

	*** date
 	gen date_ = date(date,"YMD")
 	drop date
 	rename date_ date
 	format %td date
 	gen dow = dow(date) // 0 = Sunday ... 6 = Saturday

 	*** time
	gen hour = substr(start,1,2)
	gen mins = substr(start,4,2)
	destring hour mins, replace
	gen start_ = hour + mins / 60
	gen mmd = hour*60 + mins

	*** coding
	gen delay = 1 / speed / 60 * 1000

	*** cleaning and sample
	// sum speed, d
	// drop if ~inrange(speed,0.5, 30)
	drop if length < 200
	drop if inlist(dow,0,6)

	qui merge m:1 mmd using "data/raw_other/road tech artery level/gmapi_delay_all.dta"
	assert _m!=1
	drop _m

	sort mmd

*** Plot volume and Google speed on the same direction


***** GRAPHS

	*****  DATASET
	forv j=1/2{
		
		preserve
			kdensity mmd if start!="" & `cond_`j'_`i'', bw(30) generate(volume) at(mmd) nograph
			
			gen n_trips_original = _N

			keep mmd volume delay_ipol`gm_id_`j'_`i'' n_trips_original
			duplicates drop
			gisid mmd
			rename delay_ipol`gm_id_`j'_`i'' delay_ipol
			gen odid = `gm_id_`j'_`i''
			gen direction = "`dir_`j'_`i''"

			tempfile chunk_`j'_`i'
			save 	`chunk_`j'_`i''
		restore

	}

}

*** LOAD ALL
	use `chunk_1_1', clear
	append using `chunk_2_1'

	forv i=2/`n_arteries'{
		forv j=1/2{
			append using `chunk_`j'_`i''	
		}
		
	}

	save "data/coded_road_tech/all_btm_vol_delay_density", replace

 */
*** ANALYSIS all together
	use "data/coded_road_tech/all_btm_vol_delay_density", clear

	// encode odid, gen(odid_)

	bys odid: gegen vol_tot = sum(volume)
	replace volume = volume / vol_tot 
	replace volume = volume * 1440

	bys mmd: gegen mean_vol = mean(volume)

	reg delay_ipol mean_vol, cl(odid)
	reg delay_ipol volume mean_vol, cl(odid)
	reg delay_ipol volume , cl(odid)

	reg volume mean_vol

	// fsdfds

	gen volume_shrunk = volume * 0.5 + mean_vol * 0.5
	reg delay_ipol volume_shrunk, cl(odid)
	reg delay_ipol mean_vol, cl(odid)


	// scatter delay_ipol volume
	// reghdfe delay_ipol c.volume##odid, noab
	// tsset mmd
	// xtset odid mmd
	// nl (delay_ipol = {l0=2.0} + {l1=1000.0} * (volume ^ {gamma=1.0})), vce(hac nwest 180)
	// nl (delay_ipol = {l0=2.0} + {l1=1000.0} * (volume ^ {gamma=1.0})), vce(cluster odid)

	gen hour = floor(mmd/60)
	gen is_hour = hour == mmd/60
	gen is_2hour = is_hour == 1 & mod(hour,2) == 0
	gen is_5hour = is_hour == 1 & mod(hour,5) == 0
	tab hour is_5hour

	glevelsof odid, local(all_odids)

	gen delay_city_predict = 2.09 + volume * 1.01
	sort volume, stable

	gen n_index = _n

	xtset odid mmd
	// sort odid mmd
	sort volume

	local i = 1
	foreach my_odid in `all_odids'{

		// local i=1
		// local my_odid = 31
		
		newey2 delay_ipol volume if odid == `my_odid', lag(180)
		cap drop delay_ipol_predict
		predict delay_ipol_predict if odid == `my_odid'
		cap drop delay_ipol_predict_sd
		predict delay_ipol_predict_sd, stdp

		cap drop ci_low 
		cap drop ci_hig
		gen ci_low = delay_ipol_predict - 1.96 * delay_ipol_predict_sd
		gen ci_hig = delay_ipol_predict + 1.96 * delay_ipol_predict_sd

		local ytitle = " "
		local xtitle = " "
		local ylabel = "nolabels"
		local xlabel = "nolabels"
		if (mod(`i',8) == 1){
			local ytitle = "Travel delay"
			local ylabel = ""
		}
		if (`i'>38){
			local xtitle = "Traffic density"
			local xlabel = ""
		}


		sum n_index if odid == `my_odid'
		local direction = direction[r(min)]

		twoway 	(line 	 delay_city_predict volume, lcolor(gs8) lpattern(dash)) ///
				(line 	 delay_ipol_predict volume if odid == `my_odid', lcolor(blue)) ///
				(scatter delay_ipol 		volume if odid == `my_odid' & is_hour == 1, mcolor(red) msize(small)) ///
				(scatter delay_ipol 		volume if odid == `my_odid' & is_5hour == 1, mcolor(red) mlabel(hour) msize(small) mlabcolor(red)) ///
				(line ci_low volume if odid == `my_odid', lcolor(blue) lpattern(dash)) ///
				(line ci_hig volume if odid == `my_odid', lcolor(blue) lpattern(dash)) ///
				, legend(order(1 "citywide fit" 2 "artery fit" 3 "artery") cols(2) size(tiny) symysize(2pt) symxsize(2pt) ) ///
				graphregion(color(white)) ///
				xtitle("`xtitle'") ytitle("`ytitle'") xlabel(0(2)6, grid `xlabel') ylabel(0(2)8, grid `ylabel') ///
				title("Route `my_odid' direction `direction'", position(12) ring(0) size(medium))

		graph export "`outputfolder'byodid/scatter_`my_odid'.png", replace
		graph save   "`outputfolder'byodid/scatter_`my_odid'.gph", replace

		local i = `i' + 1
	}

**** ONE GRAPH COMMON LEGEND
	cd "`outputfolder'byodid/"

	grc1leg "scatter_2" "scatter_3" "scatter_12" "scatter_13" "scatter_26" "scatter_27" ///
			"scatter_28" "scatter_29" "scatter_30" "scatter_31" "scatter_34" "scatter_35" ///
			"scatter_36" "scatter_37" "scatter_40" "scatter_41" "scatter_42" "scatter_43" ///
			"scatter_44" "scatter_45" "scatter_46" "scatter_49" "scatter_50" "scatter_51" ///
			"scatter_56" "scatter_57" "scatter_62" "scatter_63" "scatter_92" "scatter_93" ///
			"scatter_96" "scatter_97" "scatter_98" "scatter_99" "scatter_102" "scatter_103" ///
			"scatter_104" "scatter_105" "scatter_106" "scatter_107" "scatter_108" "scatter_109" ///
			"scatter_164" "scatter_165" "scatter_168" "scatter_169"  ///
			, scheme(s1color) commonscheme ///
			 legendfrom(scatter_2) ///
			col(8) /* ycommon */ xcommon imargin(zero) iscale(0.3) position(5) ring(0) 

	cd "../"
	local figure "combined_all_density"
	graph export "`figure'_1800.png", replace width(1800)
	graph export "`figure'_3600.png", replace width(3600)
	graph save   "`figure'", replace
	graph export "`figure'.eps", replace
	graph export "`figure'.pdf", replace


