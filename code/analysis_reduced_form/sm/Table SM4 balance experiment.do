
 * This do file: Randomization Balance Table 
 * This version: Oct 26 2017
 * Authors: Gabriel Kreindler

clear all
pause on
set more off

* Path based on global set in profile.do
	local path "$congestion_pricing_path"
	cd "`path'"
	adopath++ "code/ado/"

*** TABLE number
local tfolname "smtable4"
local outputfolder "paper/tables/`tfolname'/"
cap mkdir "`outputfolder'"

*** Experiment
	use "data/coded_gps_dta/treat_assign.dta", clear
	keep uidp dttreat atreat dt_high dt_low dt_info dt_ctrl strat_cell a_early a_late
	gen in_area_treat = atreat != "0 No treatment"
	assert _N==497

	rename uidp uidp_original
	tempfile treat_status
	save 	`treat_status'

*** VOT
	use "data/coded_cto/baseline_coded.dta", clear

	*** Sample
	* "Earliest survey (early wave included)"
	* "Earliest survey in late wave"
	* "Most recent survey"
	keep if baseline_sample_2 == 1
	gisid uidp_original time

	sum VOT, d
	winsor2 VOT, cuts(1 99)
	gen VOT_wl = log(VOT_w)
	sum VOT_w, d

	reg VOT_w vot_base_extra vot_slower /* am_delta_t */ , vce(cluster uidps)

	keep uidp_original time VOT_w vot_base_extra vot_slower VOT_wl 
	reshape wide VOT_w vot_base_extra vot_slower VOT_wl, i(uidp_original) j(time)

	tempfile vots
	save `vots'

*** DT SP
	use "data/coded_cto/baseline_coded.dta", replace

	*** Sample
	* "Earliest survey (early wave included)"
	* "Earliest survey in late wave"
	* "Most recent survey"
	keep if baseline_sample_2 == 1
	gisid uidp_original time

	drop if absurd_dt == 1
	drop if !inrange(dt_delta, -60,60)
	drop if dt_delta > 0 & dt_cheaper == "earlier"
	drop if dt_delta < 0 & dt_cheaper == "later"
	assert absurd_int == 0

	/* Retaining only relevant variables */
	keep dt_* uidp_original absurd_int time1 time 
	drop *_c *_c_2

	assert dt_delta == . if dt_intention ==0 // nochange is replaced from missing to 0
	replace dt_delta = 0 if dt_intention ==0 // nochange is replaced from missing to 0

	gen abs_dt_delta = abs(dt_delta)
	gen log_abs_dt_delta = log(abs_dt_delta)
	
	gen dt_slope_large = dt_p10 == 30 // 7.5-15-30 is slope "1" and 15-30-60 is slope "0"
	gen dt_cheap_early = dt_cheaper == "earlier" // Dummy for cheaper earlier

	reg abs_dt_delta dt_slope_large, r
	bys dt_cheaper time1: sum abs_dt_delta

	keep uidp_original time dt_slope_large dt_cheap_early dt_delta abs_dt_delta
	reshape wide dt_slope_large dt_cheap_early dt_delta abs_dt_delta, i(uidp_original) j(time)

	tempfile dts
	save 	`dts'


*** Get recruitment data
	use "data/coded_cto/recruitment coded.dta", clear


* merge experiment data
	merge 1:1 uidp_original using `treat_status'
	assert _m!=2
	gen in_exp = _m==3
	drop _m

* merge VOT data 
	merge 1:1 uidp_original using `vots'
	assert _m!=2
	gen in_vot = _m==3
	drop _m

* merge DT data
	merge 1:1 uidp_original using `dts'
	assert _m!=2
	gen in_dt = _m==3
	drop _m

	gen status_new = status
	replace status_new = 5 if in_exp == 1


*** Coding

	* Creating dummies for status
    gen status_0_ref  = status_new == 0
    gen status_1_inel = status_new == 1
    gen status_2_call = status_new == 2 | status_new == 3
    // gen status_3_surv = 
    gen status_4_app  = status_new == 4
    gen status_5_exp  = status_new == 5

    assert status_0_ref + status_1_inel + status_2_call + status_4_app + status_5_exp == 1

	* interview  start time (hour)
	gen start_hff = starthour + startmin / 60

	* Recruitment shift - AM/PM
	assert start_hff~=.
	gen start_shift_am = start_hff < 14 

	* recruitment time relative to shift start time (8am / 4pm)
	gen 	start_hff_rel = start_hff - 8
	replace start_hff_rel = start_hff - 16 if start_shift_am == 0

*** extremely few obs for ineligible SP, so remove:
	foreach yvar of varlist VOT_w1 VOT_w2 VOT_wl1 VOT_wl2 abs_dt_delta1 abs_dt_delta2{
		replace `yvar' = . if status_1_inel == 1
	}

	keep uidp_original km_daily stand_car new_age log_income vehicle_logprice VOT_w1 VOT_w2 VOT_wl1 VOT_wl2 abs_dt_delta1 abs_dt_delta2 ///
		 drive_days_this drive_days_other km_reading_clean start_hff start_shift_am start_hff_rel


	tempfile recruitment_vars
	save 	`recruitment_vars'


*** Pre GPS behavior vars
	*** Dummy for being at home/work
	use "data/coded_gps_dta/coded_days_15.dta", clear

	keep if post == 0
	assert inrange(dow,1,5)

	gen regular_dest = inlist(regular_commuter,2,3)

	collapse (mean) report_outstation km_daily_tot qual_gb qual_no_data dropped date_wo_trips ///
					at_h at_w at_w_all ///
					nt_ok_ km_ok_ du_ok_ ra_ok_ ch_ok_ cha_ok_ nt0_ok_ datestart ///
					regular_dest regular_commuter ///
					, by(uidp)

	isid uidp
	rename uidp uidp_original
	merge 1:1 uidp_original using `recruitment_vars'
	assert _m!=1
	keep if _m==3
	drop _m

	merge 1:1 uidp_original using `treat_status'
	assert _m==3
	drop _m	

********************************
*** ANALYSIS *** *** *** *** ***
********************************
	
file open myfile using "`outputfolder'table.tex", write replace
file write myfile /* "\resizebox{\textwidth}{!}{" _n  */ "\begin{tabular}{clcccccccccccc}" _n "\toprule" _n  ///
					  " & & \multicolumn{8}{c}{\textbf{Departure Time Treatments}} & \multicolumn{4}{c}{\textbf{Route Treatment}} \\" _n ///
					  "\addlinespace" _n ///
					  " & & \multicolumn{2}{c}{Information} & " ///
					  " \multicolumn{2}{c}{Low Rate} & " ///
					  " \multicolumn{2}{c}{High Rate} & " ///
					  " Obs. & \thead{Control \\ Mean} & " ///
					  " Route Early & & Obs. & \thead{Control \\ Mean} \\" _n ///
					  "\addlinespace" _n ///
					  " & & & (S.E.) & & (S.E.) & & (S.E.) & & & & (S.E.) \\" _n ///
					  "\addlinespace\addlinespace" _n

	label var stand_car 		"Car user"
	label var regular_dest 		"Regular destination"
	label var new_age 			"Age"
	label var vehicle_logprice 	"Log vehicle price"
	label var log_income 		"Log income"
	label var qual_gb 			"Frac days with good GPS data"
	label var at_w_all 			"Frac days present at work"
	label var nt_ok_ 			"Number of trips per day"
	label var km_ok_ 			"Total distance per day (Km.)"
	label var du_ok_ 			"Total duration per day (min)"
	label var ra_ok_ 			"Total D.T. hypothetical rate per day"
	label var cha_ok_ 			"Total Route hypothetical rate per day "

	forv i=1/12{
		local yvar: word `i' of stand_car ///
							 	regular_dest ///
							 	new_age ///
							 	vehicle_logprice ///
							 	log_income ///
								qual_gb ///
								at_w_all ///
								nt_ok_ ///
								km_ok_ ///
								du_ok_ ///
								ra_ok_ ///
								cha_ok_	

		** Departure Time
		qui reg `yvar' dt_high dt_low dt_info i.strat_cell, r

		local dt_b_info: di %4.2f _b[dt_info]	
		local dt_b_low : di %4.2f _b[dt_low]
		local dt_b_high: di %4.2f _b[dt_high]

		local dt_se_info: di %4.2f _se[dt_info]	
		local dt_se_low : di %4.2f _se[dt_low]
		local dt_se_high: di %4.2f _se[dt_high]

		forv j=1/3{
			local tr: word `j' of "info" "low" "high"

			local z_val =  abs(_b[dt_`tr'] / _se[dt_`tr'])

			if `z_val' > 2.58{
				local dt_b_`tr' = "`dt_b_`tr''$^{***}$"
			}
			else if `z_val' > 1.96{
				local dt_b_`tr' = "`dt_b_`tr''$^{**}$"
			}
			else if `z_val' > 1.63{
				local dt_b_`tr' = "`dt_b_`tr''$^{*}$"
			}
			else{
				local dt_b_`tr' = "`dt_b_`tr''"
			}
		}

		local dt_nobs = e(N)

			qui sum `yvar' if dt_ctrl == 1 & e(sample) == 1 // control mean
			local dt_control_mean: di %4.2f r(mean)

			qui estadd scalar control_mean = r(mean)
			qui estadd scalar control_sd = r(sd)
			*** Store
			qui estimates store `yvar'_dt
		
		** Area
		qui reg `yvar' a_early i.strat_cell if in_area_treat==1, r

		local area_b : di %4.2f _b[a_early]
		local area_se: di %4.2f _se[a_early]
		local area_nobs = e(N)

		local z_val =  abs(_b[a_early] / _se[a_early])

		if `z_val' > 2.58{
			local area_b = "`area_b'$^{***}$"
		}
		else if `z_val' > 1.96{
			local area_b = "`area_b'$^{**}$"
		}
		else if `z_val' > 1.63{
			local area_b = "`area_b'$^{*}$"
		}
		else{
			local area_b = "`area_b'"
		}

			qui sum `yvar' if a_late == 1 & e(sample) == 1 // control mean
			local area_control_mean: di %4.2f r(mean)
			qui estadd scalar control_mean = r(mean)
			qui estadd scalar control_sd = r(sd)
			*** Store
			qui estimates store `yvar'_area
			// qui eststo `yvar'_area
		

		** Area vs Non-Area
		// qui reg `yvar' in_area_treat, r
		// 	qui sum `yvar' if in_area_treat == 0 & e(sample) == 1 // control mean
		// 	qui estadd scalar control_mean = r(mean)
		// 	qui estadd scalar control_sd = r(sd)
		// 	*** Store
		// 	qui estimates store `yvar'_area_comp
			// qui eststo `yvar'_area_comp

		local var_label : variable label `yvar'


		file write myfile " (`i') & `var_label' & `dt_b_info' & (`dt_se_info') & " ///
												" `dt_b_low' & (`dt_se_low') & " ///
												"`dt_b_high' & (`dt_se_high') & " ///
												" `dt_nobs' & `dt_control_mean' & " ///
												" `area_b' & (`area_se') & " ///
												" `area_nobs' & `area_control_mean' \\" _n ///
												"\addlinespace" _n
	}
file close myfile

** Joint test - DT

forv i=1/5{
	local nm_test: word `i' of "dt_high" "dt_low" "dt_info" "dt_all" "area_early"
	local suffix : word `i' of "dt" "dt" "dt" "dt" "area"
	local my_test: word `i' of "(_b[dt_high] = 0)" ///
							   "(_b[dt_low] = 0)" ///
							   "(_b[dt_info] = 0)" ///
							   "(_b[dt_high] = 0) (_b[dt_low] = 0) (_b[dt_info] = 0)" ///
							   "(_b[a_early] = 0)"
	
	qui estimates restore stand_car_`suffix'
	// test (_b[dt_high] = 0) (_b[dt_low] = 0) (_b[dt_info] = 0), notest
	qui test `my_test', notest

	foreach v in  regular_dest 	new_age 	 	vehicle_logprice	log_income ///
				qual_gb 	at_w_all 	nt_ok_ 	km_ok_	du_ok_	ra_ok_{

		qui estimates restore `v'_`suffix'				
		// test (_b[dt_high] = 0) (_b[dt_low] = 0) (_b[dt_info] = 0), notest accum
		qui test `my_test', notest accum
		
	}

	qui estimates restore cha_ok__`suffix'
	// test (_b[dt_high] = 0) (_b[dt_low] = 0) (_b[dt_info] = 0), accum
	test `my_test', accum

	display _n _n _n " `nm_test' " _n _n _n
	local chi2_`nm_test': display %3.2f r(F)
	local pval_`nm_test': display %3.2f r(p)

	di `chi2_`nm_test''
	di `pval_`nm_test''
	
}

file open myfile using "`outputfolder'table.tex", write append

file write myfile	"\addlinespace\addlinespace" _n ///
					 "(13) & Joint Significance Test F-stat & " ///
						" `chi2_dt_all' & & & & & & & & `chi2_area_early' \\" _n ///
					"(14) & Joint Significance Test P-value & " ///
						" `pval_dt_all' & & & & & & & & `pval_area_early' \\ " _n ///
					"\bottomrule" _n "\end{tabular}" _n /* _n "}" */ 

file close myfile
