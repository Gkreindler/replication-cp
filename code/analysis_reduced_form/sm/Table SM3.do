
 * This do file: Selection into study - Recruitment Table
 * This version: 10 Jan 2018
 * Authors: Gabriel Kreindler

clear all
pause on
set more off

*** Table number
	local tnum "smtable3"
	local outputfolder "paper/tables/`tnum'/"
	cap mkdir "`outputfolder'"

* Path based on global set in profile.do
	local path "$congestion_pricing_path"
	cd "`path'"
	adopath ++ "code/ado/"

*** Experiment
	use "data/coded_gps_dta/treat_assign.dta", clear
	keep uidp dttreat atreat
	assert _N==497

	rename uidp uidp_original
	tempfile exp
	save 	`exp'

*** VOT
	use "data/coded_cto/baseline_coded.dta", clear

	gen am = time1 == "_am"
	
	*** VOT coding
	assert inlist(vot_slower, 10, 15)
	gen dum_slow_15 = vot_slower == 15

	assert inlist(vot_base_extra, 0, 20)
	gen dum_extra_20 = vot_base_extra == 20

	winsor2 VOT, cuts(1 99)
	gen VOT_wl = log(VOT_w)

	*** 
	reg VOT_w dum_slow_15 am, vce(cluster uidp_original)
	predict vot_r, residual
	replace vot_r = vot_r + _b[_cons]

	reg VOT_wl dum_slow_15 am, vce(cluster uidp_original)
	predict log_vot_r, residual
	replace log_vot_r = log_vot_r + _b[_cons]

	* Average value of time expressed in Rs. per hour 
	gen vot_r_perhr = VOT_w * 60 / vot_slower

	sum vot_r_perhr, d

	gcollapse (mean) vot_r log_vot_r vot_r_perhr, by(uidp_original)

	tempfile vots
	save `vots'


*********
*** DT SP
	use "data/coded_cto/baseline_coded.dta", clear

	drop if absurd_dt == 1
	drop if dt_delta > 0 & dt_cheaper == "earlier"
	drop if dt_delta < 0 & dt_cheaper == "later"

	* censor at one hour
	sum dt_delta, d
	replace dt_delta = -60 if dt_delta < - 60
	replace dt_delta =  60 if dt_delta >   60 & dt_delta != .

	* include zeros
	replace dt_delta = 0 if dt_intention == 0 // nochange is replaced from missing to 0


	/* Retaining only relevant variables */
	// keep dt_* uidp_original absurd_int time1 time 
	// drop *_c *_c_2

	gen abs_dt_delta = abs(dt_delta)
	gen log_abs_dt_delta = log(abs_dt_delta + 1)
	
	* Question randomized treatments
		assert inlist(dt_p10, 15, 30)
		gen dum_high_slope = dt_p10 == 30 // 7.5-15-30 is slope "1" and 15-30-60 is slope "0"

		assert inlist(dt_cheaper, "earlier", "later")
		gen dum_early = dt_cheaper == "earlier" // Dummy for cheaper earlier

	* time indicators
	gen dum_am = time1 == "_am"

	egen time_early = group(dum_early dum_am)

	*** residualize with respect to question randomization
	reg abs_dt_delta i.time_early, cl(uidp_original)
	predict abs_dt_r, residual
	replace abs_dt_r = abs_dt_r + _b[_cons]

	reg log_abs_dt_delta i.time_early, cl(uidp_original)
	predict log_abs_dt_r, residual
	replace log_abs_dt_r = log_abs_dt_r + _b[_cons]

	gcollapse (mean) abs_dt_r log_abs_dt_r, by(uidp_original)

	tempfile dts
	save 	`dts'







*** Get recruitment data
	use "data/coded_cto/recruitment coded.dta", clear

* merge experiment data
	merge 1:1 uidp_original using `exp'
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

*** Coding income data -- we are censoring at 
	gen income_level = exp(log_income)
	sum income_level, d
	sum income_level if in_exp == 1,d

	replace log_income =  log(87500) if (income_level > 100000) & !missing(income_level)
	replace log_income =  . if (income_level > 100000) & !missing(income_level)
	tab log_income

*** Coding
	gen self_employed = occupation == 4 if occupation != 99 & occupation != .

    * Creating dummies for status
    gen status_0_ref  = status_new == 0
    gen status_1_inel = status_new == 1
    gen status_2_call = status_new == 2 | status_new == 3
    // gen status_3_surv = 
    gen status_4_app  = status_new == 4
    gen status_5_exp  = status_new == 5

    assert status_0_ref + status_1_inel + status_2_call + status_4_app + status_5_exp == 1

    // encode occ_level_1, gen(occ_level)

*** extremely few obs for ineligible SP, so remove:
	foreach yvar of varlist vot_r vot_r_perhr log_vot_r abs_dt_r log_abs_dt_r{
		replace `yvar' = . if status_1_inel == 1
	}

*** Fraction of ineligible among non-refusals
	sum status_1_inel if status_new > 0
	local frac_inel = r(mean)

*** Weights to correct for 1/15 chance of observing ineligible and refusals
    gen has_car_info = stand_car!=.
    tab status_new has_car_info, m row
    bys status_new: gegen invw = mean(has_car_info)
    
	* among refusals, same fraction of ineligible
    gen 	obs_weight = round(1 / invw)
    replace obs_weight = round(1 / invw * (1-`frac_inel')) if status_new == 0 
    tab status_new obs_weight


*** Weights for age
    gen has_age = new_age != .
    tab status_new has_age, m row
    cap drop invw
    bys status_new: gegen invw = mean(has_age)

    * among refusals, same fraction of ineligible
    gen obs_weight_age = round(1 / invw)
    replace obs_weight_age = round(1 / invw * (1-`frac_inel')) if status_new == 0 
    tab status_new obs_weight_age

*** Other coding variables
	reg vehicle_logprice stand_car, r
	predict vehicle_logprice_r, r
	replace vehicle_logprice_r = vehicle_logprice_r + _b[_cons]

	* gender -- same missing as age
	assert missing(new_age) == missing(gender)
	gen male = gender == 0 if !missing(gender)


***************************
*** Table B : regress in_exp on other variables
***************************

*** regress on response status
	foreach yvar of varlist km_daily log_income vot_r log_vot_r abs_dt_r log_abs_dt_r vot_r_perhr {
		 qui sum `yvar' 
		 gen `yvar'_z = (`yvar' - r(mean)) / r(sd)
	}

	foreach yvar of varlist /* stand_car */ vehicle_logprice{
		 qui sum `yvar' [fweight=obs_weight]
		 gen `yvar'_z = (`yvar' - r(mean)) / r(sd)
	}

	qui sum new_age [fweight=obs_weight_age]
	gen new_age_z = (new_age - r(mean)) / r(sd)

	// replace abs_dt_r = - abs_dt_r

*** Compare all eligibles
    reg in_exp stand_car 							if status_new != 1 [fweight=obs_weight], robust
    sum in_exp if e(sample) [fweight=obs_weight]
    estadd scalar cmean = r(mean)
    estimates store c1

    reg in_exp new_age_z  							if status_new != 1 [fweight=obs_weight_age], robust
    estimates store c2

    reg in_exp vehicle_logprice_z 					if status_new != 1 [fweight=obs_weight], robust
    estimates store c3

    reg in_exp stand_car new_age_z vehicle_logprice_z 	if status_new != 1 [fweight=obs_weight], robust
    estimates store c4

    * income -- whenever there is real data, truncating above 87500
    reg in_exp log_income_z if status_new != 1, robust
    sum in_exp if e(sample)
    estadd scalar cmean = r(mean)
    estimates store c5


*** Compare to other survey + apps
    reg in_exp km_daily_z if inlist(status_new,2,3,4,5), robust
    sum in_exp if e(sample)
    estadd scalar cmean = r(mean)
    estimates store c6

    reg in_exp /* vot_r_z */ vot_r_perhr_z abs_dt_r_z 	if inlist(status_new,2,3,4,5), robust
    sum in_exp if e(sample)
    estadd scalar cmean = r(mean)
    estimates store c7


*** OUTPUT
    * file names
	cap erase "`outputfolder'panelAv1.tex"
	cap erase "`outputfolder'table_v1.tex"

	esttab c* using "`outputfolder'panelAv1.tex", ///
		append keep(stand_car new_age_z km_daily_z log_income_z vehicle_logprice_z vot_r_perhr_z abs_dt_r_z) ///
		b(%12.3f) se(%12.3f) ///
		coeflabels(	stand_car 		"Drives Car" ///
					new_age_z   	"Age (z-score)" ///
					vehicle_logprice_z "Log Vehicle Value (z-score)" ///
					log_income_z 	"Log Income (Self-Reported, z-score)" ///
					km_daily_z 		"KM Daily (Stated, z-score)" ///
					vot_r_perhr_z 		"Value of Time (Stated, z-score)" ///
					abs_dt_r_z 		"Schedule Cost (Stated, z-score)") ///
		starlevels(* .10 ** .05 *** .01) nomtitles nonumbers ///
		stats(N cmean, labels("Observations" "Fraction in Experiment") fmt("%12.0fc" "%12.2f")) booktabs ///
		substitute(_ \_ "<" "$<$" "\midrule" "\addlinespace\addlinespace") fragment


*** Write main table tex file
	file open myfile using "`outputfolder'/table_v1.tex", write replace
	file write myfile "\resizebox{1\textwidth}{!}{" _n "\begin{tabular}{lcccc@{\hspace{2em}}ccc}" _n "\toprule" _n  ///
					  " & (1) & (2) & (3) & (4) & (5) & (6) & (7) \\" _n ///
					  "\addlinespace " _n ///
					  " & \multicolumn{7}{c}{\emph{Outcome: Respondent In Experiment}} \\" _n ///
					  "\emph{Sample:} & \multicolumn{4}{c}{All Eligible Respondents} & \multicolumn{3}{c}{Survey Respondents} \\" _n ///
					  "\addlinespace " _n ///
				   	  "\ExpandableInput{tables/`tnum'/panelAv1}" _n ///
					  "\bottomrule" _n "\end{tabular}" _n "}"

	file close myfile


***************************
*** Table A : Summary Stats
***************************


forv i=1/9{
	local yvar: word `i' of male 				new_age 	stand_car 	vehicle_logprice_r 				log_income 		km_daily 	log_income 	vot_r_perhr abs_dt_r
	local name: word `i' of "Male respondent" 	"Age" 		"Car driver" "Log vehicle price (residual)" "Log Income (self-reported)" 	"Stated Daily Travel (Km/day)" 	"Log income" 	"Stated Value of Time (Rs/hr)" 	"Stated Schedule Flexibility (min)"
	local sample: word `i' of "All" "All"		"All"		"All"			 				"Survey"		"Survey"	"Survey"	"Survey" "Survey"
	local dec : word `i' of 2 	1 		 	2			1 								1				1 			2 			1 		1
	local gpre: word `i' of "0," "0," 		"0,"		"0," 							"" 				""				"" 			""			""		""
	local weig: word `i' of "[fweight=obs_weight_age]" "[fweight=obs_weight_age]" "[fweight=obs_weight]" "[fweight=obs_weight]" "" "" ""

	sum   `yvar' if in_exp == 0 & inlist(status_new,`gpre'2,3,4,5) `weig'
 	local `yvar'_mn0: di %4.`dec'f r(mean)
 	local `yvar'_sd0: di %3.`dec'f r(sd)

 	sum   `yvar' if in_exp == 1 & inlist(status_new,`gpre'2,3,4,5) `weig'
 	local `yvar'_mn1: di %4.`dec'f r(mean)
 	local `yvar'_sd1: di %3.`dec'f r(sd)

 	reg   `yvar' in_exp if inlist(status_new,`gpre'2,3,4,5) `weig', robust
 	local N: di %6.0fc e(N)
 	local ts0 = _b[in_exp] / _se[in_exp]
 	// local ts: di %4.1f `ts0'
 	local ts ""
 	if (abs(`ts0') > 1.645) {
 		local ts "`ts'*"
 	}
 	if (abs(`ts0') > 1.96 ) {
 		local ts "`ts'*"
 	}
 	if (abs(`ts0') > 2.576) {
 		local ts "`ts'*"
 	}
 	// local `yvar'_t = "`ts'"
 	loca `yvar'_fxsd: di %4.2f `=_b[in_exp] / ``yvar'_sd0''
 	local `yvar'_fxsd "``yvar'_fxsd'`ts'"

 	// local my_text_`yvar'_1 "`name' & ``yvar'_mn1'	& ``yvar'_mn0'		& ``yvar'_fxsd' & `N' \\ "
 	// local my_text_`yvar'_2 "       &[``yvar'_sd1']	& [``yvar'_sd0']	&  & \\ "
 	local my_text_`yvar'_1 "`name' & ``yvar'_mn1' & [``yvar'_sd1']	& ``yvar'_mn0' & [``yvar'_sd0']		& ``yvar'_fxsd' & `N' \\ "
 	// local my_text_`yvar'_2 "       &[``yvar'_sd1']	& [``yvar'_sd0']	&  & \\ "

}

*** OUTPUT
    * file names
	cap erase "`outputfolder'selection_summary.tex"

*** Write main table tex file

	file open myfile using "`outputfolder'/selection_summary.tex", write replace
	file write myfile "\resizebox{0.85\textwidth}{!}{" _n "\begin{tabular}{lcccccc}" _n "\toprule" _n  ///
					  " & (1) & (2) & (3) & (4) & (5) & (6) \\" _n ///
					  "\addlinespace" _n ///
					  " & \multicolumn{2}{c}{In Experiment} & \multicolumn{2}{c}{Not in Experiment} & Difference & \\" _n ///
					  "\addlinespace" _n ///
					  " & Mean & [SD] & Mean & [SD]  &  in SD units & N \\" _n ///
					  "\addlinespace" _n ///
					  "\multicolumn{4}{l}{\textit{Panel A. All Respondents Approached}} \\" _n ///
					  	"`my_text_male_1'" _n ///
					  	"\addlinespace" _n ///
					  	"`my_text_new_age_1'" _n ///
					  	"\addlinespace" _n ///
						"`my_text_stand_car_1'" _n ///
						"\addlinespace" _n ///
						"`my_text_vehicle_logprice_r_1'" _n ///
						"\addlinespace\addlinespace" _n ///
						"\multicolumn{4}{l}{\textit{Panel B. Survey Respondents}} \\" _n ///
						"`my_text_log_income_1'" _n ///
						"\addlinespace" _n ///
						"`my_text_km_daily_1'" _n ///
						"\addlinespace" _n ///
						"`my_text_vot_r_perhr_1'" _n ///
						"\addlinespace" _n ///
						"`my_text_abs_dt_r_1'" _n ///
					  "\bottomrule" _n "\end{tabular}" _n "}"

	file close myfile

***************************
*** Table C : occupations
***************************

	decode occ_coarse, gen(occ_level_1)

    tab occ_level_1 in_exp if inlist(status_new,2,3,4,5), col nof
    tab occ_level_1

    tab occ_level_1 if inlist(status_new,2,3,4) [fweight=obs_weight]
    tab occ_level_1 if inlist(status_new,5) [fweight=obs_weight]

    tab self_employed if inlist(status_new,2,3,4,5) [fweight=obs_weight]
    tab self_employed if inlist(status_new,5) [fweight=obs_weight]

    gen one = 1
    sum one if occ_level_1 != "" & inlist(status_new,2,3,4) [fweight=obs_weight]
    local n_occ_nexp = r(sum)
    sum one if occ_level_1 != "" & inlist(status_new,5) [fweight=obs_weight]
    local n_occ_exp = r(sum)

    local n_occ_nexp_f: di %5.0fc `n_occ_nexp'
    local n_occ_exp_f : di %5.0fc `n_occ_exp'

	forv i=1/9{
	    local name:    word `i' of "High level staff (C-suite, managers) and business owners(both high and low)" 	"Mobile professions" 	"Non-office/Manual labour jobs" 	"Office staff (low-medium)" 	"Other Engineers and Semi-technical" 	"Others and Retired" 	"Software and IT" 	"Specialized Professions" 	"Student" 
		local name`i': word `i' of 	"\thead{Business owner\\ or manager}" ///
									"\thead{Mobile\\ professions}" ///
									"\thead{Manual\\ jobs}" ///
									"\thead{Office\\ staff}" ///
									"\thead{Engineers,\\ Technical}" ///
									"\thead{Others,\\ Retired}" ///
									"\thead{Software\\ and IT}" ///
									"\thead{Accountant,\\ Teacher,\\ Doctor}" ///
									"Student" 
		qui sum one if occ_level_1 == "`name'" & inlist(status_new,2,3,4) [fweight=obs_weight]
		local n_nexp_`i': di %3.1f `=r(sum)/`n_occ_nexp'*100'

		qui sum one if occ_level_1 == "`name'" & inlist(status_new,5) [fweight=obs_weight]
		local n_exp_`i': di %3.1f `=r(sum)/`n_occ_exp'*100'

		di "`n_nexp_`i''"
		di "`n_exp_`i''"
	}


*** OUTPUT
    * file names
	cap erase "`outputfolder'table_occupations.tex"

*** Write main table tex file
	file open myfile using "`outputfolder'/table_occupations.tex", write replace
	file write myfile "\resizebox{\textwidth}{!}{" _n "\begin{tabular}{lcccccccccc}" _n "\toprule" _n  ///
					  " 						& (1) 		 & (2) 		  & (3) 	   & (4) 		& (5) 		 & (6) 		  & (7) 	   & (8) 		& (9) 	 	 & (10)  \\" _n ///
					  "  						& `name1'	 & `name8'	  & `name7'	   & `name5'	& `name4'	 & `name3'	  & `name2'	   & `name9'	& `name6'    & Total \\" _n ///
					  "\addlinespace " _n ///
					  "\multicolumn{4}{l}{\textit{Panel C. Survey Respondents}} \\" _n ///
					  "\addlinespace " _n ///
					  " In Experiment 		& `n_exp_1'  & `n_exp_8'  & `n_exp_7'  & `n_exp_5'  & `n_exp_4'  & `n_exp_3'  & `n_exp_2'  & `n_exp_9'  & `n_exp_6'	 & `n_occ_exp_f'  \\" _n ///
					  " Not in Experiment 	& `n_nexp_1' & `n_nexp_8' & `n_nexp_7' & `n_nexp_5' & `n_nexp_4' & `n_nexp_3' & `n_nexp_2' & `n_nexp_9' & `n_nexp_6' & `n_occ_nexp_f' \\" _n ///
					  "\bottomrule" _n "\end{tabular}" _n "}"

	file close myfile
