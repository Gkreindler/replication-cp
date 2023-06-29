
 * This do file: Table III: Road Technology Supply Estimation: Travel Delay Linear in Traffic Density
 * This version: 8 Jan 2022
 * Authors: Gabriel Kreindler
 
clear all
pause off
set more off

* Path based on global set in profile.do
	local path "$congestion_pricing_path"
	cd "`path'"

*** TABLE number
	scalar tnum="table3"
	scalar outputfolder_table="paper/tables/`=tnum'/"
	cap mkdir "`=outputfolder_table'"


***************	
** ANALYSIS Main Table (using Google Maps outcome) and appendix table (using GPS speed outcome)
***************	

**************
*** date level
**************
	use "data/coded_road_tech/density_t", clear

	merge 1:1 date using "data/coded_road_tech/google_maps_t_level"
	keep if _m==3
	drop _m

	merge 1:1 date using "data/coded_road_tech/volumes_t_level", keepusing(delay_gps)
	drop if _m==2
	drop _m

	tsset date

	*** Column 4 OLS over dates
	qui reg delay_gm density_norm, r
	local myr2=e(r2_a)

	eststo date_gm: newey2 delay_gm density_norm, lag(14)
	qui estadd scalar my_r2=`myr2'
	qui sum density_norm
	qui estadd scalar xsd = r(sd)

	*** Column 8 OLS over dates
	qui reg delay_gps density_norm, r
	local myr2=e(r2_a)

	eststo date_gps: newey2 delay_gps density_norm, lag(14)
	qui estadd scalar my_r2=`myr2'
	qui sum density_norm
	qui estadd scalar xsd = r(sd)


**************
*** hour level
**************
*** AUX
	capture program drop nl_rename_cols
	program define nl_rename_cols, eclass
		// DAYLIGHT: {ld=0.3}*daylight + 
	    nl (`1' = {l0=2.0} + {l1=1.5} * (density_norm ^ {gamma=1.0})), vce(hac nwest 3)
	    matrix A = e(b)
	    matrix colnames A = "_cons" "density_norm" "gamma"
	    matrix coleq    A = ""
	    // estimates repost b = A
	    ereturn repost b = A, rename
	end


*** Column hour level data
	use "data/coded_road_tech/density_h", clear

	merge 1:1 hour using "data/coded_road_tech/google_maps_h_level"
	keep if _m==3
	drop _m

	merge 1:1 hour using "data/coded_road_tech/volumes_h_level", keepusing(delay_gps)
	assert _m==3
	drop _m	

	tsset hour

	*** Column 1 OLS
	qui reg delay_gm density_norm, r
	local myr2=e(r2_a)

	eststo hour_gm: newey2 delay_gm density_norm, lag(3)
	qui estadd scalar my_r2=`myr2'
	qui sum density_norm
	qui estadd scalar xsd = r(sd)

	*** Column 2 Non-linear least square
	// nl (delay_gm = {l0=2.0} + {l1=1.5} * (density_norm ^ {gamma=1.0})), vce(hac nwest 6)
	// di e(r2)

	nl_rename_cols delay_gm
	qui estadd scalar my_r2=e(r2_a)
	qui sum density_norm
	qui estadd scalar xsd = r(sd)

	qui estimates store hour_gm_nl

	*** Column 5 OLS
	qui reg delay_gps density_norm, r
	local myr2=e(r2_a)

	eststo hour_gps: newey2 delay_gps density_norm, lag(3)
	qui estadd scalar my_r2=`myr2'
	qui sum density_norm
	qui estadd scalar xsd = r(sd)

	*** Column 6 Non-linear least square
	nl_rename_cols delay_gps
	qui estadd scalar my_r2=e(r2_a)
	qui sum density_norm
	qui estadd scalar xsd = r(sd)

	qui estimates store hour_gps_nl

**************
*** date x hour level -> 2SLS regressions
**************
	use "data/coded_road_tech/density_th", clear

	merge 1:1 date hour using "data/coded_road_tech/volumes_th_level", keepusing(delay_gps)
	drop if _m==2
	drop _m	

	merge 1:1 date hour using "data/coded_road_tech/google_maps_th_level"
	keep if _m==3
	drop _m

	*** SAMPLE: only weekdays
	keep if inrange(dow(date),1,5)

	tab hour, gen(ih)

	xtset date hour

	*** check that HAC Newey-West give smaller SEs
	// ivreghdfe delay_gm (density_norm = ih*), ab(date) bw(3) robust

	**** 
	sum date
	gen date_relative = date - r(min)
	gen hour_day = date_relative * 24 + hour

	sort hour_day
	// line hour hour_day

	tsset hour_day


*** 1. CLUSTERING DATE x HOUR
	eststo tsls1: ivreg2 delay_gm i.date (density_norm = ih*), cluster(date)
	weakivtest
	estadd scalar eff_f = r(F_eff)
	estadd scalar crit = r(c_TSLS_5)
	estadd local ses = "cluster date"
	qui sum density_norm
	qui estadd scalar xsd = r(sd)

	*** GPS
		eststo tsls1_gps: ivreg2 delay_gp i.date (density_norm = ih*), cluster(date)
		// weakivtest
		// estadd scalar eff_f = r(F_eff)
		// estadd scalar crit = r(c_TSLS_5)
		estadd local ses = "cluster date"
		qui sum density_norm
		qui estadd scalar xsd = r(sd)

*** 2. CLUSTERING HOUR
	eststo tsls2: ivreg2 delay_gm i.date (density_norm = ih*), cluster(hour)
	weakivtest
	estadd scalar eff_f = r(F_eff)
	estadd scalar crit = r(c_TSLS_5)
	estadd local ses = "cluster hour"
	qui sum density_norm
	qui estadd scalar xsd = r(sd)

	*** GPS
		eststo tsls2_gps: ivreg2 delay_gps i.date (density_norm = ih*), cluster(hour)
		weakivtest
		estadd scalar eff_f = r(F_eff)
		estadd scalar crit = r(c_TSLS_5)
		estadd local ses = "cluster hour"
		qui sum density_norm
		qui estadd scalar xsd = r(sd)

*** 3. HAC (144 hours lag)
	eststo tsls3: ivreg2 delay_gm i.date (density_norm = ih*), robust bw(144)
	weakivtest
	estadd scalar eff_f = r(F_eff)
	estadd scalar crit = r(c_TSLS_5)
	estadd local ses = "HAC"
	qui sum density_norm
	qui estadd scalar xsd = r(sd)

	*** GPS
		eststo tsls3_gps: ivreg2 delay_gps i.date (density_norm = ih*), robust bw(144)
		weakivtest
		estadd scalar eff_f = r(F_eff)
		estadd scalar crit = r(c_TSLS_5)
		estadd local ses = "HAC"
		qui sum density_norm
		qui estadd scalar xsd = r(sd)

*** 4. HAC (144 hours lag) restrict to 6am-10pm
	eststo tsls4: ivreg2 delay_gm i.date (density_norm = ih7-ih23) if inrange(hour,6,22), robust bw(144)
	weakivtest
	estadd scalar eff_f = r(F_eff)
	estadd scalar crit = r(c_TSLS_5)
	estadd local ses = "HAC"
	estadd local time_sample = "6am-10pm"
	qui sum density_norm
	qui estadd scalar xsd = r(sd)

	*** GPS
		eststo tsls4_gps: ivreg2 delay_gps i.date (density_norm = ih7-ih23) if inrange(hour,6,22), robust bw(144)
		weakivtest
		estadd scalar eff_f = r(F_eff)
		estadd scalar crit = r(c_TSLS_5)
		estadd local ses = "HAC"
		estadd local time_sample = "6am-10pm"
		qui sum density_norm
		qui estadd scalar xsd = r(sd)


*** 4. HAC (144 hours lag) restrict to 8am-12pm
	eststo tsls5: ivreg2 delay_gm i.date (density_norm = ih9-ih13) if inrange(hour,8,12), robust bw(144)
	weakivtest
	estadd scalar eff_f = r(F_eff)
	estadd scalar crit = r(c_TSLS_5)
	estadd local ses = "HAC"
	estadd local time_sample = "8am-12pm"
	qui sum density_norm
	qui estadd scalar xsd = r(sd)

	*** GPS
		eststo tsls5_gps: ivreg2 delay_gps i.date (density_norm = ih9-ih13) if inrange(hour,8,12), robust bw(144)
		weakivtest
		estadd scalar eff_f = r(F_eff)
		estadd scalar crit = r(c_TSLS_5)
		estadd local ses = "HAC"
		estadd local time_sample = "8am-12pm"
		qui sum density_norm
		qui estadd scalar xsd = r(sd)



*** TABLE for response to the editor 
	esttab tsls1 tsls2 tsls5, se keep(density_norm _cons) ///
		nomtitles ///
		coeflabels(density_norm "Traffic Density" _cons "Constant") ///
		stats(time_sample ses eff_f crit N, ///
		labels("Time sample"  "Standard errors" "Effective F-stat" "Critical value for 5%" "Observations" ) ///
		fmt("%s" "%s" "%12.1f" "%12.1f" "%12.0fc")) varwidth(23)


************************
*** SAVE MAIN TABLE ****
************************

	// local dec = 2
	local outputfile "`=outputfolder_table'rev2_panel_A_gm"
	esttab tsls3 tsls4 hour_gm hour_gm_nl date_gm using "`outputfile'.tex", replace ///
			b(%12.2f) se(%12.2f) nomtitles nonumbers ///
			keep(density_norm gamma _cons) ///
			coeflabels( ///
					density_norm "Traffic Density" ///
					gamma "Traffic Density Exponent \$\nu\$" ///
					_cons "Constant") ///
			starlevels(* .10 ** .05 *** .01) ///
			stats(N eff_f crit xsd my_r2, ///
					labels("Observations"  ///
						   "Effective F-Stat"  ///
						   "Critical value 5\%"  ///
						   "Traffic Density Std. Dev." ///
						   "$ Adj. R^2 $") ///
					fmt("%12.0fc" "%4.1f" "%4.1f" "%12.2f" "%12.2f")) booktabs ///
			substitute(_ \_ "<" "$<$" "\midrule" "\addlinespace\addlinespace") ///
			fragment
	
*** Write main table tex file
	file open myfile using "`=outputfolder_table'rev2_table_`=tnum'_gm.tex", write replace
	file write myfile "\begin{tabular}{lccccc}" _n "\toprule" _n  ///
					  " & (1) & (2) & (3) & (4) & (5) \\ \addlinespace" _n ///
					  "\emph{Dependent Variable:} & \multicolumn{5}{c}{Travel Delay Google Maps (min/km)} \\ \addlinespace" _n ///
					  "\emph{Sample:} & \multicolumn{2}{c}{Date $\times$ Dep Time ($ th $)} &" ///
					  		" \multicolumn{2}{c}{Dep Hour ($ h $)} & Date ($ t $) \\ \addlinespace" _n ///
					  	" & & 6am-10pm & & & \\ \addlinespace" _n ///
					  "\emph{Specification:} & 2SLS & 2SLS & OLS & NLS & OLS \\ \addlinespace" _n ///
				   	  "\ExpandableInput{tables/`=tnum'/rev2_panel_A_gm}" _n ///
					  "\bottomrule" _n "\end{tabular}" _n  
	file close myfile



**********************************
**** SAVE TABLE WITH GPS DATA ****
**********************************
	local outputfile "`=outputfolder_table'rev2_panel_A_gps"
	esttab tsls3_gps tsls4_gps hour_gps hour_gps_nl date_gps using "`outputfile'.tex", replace ///
			b(%12.2f) se(%12.2f) nomtitles nonumbers ///
			keep(density_norm gamma _cons) ///
			coeflabels( ///
					density_norm "Traffic Density" ///
					gamma "Traffic Density Exponent \$\nu\$" ///
					_cons "Constant") ///
			starlevels(* .10 ** .05 *** .01) ///
			stats(N eff_f crit xsd my_r2, ///
					labels("Observations"  ///
						   "Effective F-Stat"  ///
						   "Critical value 5\%"  ///
						   "Traffic Density Std. Dev." ///
						   "$ Adj. R^2 $") ///
					fmt("%12.0fc" "%4.1f" "%4.1f" "%12.2f" "%12.2f")) booktabs ///
			substitute(_ \_ "<" "$<$" "\midrule" "\addlinespace\addlinespace") ///
			fragment


	file open myfile using "`=outputfolder_table'rev2_table_`=tnum'_gps.tex", write replace
	file write myfile "\begin{tabular}{lccccc}" _n "\toprule" _n  ///
					  " & (1) & (2) & (3) & (4) & (5) \\ \addlinespace" _n ///
					  "\emph{Dependent Variable:} & \multicolumn{5}{c}{Travel Delay GPS Data (min/km)} \\ \addlinespace" _n ///
					  "\emph{Sample:} & \multicolumn{2}{c}{Date $\times$ Dep Time ($ th $)} &" ///
					  		" \multicolumn{2}{c}{Dep Hour ($ h $)} & Date ($ t $) \\ \addlinespace" _n ///
					  	" & & 6am-10pm & & & \\ \addlinespace" _n ///
					  "\emph{Specification:} & 2SLS & 2SLS & OLS & NLS & OLS \\ \addlinespace" _n ///
				   	  "\ExpandableInput{tables/`=tnum'/rev2_panel_A_gps}" _n ///
					  "\bottomrule" _n "\end{tabular}" _n  
	file close myfile
