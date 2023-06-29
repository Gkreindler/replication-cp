
 * This do file: GMM numerical identification
 * This version: Jan 16 2022 // Mar 19 2020
 * Authors: Gabriel Kreindler

clear all
pause on
set more off

*** TABLE number
local tname "table"
local tfolname "smtable12"
local outputfolder "paper/tables/`tfolname'/"
cap mkdir "`outputfolder'"

* Path based on global set in profile.do
	local path "$congestion_pricing_path"
	cd "`path'"

*** Load main estimation parameter values
	import delimited using "data\coded_model\simdata\all_results_fin.csv", clear
	
	sum x1
	scalar a_true = r(mean)
	
	sum x2
	scalar be_true = r(mean)

	sum x3
	scalar bl_true = r(mean)

	sum x4
	scalar gamma_true = r(mean)
	
	sum x5
	scalar sig_dt_true = r(mean)
	
	sum x6
	scalar sig_r_true = r(mean)


*** Regressions
foreach va of varlist param_1 param_4 param_5 param_2 param_6 param_7{ 
	qreg `va' x1 x4 x5 x2 x6 x7, vce(r)
	qui estimates store ests_`va'
}

*** Output
	local dec=2
	esttab ests_* using "`outputfolder'`tname'_panel_A.tex", replace ///
				b(%12.`dec'f) se(%12.`dec'f) keep(x1 x4 x5 x2 x6 x7) ///
				coeflabels(x1 "(True) Value of time $\alpha$" ///
						  x4 "(True) Penalty early $\beta_E$ " ///
						  x5 "(True) Penalty late $\beta_L$ " ///
						  x2 "(True) Switch Cost $\gamma$" ///
						  x6 "(True) Logit departure time $\sigma^{DT}$" ///
						  x7 "(True) Logit route $\sigma^{R}$") ///
				starlevels(* .10 ** .05 *** .01) nomtitles nonumbers ///
				stats(N, labels("Observations" ) fmt("%12.0fc")) booktabs ///
				substitute("<" "$<$" "\midrule" "\addlinespace") ///
				fragment 

	estimates clear

*** Write main table tex file
	file open myfile using "`outputfolder'`tname'.tex", write replace
	file write myfile "\begin{tabular}{lcccccc}" _n "\toprule" _n  ///
					  " & (1) & (2) & (3) & (4) & (5) & (6) \\" _n ///
					  "\addlinespace" _n ///
					  " & \multicolumn{6}{c}{\emph{Estimated Parameter}} \\" _n ///
					  "\addlinespace" _n ///
					  " & $\widehat{\alpha}$  & $\widehat{\beta}_E$ & $\widehat{\beta}_L$ & $\widehat{\gamma}$  & $\widehat{\sigma^{DT}}$  & $\widehat{\sigma^{R}}$ \\" _n ///
				   	  "\ExpandableInput{tables/`tfolname'/`tname'_panel_A}" _n ///
					  "\bottomrule" _n "\end{tabular}" _n ///

	file close myfile
