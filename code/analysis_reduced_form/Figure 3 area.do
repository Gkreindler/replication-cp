
 * This do file: analysis for main figure (figure 3) analysis of route charges
 * This version: Nov 23, 2021
 * Authors: Gabriel Kreindler

clear all
pause off
set more off
set matsize 10000

*** TABLE number
// local tfolname "sandbox"
// local outputfolder "$cp_git_path/paper/tables/`tfolname'/"

*** FIGURE number
	local fnum "figure3"
	local outputfolder "paper/figures/`fnum'/"
	cap mkdir "`outputfolder'"

* Path based on global set in profile.do
	local path "$congestion_pricing_path"
	cd "`path'"



*** LOAD 
	use "data\coded_gps_dta\area_coded", clear


************************
**** ANALYSIS GRAPH ****
************************

	forv i=4(-1)0{
		gen sc_`i'_early = (study_cycle == `i') * a_early
		gen sc_`i'_late  = (study_cycle == `i') * a_late
	}


/*
	bys uidp: gegen has_0=max(study_cycle==0)
	bys uidp: gegen has_1=max(study_cycle==1)
	bys uidp: gegen has_2=max(study_cycle==2)
	bys uidp: gegen has_3=max(study_cycle==3)

	gen mysample = has_0 + has_1 + has_2 + has_3
	tab mysample u1

	keep if mysample == 4 
*/
	
	gunique uidp 

*** collapse to mean detour route rate per uidp -- 
	gcollapse (mean) is_long_route, by(uidp study_cycle sc_* a_early a_late)


*** run regressions with commuter FEs
	eststo r_early: reghdfe is_long_route sc_?_early if a_early == 1, ab(uidp) vce(cl uidp)
	sum is_long_route if sc_0_early == 1 
	scalar control_mean_early = r(mean)

	eststo r_late: reghdfe is_long_route sc_?_late if a_late == 1, ab(uidp) vce(cl uidp)
	sum is_long_route if sc_0_late == 1 
	scalar control_mean_late = r(mean)


*** Plot!
	coefplot (r_early, omit drop(_cons) offset(-0.05) ///
						transform(* = (@+`=control_mean_early')) ///
						mcolor(blue) msize(medlarge)  ciopts(color(gray) lwidth(medthick))) ///
			 (r_late , omit drop(_cons) offset(+0.05) ///
			 			transform(* = (@+`=control_mean_late')) ///
			 			mcolor(red) msize(medlarge) msymbol(square) ciopts(color(gray) lwidth(medthick))) ///
			, vertical  ///
			order(sc_0_early sc_0_late ///
					 sc_1_early sc_1_late ///
					 sc_2_early sc_2_late ///
					 sc_3_early sc_3_late ///
					 sc_4_early sc_4_late) ///
			coeflabels(sc_0_* = `" "before" "experiment" "' ///
					 sc_1_* = "week 1" ///
					 sc_2_* = "week 2" ///
					 sc_3_* = "week 3" ///
					 sc_4_* = "week 4") 	///
			relocate(sc_0_early=0 ///
					 sc_1_early=1 ///
					 sc_2_early=2 ///
					 sc_3_early=3 ///
					 sc_4_early=4 ///
					 sc_0_late =0 ///
					 sc_1_late =1 ///
					 sc_2_late =2 ///
					 sc_3_late =3 ///
					 sc_4_late =4) ///
			legend(order(2 "Charges in week 1" 4 "Charges in week 4") pos(6) cols(2)) ///
			graphregion(color(white)) ///
			ytitle("Fraction" "detour" "route", orientation(horizontal)) ///
			ylabel(0(0.1)0.5, angle(0)) ///
			yline(`=control_mean_early', lcolor(blue%40) lpattern(dash) lwidth(0.5)) ///
			yline(`=control_mean_late', lcolor(red%40) lpattern(dash) lwidth(0.5)) 


*** Save
	local outputfile "`outputfolder'figure_vott_byweek"
	graph export "`outputfile'.eps", replace
	graph export "`outputfile'.png", replace
	graph export "`outputfile'.pdf", replace
