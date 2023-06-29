
 * This do file: Prepares for demand estimation
 * This version: Oct 11, 2017
 * Changes: proper day level number of trips per day by bin
 * Authors: Gabriel Kreindler

** RESTRICTED TO AM **

 /* Plan:
at trip level (?)
- departure time  			OK
- intersect with area 		OK

at ID level
- treatment status  		OK
-- DT treat 				OK
-- Area treat 				OK
-- Area control 			OK
- peak value 				OK	
- area detour 				OK
- typical KM 				OK
- normal mu and sigma for h_arrival

at Google Map level
- travel time at each HD, one way

  */

clear all
pause on
set more off

* Path based on global set in profile.do
	local path "$congestion_pricing_path"
	cd "`path'"
	adopath ++ "code/ado/"

************************
*** TRIPS micro data ***
************************
*** Visit H/W/W2
	use "data/coded_gps_dta/coded_days_15.dta", clear

	drop if inlist(study_cycle,9,99) // drop post-experiment
	keep if qual_gb == 1

* checks
	assert sample_analysis == qual_gb
	assert date != mdy(4,26,2017)
	assert qual_good + qual_bad + qual_halfday + qual_no_data == 1
	isid uidp date
	sort uidp date

************************
**** OVERALL SAMPLE ****
************************
	gunique uidp
*** Regular commuters only
	keep if inrange(regular_commuter,2,3)
	gunique uidp

*** days with trips only
	local variant = "v1"

	* V1 - ever at w
	*    - use cha0_all_
	if ("`variant'" == "v1") keep if at_w_all == 1	

	* V2 - nt0_all_hwp_am > 0
	*    - use cha0_all_am
	if ("`variant'" == "v2") keep if nt0_all_hwp_am > 0

*** AREA ONLY
	keep if in_area_treatment == 1

*** SAMPLE drop 11 respondents with faulty AREA locations
	// gunique uidp if sample_area_drop == 1 
	drop if sample_area_drop == 1 

*** AM
	// keep if am_sample == 1 
	// unique uidp

*** H-W-P trips
	// keep if trip_hwp == 1

****************
**** CODING ****
****************
	* V1 - ever at w
	*    - use cha0_all_
	if ("`variant'" == "v1"){
		tab cha0_all_
		local myvar cha0_all_
	}

	* V2 - nt0_all_hwp_am > 0
	*    - use cha0_all_am
	if ("`variant'" == "v2"){
		tab cha0_all_am
		local myvar cha0_all_am
	}

	* definition of the area control group using weeks 1 and 4 of the post period
	gen area_control = inlist(study_cycle,1,4) & area_treat == 0

	sort uidp
	* count number of post
	by uidp: egen npost=sum(post)

	* count number of treatment
	by uidp: egen ntr=sum(area_treat)
	
	* count number of control in weeks 1 and 4 (post)
	by uidp: egen nct=sum(area_control)

	* count number of control in non-treatment weeks 1 - 4 (so includes week 2 and 3 always)
	by uidp: egen nct23=sum((1 - area_treat)*post)

	* dummies for positive
	gen npost0 = npost > 0
	gen ntr0   = ntr > 0
	gen nct0   = nct > 0
	gen nct230 = nct23 > 0

	gen n_ct_tr = nct0 & ntr0
	gen n_ct23_tr = nct230 & ntr0


* average area intersect PRE
	by uidp: egen cha_pre = mean(`myvar') if post == 0

* average area intersect POST non-treated (include weeks 2,3)
	by uidp: egen cha_ct23 = mean(`myvar') if post == 1 & area_treat == 0 

* average area intersect POST control (only control week 1,4)
	by uidp: egen cha_ct = mean(`myvar') if post == 1 & area_control == 1

* average area intersect POST treated
	by uidp: egen cha_tr = mean(`myvar') if post == 1 & area_treat == 1

	qui fill_by_group cha_pre cha_ct23 cha_ct cha_tr, fillby(uidp) replace 

/* *** Interlude: check that most weekly observations are all main or all detour (65%)
	keep if inrange(study_cycle,1,4)
	gen nobs=1
	gcollapse (mean) cha0_all_am (sum) nobs, by(study_cycle uidp)
	tab nobs
	tab cha0_all_am
***  */

	cap drop o1
	bys uidp: gen o1=_n==1
	sum npost nct23 nct ntr if o1==1
	sum npost0 nct230 nct0 ntr0 if o1==1
	sum cha_pre cha_ct23 cha_ct cha_tr if o1==1

	sum nct23 if cha_ct23 ==1 & cha_tr==1 & o1==1
	sum   ntr if cha_ct23 ==1 & cha_tr==1 & o1==1

	tab n_ct_tr if o1==1
	tab n_ct23_tr if o1==1
	tab nct230 ntr0 if o1==1

	keep if o1==1
	keep  uidp nct23 nct ntr cha_ct23 cha_ct cha_tr
	order uidp nct23 nct ntr cha_ct23 cha_ct cha_tr

	*** Drop 21 observations without any post data
	drop if nct23 + nct + ntr == 0

	gisid uidp
	sort uidp
	save  "data/coded_model/area_het.dta", replace

	// twoway	(hist cha_tr  , start(0) w(0.1) fcolor(gs12) lcolor(black) lwidth(medthick)) ///
	// 		(hist cha_ct23, start(0) w(0.1) fcolor(none) lcolor(red) lwidth(medthick)) ///
	// 		, legend(order(1 "treat" 2 "ctrl"))

