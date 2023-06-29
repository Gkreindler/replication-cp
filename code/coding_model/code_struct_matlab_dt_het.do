
 * This do file: Create multiple samples (including for structural) and compare RF results 
 * This version: June 12, 2019
 * Authors: Gabriel Kreindler

clear all
pause off
set more off
set matsize 10000


* Path based on global set in profile.do
	local path "$congestion_pricing_path"
	cd "`path'"
	adopath ++ "code/ado/"

*** LOAD 
	use "data/coded_gps_dta/coded_trips_15.dta", clear

****************
**** SAMPLE ****
****************

	*** Drop dates after the study ends (=9) and during the study but not in the experiment (=99)
	assert inlist(study_cycle,0,1,2,3,4,5,9,99)
	drop if inlist(study_cycle,9,99) // drop post-experiment

	* Drop DAYS without trips
	assert trip == 1 - date_wo_trips
	keep if trip == 1

	* sample if all good+bad days
	assert sample_analysis == qual_gb

	* Check April 26 is dropped, which was a holiday and wrongly included in participants' timelines
	assert date != mdy(4,26,2017)

	* check gunique ids
	gisid uidp date chain
	sort uidp date chain

	*** CODING 
	assert qual_good + qual_bad + qual_halfday + qual_no_data == 1
	assert regular_commuter != .

	* AM early:
	// gen am_e_sample = am_pk1_sample == 1 & t_start < pk_am
	// gen pm_l_sample = pm_pk1_sample == 1 & t_start > pk_pm

********************
**** DT SAMPLES ****
********************

	* sample_trip_ok is a definiton at the trip level for a reliably measured trip 
	* sample_analysis == qual_gb is a day level indicator for good or "bad" (=medium)

		gen dt_sample_full = 	sample_analysis == 1 & ///
								sample_trip_ok == 1 & ///
								am_sample == 1 & ///
								sample_3w == 1

	*** Table 3 regular commuters sample
		gen dt_sample_hwp = 	dt_sample_full & ///	
								inlist(regular_commuter,2,3) & ///
								trip_hwp == 1

	*** Sample for structural estimation
		gen dt_sample_struct = 	dt_sample_full & ///	
								inlist(regular_commuter,2,3) & ///
								trip_hwp == 1 & (oh==1 & dw==1) & ///
								inrange(t_rel, -2.5, 2.5)

		gunique uidp if dt_sample_struct == 1

		by uidp: gegen n_pos = sum(post * dt_sample_struct)
		by uidp: gegen n_pre = sum((1-post) * dt_sample_struct)

		// replace dt_sample_struct = dt_sample_struct & ///
		// 						   n_pos > 1 & n_pre > 1

		replace dt_sample_struct = dt_sample_struct & n_pos > 0 & n_pre > 0

		gunique uidp if dt_sample_struct == 1

		drop n_pre n_pos

	*** Sample for structural estimation -- Pre-peak sample (used for heterogeneity)
		gen am_pre_sample = am_sample == 1 & inrange(t_rel, -2, 0)

		cap drop temp
		bys uidp: gegen temp = sum((1-post) * dt_sample_struct * am_pre_sample)
		bys uidp: gegen idfxnpre = mean(temp)
		drop temp
		bys uidp: gegen temp = sum(   post  * dt_sample_struct * am_pre_sample)
		bys uidp: gegen idfxnpos = mean(temp)
		drop temp am_pre_sample

		* ids with at least one pre-peak observation both pre- and post- experiment
		gen dt_sample_struct_pre = dt_sample_struct & idfxnpre >=1 & idfxnpos >= 1
		drop idfxnp*

	* sorting again
		sort uidp date chain

		assert date_wo_trips == 0

		tempfile alldata
		save 	`alldata'

*************************************************
*** Code and save for struct
*************************************************

	*** SAMPLE
		keep if dt_sample_struct == 1

	*** Compute departure time mean and SD in pre period
		by uidp: gegen n_pos = sum(post)
		by uidp: gegen n_pre = sum(1-post)

		cap drop temp
		by uidp: egen temp = mean(t_rel) if post == 0
		by uidp: egen norm_mean = mean(temp) // extend to all obs

		cap drop temp
		by uidp: egen temp = sd(t_rel) if post == 0
		by uidp: egen norm_sd = mean(temp) // extend to all obs
		drop temp

	*** fill in SD with median if missing
		by uidp: gen o1=_n==1
		sum norm_sd if o1==1, d
		replace norm_sd = r(p50) if norm_sd == .
		sum norm_mean if o1==1, d 
		replace norm_mean = r(p50) if norm_mean == .
		drop o1


	********************************************************************
	*** Compute average daily number of trips in each bin, pre/post  ***
	********************************************************************
		isid uidp date chain
		sort uidp date chain

	*** DT by departure time bin
		scalar nbins = 60
		scalar h_min = -2.5
		scalar h_max = +2.5

		scalar binwidth = (h_max - h_min) / nbins
		scalar binwidth_half = binwidth / 2

	*** Generate number of trips by date and bin
		sort uidp date chain
		forv i=1/`=nbins+1'{
			local trel_min = h_min + (`i'-1) * binwidth - binwidth_half
			local trel_max = h_min +  `i'    * binwidth - binwidth_half

			gen bin`i' = `trel_min' <= t_rel & t_rel < `trel_max'
			qui by uidp date: egen ntrips`i' = sum(bin`i')
			replace ntrips`i' = 100 * ntrips`i'
			label var ntrips`i' "Number of trips in bin `i' (`trel_min'-`trel_max') today"
			
		} // t_rel
		
	*** Average by pre/post for each person
	preserve
		sort uidp date chain
		by uidp date: gen od1 = _n==1
		keep if od1==1
		
		forv i=1/`=nbins+1'{
			qui by uidp: egen mean_dt_pre`i' = mean(ntrips`i') if post == 0
			qui by uidp: egen mean_dt_pos`i' = mean(ntrips`i') if post == 1
			label var mean_dt_pre`i' "Mean number of trips in bin `i' per day - PRE"
			label var mean_dt_pos`i' "Mean number of trips in bin `i' per day - POST"

			qui fill_by_group mean_dt_pre`i' mean_dt_pos`i', fillby(uidp) replace
		}

		keep uidp mean_dt_pre* mean_dt_pos*
		duplicates drop

		tempfile dt_ntrips_prepost
		save 	`dt_ntrips_prepost'
	restore

		merge m:1 uidp using `dt_ntrips_prepost'
		assert _m==3
		drop _m


	*************************************
	*** Individual changes in charges *** (AM, 2h pre peak sample)  N = 166
	*************************************

	preserve
		gen am_pre_sample = am_sample == 1 & inrange(t_rel, -2, 0)
		keep if am_pre_sample == 1

		bys uidp: gegen idfxnpre = sum(1-post)
		bys uidp: gegen idfxnpos = sum(post)

		gunique uidp if idfxnpre >= 1 & idfxnpos >= 1
		keep if idfxnpre >= 1 & idfxnpos >= 1
			
		statsby _b _se n=e(N), by(uidp dt_wcha idfxnpre idfxnpos) clear: reg c_dt_rate post, r
		
		* drop uidpns where the data is constant (no real SE)
		drop if _se_post == 0

		mdesc

		// drop if _b_post < -90 // one complete outlier "L0221174033" (treatment group, so this is conservative)
		assert _b_post != .
		keep if _b_post != .

		assert _eq2_n == idfxnpre + idfxnpos

		rename _b_post _b_post_old
		// ebayes _b_post_old _se_post, gen(_b_post) by(dt_wcha) tol(0.00001)
		ebayes _b_post_old _se_post dt_wcha, gen(_b_post) tol(0.00001)

		reg _b_post dt_wcha, r
		reg _b_post_old dt_wcha, r

		sum _b_post if dt_wcha==0
		sum _b_post if dt_wcha==1

		sum _b_post_old if dt_wcha==0
		sum _b_post_old if dt_wcha==1

		/* twoway 	(kdensity _b_post_old if dt_wcha == 0, lcolor(gs10)  lwidth(medthick) bw(6)) ///
				(kdensity _b_post_old if dt_wcha == 1, lcolor(black) lwidth(medthick) bw(6)) ///
				, legend(off) ///
				graphregion(color(white)) xtitle("Individual Change in Charges (Post minus Pre)") xlabel(-100(25)50, grid) ///
				ytitle("") ylabel(none)   scale(1.2) name(var_old, replace)

		twoway 	(kdensity _b_post if dt_wcha == 0, lcolor(gs10)  lwidth(medthick) bw(6)) ///
				(kdensity _b_post if dt_wcha == 1, lcolor(black) lwidth(medthick) bw(6)) ///
				, legend(off) ///
				graphregion(color(white)) xtitle("Individual Change in Charges (Post minus Pre)") xlabel(-100(25)50, grid) ///
				ytitle("") ylabel(none)   scale(1.2) name(var_eb, replace)

		gr combine var_old var_eb, col(1) iscale(1)  */

		mdesc
		keep uidp _b_post _b_post_old idfxnpre idfxnpos
		rename _b_post idfxbe
		rename _b_post_old idfx

		tempfile idfx
		save 	`idfx'

		gunique uidp

	restore

		merge m:1 uidp using `idfx'
		assert _m!=2
		drop _m
		gunique uidp

		recode idfxbe idfx idfxnpre idfxnpos (.=0)

	********************
	*** Final coding ***
	********************
		keeporder uidp norm_mean norm_sd mean_dt_pre* mean_dt_pos* idfxbe idfx idfxnpre idfxnpos n_pre n_pos
		gduplicates drop
		gisid uidp

		mdesc mean_dt_*
		gen sample_dt_pre = mean_dt_pre1 != .
		gen sample_dt_pos = mean_dt_pos1 != .
		forv i=1/`=nbins'{
			assert mean_dt_pre`i' != . if sample_dt_pre==1
			assert mean_dt_pos`i' != . if sample_dt_pos==1
		}

		save "data/coded_model/dt_het.dta", replace

	clear
	// */ 
