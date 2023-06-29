
 * This do file: Loads and code area candidates (congesiton areas that were not chosen ultimately)
 * This version: 8 August, 2017
 * Authors: Gabriel Kreindler

clear all
pause off
set more off

* Path based on global set in profile.do

	local path "$congestion_pricing_path"
	cd "`path'"


*** Get final predicted detour, base, etc
	use "data/treatment/treatment roster noPII.dta"
	keep if meeting == "done"
	keep if atreat != "0 No treatment"
	keep uidp adetourtime aradius 
	rename adetourtime 	ac_detour_chosen
	rename aradius 		ac_area_radius
	
	tempfile chosen
	save 	`chosen'


**************
*** Load trips
	import delimited using "data/treatment/all_area_candidates.csv", clear
	isid uidp area_radius idx

	* describe missing values
	mdesc

	* a bit of cleaning
	count if valid == 1 & valid_long == 0 & valid_short == 0
	assert r(N) == 7
	replace valid_short = 1 if valid == 1 & valid_long == 0 & valid_short == 0 & inrange(detour_time,1,3)
	replace valid_long  = 1 if valid == 1 & valid_long == 0 & valid_short == 0 & inrange(detour_time,14,15)
	assert valid_long + valid_short == 1 if valid == 1

	* check
	gen diff = (avoid_dur - base_dur) / 60 - detour_time
	sum diff if detour_time != -1, d
	assert inrange(diff,-0.1,0.1) if detour_time != -1

	foreach va of varlist base_dur avoid_dur detour_time detour_dist{
		replace `va' = . if valid == 0
		gen `va'_short = `va' if valid_short == 1
		gen `va'_long  = `va' if valid_long == 1
	}

	gcollapse (sum) valid_n=valid valid_short_n=valid_short valid_long_n=valid_long ///
			 (mean) valid_frac=valid ///
			 		base_t_mean=base_dur 		base_t_mean_sh=base_dur_short 		base_t_mean_lg=base_dur_long ///
			 		avoid_t_mean=avoid_dur 		avoid_t_mean_sh=avoid_dur_short 	avoid_t_mean_lg=avoid_dur_long ///
			 		detour_t_mean=detour_time 	detour_t_mean_sh=detour_time_short	detour_t_mean_lg=detour_time_long ///
			 		detour_d_mean=detour_dist 	detour_d_mean_sh=detour_dist_short	detour_d_mean_lg=detour_dist_long ///
			 (sd) 	base_t_sd=base_dur 			base_t_sd_sh=base_dur_short 		base_t_sd_lg=base_dur_long ///
			 		avoid_t_sd=avoid_dur 		avoid_t_sd_sh=avoid_dur_short 		avoid_t_sd_lg=avoid_dur_long ///
			 		detour_t_sd=detour_time 	detour_t_sd_sh=detour_time_short	detour_t_sd_lg=detour_time_long ///
			 		detour_d_sd=detour_dist 	detour_d_sd_sh=detour_dist_short	detour_d_sd_lg=detour_dist_long ///
			 , by(uidp)

	* indicators for at least one candidate in short
	gen valid_short = valid_short_n > 0
	gen valid_long  = valid_long_n > 0

	foreach va of varlist _all{
		if "`va'" != "uidp"{
			rename `va' ac_`va'
		}
	}

* merge chosen values
	merge 1:1 uidp using `chosen'
	assert _m==3
	drop _m

* save
	save "data/coded_gps_dta/area_candidates.dta", replace
