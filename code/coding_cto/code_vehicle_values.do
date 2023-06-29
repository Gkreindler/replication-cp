* This do file: Clean vehicle values
* This version: 15 August 2018
* Authors: Gabriel Kreindler

clear all
pause off
set more off

* Path based on global set in C:\ado\profile.do
	local path "$bang_launch1"
	cd "`path'"

** Load recruitment
	use "data/raw_other/vehicle_values/recruit_with_dict.dta", clear


** Merge using 3 vars
	merge m:1 standvehtype standbrand standmod using "data/raw_other/vehicle_values/olx_stats_brand_model.dta"
	// br uidp_original standvehtype standbrand standmod this_car_or_moto this_vehicle_model if _m==1 & standbrand!="999" & standmod!="999"
	drop if _m==2
	count if _m==1 &  standbrand!="999" & standmod!="999"
	assert r(N)==296
	drop _m

** Merge using 2 vars
	merge m:1 standvehtype standbrand using "data/raw_other/vehicle_values/olx_stats_brand.dta"
	// br uidp_original standvehtype standbrand standmod this_car_or_moto this_vehicle_model if _m==1 & standbrand!="999"
	drop if _m==2
	count if _m==1 &  standbrand!="999"
	assert r(N)==34
	drop _m

*** Generate Vehicle Value prediction
	gen inv_std_bm = 1/sdlp_bm
	gen inv_std_b  = 1/sdlp_b
	recode inv_std_bm inv_std_b (.=0)

	gen 	logprice = (avglp_bm * inv_std_bm + avglp_b * inv_std_b) / (inv_std_bm + inv_std_b)

	* looks good, except a suspicious overlap between moto and scooter
	/* twoway 	(kdensity logprice if standvehtype == "Car (4 wheeler)", lcolor(red)) ///
			(kdensity logprice if standvehtype == "Moto (2 wheeler)", lcolor(blue)) ///
			(kdensity logprice if standvehtype == "Scooter (2 wheeler)", lcolor(black)) ///
			, legend(order(1 " car" 2 "moto" 3 "scooter")) */
    
*** reliability 

*** save
	keep uidp_original logprice  standvehtype standbrand standmod
	save "data/coded_cto/vehicle_values_prices.dta", replace

