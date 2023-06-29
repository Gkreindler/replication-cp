* This do file: Codes and merges auxiliary data to recruitment data
 * This version: Aug 15 2018
 * Authors: Keerthana Jagadeesh, Gabriel Kreindler

clear all
pause off
set more off

* Path based on global set in C:\ado\profile.do
	local path "$congestion_pricing_path"
	cd "`path'"

	
* Load the output files from the above do files
   * occupation output
   // use "data/recruitment_temp/recruit_occupation.dta", clear
   use "data/coded_cto/occupations_clean.dta", clear 
   tempfile occupation
   save `occupation'
   
   * petrol_pump output
   // use "data/recruitment_temp/recruit_petrol pump name.dta"
   use "data/coded_cto/pumps_clean.dta", clear
   tempfile petrol_pump_name
   save `petrol_pump_name'
   
   * vehicle_prices output
   use "data/coded_cto/vehicle_values_prices.dta", clear
   tempfile vehicle_prices
   save `vehicle_prices'
   


**** MERGE ALL TOGETHER + final coding
* Load recruitment pooled data
   use "data/coded_cto/recruitment pooled.dta"
   
* Merge recruitment pooled with the above saved datasets to add the new columns
   
   merge m:1 uidp_original using `occupation'
   assert _m == 3
   drop _m
   
   merge m:1 uidp_original using `petrol_pump_name'
   assert _m == 3
   drop _m
   
   merge m:1 uidp_original using `vehicle_prices'
   //assert _m == 3 // 11,450 contradictions because there is no info for this_car_or_moto this_vehicle_brand this_vehicle_brand
   count if _m != 3 & this_car_or_moto != .
   assert r(N) == 10

   count if _m==2
   assert r(N) == 3
   drop if _m==2

   drop _m
   rename logprice vehicle_logprice

* Coding
    *replacing "." for 99 and 88
   destring income_keep, generate(income_keep2)
   recode income_keep2  (99=.) (88=.)
   
   replace income_keep2 =   5000 if income_keep2 == 0 
   replace income_keep2 =  12500 if income_keep2 == 1 
   replace income_keep2 =  17500 if income_keep2 == 1.5 
   replace income_keep2 =  22500 if income_keep2 == 2 
   replace income_keep2 =  27500 if income_keep2 == 2.5 
   replace income_keep2 =  35000 if income_keep2 == 3 
   replace income_keep2 =  45000 if income_keep2 == 4
   replace income_keep2 =  62500 if income_keep2 == 5
   replace income_keep2 =  87500 if income_keep2 == 7.5
   replace income_keep2 = 175000 if income_keep2 == 10
   replace income_keep2 = 375000 if income_keep2 == 25
   replace income_keep2 = 500000 if income_keep2 == 50

   * Logging income
   generate log_income = ln(income_keep2)

   * Car dummy
   encode   standvehtype, gen(standvehtype_e)
   gen   stand_car = standvehtype == "Car (4 wheeler)" if standvehtype != "999"
   replace stand_car = 1 if stand_car == . & this_car_or_moto == 1
   replace stand_car = 1 if stand_car == . & inlist(this_car_or_moto,2,3)
   replace stand_car = . if this_car_or_moto == .

   * Age
    gen new_age = 24 if agebracket == 1
    replace new_age = 35 if agebracket == 2
    replace new_age = 45 if agebracket == 3
    replace new_age = 55 if agebracket == 4
    replace new_age = 60 if agebracket == 5

    * KM reading - car age
   gen km_reading_clean = km_reading
   sum km_reading if km_reading > 99, d
   replace km_reading_clean = . if ~inrange(km_reading_clean, r(p1), r(p99))

**** NEW (3/10/2020) clean and extrapolate income
  gen income_outlier = !missing(income_keep2) & income_keep2 > 100000

  * Regression model 1 -- no vehicle value
  preserve 
    reghdfe log_income agebracket if income_outlier == 0, ab(occ_medium, savefe) vce(r) resid
    bys occ_medium: gegen temp = mean(__hdfe1__)
    replace __hdfe1__ = temp
    predict log_income_predicted_1, xb
    replace log_income_predicted_1 = log_income_predicted_1 + __hdfe1__
    keep uidp_original log_income_predicted_1
    tempfile fit_1
    save    `fit_1'
  restore

  * Regression model 2 -- no occupation
  preserve 
    reghdfe log_income vehicle_logprice agebracket if income_outlier == 0, ab(income_outlier) vce(r) resid
    predict log_income_predicted_2, xb
    keep uidp_original log_income_predicted_2
    tempfile fit_2
    save    `fit_2'
  restore

  merge 1:1 uidp_original using `fit_1'
  assert _m==3
  drop _m

  merge 1:1 uidp_original using `fit_2'
  assert _m==3
  drop _m

* Regression model 2
  reghdfe log_income vehicle_logprice agebracket if income_outlier == 0, ab(occ_medium, savefe) vce(r) resid
  bys occ_medium: gegen temp = mean(__hdfe1__)
  replace __hdfe1__ = temp
  predict log_income_predicted, xb

  replace log_income_predicted = log_income_predicted + __hdfe1__
  replace log_income_predicted = log_income_predicted_1 if missing(log_income_predicted) // 11 obs
  replace log_income_predicted = log_income_predicted_2 if missing(log_income_predicted) // 10 obs

  * final version -- mostly overlap with data, fill in where missing or outlier
  gen     log_income_with_pred = log_income_predicted
  replace log_income_with_pred = log_income if income_outlier == 0 & !missing(income_keep2)

  drop _reghdfe_resid __hdfe1__ log_income_predicted_1 log_income_predicted_2 temp

  label variable income_outlier          "Self-reported income above 100,000 INR"
  label variable log_income              "Self-reported income (log)"
  label variable log_income_with_pred    "Self-reported income (log) with fitted vals for outliers and missing"
  label variable log_income_predicted    "Self-reported income (log) fitted vals"

   * save new recruitment coded dta that contains the new variables from the above do files
   compress
   save "data/coded_cto/recruitment coded.dta", replace
