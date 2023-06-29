* This do file: Merge Clean Occupation Data
* This version: 21 Sept, 2017
* Authors: Keerthana Jagadeesh

clear all
pause off
set more off

* Path based on global set in C:\ado\profile.do
  local path "$congestion_pricing_path"
	cd "`path'"


//---------- Occupation data cleaning --------- \\

* Loading treatment csv (from google spreadsheet)
   * the latest file is on August 31st
   import delimited "data/raw_other/occupations_clean_treatment.csv", varname(1) clear
   rename broadcategories t_broadcategories
   rename indetailcategories t_indetailcategories
   drop correctoccupation
   tempfile treatment_data
   save `treatment_data'
	
* Loading occupation csv data (from google spreadsheet)
   * the latest file is on August 31st
   import delimited "data/raw_other/occupations_clean.csv", varname(1) clear
   rename uidp_original uidp
   rename broadcategory all_broadcategory
   rename indetailcategoryedited all_indetailcategory
   
* Merge treatment occupation data to entire recruitement occupation data
   merge m:1 uidp using `treatment_data'
   assert _m!=2
    
   ** Replacing corrected occupations for treatment respondents
   replace all_broadcategory = t_broadcategories if _m ==3 
   replace all_indetailcategory = t_indetailcategories if _m ==3 
   keep uidp all_indetailcategory all_broadcategory
   rename uidp uidp_original
   tempfile merged_occupation
   save `merged_occupation'
   
* Load recruitment pooled data and merge with the occupation data
   use "data/coded_cto/recruitment pooled.dta", clear
   merge m:1 uidp_original using `merged_occupation'
   count if _m==2
   assert r(N) == 2
   drop if _m==2
   
   assert _m!=2
   drop _merge
   rename all_indetailcategory occ_indetail
   rename all_broadcategory occ_broad
   
* Extracting one new occupation column from occ_indetail and occ_broad
   gen occ_level_1 = ""
   
   * occ_level_1 grouping eliminates all the smaller frequency occupations in occ_indetail 
   replace occ_level_1 = "Specialized Professions" if inlist(occ_indetail, "Lawyer", "Accountant", "Doctor", "Teacher/Lecturer")
   replace occ_level_1 = "Mobile professions" if inlist(occ_indetail,"Driver/Delivery", "Field worker", "Real Estate", "Sales person")
   replace occ_level_1 = "Software and IT" if inlist(occ_indetail,"Software Engineers and developers", "IT Company professional")
   replace occ_level_1 = "Other Engineers and Semi-technical" if inlist(occ_indetail,"Engineer (vague)", "Semi-Technical staff", "Specialized engineers")
   replace occ_level_1 = "Office staff (low-medium)" if inlist(occ_indetail,"Bank worker", "Cashier", "Office staff (low- medium)", "Supervisor", "Government employee")
   replace occ_level_1 = "High level staff (C-suite, managers) and business owners(both high and low)" if inlist(occ_indetail,"Business Owner", "C-level staff", "Manager", "Shop owner: medical, mobile, vegetable, flower, clothes, others")
   replace occ_level_1 = "Non-office/Manual labour jobs" if inlist(occ_indetail, "Building/Construction", "Contract workers", "Contractor", "Hotel/food Industry", "Manual/Non-Office workers")
   replace occ_level_1 = "Student" if inlist(occ_indetail, "Student")
   replace occ_level_1 = "Others and Retired" if inlist(occ_indetail, "Other", "Retired")
   
   //We decided that occ_level_2 - based on travel was not necessary so I have left it out

*** Clean, save
   encode occ_indetail, gen(occ_detail)
   encode occ_broad   , gen(occ_medium)
   encode occ_level_1 , gen(occ_coarse)

   *** save hierarchy
   preserve
      drop if occ_detail == .
      keep occ_detail occ_coarse
      duplicates drop
      isid occ_detail
      order occ_coarse
      sort occ_coarse occ_detail
      export delimited "data/coded_cto/occupations_hierarchy_1.csv", replace
   restore

   preserve
      drop if occ_detail == .
      keep occ_detail occ_medium
      duplicates drop
      order occ_medium
      sort occ_medium occ_detail
      export delimited "data/coded_cto/occupations_hierarchy_2.csv", replace
   restore

   keep uidp_original occ_detail occ_medium occ_coarse

   save "data/coded_cto/occupations_clean.dta", replace


   
