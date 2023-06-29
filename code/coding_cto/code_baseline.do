
 * This do file: appends and codes all baseline (stated preference survey) data, calculates VOT 
 * This version: 14 Feb, 2017, 22 June, 2017 (VOT), 25 October 2018 (revision)
 * Authors: Gabriel Kreindler, Ashwin MB	
 
clear all
pause off
set more off

* Path based on global set in C:\ado\profile.do
	local path "$congestion_pricing_path"
	local pathcode "$cp_code_path"
	cd "`path'"

* Load recruitment data
	use "data/coded_cto/recruitment pooled.dta", clear
	keep uidp_original uidps_original /* phone_no add_phone__1 name */
	rename uidps_original uidps
	// tostring add_phone__1 name, replace
	// tostring phone_no, replace
	// tostring name, replace
	tostring uidp_original, replace
	tostring uidps, replace
	tempfile recruit
	save 	`recruit'

	
*** Refresh dta files:
	* Note: scripts not included because the replication package only contains the output
	// do "`pathcode'clean_cto/import_baseline_v1.do"
	// do "`pathcode'clean_cto/import_baseline_v2.do"
	// do "`pathcode'clean_cto/import_baseline_v3.do"
	// do "`pathcode'clean_cto/import_baseline_v4.do"


**********************************************
*** Combine different versions of raw data ***
**********************************************
	// use "data/raw_cto/recruit/Baseline v1.dta", clear
	// gen form=1

	use "data/raw_cto/baseline/Baseline v2 noPII.dta", clear
	gen form=2

	// fixes
	tostring dest_flex, gen(temp)
	drop dest_flex
	rename temp dest_flex

	destring *_dep_total_*, replace

	append using "data/raw_cto/baseline/Baseline v3 noPII.dta", 
	recode form (.=3)

	append using "data/raw_cto/baseline/Baseline v4 noPII.dta", 
	recode form (.=4)

	rename uidp uidp_office
	replace uidp_office = upper(uidp_office)

	// list if key == "uuid:11e917d3-afc7-4e13-9b90-911e28642dbc"

* sort 
	sort ru_section_order

* clean - rename long variables
	rename am_dt_intention_noresponse_chec am_dt_noresp_check
	rename pm_dt_intention_noresponse_chec pm_dt_noresp_check
	rename am2_dt_intention_noresponse_chec am2_dt_noresp_check

* clean - delete useless variables
	assert short_uidp == short_uidp_confirm
	drop short_uidp_confirm
	drop am_dt_cheaper_k am_dt_more_exp_k pm_dt_cheaper_k pm_dt_more_exp_k
	drop *_dt_intention_direction *_dt_intention_direction_k
	drop *dep_time_0_mm_string *dep_time_0_hhmm *dep_time_1_mm_string *dep_time_1_hhmm *dep_time_mean_tot_mm *dep_time_mean_mm_string *dep_time_mean_hhmm
	drop *_dep_delta_*_frac /**_dep_total_**/
	drop am_dt_toll_10 am_dt_toll_20 pm_dt_toll_10 pm_dt_toll_20

	foreach va of varlist *confirm_large*{
		assert inlist(`va',.,1)
		drop `va'
	}

* destring numeric variables	
	destring ru*, replace
	destring cell_* *_delta_t *vot_base_extra *vot_slower *dt_base_toll *dt_p10 *dt_p20 , replace
	destring *_dep_time_*_mm *_dep_time_*_hh, replace
	destring *mean_dur, replace
	destring *vot_dur *vot_dur_slower, replace

* move all am2_ variables back to am_
	gen section_order_AM = ru_section_order < 0.5
	foreach va of varlist *am2*{
		local `va'no2 = subinstr("`va'", "am2_", "am_", .)
		replace ``va'no2' = `va' if section_order_AM == 0 & partial_survey == 2
		drop `va'
	}
	
* order preloaded variables at the correct spot
	order am_delta_t, before(am_mean_dur)
	order pm_delta_t, before(pm_mean_dur)

	order am_vot_base_extra am_vot_slower, before(am_vot_dur)
	order pm_vot_base_extra pm_vot_slower, before(pm_vot_dur)

	order am_dt_base_toll - am_dt_p20, before(am_dt_intention)
	order pm_dt_base_toll - pm_dt_p20, before(pm_dt_intention)

	* cosmetic
	* remove ugly html labels
	foreach var of varlist _all {
		label var `var' ""
	}

*********************************************************
*** Date variable, Typo Cleaning, Merge Recruit data  ***
*********************************************************

*** basic date and time coding
	* dates
	gen startdate = dofc(starttime)
	format startdate %td
	gen startmonth=month(startdate)
	format startmonth %tm
	gen startday=day(startdate)
	gen starthour = hh(starttime)
	gen startmin  = mm(starttime)
	gen startsec  = ss(starttime)

	gen enddate = dofc(endtime)
	format enddate %td
	gen endday=day(enddate)
	gen endhour = hh(endtime)
	gen endmin  = mm(endtime)
	
	gen subdate = dofc(submissiondate)
	format subdate %td
	gen subday=day(subdate)
	gen subhour = hh(submissiondate)
	gen submin  = mm(submissiondate)

	* generate other times
	* TODO: minutes relative to starttime
	foreach va of varlist time_*{
		cap drop temp
		gen double temp = clock(`va',"YMD hms")
		drop `va'
		rename temp `va'
		format %tc `va' 
	}

	replace short_uidp = upper(short_uidp)
	rename short_uidp uidps

*** corrections
	drop if key == "uuid:8cbd8e50-107e-4255-8e94-ed47ca176ec5" // Sharath suspicious survey
    
	replace uidps = "S0188" if key == "uuid:f3bda184-d005-486b-bd7f-3fef30520ea1"
	replace uidps = "S0277" if key == "uuid:21f67884-76a2-4fc2-aa92-bdb209643bba"
	replace uidps = "S0277" if key == "uuid:6e558cbc-0759-4ae4-8547-c1bf12a4f9e0"
	replace uidps = "S0434" if key == "uuid:da3ed8e3-c69e-4514-ae21-49b806c2af4a"
	replace uidps = "S0287" if key == "uuid:925437c6-33af-462f-9f5b-fbe7901b4c6e"
	replace uidps = "S2416" if key == "uuid:9ac8ed1f-19ea-4a18-b5f1-68b0b9fc31a3" // original seems correct but changed in recruitment pooled
	replace uidps = "S0670" if key == "uuid:698e69d2-95be-4613-95d9-440bd9406901"
	replace uidps = "S0265" if key == "uuid:d8fdff51-22e3-4ace-87ca-b4b79631b309"
	replace uidps = "S4099" if key == "uuid:71ce07da-329f-424a-b3f4-97699bc9b192"
	
	replace uidps = "S6052" 			if key == "uuid:15a81375-b04e-4382-8f51-764aa48000b2" // name and phone match this uidps
	replace uidp_office = "L0408105809" if key == "uuid:15a81375-b04e-4382-8f51-764aa48000b2" // name and phone match this uidps

	replace uidp_office = "L0424170750_63AE" if key == "uuid:53dfc734-ec41-4b74-8dce-9d35ccf8267f"
	replace uidp_office = "L0424173738_33F9" if key == "uuid:a8acc459-d0d5-4b3a-bd1d-ef30e06246aa"
	replace uidp_office = "L0425090154_54AB" if key == "uuid:54899885-254a-47bb-ba7a-11f70eee2005"
	replace uidp_office = "L0425095940_E105" if key == "uuid:3ba48d96-d9e5-48e0-98e8-b5ca210131ab"
	replace uidp_office = "L0425102633_A999" if key == "uuid:210a40a0-b95c-4a6c-acbe-343cca7c908f"
	replace uidp_office = "L0425104320_33B2" if key == "uuid:2065358a-7aba-413e-b631-331f38d1f761"
	replace uidp_office = "L0425150207_BE89" if key == "uuid:948327c5-8a09-4c4b-b0cc-aef79595cd77"
	replace uidp_office = "L0303105651" if key == "uuid:05f6d775-853c-47a6-ba38-631338c2068a"
	replace uidp_office = "L0308151018" if key == "uuid:87e42ffd-0789-4fd1-b251-72c494ef20ca"
	replace uidp_office = "L0413095839" if key == "uuid:9e19a4bb-a996-4ca1-b174-6a892e4dc3a5"
	replace uidp_office = "L0419084801" if key == "uuid:42527980-02b5-45e7-84a7-98dc314eac7d"
	replace uidp_office = "L0422183529_EB4D" if key == "uuid:d36fe78a-ec66-4b38-b2a7-1f869c036e91"
	replace uidp_office = "L0424100846_AFE3" if key == "uuid:821ce59a-4adf-4c2b-9615-b333cd070674"

	replace uidp_office = "L0318131309" if key == "uuid:fe3a2c6d-59af-41b5-a39e-37aa31921f4e"
	replace uidp_office = "L0413113427" if key == "uuid:f2121f63-1c68-4041-b2d3-f484977855d4"
	replace uidp_office = "L0411193019" if key == "uuid:281e5030-a4a7-429b-877b-6d85b9834e41"
	replace uidp_office = "L0425130906_3300" if key == "uuid:1dfabe73-7a26-4836-955f-a5b0736ee959"
	replace uidp_office = "L0420140533" if key == "uuid:46815a48-8373-48ce-9c1f-53d3d1226911"
	replace uidp_office = "L0419174134" if key == "uuid:94529d2d-e305-4af2-8c22-9cd46f261c67"

    replace uidp_office = "L0516184700_C5E7" if key == "uuid:4727c1a8-b875-4920-8852-62a3c8855870"
    replace uidp_office = "L0513185027_AB24" if key == "uuid:f0c52aff-a4a1-4fc9-91b8-4fd40eae5a92"
    replace uidp_office = "L0504193107_4BC7" if key == "uuid:d235c173-520c-440d-ad28-47c65c9b3fa8"
    replace uidp_office = "L0504135613_536B" if key == "uuid:76378aca-0c8b-45b3-82a4-2a2a87d1971e"
    replace uidp_office = "L0427155549_AEA0" if key == "uuid:c91cfe44-dbc8-4fc5-8754-423a03da7e64"
	replace uidp_office = "L0425194311_634C" if key == "uuid:d40a828d-804a-499e-833e-c929ec02cb64"
    replace uidp_office = "L0411103832" if key == "uuid:d7c49d2a-1902-450e-aa97-eeeba5966df9"
	replace uidp_office = "L0505180850_130F" if key == "uuid:50167f92-f393-4fb7-9003-81f4cd0c97b8"
	replace uidp_office = "L0525104002_B610" if key == "uuid:42364e2c-71b4-4a68-89f4-88c2f0efa4fb"
	replace uidp_office = "L0530165130_DC0A" if key == "uuid:96fd91fa-7890-4481-a315-88f195d70683"
	replace uidp_office = "L0511153309_39CA" if key == "uuid:e8493cb0-9845-46b5-9c56-ae644e5dab6c"
	replace uidp_office = "L0515122455_9F89" if key == "uuid:a2e3a21f-c733-4fbf-b02a-e911715d359f"
	replace uidp_office = "L0410174301" if key == "uuid:54a7e00d-6dd3-4770-8e9e-d77c74c575b3"
	replace uidp_office = "L0413112809" if key == "uuid:9feb73f8-3089-4487-a586-ed48e045fc7c"
    replace uidp_office = "L0603154038_4C0B" if key == "uuid:8c709c9b-cc7c-4719-804a-d7b0422e04be"
	replace uidp_office = "L0422151139_3F77" if key == "uuid:5abe6ce6-9d89-4e30-b676-fea816779a26"
    replace uidp_office = "L0421201319" if key == "uuid:6ba76e30-bd62-4fb0-8d67-0c54a9f6598d"
    replace uidp_office = "L0527184949_2818" if key == "uuid:186f9617-f700-47c1-bd59-6c60ac9d7872"
    replace uidp_office = "L0425195539_D949" if key == "uuid:5b6275cc-a818-4263-a285-d0e4c3798ca8"
    replace uidp_office = "L0428180732_8E18" if key == "uuid:adab708a-a6a5-49f1-b862-fc95bc108b9b"
    replace uidp_office = "L0422121301_436E" if key == "uuid:2d85924f-5699-4391-8859-7e4f5ba09524"

    * assert errors on 6/10/2019. Names and phone numbers look ok
    replace uidp_office = "L0410113030" if key =="uuid:728d5b78-523b-42e8-bd58-28eeb74179d3"
    replace uidp_office = "L0410113028" if key =="uuid:9029ae66-9111-4268-b4ea-500895e75b7f"
	replace uidp_office = "L0410113031" if key == "uuid:e0d16c8c-3e3d-43f8-bbc1-fd3c98206204"
	replace uidp_office = "L0410113042" if key == "uuid:618a30e5-5233-4816-a701-99c996193629"
	

*** clean uidps and merge with recruitment data
	merge m:1 uidps using `recruit'
	assert _m!=1
	// cap assert _m!=1
	// if _rc!=0{
	// 	br uidps name_office phone_no_office key if _m==1
	// 	STOPHERE
	// }
	keep if _m==3
	drop _m

	// assert (phone_no == phone_no_office | add_phone__1 == phone_no_office) & (uidp_office == uidp_original | form <=3)
	// if _rc!=0{
	// 	br uidps uidp_office uidp_original name_office name phone_no_office phone_no add_phone__1 key ///
	// 			if !((phone_no == phone_no_office | add_phone__1 == phone_no_office) & (uidp_office == uidp_original | form <=3))
	// 	STOPHERE
	// }

*** Combine AM and PM if they are in different surveys
	gen has_am = am_dep_time_0 !=.
	gen has_pm = pm_dep_time_0 !=.

	tab has_am partial_survey
	tab has_pm partial_survey, m

	gen 	partial_actual = 0
	replace partial_actual = 1 if has_am==1 & has_pm==0
	replace partial_actual = 2 if has_am==0 & has_pm==1
	replace partial_actual = 3 if has_am==1 & has_pm==1

	assert inlist(partial_actual,1,2,3) 

*** SAVE 
	save "data/coded_cto/baseline_pooled.dta", replace

