
 * This do file: Code survey to record device id for smartphone apps installed while on the phone (recruitment done previously)
 * This version: 6 August. 2018 (2 February, 2017)
 * Authors: Gabriel Kreindler
 
clear all
pause off
set more off

* Path based on global set in C:\ado\profile.do
	local path "$congestion_pricing_path"
	cd "`path'"
	
*** These surveys were submitted in error and should be dropped
	import delimited using "data/raw_cto/corrections_drop.csv", clear varnames(1)  stringc(1)
	tempfile key_to_drop
	save 	`key_to_drop'

*** Refresh dta files:
	* Note: script not included because the replication package only contains the output
	// do "code/clean_cto/import_devid_v1.do"

**********************************************
*** Combine different versions of raw data ***
**********************************************
	use "data/raw_cto/devid/Devid v1 noPII.dta", clear
	gen form=1

*** Drop duplicate and erroneous surveys 
	merge 1:1 key using `key_to_drop'
	assert _m!=2
	drop if _m==3
	drop _m

	drop if key == "uuid:c7df46d3-798e-40a8-a9ae-b67780ab448e"
	drop if key == "uuid:e0ac6acd-1bdd-4bc2-b976-80e3118e15a3"
	drop if key == "uuid:d279e340-eb05-4a0f-9b9a-795700599f74" 
	drop if key == "uuid:f15bfe2b-5590-4766-a9b2-532aae618a71" 
	drop if key == "uuid:728e6f98-8dc1-4d4a-ba3c-0e545408d620"
	drop if key == "uuid:49a66bf6-860c-4fdc-b019-3c6073a126a0"
	drop if key == "uuid:1ed4fd31-6038-4b61-a4d6-0f844ceb0d1a"
	drop if key == "uuid:1af2ee07-78a2-481d-9492-96d09518fd01"
	drop if key == "uuid:acfeaad2-2f24-4bdb-b008-9859ab822719"
	drop if key == "uuid:a98ff582-7e75-4ef0-9104-3e921618345a"

	* 
	gen startdate = dofc(starttime)
	format %td startdate
		
*** Fix Typos
	replace uidp = "L0426103143_A9E9" if key == "uuid:c7a06e06-94e6-4b62-a4f5-4b27f1c5d077"
	replace uidp = "L0426100005_B2C4" if key == "uuid:6207d9be-0f9c-4678-80a7-0734f4526cc0"
	replace uidp = "L0426171152_C3F1" if key == "uuid:e6a15377-9b65-490e-b316-cc8657e8f63f"
	replace uidp = "L0426104654_AC6E" if key == "uuid:f5c2ef2a-6209-41e8-bfb9-4e6b64c58e44" 
	replace short_uidp = 	 "S11606" if key == "uuid:f5c2ef2a-6209-41e8-bfb9-4e6b64c58e44" 
	replace uidp = "L0422165605_6378" if key == "uuid:fa8a536e-559b-4d1e-be1c-be00bc7584c8"
	replace uidp = "L0422184004_B71F" if key == "uuid:0ba39193-4daa-496a-bd96-685c2755ad3b"
	replace uidp = "L0422100656_08DF" if key == "uuid:cd91dabc-2a6a-450d-b2f4-e10e5236d982"

	replace deviceid_dcs = "9F6GA5JMB" if deviceid_dcs == "DASIDJASD" 
	replace deviceid_dcs = "KPBDUZPTP" if key == "uuid:21d3c2d6-6b02-46f7-bbd3-f98523d58b72"

	assert uidp == "L0303185408" if key == "uuid:2d6b9ab8-f92b-4f0a-b711-e77b8366c891"
	replace deviceid_dcs = "VGN2LMWXH" if key == "uuid:2d6b9ab8-f92b-4f0a-b711-e77b8366c891"

	assert uidp == "L0301084231" if key == "uuid:e2d34ed7-55bf-4cb0-a96f-515164be37bb"
	replace deviceid_dcs = "Z244UPJRT" if key == "uuid:e2d34ed7-55bf-4cb0-a96f-515164be37bb"

	assert uidp == "L0304173742" if key == "uuid:dbd02b01-9db4-4631-84ed-601ee6a9e4f5"
	replace deviceid_dcs = "NJFPX5KFW" if key == "uuid:dbd02b01-9db4-4631-84ed-601ee6a9e4f5"
	
	assert short_uidp == "s6505" if key == "uuid:11ba92d5-733f-43c6-9cb7-6dbab6f5583c"
	replace short_uidp = "S6377" if key == "uuid:11ba92d5-733f-43c6-9cb7-6dbab6f5583c"

	assert uidp == "L0410162049" if key == "uuid:ff3a9fe1-a190-4ced-a586-c9c49ebcb0cc" 
	replace short_uidp = "S6669" if key == "uuid:ff3a9fe1-a190-4ced-a586-c9c49ebcb0cc"

*** clean and code
	keep flyer_code short_uidp surveyor_code deviceid_dcs uidp /* name_office */ ///
		 submissiondate has_smart_brand has_3g_active wants_rech ///
		 /* phone_no_recharge */ phone_carrier phone_prepost app_install_rech_c ///
		 app_install_2_c activated_autostart app_install_3_c ///
		 key starttime

	// rename name_office name
	rename deviceid_dcs deviceid

	foreach va of varlist _all{
		rename `va' devid_`va'
	}

	* make this choice zero
	gen devid_confirm_rech_main = 0
	// rename devid_phone_no_recharge devid_phone_no_alternate

	rename devid_short_uidp uidps

	replace devid_uidp = upper(devid_uidp)
	replace uidps = upper(uidps)
	replace devid_deviceid = upper(devid_deviceid)

	* correct more erronous uidp
	assert devid_uidp == "L0209112434" if uidps == "S0402"
	replace devid_uidp = "L0209105424" if uidps == "S0402"
	
* save
	drop devid_key
	isid devid_deviceid
	isid uidps

	save "data/coded_cto/devid_short.dta", replace
