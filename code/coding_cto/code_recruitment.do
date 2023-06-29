 
 * This do file: Cleans, appends and codes all (Survey CTO) recruitment data
 * This version: 6 Aug, 2018 (6 Feb, 2017)
 * Authors: Anupriya Khemka, Gabriel Kreindler

clear all
pause off
set more off
// cap log close

* Path based on global set in C:\ado\profile.do
	local path "$congestion_pricing_path"
	cd "`path'"

* Corrections to the deviceid_dcs
	import excel using "data/raw_cto/corrections_recruit_v1.xlsx", clear firstrow
	
	keep if forreal==1
	save "data/temp/corrections_recruit_v1", replace
	
*** Refresh dta files:
	* Note: script not included because the replication package only contains the output
	// do "code/clean_cto/import_recruit_v1.do"
	// do "code/clean_cto/import_recruit_v2.do"
	// do "code/clean_cto/import_recruit_v3.do"
	// do "code/clean_cto/import_recruit_v4.do"
	// do "code/clean_cto/import_recruit_v5.do"
	// do "code/clean_cto/import_recruit_v6.do"

*** check the crosswalk exist - otherwise, create (strings)
	capture confirm file "data/raw_cto/crosswalk_uidp_suidp.dta"
  	if _rc==0 {
  		display "The crosswalk file exists."
  	}
  	else{
		clear all
		gen key = ""
		gen uidp = ""
		gen uidps = ""
		save "data/raw_cto/crosswalk_uidp_suidp.dta", replace
		display "The crosswalk file does not exist - created and saved."
	}


**********************************************
*** Combine different versions of raw data ***
**********************************************
	use "data/raw_cto/recruit/Recruitment v1 noPII.dta", clear
	gen form=1
	
	append using "data/raw_cto/recruit/Recruitment v2 noPII.dta"
	replace form=2 if form==.
	// order phone_carrier_0 phone_prepost_0, after(phone_no)

	append using "data/raw_cto/recruit/Recruitment v3 noPII.dta"
	replace form=3 if form==.
	// order phone_carrier_0 phone_prepost_0, after(phone_no)

	label drop surveyor_code phone_carrier has_smart_brand
	append using "data/raw_cto/recruit/Recruitment v4 noPII.dta"
	replace form=4 if form==.

	append using "data/raw_cto/recruit/Recruitment v5 noPII.dta"
	replace form=5 if form==.

	append using "data/raw_cto/recruit/Recruitment v6 noPII.dta"
	replace form=6 if form==.

* cosmetic
	* remove html labels
	foreach var of varlist _all {
		label var `var' ""
	}

*** Drop test and training surveys
	drop if key == "uuid:29e71409-797e-415f-9260-16dd9c7acd11"
	drop if key == "uuid:6cfa15eb-c02f-4daa-a71a-48f1c6ff8a9b"
	drop if key == "uuid:1ec03b88-a70f-45d9-a66b-5bd55e13aded"
	drop if key == "uuid:60c1ee01-f654-4e3f-b31e-0bc51f5b3b21"
	drop if key == "uuid:000823e7-f220-4431-b08c-f1258df494db"

*** Drop duplicate surveys (some due to manual restore on tablets)
	drop if key == "uuid:be9e7665-ba55-43fa-973b-ba89fb8ead3f"
	drop if key == "uuid:db1da38a-8665-4586-bb21-106f3e8b89a5" 
	drop if key == "uuid:34153670-1368-40c4-942d-a1f192a251da"
	drop if key == "uuid:dac90fff-b062-4a1b-87b6-f98a1463b0cb"
	drop if key == "uuid:efb57ff9-d326-4ed6-bfb4-6a548aeb37d1"
	drop if key == "uuid:e6801e81-0327-49ef-a2d6-705cf0a7bb0b"
	drop if key == "uuid:b663efc8-ab87-457b-8651-a7c143ed7175"
	drop if key == "uuid:5d8359ad-5d33-4af0-88b3-e964a8b7e0f9"
	drop if key == "uuid:63f35263-4aa1-4e98-8a4a-5b360ee03eb4"
	drop if key == "uuid:e19975b4-fffe-4624-a6da-ca2f472fc493"
	drop if key == "uuid:b415c0b5-160a-457f-92f4-cd7d8ef496e2"
	drop if key == "uuid:a730bd2c-4b44-4016-a216-32c2c2191c54"
	drop if key == "uuid:215a756c-ad3e-40e4-a38e-1f4ac7ee0a90"
	drop if key == "uuid:917e1711-d4c9-4412-99b1-ef8baf4a0b8c"
	drop if key == "uuid:c1033f37-3bbf-4ecf-aa50-b3215b91f9c9"
	drop if key == "uuid:1427a413-6aa3-4b51-ad59-10ace072a5f2"
	drop if key == "uuid:0977196c-2dd6-48b9-8dd5-f28b83ecc90b"
	drop if key == "uuid:ff0eb94c-a87e-48e3-9d93-7b9bf5be6d05"
	drop if key == "uuid:64110826-714d-4dd0-962b-53508c556a35"
	drop if key == "uuid:4b946394-da21-42b5-810a-48e6cb574f77"
	drop if key == "uuid:4c0bf201-476b-409b-a4f0-e333e801a6a1"
	drop if key == "uuid:a8d695f5-8dca-42d8-850f-feb5d794f5f4"

*** Pure "call later" from the office should not be logged as suerveys:
	drop if key == "uuid:4ded902a-e495-4ff5-bde4-050ef7a40612"
	drop if key == "uuid:ccf0e432-1d90-4fc1-8aa7-a28a6e495851"
	drop if key == "uuid:94dfd94e-8da5-4622-ac5e-438828c14060"
	drop if key == "uuid:1093d782-b7a9-4649-93d7-cad7a83ede2b"

	* this survey was superseeded by a second "from office" (3rd March,2017):
	drop if key == "uuid:f72a025d-2958-4bed-8b14-035bf9a658be"

	* this survey was superseeded by a second "from office":
	drop if key == "uuid:96aec724-f4d5-4a50-876f-fa55c012032f"

	* this survey was done even though the original recruitment survey was complete
	drop if key == "uuid:ebe5eef0-99bd-4f44-96fb-4adb53526a87"

	* this survey was done even though the original recruitment survey was complete
	drop if key == "uuid:c4de6417-6c92-48d0-97ce-298442581103" 

	* these surveys were done twice for same person (interviewed twice)
	drop if key == "uuid:79ebe5e3-fe37-48cc-bd4d-e2519ed64c05"
	drop if key == "uuid:cfb43014-65dc-4582-99fd-265e911e5260"


*** Drop surveys automatically and erronously submited twice
	* ru1 is a random number with very high precision. If two surveys have same ru1 and starttime they are identical surveys (checked manually)
	duplicates tag starttime ru1 /* name phone_no */, gen(dups)
	sort starttime ru1 /* name phone_no */ key, stable
	by starttime ru1 /* name phone_no */ (key): gen to_drop = dups * (_n==1)
	count if to_drop == 1
	assert r(N) == 4
	drop if to_drop == 1
	isid ru1

*** check for legitimate duplicates in starttime
	* some surveys start at the same second. Add one second to make starttime a unique ID
	* isid starttime fails because after April 24th the unique id also includes last 4 digits in key
	replace starttime = starttime + 1000 if key == "uuid:4c4a1985-4efa-4c10-a0fa-caa0c3da7210"
	replace starttime = starttime + 1000 if key == "uuid:237554f6-9bd2-4091-b655-21843089c52c"
	replace starttime = starttime + 1000 if key == "uuid:b9a893b5-a8ed-4143-8c75-4f2f4ac41330"
	replace starttime = starttime + 1000 if key == "uuid:4cc4d277-59af-4c76-92e8-241eef034e3e"
	replace starttime = starttime + 1000 if key == "uuid:80109e0c-5908-421b-bd4f-8f27c9f8b9a6"
	replace starttime = starttime + 1000 if key == "uuid:9738c599-e095-4c7b-a06e-84997212dc3d"
	replace starttime = starttime + 1000 if key == "uuid:a4db21dd-fa9e-47a6-8048-45d2b6bd7ea7"
	replace starttime = starttime + 1000 if key == "uuid:933eace8-e250-479c-9d57-e25f00ca1966"
	replace starttime = starttime + 1000 if key == "uuid:dc12afc9-fc28-4cef-9af4-ae6579950460"
	replace starttime = starttime + 1000 if key == "uuid:9a966e19-d16d-45da-a10b-bdff10f466d9"
	replace starttime = starttime + 1000 if key == "uuid:cdaf9cc7-fec4-41c3-845f-c5162f7d2413"
	replace starttime = starttime + 1000 if key == "uuid:7c6590b0-854a-4e6c-b596-c668b6549bea"
	replace starttime = starttime + 1000 if key == "uuid:80a1f937-440e-4d8c-9512-126d84ebffa0"
	replace starttime = starttime + 1000 if key == "uuid:d8a47912-86c0-43cd-9b55-cc11569541b1"
	replace starttime = starttime + 1000 if key == "uuid:93e1df0d-f351-4200-9a4c-59f2a96f57e9"
	replace starttime = starttime + 1000 if key == "uuid:d0f0c869-0bb8-4ae0-9e13-12a722f9a75b"
	replace starttime = starttime + 1000 if key == "uuid:6858a05b-3c8c-4aa3-b577-6a5ed639e9de"
	replace starttime = starttime + 1000 if key == "uuid:3b51ee4c-89bd-4d51-a021-da738d138da9"
	replace starttime = starttime + 1000 if key == "uuid:7b684480-b5cc-4e84-bac6-9f8449600940"
	replace starttime = starttime + 1000 if key == "uuid:ef9b9200-4c8a-486b-a738-16c92d302310"
 	replace starttime = starttime + 1000 if key == "uuid:3192449d-b79d-4bb4-9bd1-04a648fe38a0"
 	replace starttime = starttime + 1000 if key == "uuid:b4ee00c2-0307-4007-9144-caa539d96a4c"
 	replace starttime = starttime + 1000 if key == "uuid:cbed9a7c-22ac-402e-a69d-15b2741f9136"
 	replace starttime = starttime + 1000 if key == "uuid:d0992ce4-037c-41e5-9977-651d8696f24b"
	replace starttime = starttime + 1000 if key == "uuid:eb2b968c-4e59-437f-be51-101871477f78"
	replace starttime = starttime + 1000 if key == "uuid:f39c80a5-06e6-4ad3-85a2-f20d5cb6462e"
	replace starttime = starttime + 1000 if key == "uuid:bdc7efc3-65e7-42e8-a7ac-0711165a70e4"
	replace starttime = starttime + 1000 if key == "uuid:4afb1f4c-70ed-46f5-93fa-60d1e3da12a5"
	replace starttime = starttime + 1000 if key == "uuid:a5246861-c0cb-40ef-b1ad-d50bac992852"


*** basic coding corrections
	replace gender = 1 if key == "uuid:63bfb427-27a3-4762-a197-4d7ec5243f36" 

*** generate other times
	foreach va of varlist time_*{
		cap drop temp
		gen double temp = clock(`va',"YMD hms")
		drop `va'
		rename temp `va'
		format %tc `va' 
	}

*** Fix wrong tablet time for three entries
	gen n=_n
	local i = 1
	foreach key in "uuid:919d3ffb-6775-448b-83ac-2f70d7a3a5e8" "uuid:9148937c-bee8-4b63-bc05-2f9306a06ac7" "uuid:e081ac14-dcd5-406a-9175-288e288cc82b" ///
				   "uuid:a9baca19-da62-4929-baab-d5f273506eb4" "uuid:0eee627b-dc98-4c30-9f8b-6070e5262d6e" "uuid:45c3500b-e557-43b1-81be-3a7ba09b0a52" "uuid:8ee19366-cb33-4fb6-b007-0d06cc79756c" ///
				"uuid:79f46298-567c-495c-aaa0-e53a403456ff" "uuid:d54e25c7-3b68-4c7b-9e19-63e38699fc7d" "uuid:badfcbed-8261-4080-9e45-b10c6d62acf1" "uuid:0cc74d70-8844-43c5-b4a9-d60ead93caaf" ///
				"uuid:a441f7fe-a7c2-4c25-819d-fcd8e4dbb891" "uuid:0e2c6a5f-3967-41ba-b831-fb30a45e8923" "uuid:b2cc5745-ddb5-4156-847e-2cb25f63fec5" "uuid:b54e8813-7bcd-41b3-ba7c-8599e24ef7f0" ///
				"uuid:7c9b5eb5-149a-45bd-9d94-659a692e6e98" "uuid:2b467188-df11-4228-9a18-a9008ef732ab" "uuid:80695634-38d3-4bec-ac7a-59ad586fe34a" "uuid:1698c6cc-9cb0-4db4-9354-07c42352d3af" ///
				"uuid:a5707ef7-2177-457a-9009-9cb5ff2b93bb" "uuid:cc879e97-94f5-4373-9d84-09c17f9ea1a3" "uuid:f1a0641f-b8f3-4bf9-ac2c-190392f9fc7f" "uuid:0613f470-3090-431e-abf9-8a12d3351584" ///
				"uuid:26f66464-be29-43e6-a799-6ffff8ff6c8c" "uuid:9415cc13-869d-4ff2-8851-54031f65c891" "uuid:29aee240-1569-4905-b04e-e43530de3b68" "uuid:776ad6e4-86fc-419f-864c-a7f88f1be50b" ///
				"uuid:7f49c657-13fa-403b-93b5-f2df309db39c" "uuid:f194ef6f-7846-4736-b08f-4f3ff193c021" "uuid:e6b68646-cde0-48b7-b730-0238e0d8b63d" "uuid:f0e68e23-1fe9-4bc6-adfc-69492e318554" ///
				"uuid:0173c1b9-a168-4502-af97-cf4b0ce06588" "uuid:972ad5d3-3027-48c6-a469-182b1b4f596a" "uuid:0a2114db-3e3b-4c73-a237-a6b30e8f9d34" "uuid:df9a9826-e584-48bf-904c-fd0d8b7861d2" ///
				"uuid:c7ff9762-c6be-4e8e-bcb8-d8ca4a3aa897" "uuid:d9aff1f1-5e69-4f0b-8ebb-811c1d4ed409" "uuid:fcdd1648-085f-4310-8203-46c8e4cf8850" "uuid:8f2a4442-2ee9-48d9-895e-13a7bceb6016" {
		sum n if key == "`key'"
		local index = r(mean)
		local starttime = starttime[`index']
		local time_diff = Cmdyhms(4,10,2017,11,30,`i') - `starttime'
		foreach va of varlist time_elig time_basic_info time_app_consent time_app_1 time_app_2 time_app_3 time_consent_signature time_occ time_km time_main time_obs_info time_end starttime endtime{
			replace `va' = `va' + `time_diff' if key == "`key'"
		}
		local i = `i' + 1
	}
	drop n

***********************************
*** Basic date and time coding ****
***********************************

	gen startdate = dofc(starttime)
	format startdate %td
	gen startmonth=month(startdate)
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


***********************************
*** UNIQUE ID AND SHORT VERSION ***
***********************************
*** generate unique id UIDP

	gen key_end = upper(substr(key,-4,4))
	isid ru1
	isid starttime key_end
	sort starttime key_end
	gen uidp = "L" + string(startmonth,"%02.0f") + string(startday,"%02.0f") +  ///
					  string(starthour,"%02.0f") + string(startmin,"%02.0f") + string(startsec,"%02.0f")

	gen uidp_simple = uidp

	* make uidp more complex starting on April 24th
	replace uidp = uidp + "_" + key_end if startdate >= mdy(4,22,2017)
	isid uidp

*** merge existing correspondence (crosswalk) between uidp and uidps (short uidp)
	
	* Existing correspondence between UIDP and short uidp (UIDPS)
	merge 1:1 key using "data/raw_cto/crosswalk_uidp_suidp.dta"
	assert _m!=2
	assert uidps!="" if _m==3
	assert uidps=="" if _m!=3
	drop _m

	*** all short uidps exist
	assert uidps != ""
	isid uidps


*** short uidp inputed manually for OFFICE shifts
	replace short_uidp = (length(short_uidp) < 5) * "S" + (4 - length(short_uidp)) * "0" + short_uidp if survey_pump_type == 50
	replace short_uidp = upper(short_uidp)

	* fix office call after office call
	* point it to *original* field survey
	replace short_uidp = "S0112" if uidp == "L0211104552"

*** Another round of dropping surveys:
**** DROP erronous surveys
	drop if key == "uuid:0a3b157d-2429-4d5d-818b-f471791345e7"
	drop if key == "uuid:f4e52d4e-edc4-4c34-af27-9d82b68b3321"
	
	* Suspicious survey from (two consecutive surveys with same name)
	drop if key == "uuid:80ca836f-32d8-4719-8a9e-5de6a5d98c8f" 

	* Possibly fake surveys
	drop if inlist(uidp, "L0410164011", "L0418155849", "L0419165605")

	* Same device ID 
	drop if inlist(uidp, "L0418121445", "L0419095142")

	* Surveyed twice (L0308090144, L0311111158), drop the latter (L0311111158)
	drop if key == "uuid:da44bef9-2b11-449e-b4f6-4bf44286d64c"

	* Other
	drop if uidp == "L0224194912"

	* Drop training day April 5th
	drop if startmonth == 4 & startday == 5

*************************
*** CLEAN APP STATUS ****
*************************

* No app was installed for these surveys (nonsensical deviceid)
	rename deviceid tabletid
	rename deviceid_dcs deviceid
	replace deviceid = upper(deviceid)

	replace app_consent 		= 0 	if key == "uuid:6aac91fa-1529-4470-a6ba-d537d4e5c655"
	replace app_install_success = 0 	if key == "uuid:6aac91fa-1529-4470-a6ba-d537d4e5c655"
	replace deviceid 	 		= "" 	if key == "uuid:6aac91fa-1529-4470-a6ba-d537d4e5c655"

	replace app_consent 		= 0 	if uidp == "L0221183146"
	replace app_install_success = 0 	if uidp == "L0221183146"
	replace deviceid 	 		= "" 	if uidp == "L0221183146"

	replace app_consent 		= 0 	if uidp == "L0310163722"
	replace app_install_success = 0 	if uidp == "L0310163722"
	replace deviceid 	 		= "" 	if uidp == "L0310163722"

	replace app_consent 		= 0 	if uidp == "L0412202048"
	replace app_install_success = 0 	if uidp == "L0412202048"
	replace deviceid 	 		= "" 	if uidp == "L0412202048"


********************************************************
*** Merge app device ID from separate (phone) survey ***
********************************************************

	merge 1:1 uidps using "data/coded_cto/devid_short.dta", keepusing()
	drop if devid_uidp == "L0223182039" // dropped above
	assert _m!=2
	gen updated_devid = _m==3

	assert devid_uidp == uidp if _m==3
	drop devid_flyer_code devid_uidp

	cap drop devid_name

	// devid_starttime devid_submissiondate devid_surveyor_code

	foreach va of varlist devid_*{
		local root_va = substr("`va'",7,30)
		gen `root_va'_devid_survey = `va'
	}

	drop devid_starttime devid_submissiondate devid_surveyor_code

	foreach va of varlist devid_*{
		local root_va = substr("`va'",7,30)
		replace `root_va' = `va' if _m==3
		drop `va'
	}

	* change to agreed for app
	replace app_consent = 1 if _m==3
	replace app_install_success = 1 if _m==3

	* keep track of devids from later survey
	gen devid_later_survey = _m==3
	drop _m

*** Device id corrections (sometimes surveyors entered the device ID incorrectly)
	merge 1:1 uidp using "data/temp/corrections_recruit_v1.dta", keepusing(deviceid*)
	assert  _m != 2
	assert deviceid == deviceid_wrong if _m == 3

	tab surveyor_code if _m==3

	* replace with correct deviceid version
	replace deviceid = deviceid_correct if _m==3
	drop _m
	drop deviceid_wrong deviceid_correct


****************************
*** Cleaning and Coding ****
****************************
	* destring
	destring elig_calc inel_later do_survey ru1 ru2 ru3 formdef_version, replace

	* recharge
	cap gen recharge_phone = phone_no
	rename phone_carrier_0 recharge_carrier
	rename phone_prepost_0 recharge_prepost

	cap replace recharge_phone = phone_no_alternate if confirm_rech_main == 0
	replace recharge_carrier = phone_carrier 	if confirm_rech_main == 0
	replace recharge_prepost = phone_prepost  	if confirm_rech_main == 0

	cap drop phone_no_alternate 
	drop confirm_rech_main


/* code STATUS reminder:
	consent_1:
		consent_1	1	Yes
		consent_1	2	Call later (NO TIME for eligibility)
		consent_1	4	Complete Refusal

	consent_2:		
		consent_2	1	Yes
		consent_2	2	Call later
		consent_2	4	I do not want to install the app.

	app_consent:		
		consent_app	1	Yes
		consent_app	0	No
		// consent_app	2	Later - no time
		// consent_app	3	Later - install failed

	app_install_success:
		consent_app	1	Yes
		consent_app	2	Later - no time
		consent_app	3	Later - install failed
*/

* fix early form version where app_install_success not there
	cap gen app_install_success = .
	replace app_install_success = 1 if formdef_version <= 1702012109 & deviceid != ""
	replace app_install_success = 2 if formdef_version <= 1702012109 & app_consent == 2
	replace app_install_success = 3 if formdef_version <= 1702012109 & app_consent == 3
	replace app_consent = 1 		if formdef_version <= 1702012109 & inlist(app_consent,2,3)

	assert inlist(app_consent,.,0,1)

*** Code detailed survey status
	gen 	status_d = 0	if consent_1 == 4
	replace status_d = 1	if consent_1 == 2
	replace status_d = 2	if consent_1 == 1 & elig_calc == 0 & inel_later == 0
	replace status_d = 3	if consent_1 == 1 & elig_calc == 0 & inel_later == 1 & call_inel_consent == 0
	replace status_d = 4	if consent_1 == 1 & elig_calc == 0 & inel_later == 1 & call_inel_consent == 1
	replace status_d = 5	if consent_1 == 1 & elig_calc == 1 & consent_2 == 4 & consent_3 == 0
	replace status_d = 6	if consent_1 == 1 & elig_calc == 1 & consent_2 == 4 & consent_3 == 1
	replace status_d = 7	if consent_1 == 1 & elig_calc == 1 & consent_2 == 2
	replace status_d = 8	if consent_1 == 1 & elig_calc == 1 & consent_2 == 1 & app_consent == 0
	replace status_d = 9	if consent_1 == 1 & elig_calc == 1 & consent_2 == 1 & inlist(app_install_success,2,3)
	replace status_d = 10	if consent_1 == 1 & elig_calc == 1 & consent_2 == 1 & app_install_success == 1

	label var status_d "Detailed Status after Recruitment Survey"

	label define status_d ///
	0 "0 Refusal" ///
	1 "1 Call Later" ///
	2 "2 Inel, Drop" ///
	3 "3 Inel, Refusal" ///
	4 "4 Inel, Call later" ///
	5 "5 Elig, No app, Refusal" ///
	6 "6 Elig, No app, Call later" ///
	7 "7 Elig, Call later" ///
	8 "8 Elig, OK, No app" ///
	9 "9 Elig, OK, App later" ///
	10 "10 Elig, OK, App OK" 

	label values status_d status_d
	tab status_d, m

	gen status_d_code = status_d

* Summary status
	gen 	status = 0 if inlist(status_d,0,5)
	replace status = 1 if inlist(status_d,2,3)
	replace status = 2 if inlist(status_d,1,4,6,7)
	replace status = 3 if inlist(status_d,8,9)
	replace status = 4 if inlist(status_d,10)

	label var status "Summary Status after Recruitment Survey"

	label define status ///
	0 "0 Refusal" ///
	1 "1 Ineligible" ///
	2 "2 Call Later" ///
	3 "3 Survey" ///
	4 "4 App"

	label values status status
	tab status, m

* eligible / Ineligible
	gen 	eligible = elig_calc == 1
	replace eligible = . if inlist(status_d,0,1)
	label var eligible "Detailed Status after Recruitment Survey"
	tab eligible, m

* in study
	gen instudy = inlist(status_d, 1, 4, 6, 7, 8, 9, 10)
	label var instudy "In Study (excludes refusal and some ineligible)"

* recruitment survey done
	gen 	recruit_done = inlist(status_d, 4, 8, 9, 10) if instudy == 1
	label var recruit_done "Recruitment Survey done (0=pending,.=no action)"

* todo app
	gen 	app_done = status_d == 10 if instudy == 1 & status_d != 4 // 4 = inel call later - no app
	label var app_done "Basline Survey done (0=pending,.=no action)"

* language
	decode language_contact, gen(temp)
	drop language_contact
	rename temp language_contact
	replace language_contact_all = subinstr(language_contact_all,"30","Other",.)
	replace language_contact_all = subinstr(language_contact_all,"60","NotAsked",.)
	replace language_contact_all = subinstr(language_contact_all,"99","",.)
	replace language_contact_all = subinstr(language_contact_all,"3","E",.)
	replace language_contact_all = subinstr(language_contact_all,"2","H",.)
	replace language_contact_all = subinstr(language_contact_all,"1","K",.)
	replace language_contact_all = subinstr(language_contact_all,"4","Tam",.)
	replace language_contact_all = subinstr(language_contact_all,"5","Tel",.)

	gen language = language_contact + "; " + language_contact_all + " " + language_c if language_contact!="" & language_contact_all!=""
	drop language_contact* language_c


****************************************
**** Merge "Call Later"'s to Office ****
****************************************
	* At recruitment, some surveys are "call later"
	* For some of them, there is a corresponding survey from the office (survey_pump_type == 50)
	* We keep the office survey, but we want to merge some information from the original (call later) survey
	* Finally, we drop the original survey 
	* So, we always keep the later survey

	*** Pure call laters from the office should not be logged as suerveys:
	assert !(inlist(status_d, 1, 7) & survey_pump_type == 50)

	* 
	tab survey_pump_type

	* save call later details from the field (to merge with office calls)
	preserve
		* call later (except those from office)
		keep if status == 2 & survey_pump_type != 50
		keep uidp uidps /* phone_no */ this_car_or_moto this_vehicle_brand this_vehicle_model ///
					other_notes number_passengers2 gender agebracket check_flyer recharge_*
		rename uidps short_uidp
		rename uidp uidp_original

		tempfile call_later_obs_data
		save 	`call_later_obs_data'
	restore

	* the original (short) id from surveys done in the office (this is the id of the original "call later" survey) 
	preserve
		keep if survey_pump_type == 50 
		keep short_uidp
		rename short_uidp uidps

		tempfile from_office_to_drop
		save 	`from_office_to_drop'
	restore

	* update in the office survey
	merge m:1 short_uidp using `call_later_obs_data', update
	drop if _m==2
	assert _m!=5
	drop _m

	* drop original call later survey
	merge m:1 uidps using `from_office_to_drop'
	assert _m != 2
	assert status == 2 if _m==3
	drop if _m==3
	drop _m

	rename short_uidp uidps_original
	replace uidp_original  = uidp  if uidp_original  == ""
	replace uidps_original = uidps if uidps_original == ""
	gen in_office = survey_pump_type == 50 

*************************************************
**** EXTRA MINOR CODING 
*************************************************

* title
	// replace name = name_office if survey_pump_type == 50
	// replace name = gender * "Ms " + (1-gender) * "Mr " + trim(proper(name)) if name!=""

* clean leftover office variables
	// di "Phone number poorly entered in the office survey:"
	// list phone_no* if phone_no != phone_no_office & survey_pump_type == 50
	drop /* phone_no_office name_office */ short_uidp short_uidp_confirm

*************************************************
**** OUTPUT
*************************************************
* cosmetic
	* order
	order uidp* 
	cap order recharge_phone, before(recharge_carrier)
	order devid_later_survey, after(deviceid)

* final check
	isid uidps 
	isid uidp
	isid uidps_original
	isid uidp_original

	tab status

*** Save pooled data
	save "data\coded_cto\recruitment pooled.dta", replace
	export delimited using "data\coded_cto\recruitment_pooled.csv", replace

	keep if status==4
	export delimited using "data\coded_cto\recruitment_pooled_app.csv", replace	
