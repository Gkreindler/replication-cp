
 * This do file: Create simple crosswalk between survey id (UIDP) and app id (deviceid)
 * This version: 6 August, 2018
 * Authors: Gabriel Kreindler

clear all
set more off

*Path based on global set in C:\ado\profile.do
	local path "$congestion_pricing_path"
	cd "`path'"

*** List of devices from the baseline survey * save deviceid - uidp correspondence	
	use "data\coded_cto\recruitment pooled.dta", clear
	keep if deviceid != ""
	keep deviceid uidp_original /* name phone_no */
	isid deviceid
	rename uidp_original uidp
	export delimited using "data/coded_cto/crosswalk_uidp_deviceid.csv", replace

