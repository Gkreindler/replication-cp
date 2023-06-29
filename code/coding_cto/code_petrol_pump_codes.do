* This do file: Clean petrol pump (gas station) unique IDs
* This version: 15 August 2018
* Authors: Keerthana Jagadeesh, Gabriel Kreindler

clear all
pause off
set more off

cap ssc install  sdecode

* Path based on global set in C:\ado\profile.do
    local path "$congestion_pricing_path"
	cd "`path'"


* Load recruitment data
   use "data/coded_cto/recruitment pooled.dta", clear
 
* Creating a new var called pump_type that will be used for combining pump_type & pump_code
   gen pump_type = ""
   replace pump_type = "IOC" if survey_pump_type == 1
   replace pump_type = "BP" if survey_pump_type == 2
   replace pump_type = "HP" if survey_pump_type == 3
   replace pump_type = "SHELL" if survey_pump_type == 4
   replace pump_type = "Other/Office" if survey_pump_type == 30 | survey_pump_type == 50
 
   * decodes pump_code from numeric to string for combining pump_type & pump_code
   sdecode survey_pump_code, gen(pump_code)
   
   * Cleaning so that we generate a pump_name that is similar to how the pump_name is in the combined phase1 and phase2 file
   
   replace pump_code = "003" if survey_pump_code == 3
   replace pump_code = "004" if survey_pump_code == 4
   replace pump_code = "006" if survey_pump_code == 6
   replace pump_code = "008" if survey_pump_code == 8
   replace pump_code = "010" if survey_pump_code == 10
   replace pump_code = "013" if survey_pump_code == 13
   replace pump_code = "014" if survey_pump_code == 14
   replace pump_code = "035" if survey_pump_code == 35
   replace pump_code = "007" if survey_pump_code == 70
   replace pump_code = "086" if survey_pump_code == 86
   
   gen pump_name = pump_type + pump_code
   
   //cleaning up data entry errors by surveyors
   replace pump_name = "BP003" if pump_name == "BP0"
   replace pump_name = "BP004" if pump_name == "BP1"
   replace pump_name = "BP010" if pump_name == "BP101"
   replace pump_name = "HP104" if pump_name == "BP104"
   replace pump_name = "HP123" if pump_name == "BP123"
   replace pump_name = "BP013" if pump_name == "BP23"
   replace pump_name = "BP013" if pump_name == "BP16"
   replace pump_name = "BP013" if pump_name == "BP31"
   replace pump_name = "BP004" if pump_name == "BP44"
   replace pump_name = "BP006" if pump_name == "BP66"
   replace pump_name = "BP086" if pump_name == "BP83"
   replace pump_name = "BP003" if pump_name == "BP83"
   replace pump_name = "BP008" if pump_name == "BP88"
   replace pump_name = "IOC208" if pump_name == "BP208"
   replace pump_name = "HP104" if pump_name == "HP003"
   replace pump_name = "BP035" if pump_name == "HP035"
   replace pump_name = "HP104" if pump_name == "HP101"
   replace pump_name = "HP104" if pump_name == "HP105"
   replace pump_name = "HP115" if pump_name == "HP11"
   replace pump_name = "HP123" if pump_name == "HP12350"
   replace pump_name = "HP123" if pump_name == "HP133"
   replace pump_name = "HP115" if pump_name == "HP155"
   replace pump_name = "BP003" if pump_name == "IOC003"
   replace pump_name = "BP004" if pump_name == "IOC004"
   replace pump_name = "BP013" if pump_name == "IOC013"
   replace pump_name = "BP086" if pump_name == "IOC086"
   replace pump_name = "IOC204" if pump_name == "IOC104"
   //replace pump_name = "HP104" if pump_name == "IOC104"
   replace pump_name = "IOC185" if pump_name == "IOC105"
   replace pump_name = "IOC165" if pump_name == "IOC115"
   replace pump_name = "IOC185" if pump_name == "IOC15"
   replace pump_name = "IOC185" if pump_name == "IOC182"
   replace pump_name = "IOC185" if pump_name == "IOC189"
   replace pump_name = "IOC187" if pump_name == "BP187"
   replace pump_name = "IOC298" if pump_name == "IOC198"
   replace pump_name = "IOC223" if pump_name == "IOC2"
   replace pump_name = "IOC204" if pump_name == "IOC201"
   replace pump_name = "IOC208" if pump_name == "IOC205"
   replace pump_name = "IOC208" if pump_name == "IOC2080"
   replace pump_name = "SHELL337" if pump_name == "IOC347"
   replace pump_name = "BP006" if pump_name == "SHELL006"
   replace pump_name = "BP010" if pump_name == "SHELL010"
   replace pump_name = "BP086" if pump_name == "SHELL086"
   replace pump_name = "SHELL337" if pump_name == "SHELL334"
   replace pump_name = "IOC204" if pump_name == "HP999" //Other - looked at the submission date to determine
   replace pump_name = "IOC204" if pump_name == "HP555" //Other - looked at the submission date to determine
   replace pump_name = "IOC204" if pump_name == "BP6666" //Other - looked at the submission date to determine
   replace pump_name = "HP123" if pump_name == "HP2" //Training day petrol pumps - feb 2nd, 4th and till Feb 10th
   replace pump_name = "HP113" if pump_name == "HP1" //Training day petrol pumps - feb 2nd, 4th and till Feb 10th
   replace pump_name = "IOC205" if pump_name == "IOC1" //Training day petrol pumps - feb 2nd, 4th and till Feb 10th
   
   *surveys submitted from office or others mapped to the petrol pump
   replace pump_name = "BP003" if pump_name == "Other/Office003"
   replace pump_name = "BP004" if pump_name == "Other/Office004"
   replace pump_name = "BP086" if pump_name == "Other/Office086"
   replace pump_name = "IOC223" if pump_name == "Other/Office223"   
   
   replace pump_name = "BP070" if pump_name == "BP007"
   
   *cleaning up the strings
   replace pump_name = trim(pump_name)

   drop pump_type
   drop pump_code
  
   keep uidp_original pump_name
   
   save "data/coded_cto/pumps_clean.dta", replace 
   
   
