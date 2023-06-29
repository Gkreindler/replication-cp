cap program drop fill_by_group
program define fill_by_group
	syntax varlist, fillby(varname) [replace]
	sort `fillby'
	foreach va of varlist `varlist'{
		by `fillby': egen `va'__ = mean(`va')
		if "`replace'" != ""{
			replace `va' = `va'__
		}
		else{
			gen `va'_full = `va'__
		}
		drop `va'__
	}
end
