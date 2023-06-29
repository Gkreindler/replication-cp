cap program drop weekdays_bw
program define weekdays_bw
	syntax varlist(min=2 max=2), GENerate(name)
	gen __temp_1 = `1'
	gen __temp_2 = `2'

	* define Monday date in 1's week
	* define Monday date in 2's week
	gen __start_1 = dofw(wofd(__temp_1)) // __start_1 - dow(__start_1)
	gen __start_2 = dofw(wofd(__temp_2))
	gen `generate' = (__start_2 - __start_1) * 5 / 7 + dow(__temp_2) - dow(__temp_1) + 1
	drop __*
end
