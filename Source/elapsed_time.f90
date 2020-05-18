	subroutine elapsed_time(values_s,values_e)

		use param
		implicit none
		integer, dimension(8) :: values_s,values_e
		integer :: days, hours, minutes, seconds
		integer, dimension(12) :: days_month = (/31, 28, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31/)

		if (mod(values_e(1),4) == 0) days_month(2) = 29
		days = values_e(3) - values_s(3)
		if (values_e(1) > values_s(1)) then
			days = days + 31
		else if (values_e(2) > values_s(2)) then
			days = days + days_month(values_s(2))
		end if
		hours = 24*days + values_e(5) - values_s(5)
		minutes = values_e(6) - values_s(6)
		if (minutes < 0) then
			hours = hours - 1
			minutes = minutes + 60
		end if
		seconds = values_e(7) - values_s(7)
		if (seconds < 0) then
			minutes = minutes  - 1
			seconds = seconds + 60
		end if

		open (unit=6,file="farprt",status="old",POSITION="APPEND")		
		write(6,'(" elapsed time:",i4," hours",i3," minutes",i3," seconds")') hours,minutes,seconds
		close(6)

	end subroutine elapsed_time
