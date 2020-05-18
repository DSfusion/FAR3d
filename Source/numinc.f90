	subroutine numinc

		use param
		use cotrol
		implicit none

		integer :: i
		character(len=2), dimension(36) :: charlist=(/"0","1","2","3","4","5","6","7","8","9","a","b","c","d","e","f","g","h","i", &
																	 "j","k","l","m","n","o","p","q","r","s","t","u","v","w","x","y","z"/)

		do i=1,3
			numruno(i)=numrun(i)
		end do
		do i=1,36
			if (numrun(3) == charlist(i)) exit
		end do
		if (i < 36) then
			numrun(3)=charlist(i+1)
		else
			open (unit=6,file="farprt",status="old",POSITION="APPEND")	
			write(6,'(" error stop in numinc. numrun(3) got to ""z"".")')
			close(6)
			stop
		end if

	end subroutine numinc