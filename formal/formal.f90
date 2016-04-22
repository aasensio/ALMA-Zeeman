module formal
use iso_c_binding, only : c_double, c_int
implicit none
	real(c_double), parameter, dimension(4,4) :: identity = reshape((/1.d0,0.d0,0.d0,0.d0,0.d0,1.d0,0.d0,0.d0,0.d0,0.d0,1.d0,0.d0,0.d0,0.d0,0.d0,1.d0/), (/4,4/))
	
contains
! ---------------------------------------------------------
! Performs the formal solution of the RT equation for a plane-parallel magnetized atmosphere
! for a given source function and opacity
! ---------------------------------------------------------
	subroutine c_delopar_source(n, height, opacity, source, boundary, stokesOut) bind(c)
	integer(c_int), intent(in) :: n
	real(c_double), intent(in), dimension(n) :: height
	real(c_double), intent(in), dimension(7,n) :: opacity
	real(c_double), intent(in), dimension(n) :: source
	real(c_double), intent(in), dimension(4) :: boundary
	real(c_double), intent(out), dimension(4) :: stokesOut
	
	real(c_double) :: formal_sol_polarized(4), Inten(4)
	integer(c_int) :: k, km, kp, k0, kf
	real(c_double) :: chim, chi0, chip, dtp, dtm, exu
	real(c_double) :: psim, psi0, psip, psim_lin, psi0_lin, dm, dp

	integer(c_int) :: i, j
	real(c_double), allocatable :: ab_matrix(:,:,:), source_vector(:,:)
	real(c_double) :: sm(4), s0(4), sp(4), mat1(4,4), mat2(4,4), mat3(4,4)

		allocate(ab_matrix(4,4,n))
		allocate(source_vector(4,n))
		
		call fill_matrix(opacity,ab_matrix)

! Transform K into K* and then into K'
		do i = 1, 4
			do j = 1, 4
				ab_matrix(i,j,:) = ab_matrix(i,j,:) / opacity(1,:)
			enddo
			ab_matrix(i,i,:) = ab_matrix(i,i,:) - 1.d0
			source_vector(i,:) = opacity(i,:) * source / opacity(1,:)			
		enddo

! Boundary condition
		k0 = 2
		kf = n
		Inten = boundary

		do k = k0, kf

! Parabolic short-characteristics
			if (k /= kf) then
				km = k - 1
				kp = k + 1
				chim = opacity(1,km)
				chi0 = opacity(1,k)
				chip = opacity(1,kp)
				sm = source_vector(:,km)
				s0 = source_vector(:,k)
				sp = source_vector(:,kp)
				dm = dabs((height(k) - height(km)))
				dp = dabs((height(kp) - height(k)))
			else
! Linear short-characteristics
				km = k - 1
				chim = opacity(1,km)
				chi0 = opacity(1,k)
				chip = 0.d0
				sm = source_vector(:,km)
				s0 = source_vector(:,k)
				sp = 0.d0
				dm = dabs((height(k) - height(km)))
				dp = 0.d0
			endif

			dtm = 0.5d0 * (chim + chi0) * dm
			dtp = 0.5d0 * (chi0 + chip) * dp

			if (dtm >= 1.d-4) then
				exu = dexp(-dtm)
			else
				exu = 1.d0 - dtm + 0.5d0 * dtm**2.d0
			endif

			call lin_sc(dtm,psim_lin,psi0_lin)
			mat1 = exu * identity - psim_lin*ab_matrix(:,:,km)
			mat2 = identity + psi0_lin * ab_matrix(:,:,k)
			call invert(mat2)

			if (k /= kf) then
				call par_sc(dtm,dtp,psim,psi0,psip)
				Inten = matmul(mat2,matmul(mat1,Inten) + psim*sm + psi0*s0 + psip*sp)
			else
				call lin_sc(dtm,psim,psi0)
				Inten = matmul(mat2,matmul(mat1,Inten) + psim*sm + psi0*s0)
			endif						

		enddo

		deallocate(ab_matrix)
		deallocate(source_vector)

		stokesOut = Inten

	end subroutine c_delopar_source
	
	
! ---------------------------------------------------------
! Performs the formal solution of the RT equation for a plane-parallel magnetized atmosphere
! for a given source function and opacity
! ---------------------------------------------------------
	subroutine c_delopar(n, height, opacity, emissivity, boundary, stokesOut) bind(c)
	integer(c_int), intent(in) :: n
	real(c_double), intent(in), dimension(n) :: height
	real(c_double), intent(in), dimension(7,n) :: opacity
	real(c_double), intent(in), dimension(4,n) :: emissivity
	real(c_double), intent(in), dimension(4) :: boundary
	real(c_double), intent(out), dimension(4) :: stokesOut
	
	real(c_double) :: formal_sol_polarized(4), Inten(4)
	integer(c_int) :: k, km, kp, k0, kf
	real(c_double) :: chim, chi0, chip, dtp, dtm, exu
	real(c_double) :: psim, psi0, psip, psim_lin, psi0_lin, dm, dp

	integer(c_int) :: i, j
	real(c_double), allocatable :: ab_matrix(:,:,:), source_vector(:,:)
	real(c_double) :: sm(4), s0(4), sp(4), mat1(4,4), mat2(4,4), mat3(4,4)

		allocate(ab_matrix(4,4,n))
		allocate(source_vector(4,n))
		
		call fill_matrix(opacity,ab_matrix)

! Transform K into K* and then into K'
		do i = 1, 4
			do j = 1, 4
				ab_matrix(i,j,:) = ab_matrix(i,j,:) / opacity(1,:)
			enddo
			ab_matrix(i,i,:) = ab_matrix(i,i,:) - 1.d0
			source_vector(i,:) = emissivity(i,:) / opacity(1,:)			
		enddo

! Boundary condition
		k0 = 2
		kf = n
		Inten = boundary

		do k = k0, kf

! Parabolic short-characteristics
			if (k /= kf) then
				km = k - 1
				kp = k + 1
				chim = opacity(1,km)
				chi0 = opacity(1,k)
				chip = opacity(1,kp)
				sm = source_vector(:,km)
				s0 = source_vector(:,k)
				sp = source_vector(:,kp)
				dm = dabs((height(k) - height(km)))
				dp = dabs((height(kp) - height(k)))
			else
! Linear short-characteristics
				km = k - 1
				chim = opacity(1,km)
				chi0 = opacity(1,k)
				chip = 0.d0
				sm = source_vector(:,km)
				s0 = source_vector(:,k)
				sp = 0.d0
				dm = dabs((height(k) - height(km)))
				dp = 0.d0
			endif

			dtm = 0.5d0 * (chim + chi0) * dm
			dtp = 0.5d0 * (chi0 + chip) * dp

			if (dtm >= 1.d-4) then
				exu = dexp(-dtm)
			else
				exu = 1.d0 - dtm + 0.5d0 * dtm**2.d0
			endif

			call lin_sc(dtm,psim_lin,psi0_lin)
			mat1 = exu * identity - psim_lin*ab_matrix(:,:,km)
			mat2 = identity + psi0_lin * ab_matrix(:,:,k)
			call invert(mat2)

			if (k /= kf) then
				call par_sc(dtm,dtp,psim,psi0,psip)
				Inten = matmul(mat2,matmul(mat1,Inten) + psim*sm + psi0*s0 + psip*sp)
			else
				call lin_sc(dtm,psim,psi0)
				Inten = matmul(mat2,matmul(mat1,Inten) + psim*sm + psi0*s0)
			endif						

		enddo

		deallocate(ab_matrix)
		deallocate(source_vector)

		stokesOut = Inten

	end subroutine c_delopar

!--------------------------------------------------------------
! Inversion of a 4x4 matrix
!--------------------------------------------------------------
	subroutine invert(a)
	real(c_double) :: a(4,4)
	real(c_double) :: b(4,4), det, maxim, fabsmax
	! First some tests of singularity
		b = dabs(a)
		maxim = maxval(b)
		fabsmax = 1.d0 / maxim
		if (maxim == 0.d0) then
			print *, 'Singularity in the inversion'
!			stop
		endif

		a = a * fabsmax

   	b(1,1) = a(2,2) * a(3,3) * a(4,4) + a(2,3) * a(3,4) * a(4,2)&
      	+ a(2,4) * a(3,2) * a(4,3) - a(2,2) * a(3,4) * a(4,3)&
	   	- a(2,3) * a(3,2) * a(4,4) - a(2,4) * a(3,3) * a(4,2)
   	b(2,1) = a(2,3) * a(3,1) * a(4,4) + a(2,4) * a(3,3) * a(4,1)&
	   	+ a(2,1) * a(3,4) * a(4,3) - a(2,3) * a(3,4) * a(4,1)&
   		- a(2,4) * a(3,1) * a(4,3) - a(2,1) * a(3,3) * a(4,4)
   	b(3,1) = a(2,4) * a(3,1) * a(4,2) + a(2,1) * a(3,2) * a(4,4)&
	   	+ a(2,2) * a(3,4) * a(4,1) - a(2,4) * a(3,2) * a(4,1)&
   		- a(2,1) * a(3,4) * a(4,2) - a(2,2) * a(3,1) * a(4,4)
   	b(4,1) = a(2,1) * a(3,3) * a(4,2) + a(2,2) * a(3,1) * a(4,3)&
	   	+ a(2,3) * a(3,2) * a(4,1) - a(2,1) * a(3,2) * a(4,3)&
   		- a(2,2) * a(3,3) * a(4,1) - a(2,3) * a(3,1) * a(4,2)
   	b(1,2) = a(3,2) * a(4,4) * a(1,3) + a(3,3) * a(4,2) * a(1,4)&
	   	+ a(3,4) * a(4,3) * a(1,2) - a(3,2) * a(4,3) * a(1,4)&
   		- a(3,3) * a(4,4) * a(1,2) - a(3,4) * a(4,2) * a(1,3)
   	b(2,2) = a(3,3) * a(4,4) * a(1,1) + a(3,4) * a(4,1) * a(1,3)&
	   	+ a(3,1) * a(4,3) * a(1,4) - a(3,3) * a(4,1) * a(1,4)&
   		- a(3,4) * a(4,3) * a(1,1) - a(3,1) * a(4,4) * a(1,3)
   	b(3,2) = a(3,4) * a(4,2) * a(1,1) + a(3,1) * a(4,4) * a(1,2)&
	   	+ a(3,2) * a(4,1) * a(1,4) - a(3,4) * a(4,1) * a(1,2)&
   		- a(3,1) * a(4,2) * a(1,4) - a(3,2) * a(4,4) * a(1,1)
   	b(4,2) = a(3,1) * a(4,2) * a(1,3) + a(3,2) * a(4,3) * a(1,1)&
	   	+ a(3,3) * a(4,1) * a(1,2) - a(3,1) * a(4,3) * a(1,2)&
   		- a(3,2) * a(4,1) * a(1,3) - a(3,3) * a(4,2) * a(1,1)
   	b(1,3) = a(4,2) * a(1,3) * a(2,4) + a(4,3) * a(1,4) * a(2,2)&
   		+ a(4,4) * a(1,2) * a(2,3) - a(4,2) * a(1,4) * a(2,3)&
	   	- a(4,3) * a(1,2) * a(2,4) - a(4,4) * a(1,3) * a(2,2)
   	b(2,3) = a(4,3) * a(1,1) * a(2,4) + a(4,4) * a(1,3) * a(2,1)&
	   	+ a(4,1) * a(1,4) * a(2,3) - a(4,3) * a(1,4) * a(2,1)&
   		- a(4,4) * a(1,1) * a(2,3) - a(4,1) * a(1,3) * a(2,4)
   	b(3,3) = a(4,4) * a(1,1) * a(2,2) + a(4,1) * a(1,2) * a(2,4)&
	   	+ a(4,2) * a(1,4) * a(2,1) - a(4,4) * a(1,2) * a(2,1)&
   		- a(4,1) * a(1,4) * a(2,2) - a(4,2) * a(1,1) * a(2,4)
   	b(4,3) = a(4,1) * a(1,3) * a(2,2) + a(4,2) * a(1,1) * a(2,3)&
	   	+ a(4,3) * a(1,2) * a(2,1) - a(4,1) * a(1,2) * a(2,3)&
   		- a(4,2) * a(1,3) * a(2,1) - a(4,3) * a(1,1) * a(2,2)
   	b(1,4) = a(1,2) * a(2,4) * a(3,3) + a(1,3) * a(2,2) * a(3,4)&
	   	+ a(1,4) * a(2,3) * a(3,2) - a(1,2) * a(2,3) * a(3,4)&
   		- a(1,3) * a(2,4) * a(3,2) - a(1,4) * a(2,2) * a(3,3)
   	b(2,4) = a(1,3) * a(2,4) * a(3,1) + a(1,4) * a(2,1) * a(3,3)&
	   	+ a(1,1) * a(2,3) * a(3,4) - a(1,3) * a(2,1) * a(3,4)&
   		- a(1,4) * a(2,3) * a(3,1) - a(1,1) * a(2,4) * a(3,3)
   	b(3,4) = a(1,4) * a(2,2) * a(3,1) + a(1,1) * a(2,4) * a(3,2)&
	   	+ a(1,2) * a(2,1) * a(3,4) - a(1,4) * a(2,1) * a(3,2)&
   		- a(1,1) * a(2,2) * a(3,4) - a(1,2) * a(2,4) * a(3,1)
   	b(4,4) = a(1,1) * a(2,2) * a(3,3) + a(1,2) * a(2,3) * a(3,1)&
	   	+ a(1,3) * a(2,1) * a(3,2) - a(1,1) * a(2,3) * a(3,2)&
   		- a(1,2) * a(2,1) * a(3,3) - a(1,3) * a(2,2) * a(3,1)

		det = a(1,1) * b(1,1) + a(1,2) * b(2,1) + a(1,3) * b(3,1) + a(1,4) * b(4,1)

		a = b * (fabsmax / det)

	end subroutine invert

!----------------------------------------------------------------
! Given the values of the Zeeman opacities, fill the absorption matrix
! for each depth point
!----------------------------------------------------------------
		subroutine fill_matrix(opacity,ab_matrix)
		real(c_double) :: opacity(:,:), ab_matrix(:,:,:)

			ab_matrix(1,1,:) = opacity(1,:)   !eta_I
			ab_matrix(2,2,:) = opacity(1,:)   !eta_I
			ab_matrix(3,3,:) = opacity(1,:)   !eta_I
			ab_matrix(4,4,:) = opacity(1,:)   !eta_I
			ab_matrix(1,2,:) = opacity(2,:)   !eta_Q
			ab_matrix(2,1,:) = opacity(2,:)   !eta_Q
			ab_matrix(1,3,:) = opacity(3,:)   !eta_U
			ab_matrix(3,1,:) = opacity(3,:)   !eta_U
			ab_matrix(1,4,:) = opacity(4,:)   !eta_V
			ab_matrix(4,1,:) = opacity(4,:)   !eta_V
			ab_matrix(2,3,:) = opacity(7,:)   !rho_V
			ab_matrix(3,2,:) = -opacity(7,:)  !-rho_V
			ab_matrix(2,4,:) = -opacity(6,:)  !-rho_U
			ab_matrix(4,2,:) = opacity(6,:)   !rho_U
			ab_matrix(3,4,:) = opacity(5,:)   !rho_Q
			ab_matrix(4,3,:) = -opacity(5,:)  !-rho_Q

		end subroutine fill_matrix

! ---------------------------------------------------------
! Solves RTE using short-characteristics. Source function is supposed to vary parabolically
! between points M, O and P and RTE is integrated analitically, resulting in
!            I(O)=I(M)*exp(-tau(MO))+psi(M)*S(M)+psi(O)*S(O)+psi(P)*S(P)
! where psi are functions of the optical distance tau(MO) and tau(OP)
! It returns the value of the three psi coefficients
! ---------------------------------------------------------
   subroutine par_sc(dtm, dtp, psim, psi0, psip)
   real(c_double) :: short_car_parab
   real(c_double), INTENT(IN) :: dtm, dtp
	real(c_double), INTENT(INOUT) :: psim, psi0, psip
   real(c_double) :: exu, u0, u1 ,u2, d2, d3, d4

         if (dtm.ge.1.d-4) then
            exu=dexp(-dtm)
            u0=1.d0-exu
            u1=dtm-1.d0+exu
            u2=dtm**2-2.d0*dtm+2.d0-2.d0*exu
         else
            d2=dtm**2
            d3=dtm**3
            d4=dtm**4
            u0=dtm-(d2/2.d0)
            u1=(d2/2.d0)-(d3/6.d0)
            u2=(d3/3.d0)-(d4/12.d0)
        endif

        if (dtm * dtp /= 0.d0 .and. dtm**2 /= 0.d0 .and. dtp**2 /= 0.d0) then
			  psim=(u2-u1*(dtp+2.d0*dtm))/(dtm*(dtm+dtp))+u0
      	  psi0=(u1*(dtm+dtp)-u2)/(dtm*dtp)
      	  psip=(u2-dtm*u1)/(dtp*(dtm+dtp))
	 	  else
		  	  psim = 0.d0
			  psi0 = 0.d0
			  psip = 0.d0
		  endif

   end subroutine par_sc

  ! ---------------------------------------------------------
! Solves RTE using short-characteristics. Source function is supposed to vary linearly
! between points M and O and RTE is integrated analitically, resulting in
!            I(O)=I(M)*exp(-tau(MO))+psi(M)*S(M)+psi(O)*S(O)
! where psi are functions of the optical distance tau(MO)
! It returns the value of the two psi coefficients
! ---------------------------------------------------------
   subroutine lin_sc(dtm, psim, psi0)
   real(c_double) :: short_car_linea
   real(c_double), INTENT(IN) :: dtm
	real(c_double), INTENT(INOUT) :: psim, psi0
   real(c_double) :: exu, u0, u1, c0, cm, d2

      if (dtm.ge.1.d-4) then
         exu=dexp(-dtm)
         u0=1.d0-exu
         u1=dtm-1.d0+exu

         c0=u1/dtm
         cm=u0-c0
      else
         d2=dtm**2.d0
         c0=(dtm/2.d0)-(d2/6.d0)
         cm=(dtm/2.d0)-(d2/3.d0)
      endif
		psi0 = c0
		psim = cm

   end subroutine lin_sc

end module formal