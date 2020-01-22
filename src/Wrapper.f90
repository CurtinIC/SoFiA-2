subroutine boxcar_x(a,x,y,z,filter_radius)
use omp_lib
integer:: i,j,k,x,y,z,niter,filter_radius,filter_size
real*4,dimension(x,y,z) :: a
real*4,allocatable :: tmp(:,:)
real*4::inv_filter_size

filter_size=2*filter_radius+1
inv_filter_size=1/filter_size
allocate(tmp(x+filter_size,y))

!$OMP PARALLEL shared(a,z) private(i,j,k,tmp,n_j)
!$OMP DO
do k=1,z
        call compute_along_x(tmp,a(:,:,k),x,y,filter_radius)
enddo
!$OMP END DO NOWAIT
!$OMP END PARALLEL


end subroutine

subroutine boxcar_y(a,x,y,z,filter_radius)
use omp_lib
integer:: i,j,k,x,y,z,niter,filter_radius,filter_size
real*4,dimension(x,y,z) :: a
real*4,allocatable :: tmp(:,:)
real*4::inv_filter_size

filter_size=2*filter_radius+1
inv_filter_size=1/filter_size
allocate(tmp(x,y+filter_size))

!$OMP PARALLEL shared(a,z) private(i,j,k,tmp,n_j) 
!$OMP DO
do k=1,z
        call compute_along_y(tmp,a(:,:,k),x,y,filter_radius)
enddo
!$OMP END DO NOWAIT
!$OMP END PARALLEL

end subroutine

subroutine boxcar_z(a,x,y,z,filter_radius)
use omp_lib
integer:: i,j,k,x,y,z,niter,filter_radius,filter_size
real*4,dimension(x,y,z) :: a
real*4,allocatable :: tmp(:,:)

filter_size=2*filter_radius+1
allocate(tmp(x+filter_size,y+filter_size))
!Check and see if you'll need ordered
!$OMP PARALLEL shared(a,z) private(i,j,k,tmp,n_j) 
!$OMP DO
do j=y,1,-1
        call compute_along_z(tmp,a(:,j,:),x,z,filter_radius)
enddo
!$OMP END DO NOWAIT
!$OMP END PARALLEL
end subroutine


subroutine fortran_copy_data(source,dest,x,y,z)
use omp_lib
integer:: i,j,k,x,y,z
real*4,dimension(x,y,z) :: source,dest

!$OMP PARALLEL shared(source,dest,x,y,z) private(i,j,k)
!$OMP DO 
do k=1,z
	do j=1,y
		do i=1,x
			dest(i,j,k)=source(i,j,k)
		enddo
	enddo
enddo
!$OMP END DO NOWAIT
!$OMP END PARALLEL
end subroutine


subroutine fortran_sign_copy(mask,dest,x,y,z,val)
use omp_lib
integer:: i,j,k,x,y,z
real*4 :: val
integer*4,dimension(x,y,z):: mask
real*4,dimension(x,y,z) :: dest

!$OMP PARALLEL shared(mask,dest,x,y,z,val) private(i,j,k)
!$OMP DO 
do k=1,z
        do j=1,y
                do i=1,x
			if (mask(i,j,k) .ne. 0) then
				dest(i,j,k)=SIGN(val,dest(i,j,k))
			endif
			!A nice simd if mask is binary dest(i,j,k)=SIGN(val,dest(i,j,k))*mask(i,j,k)
                enddo
        enddo
enddo
!$OMP END DO NOWAIT
!$OMP END PARALLEL
end subroutine

subroutine fortran_copy_blank(source,dest,x,y,z)
use omp_lib
use ieee_arithmetic,only : IEEE_Value, IEEE_QUIET_NAN
integer:: i,j,k,x,y,z
real*4,dimension(x,y,z) :: source,dest

!$OMP PARALLEL shared(source,dest,x,y,z) private(i,j,k)
!$OMP DO 
do k=1,z
        do j=1,y
                do i=1,x
			if (isnan(source(i,j,k))) then
				dest(i,j,k)=IEEE_VALUE(dest(i,j,k),IEEE_QUIET_NAN)
			endif
                enddo
        enddo
enddo
!$OMP END DO NOWAIT
!$OMP END PARALLEL
end subroutine

subroutine fortran_mask32(source,mask,x,y,z,thresh,val)
use omp_lib
integer:: i,j,k,x,y,z
integer*4 :: val
real*4 :: thresh
integer*4,dimension(x,y,z):: mask
real*4,dimension(x,y,z) :: source

!$OMP PARALLEL shared(mask,dest,x,y,z,val) private(i,j,k)
!$OMP DO 
do k=1,z
        do j=1,y
                do i=1,x
			if(source(i,j,k)>thresh .or. source(i,j,k)<-thresh) then
                        	mask(i,j,k)=val
			endif
                enddo
        enddo
enddo
!$OMP END DO NOWAIT
!$OMP END PARALLEL
end subroutine


subroutine compute_along_x(tmp,a,x,y,radius)
integer::x,y,radius
real*4,dimension(x+2*radius+1,y)::tmp
real*4,dimension(x,y)::a
integer :: i,j,tmp_cnt
integer :: filter_size
real*4 :: inv_filter_size

filter_size=2*radius+1
inv_filter_size=1/real(filter_size)
tmp(:,:)=0
tmp(1+radius:x+radius,:)=a(:,:)
do j=1,y
        do i=1,x
                if(isnan(tmp(radius+i,j))) then
                        tmp(radius+i,j)=0
                endif
        enddo
enddo

!Boundary values across Y
do j=1,y
        a(x,j)=0
        do tmp_cnt=0,filter_size-1
                a(x,j)=a(x,j)+tmp(x+tmp_cnt,j)
        enddo
        a(x,j)=a(x,j)*inv_filter_size
enddo

!Convolution
do j=1,y
        do i=x-1,1,-1
                a(i,j)=a(i+1,j)+((tmp(i,j)-tmp(i+filter_size,j))*inv_filter_size)
        end do
end do
!a(x,:)=0
end subroutine


subroutine compute_along_y(tmp,a,x,y,radius)
integer::x,y,radius
real*4,dimension(x,y+2*radius+1)::tmp
real*4,dimension(x,y)::a
integer :: i,j,tmp_cnt
integer :: filter_size
real*4 :: inv_filter_size

filter_size=2*radius+1
inv_filter_size=1/real(filter_size)
tmp(:,:)=0
tmp(:,1+radius:y+radius)=a(:,:)
do j=1,y
        do i=1,x
                if(isnan(tmp(i,j+radius))) then
                        tmp(i,j+radius)=0
                endif
        enddo
enddo

!Boundary values across X
do i=1,x
	a(i,y)=0
        do tmp_cnt=0,filter_size-1
                a(i,y)=a(i,y)+tmp(i,y+tmp_cnt)
        enddo
	a(i,y)=a(i,y)*inv_filter_size
enddo

!Convolution
do j=y-1,1,-1
        do i=1,x
                a(i,j)=a(i,j+1)+((tmp(i,j)-tmp(i,j+filter_size))*inv_filter_size)
        end do
end do
!a(:,y)=0
end subroutine


subroutine compute_along_z(tmp,a,x,z,radius)
integer::x,z,radius
real*4,dimension(x,z+2*radius+1)::tmp
real*4,dimension(x,z),intent(out)::a
real*4::inv_filter_size
integer :: i,j,tmp_cnt
integer :: filter_size

filter_size=2*radius+1
inv_filter_size=1/real(filter_size)
tmp(:,:)=0
tmp(:,1+radius:z+radius)=a(:,:)
do j=1,z
        do i=1,x
                if(isnan(tmp(i,j+radius))) then
                        tmp(i,j+radius)=0
                endif
        enddo
enddo

!Boundary values across X
do i=1,x
        a(i,z)=0
        do tmp_cnt=0,filter_size-1
                a(i,z)=a(i,z)+tmp(i,z+tmp_cnt)
        enddo
        a(i,z)=a(i,z)*inv_filter_size
enddo

!Convolution
do j=z-1,1,-1
        do i=1,x
                a(i,j)=a(i,j+1)+((tmp(i,j)-tmp(i,j+filter_size))*inv_filter_size)
        end do
end do

!a(:,z)=0

end subroutine

