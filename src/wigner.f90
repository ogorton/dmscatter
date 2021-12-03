module wigner

    implicit none
    real(kind=8), allocatable :: threej_table(:,:,:,:,:,:)
    real(kind=8), allocatable :: sixj_table(:,:,:,:,:,:)
    integer :: tablemin2j = 0
    integer :: tablemax2j = 12
    integer :: tablemin_used
    integer :: tablemax_used

    contains

        function logfac(n)
            ! Computes log(n!)
            implicit none
            integer :: n
            real*8 :: logfac
            integer :: i
            if (n<=0) then
                logfac = 0d0
                return
            end if   
            logfac = 1d0
            do i = 1, n
                logfac = logfac * dble(i)
            end do
            logfac = log(logfac)
            return
        end function logfac

        function logdoublefac(n)

            implicit none
            integer :: n
            real(kind=8) :: logdoublefac
            integer :: i, imin

            logdoublefac = 0d0

            if (n <= 0) return

            if (mod(n,2) == 0) then
                ! n is even
                imin = 2
            else
                ! n is odd
                imin = 1
            end if

            logdoublefac = 1d0
            do i = imin, n, 2
                logdoublefac = logdoublefac * dble(i)
            end do

            logdoublefac = log(logdoublefac)

            return

        end function logdoublefac

        function triangle(two_j1, two_j2, two_j3) result(delta)

            implicit none
            integer :: two_j1, two_j2, two_j3
            integer :: c1, c2, c3, c4
            real(kind=8) :: delta

            delta = 0.0D0
            c1 = two_j1 + two_j2 - two_j3
            c2 = two_j1 - two_j2 + two_j3
            c3 =-two_j1 + two_j2 + two_j3
            c4 = two_j1 + two_j2 + two_j3

            if (c1<0) return
            if (c2<0) return
            if (c3<0) return
            if (c4+1<0) return

            if (mod(c1,2)/=0) return
            if (mod(c2,2)/=0) return
            if (mod(c3,2)/=0) return
            if (mod(c4,2)/=0) return

            delta = exp(0.5d0*(logfac(c1/2)+logfac(c2/2)&
                +logfac(c3/2)-logfac(c4/2+1)))

            return
        end function triangle

        function vector_couple(j1, m1, j2, m2, jc, mc) result(cg)

            ! Computes the clebsh-gordon vector-coupling coefficient
            !  (j1/2 m1/2 j2/2 m2/2 | j1/2 j2/2 jc/2 mc/2)
            ! using algebraic expressions in factorials from Edmonds.

            implicit none
            integer :: j1,j2,jc,m1,m2,mc
            real(kind=8) :: cg, fac_prod, fac_sum, den
            integer :: t1, t2, t3, t4
            integer :: d1, d2, d3, d4, d5, d6
            integer :: dd1, dd2
            integer :: zmin, zmax, z

            cg = 0d0
            if (m1+m2 /= mc) return
            if (j1<0) return
            if (j2<0) return
            if (jc<0) return
            if (abs(m1)>j1) return
            if (abs(m2)>j2) return
            if (abs(jc)>jc) return
            if (mod(j1+m1,2) /= 0) return
            if (mod(j2+m2,2) /= 0) return
            if (mod(jc+mc,2) /= 0) return            

            t1 = ( j1 + j2 - jc)/2
            t2 = ( j1 - j2 + jc)/2
            t3 = (-j1 + j2 + jc)/2
            t4 = ( j1 + j2 + jc)/2
            if (t1<0) return
            if (t2<0) return
            if (t3<0) return

            d1 = (j1 + m1)/2 
            d2 = (j1 - m1)/2  
            d3 = (j2 + m2)/2 
            d4 = (j2 - m2)/2 
            d5 = (jc + mc)/2 
            d6 = (jc - mc)/2 

            dd1 = (jc - j2 + m1)/2
            dd2 = (jc - j1 - m2)/2

            fac_prod = sqrt(dble(jc)+1d0) &
                * triangle(j1,j2,jc) &
                * exp(0.5d0 * (logfac(d1) + logfac(d2) &
                              +logfac(d3) + logfac(d4) &
                              +logfac(d5) + logfac(d6)))

            zmin = max(0, -dd1, -dd2)
            zmax = min(t1, d2, d3) 

            fac_sum = 0d0
            do z = zmin, zmax
                den =- (logfac(z) + logfac(t1-z) + logfac(d2-z) &
                      + logfac(d3-z) + logfac(dd1+z) + logfac(dd2+z))

                fac_sum = fac_sum + (-1) ** (z) * exp(den)
            end do

            cg = fac_prod * fac_sum

            return
        end function vector_couple

        function threej(j1, j2, j3, m1, m2, m3) result(tj)

            ! Computes the Wigner 3-J symbol with arguments
            !   j1/2 j1/2 j3/2
            !   m1/2 m2/2 m3/2
            ! using clebsh-gordon vector-coupling coefficients

            implicit none
            integer :: j1,j2,j3,m1,m2,m3
            real(kind=8) :: tj

            tj = (-1) **((j1 - j2 - m3)/2)/sqrt(dble(j3)+1d0) &
                * vector_couple(j1, m1, j2, m2, j3, -m3)

        end function threej

        subroutine threej_table_init(min2j, max2j)

            implicit none
            integer, optional :: min2j, max2j
            integer :: a, b
            integer :: i, j, k, l, m, n

            integer(kind=8) :: ti, tf, clock_rate
            real(kind=8) :: dummy_real

            call system_clock(count_rate = clock_rate)
            call system_clock(count = ti)

            print*,"Initializing three-j symbol table..."
            if (.not. present(min2j)) then
                a = tablemin2j
            else
                a = min2j
                tablemin2j = a
            end if
            if (.not. present(max2j)) then
                b = tablemax2j
            else
                b = max2j
                tablemax2j = b
            end if
            print*,"Table min. 2J:",a
            print*,"Table max. 2J:",b
            print*,'Memory required (MB):',real(sizeof(dummy_real) &
                                    * (a-b+1)**6,4)*10d0**(-6)

            if (allocated(threej_table)) deallocate(threej_table)
            allocate(threej_table(a:b,a:b,a:b,a:b,a:b,a:b))
            threej_table = 0d0

            !$OMP PARALLEL DO private(n,m,l,k,j,i) schedule(dynamic, 1)
            do n = a, b
              do m = a, b
                do l = a, b
                  do k = a, b
                    do j = a, b
                      do i = a, b
                        threej_table(i,j,k,l,m,n) = threej(i,j,k,l,m,n)
                      end do
                    end do
                  end do
                end do
              end do
            end do

            call system_clock(count = tf)
            print*,'Table has been saved to memory.'
            print*,'Seconds to initialize:',real((tf-ti))/real(clock_rate)
        end subroutine threej_table_init 

        function threej_lookup(j1,j2,j3,l1,l2,l3) result(tj)

            implicit none
            integer :: j1,j2,j3,l1,l2,l3
            real(kind=8) :: tj
            integer :: minj, maxj

            minj = min(j1,j2,j3,l1,l2,l3)
            maxj = max(j1,j2,j3,l1,l2,l3)
            tablemin_used = min(tablemin_used, minj)
            tablemax_used = max(tablemax_used, maxj)

            if (minj < tablemin2j .or. maxj > tablemax2j) then
                tj = threej(j1,j2,j3,l1,l2,l3)
            else
                tj = threej_table(j1,j2,j3,l1,l2,l3)
            end if

            return
        end function threej_lookup        

        function sixj(j1,j2,j3,l1,l2,l3) result(sj)

            ! Computes the wigner six-j symbol with arguments
            !    j1/2 j2/2 j3/2
            !    l1/2 l2/2 l3/2
            ! using explicit algebraic expressions from Edmonds (1955/7)

            implicit none
            integer :: j1,j2,j3,l1,l2,l3
            integer :: z, zmin, zmax
            integer  :: n2, n3, n4
            integer  :: d1, d2, d3, d4
            real(kind=8) :: sj
            real(kind=8) :: fac_sum, triangle_prod, num, den

            sj = 0d0
            if ( j1 < 0) return
            if ( j2 < 0) return
            if ( j3 < 0) return
            if ( l1 < 0) return
            if ( l2 < 0) return
            if ( l3 < 0) return
            if ( j1 < abs(j2-j3) .OR. j1 > j2+j3 ) return
            if ( j1 < abs(l2-l3) .OR. j1 > l2+l3 ) return
            if ( l1 < abs(j2-l3) .OR. l1 > j2+l3 ) return
            if ( l1 < abs(l2-j3) .OR. l1 > l2+j3 ) return            

            n2 = (j1+j2+l1+l2)/2
            n3 = (j2+j3+l2+l3)/2
            n4 = (j3+j1+l3+l1)/2

            d1 = (j1+j2+j3)/2
            d2 = (j1+l2+l3)/2
            d3 = (l1+j2+l3)/2
            d4 = (l1+l2+j3)/2

            triangle_prod = triangle(j1,j2,j3) &
                           *triangle(j1,l2,l3) &
                           *triangle(l1,j2,l3) &
                           *triangle(l1,l2,j3)

            zmin = max(d1, d2, d3, d4)
            zmax = min(n2, n3, n4)

            fac_sum = 0d0

            do z = zmin, zmax
                num = logfac(z+1)
                den = logfac(n2-z) + logfac(n3-z) + logfac(n4-z) &
                    + logfac(z-d1) + logfac(z-d2) + logfac(z-d3) + logfac(z-d4)
                fac_sum = fac_sum + (-1) ** (z) * exp(num - den)
            end do

            sj = triangle_prod * fac_sum

           return

        end function sixj

        subroutine sixj_table_init(min2j, max2j)

            implicit none
            integer, optional :: min2j, max2j
            integer :: a, b
            integer :: i, j, k, l, m, n

            integer(kind=8) :: ti, tf, clock_rate
            real(kind=8) :: dummy_real

            call system_clock(count_rate = clock_rate)
            call system_clock(count = ti)

            print*,"Initializing six-j symbol table..."
            if (.not. present(min2j)) then
                a = tablemin2j
            else
                a = min2j
                tablemin2j = a
            end if
            if (.not. present(max2j)) then
                b = tablemax2j 
            else
                b = max2j
                tablemax2j = b
            end if
            print*,"Table min. 2J:",a
            print*,"Table max. 2J:",b
            print*,'Memory required (MB):',real(sizeof(dummy_real) &
                                    * (a-b+1)**6,4)*10d0**(-6)

            if (allocated(sixj_table)) deallocate(sixj_table) 
            allocate(sixj_table(a:b,a:b,a:b,a:b,a:b,a:b))
            sixj_table = 0d0

            !$OMP PARALLEL DO private(n,m,l,k,j,i) schedule(dynamic,1)
            do n = a, b
              do m = a, b
                do l = a, b
                  do k = a, b
                    do j = a, b
                      do i = a, b
                        sixj_table(i,j,k,l,m,n) = sixj(i,j,k,l,m,n)
                      end do
                    end do
                  end do
                end do
              end do
            end do

            call system_clock(count = tf)
            print*,'Table has been saved to memory.'
            print*,'Seconds to initialize:',real((tf-ti))/real(clock_rate)
        end subroutine sixj_table_init                                

        function sixj_lookup(j1,j2,j3,l1,l2,l3) result(sj)
            
            implicit none
            integer :: j1,j2,j3,l1,l2,l3
            real(kind=8) :: sj
            integer :: minj, maxj

            minj = min(j1,j2,j3,l1,l2,l3)
            maxj = max(j1,j2,j3,l1,l2,l3)
            tablemin_used = min(tablemin_used, minj)
            tablemax_used = max(tablemax_used, maxj)

            if (minj < tablemin2j .or. maxj > tablemax2j) then
                sj = sixj(j1,j2,j3,l1,l2,l3)
            else
                sj = sixj_table(j1,j2,j3,l1,l2,l3)
            end if

            return
        end function sixj_lookup

        function ninej(j1,j2,j3,j4,j5,j6,j7,j8,j9) result(nj)
            
            implicit none
            integer :: j1,j2,j3,j4,j5,j6,j7,j8,j9
            integer :: z, zmin, zmax
            real(kind=8) :: nj

            nj = 0d0

            if (j1 < 0) return
            if (j2 < 0) return
            if (j3 < 0) return
            if (j4 < 0) return
            if (j5 < 0) return
            if (j6 < 0) return
            if (j7 < 0) return
            if (j8 < 0) return
            if (j9 < 0) return

            zmin = max(abs(j1-j9), abs(j4-j8), abs(j2-j6))
            zmax = min(j1+j9, j4+j8, j2+j6)

            do z = zmin, zmax
                nj = nj + (-1) ** z * (z + 1d0) &
                    * sixj_lookup(j1,j4,j7,j8,j9, z) &
                    * sixj_lookup(j2,j5,j8,j4, z,j6) &
                    * sixj_lookup(j3,j6,j9, z,j1,j2)
            end do

            return

        end function ninej

end module wigner        
