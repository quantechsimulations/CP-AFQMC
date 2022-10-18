Module ONE_BODY
use SYSTEM
Implicit none

Contains 

subroutine free_hamiltonian(H_K)
Real(8), Dimension(:,:),Allocatable :: H_K
Real(8) :: dx, dy, x, y, xij, yij, r2, rLx, rLy, rLz
Integer :: Nsites
Integer :: i, j, q

Nsites = Lx*Ly*Lz
Allocate(H_K(1:2*Nsites,1:2*Nsites))
Allocate(sites(1:Nsites,1:2))

H_K=zero;

        dx=0.5_8
        dy=sqrt(3._8)/2._8
        rLx = real(Lx,8)
        rLy = real(Ly,8)
        rLz = real(Lz,8)
        
        
        q = 1
        do i = 1, Ly
                do j = 1, Lx
                        if (mod(i,2) .eq. 0) then
                                x = j-1
                        else
                                x = j-1 + dx
                        end if
                        y = (1-i)*dy
                        sites(q,1)=x
                        sites(q,2)=y
                        q=q+1
                end do
        end do

        
        do i=1,Nsites
                do j=1,Nsites

                        xij = abs(sites(i,1)-sites(j,1))
                        xij = min(xij,rLx-xij)
                        yij = abs(sites(i,2)-sites(j,2))
                        yij = min(yij,rLy*dy-yij)

                        r2=xij**2+yij**2
                        if (r2 .le. 1.1) then
                                if (r2.gt. 0.9) then
                                        H_K(i,j)=H_K(i,j)-1
                                        H_K(i+Nsites,j+Nsites)=H_K(i+Nsites,j+Nsites)-1
                                end if
                        end if
                end do
        end do

end

end Module
