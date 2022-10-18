Module INITIALIZATION
use SYSTEM
use ONE_BODY
Implicit none
Real(8), Allocatable:: psiu(:,:), nu(:), psid(:,:), nd(:)
Integer, Dimension(4), save :: mm=(/ 502,1521,4071,2107/), &
                             & ll=(/   0,   0,   0,1/)
contains

subroutine init
Implicit none
Real(8), parameter:: epse=10._8**(-8), epsn=10._8**(-4)
Real(8), Dimension(:,:), Allocatable:: H_0c, rho_0c
Real(8), Dimension(:,:), Allocatable:: rvec,lvec,irvec
Real(8), Dimension(:,:), Allocatable::psith, g, HG
Real(8), Dimension(:), Allocatable:: sumHG
Real(8), Dimension(:), Allocatable:: work, eval, evalz
Real(8) :: E_K, E_V, detr
Real(8) :: n_int
Integer, Dimension(:), Allocatable:: ipiv
Integer :: Nsites, Npart
Integer :: lwork, info, i,j,k

Nsites = Lx*Ly*Lz
Npart = N_u + N_d

! constructing the free-particle Hamiltonian
call free_hamiltonian(H_0)
Allocate(H_0c(1:2*Nsites,1:2*Nsites))
Allocate(rho_0(1:2*Nsites,1:2*Nsites))
Allocate(rho_0c(1:2*Nsites,1:2*Nsites))
Allocate(lvec(1:2*Nsites,1:2*Nsites))
Allocate(rvec(1:2*Nsites,1:2*Nsites))
Allocate(irvec(1:2*Nsites,1:2*Nsites))
Allocate(eval(1:2*Nsites))
Allocate(evalz(1:2*Nsites))

H_0c = H_0 ! copy of H_0
H_0 = -half*tau*H_0

! finding optimal lwork
lwork = 4*Nsites
Allocate(work(1:lwork))
call dgeev('N','V',2*Nsites,H_0,2*Nsites,eval,evalz,lvec,2*Nsites,rvec,2*Nsites, &
        &  work,-1,info)
lwork = int(work(1))
deallocate(work)
!...
Allocate(work(1:lwork))
Allocate(ipiv(1:2*Nsites))

! computing exp(-0.5*tau*H_0) using diagonalization: if A = VDV**(-1)
! then exp(A)=Vexp(D)V**{-1}; V is the enigen vector's matrix
! and D is the eigenvalue's one.

H_0 = H_0c
H_0 = -half*tau*H_0

!Finding the eigenvectors and eigenvalues

call dgeev('N','V',2*Nsites,H_0,2*Nsites,eval,evalz,Lvec,2*Nsites,rvec,2*Nsites, &
        &  work,lwork,info)

! constructing the diagonal matrix exp(D)

rho_0 = zero
do i=1, 2*Nsites
rho_0(i,i) = exp(eval(i)) 
end do
rho_0c = rho_0 ! copy of rho_0
!performing the product V*exp(D)
call dgemm('N','N',2*Nsites,2*Nsites,2*Nsites,one,rvec,2*Nsites,rho_0c, &
        &  2*Nsites,zero,rho_0,2*Nsites)
rho_0c=rho_0
! inversion of rvec
irvec = rvec
call matinv(irvec,2*Nsites,detr) 

! calculating the product (Vexp(D))*V**(-1)

call dgemm('N','N',2*Nsites,2*Nsites,2*Nsites,one,rho_0c,2*Nsites,irvec, &
        &  2*Nsites,zero,rho_0,2*Nsites)

! Computing the autovectors of H_0 to construct Psi_T

H_0=H_0c
! read psi_T from file psi.in

open(13,file='psi.in')
!open(13,file='psi_f.in')


Allocate(Psi_T(1:2*Nsites,1:Npart))

do i=1, Npart
do j=1, 2*Nsites
 read(13,*) Psi_T(j,i) !spin-up section
end do
end do

Allocate(psith(1:N_u+N_d,1:2*Nsites))


do i=1,2*Nsites
do j=1,N_u+N_d
psith(j,i)=Psi_T(i,j)
end do
end do

Allocate(g(1:2*Nsites,1:2*Nsites))
Allocate(HG(1:2*Nsites,1:2*Nsites))
Allocate(sumHG(1:2*Nsites))

call dgemm('N','N', 2*Nsites, 2*Nsites, Npart, one, Psi_T, 2*Nsites, psith, &
        &  Npart,zero,g,2*Nsites)

! potential energy
n_int = zero
do i=1,Nsites
   n_int = n_int + g(i,i)*g(i+Nsites,i+Nsites);
   n_int = n_int - g(i,i+Nsites)*G(i+Nsites,i);
end do
E_V = U*n_int

! kinetic energy

do i=1,2*Nsites
do j=1,2*Nsites
HG(i,j) = H_0(i,j)*g(i,j)
end do
end do
sumHG=zero
do i=1,2*Nsites
do j=1,2*Nsites
sumHG(i)=sumHG(i)+HG(i,j)
end do
end do
E_K=0
do i=1,2*Nsites
E_K=E_K+sumHG(i)
end do

!total energy


E_T=E_V+E_K

! Initial configuration of walkers

Allocate(Psi_wlk(1:2*Nsites,1:Npart,1:N_wlk))

do i=1, N_wlk

Psi_wlk(:,:,i)=Psi_T(:,:)

end do

! Initial weigth's

Allocate(w(1:N_wlk))

w = one

! Initial overlap's with Psi_T

Allocate(O_T(1:N_wlk))

O_T = one

! exponent of the prefactor exp(-deltau*(-E_T)) in the ground state projector 
! fac_norm also include -0.5*U*(N_up+N_dn), the exponent of the prefactor in the
!Hirsch transformation
fac_norm=(E_T-half*U*Npart)*tau
!hgamma in Hirsch's transformation
hgamma=acosh(exp(half*tau*U))
! aux_fld is the 2x2 matrix containing all the possible values of the quantity
!exp(-gamma*s(sigma)*x_i)
aux_fld=zero
! The first index corresponds to spin up or down
! The second index corresponds to the auxiliary field x_i=1 or x_i=-1
 do i=1,2
    do j=1,2
        aux_fld(i,j)=exp(hgamma*(-1)**(i+j))
    end do
end do

Allocate (afs(1:N_bp,1:Nsites,1:N_wlk))
Allocate (parents(1:2*Nsites,1:Npart,1:N_wlk))

!deallocate arrays

deallocate(H_0c)
deallocate(rho_0c)
deallocate(lvec)
deallocate(rvec)
deallocate(irvec)
deallocate(eval)
deallocate(evalz)
deallocate(g)
deallocate(HG)
deallocate(sumHG)

!Pauli matrices

sgmx=zero
sgmx(1,2) = one
sgmx(2,1) = one

sgmy=zero
sgmy(1,2) = cmplx(zero,-one)
sgmy(2,1) = cmplx(zero,one)

sgmz=zero
sgmz(1,1) = one
sgmz(2,2) = -one

end

      subroutine matinv(aa,nsub,det)
      use SYSTEM
      Implicit none
      Integer, parameter :: nmax=65536
      Integer, parameter :: nbmx=256
      Integer, parameter :: nbmxs=nbmx*nbmx

!
! routine to calculate inverse and determinant of matrix a
! assumed to be dimensioned a(nsub,nsub).
! the matrix a is replaced by its inverse.
!
      Real(8) :: aa(nbmxs), atemp(nmax), adiag, adiagi
      Real(8) :: t, determ, det
      Integer :: ipivot(nmax), idiag
      Integer :: n, nsub, iclm, itemp
      Integer :: i, j, k, jn

      n=nsub
      do 5 i=1,n
    5 ipivot(i)=i
!
! initialize determinant
!
      determ=one
!
! loop through columns
!
      iclm=-n
      do 10 i=1,n
      iclm=iclm+n
!
! loop through rows and select row with largest element
!
      adiag=aa(ipivot(i)+iclm)
      idiag=i
      do 15 k=i+1,n
      if (abs(aa(ipivot(k)+iclm)).gt.abs(adiag)) then
                                          adiag=aa(ipivot(k)+iclm)
                                          idiag=k
                                          endif
   15 continue
!
! interchange pointers if row different from
! original is selected and change sign of determinant because
! of interchange
!
      if (idiag.ne.i) then
                      determ=-determ
                      itemp=ipivot(i)
                      ipivot(i)=ipivot(idiag)
                      ipivot(idiag)=itemp
                      endif
!
! update determinant
!
      determ=adiag*determ
      aa(ipivot(i)+iclm)=one
      adiagi=one/adiag
!
! scale row by inverse diagonal
!
      call sscal(n,adiagi,aa(ipivot(i)),n)
!
! loop through other rows
! if not current row, then row reduce
!
      do 20 j=1,n
      if (j.ne.ipivot(i)) then
                  t=-aa(j+iclm)
                  aa(j+iclm)=zero
                  call thissaxpy(n,t,aa(ipivot(i)),n,aa(j),n)
                  endif
   20 continue
   10 continue
!
! interchange elements to unpivot inverse matrix
! the following is equivalent to:
!      anew(i,ipivot(j))=aold(ipivot(i),j)
!
      jn=-n
      do 30 j=1,n
      jn=jn+n
      do 40 i=1,n
   40 atemp(i)=aa(i+jn)
      do 50 i=1,n
   50 aa(i+jn)=atemp(ipivot(i))
   30 continue
      do 55 j=1,n
   55 ipivot(j)=(ipivot(j)-1)*n
      do 60 i=1,n
      jn=-n
      do 70 j=1,n
      jn=jn+n
   70 atemp(j)=aa(i+jn)
      do 80 j=1,n
   80 aa(i+ipivot(j))=atemp(j)
   60 continue
      det=determ
      return
      end
      subroutine sscal(n,scalor,x,nskip)
! routine to scale x by scalor
! n=number of elements in x to be scaled
! nskip=stride

      Implicit none
      !Integer, parameter :: nbmxs=nbmx*nbmx
      Real(8) :: x(1), scalor
      Integer :: n, nskip
      Integer :: i, ix

      ix=1-nskip
      do 10 i=1,n
      ix=ix+nskip
   10 x(ix)=scalor*x(ix)
      return
      end

      subroutine thissaxpy(n,aa,x,nxskip,y,nyskip)
!
! routine to calculate y=a*x+y where
! x and y are arrays, and a is a scalar
! n=number of elements in x and y to calculate
! nxskip=x stride
! nyskip=y stride
!
      Implicit none
      !Integer, parameter :: nbmxs=nbmx*nbmx
      Real(8) :: x(1),y(1),aa
      Integer :: n, nxskip, nyskip
      Integer :: i, ix, iy

      ix=1-nxskip
      iy=1-nyskip
      do 10 i=1,n
      ix=ix+nxskip
      iy=iy+nyskip
   10 y(iy)=aa*x(ix)+y(iy)
      return
      end

        Real(8) function rannyu()
        Real(8), Parameter :: ooto12=1.0_8/4096.0_8
        Integer, Parameter :: itwo12=4096
        Integer :: i1, i2, i3, i4
        i1=ll(1)*mm(4)+ll(2)*mm(3)+ll(3)*mm(2)+ll(4)*mm(1)
        i2=ll(2)*mm(4)+ll(3)*mm(3)+ll(4)*mm(2)
        i3=ll(3)*mm(4)+ll(4)*mm(3)
        i4=ll(4)*mm(4)
        ll(4)=mod(i4,itwo12)
        i3=i3+i4/itwo12
        ll(3)=mod(i3,itwo12)
        i2=i2+i3/itwo12
        ll(2)=mod(i2,itwo12)
        ll(1)=mod(i1+i2/itwo12,itwo12)
        rannyu=ooto12*(float(ll(1)) + ooto12*(float(ll(2)) &
       & +ooto12*(float(ll(3)) + ooto12*(float(ll(4))))))
      end function rannyu

      subroutine setrn(iseed)
        Integer, Dimension(4) :: iseed
        Integer :: i, j
        Integer :: isn, ipe, ipd, id
        do j = 1, 4
          isn = 0
          do i = 1, 4
            ipe = 4 - i
            ipd = 10 ** ipe
            id = iseed(j) / ipd
            isn = isn + id * 8 ** ipe
            iseed(j) = iseed(j) - id * ipd
          end do
          iseed(j) = isn
        end do
        ll=iseed
        ll(4)=2*(ll(4)/2)+1
      end subroutine setrn

      subroutine savern(iseed)
        Integer, Dimension(4) :: iseed
        Integer :: i, j
        Integer :: ipe, ipo, id, isn
        iseed=ll
        do j = 1,4
          isn = 0
          do i = 1,4
            ipe = 4 - i
            ipo = 8 ** ipe
            id = iseed(j) / ipo
            isn = isn + id * 10 ** ipe
            iseed(j) = iseed(j) - ipo * id
          end do
          iseed(j) = isn
       end do
      end subroutine savern


end Module
