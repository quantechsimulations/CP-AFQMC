Module WALK
use SYSTEM
use ONE_BODY
use INITIALIZATION
Implicit none
Real(8), Dimension(:,:), Allocatable:: Psi, Psi_c
Real(8), Dimension(:,:), Allocatable:: iO, iOc, eta
Real(8) :: En, Enp, Enk
Complex(8) :: chi
Real(8), Dimension(ncmx):: spc, tpc
Real(8), Allocatable:: ssc(:), nnc(:), ssij(:,:) , nnij(:,:), nij(:,:)
Integer:: Nsites, Npart, imeasure_mix, imeasure_bp, ipop_cntrl, istblz, ibp
Real(8):: dx,dy,rLx,rLy
contains

subroutine stpwlk(El,Epl,Ekl,chil,spcl,tpcl,wsum)
Implicit none
Real(8):: El,Epl,Ekl
Real(8):: chil
Real(8):: wsum
Real(8), Dimension(ncmx):: spcl, tpcl
Integer:: iwlk, isite, i, icor, jcor

Nsites=Lx*Ly*Lz
Npart=N_u+N_d
rLx = real(Lx,8)
rLy = real(Ly,8)

Allocate(Psi(1:2*Nsites,1:Npart))
if (imeasure_bp.eq.1) then
Allocate(eta(1:2*Nsites,1:Npart))
end if

El=zero
Epl=zero
Ekl=zero

dx = 0.5_8
dy = dx*sqrt(3._8)

ssloc = zero
nnloc = zero
sscor = zero
nncor = zero
ncor = zero
chil = zero
spcl = zero
tpcl = zero

Allocate(iO(1:Npart,1:Npart))
Allocate(iOc(1:Npart,1:Npart))
do iwlk = 1, N_wlk
        Psi(:,:)=Psi_wlk(:,:,iwlk)
        if (w(iwlk).gt.zero) then !1
                ! Multiply the weigth by the pre-factor correct by the Hirsch transformation
                w(iwlk)=w(iwlk)*exp(fac_norm)
                !kinetic propagation
                call k_prop(iwlk)
        end if !iwalk1
        !potential propagation site by site
        if (w(iwlk).gt.zero) then !2
                !loop over lattice sites
                do isite=1,Nsites
                        if (w(iwlk).gt.zero) call v_prop(iwlk,isite) !3
                end do ! isite

        end if !iwalk2

        !the other kinetic propagation
        if (w(iwlk).gt.zero) then !4
                call k_prop(iwlk)
        end if !iwlk4

        ! measure the energy if imeasure=1
        if (imeasure_bp.eq.1) then
                if (w(iwlk).gt.zero) then !5
                        call back_prop(iwlk)
                        Allocate(ssij(1:4*Lx+1,1:2*Ly+1))
                        Allocate(nnij(1:4*Lx+1,1:2*Ly+1))
                        Allocate(nij(1:4*Lx+1,1:2*Ly+1))
                        call measure_bp(iwlk)
                        Epl = Epl + w(iwlk)*Enp
                        Ekl = Ekl + w(iwlk)*Enk
                        chil = chil + w(iwlk)*sqrt(chi*CONJG(chi))
                        do jcor = 1, 2*Ly + 1
                        do icor = 1, 4*Lx + 1
                        sscor(icor,jcor) = sscor(icor,jcor) + w(iwlk)*ssij(icor,jcor)
                        nncor(icor,jcor) = nncor(icor,jcor) + w(iwlk)*nnij(icor,jcor)
                        ncor(icor,jcor) = ncor(icor,jcor) + w(iwlk)*nij(icor,jcor)
                        end do
                        end do
                        do icor = 1, ncmx
                        spcl(icor) = spcl(icor) + w(iwlk)*spc(icor)
                        tpcl(icor) = tpcl(icor) + w(iwlk)*tpc(icor)
                        end do
                        deallocate(ssij)
                        deallocate(nnij)
                        deallocate(nij)
                end if !iwlk5
        end if !imeasure_bp
        if (imeasure_mix.eq.1) then
                if (w(iwlk).gt.zero) then
                        call measure_mix
                        if (w(iwlk).gt.zero) then
                        El=El +w(iwlk)*En
                        end if
                end if
        end if

        if (w(iwlk).gt.zero) Psi_wlk(:,:,iwlk)=Psi(:,:)
end do !iwalk

wsum=zero
do iwlk=1,N_wlk
if (w(iwlk).gt.zero) wsum=wsum+w(iwlk)
end do


El=El/wsum

Epl = Epl/wsum
Ekl = Ekl/wsum

sscor = sscor/wsum
nncor = nncor/wsum
ncor = ncor/wsum
chil = chil/wsum
spcl = spcl/wsum
tpcl = tpcl/wsum

deallocate(Psi)
deallocate(iO)
deallocate(iOc)
if (imeasure_bp.eq.1) then
deallocate(eta)
end if


if (istblz.eq.1) call stblz
if (ipop_cntrl.eq.1) call pop_cntrl

end

subroutine k_prop(jwlk)
Implicit none
Real(8) :: o_new, o_ratio, det
Integer :: jwlk, i

Allocate(Psi_c(1:2*Nsites,1:Npart))

Psi_c = Psi
call dgemm('N','N',2*Nsites,Npart,2*Nsites,one,rho_0,2*Nsites,Psi_c, &
        &  2*Nsites,zero,Psi,2*Nsites)


call dgemm('C','N',Npart,Npart,2*Nsites,one,Psi_T,2*Nsites,Psi, &
        &  2*Nsites,zero,iO,Npart)


call matinv(iO,Npart,det)
iOc=iO
call matinv(iOc,Npart,det)

o_new=one/(det)
o_ratio=o_new/O_T(jwlk)

if (o_ratio.gt.zero) then
O_T(jwlk)=o_new
w(jwlk)=w(jwlk)*o_ratio
else
w(jwlk)=zero
end if
deallocate(Psi_c)
end

subroutine v_prop(kwlk,ksite)
Implicit none
Real(8):: Gii(2)
Real(8), Dimension(:), Allocatable:: Psi_s, PsiT_s
Real(8), Dimension(:,:), Allocatable:: Psi_2s, PsiT_2s
Real(8), Dimension(:,:,:), Allocatable:: Psin
Real(8), Dimension(:), Allocatable:: temp1
Real(8), Dimension(:,:), Allocatable:: temp2,temp3,temp5,temp6
Real(8):: o_ratio_temp(2)
Real(8):: o_ratio_real(2), sum_o
Real(8):: id(2,2), axmat(2,2), temp4(2,2)
Real(8):: o_new, det, o_ratio, at,bt,ct,dt,dets
Integer:: xspin, kwlk, ksite, k, l

Gii=zero

Allocate(temp1(1:Npart))
Allocate(Psi_2s(1:2,1:Npart))
Allocate(PsiT_2s(1:Npart,1:2))
Allocate(temp2(1:2,1:Npart))
Allocate(temp3(1:Npart,1:2))
Allocate(temp5(1:Npart,1:2))
Allocate(temp6(1:Npart,1:Npart))


do k=1, Npart
Psi_2s(1,k) = Psi(ksite,k)
PsiT_2s(k,1) = Psi_T(ksite,k)
Psi_2s(2,k) = Psi(ksite+Nsites,k)
PsiT_2s(k,2) = Psi_T(ksite+Nsites,k)
end do



!calculate the overlaps
! calculate the probability of each aux field
do k=1,2
axmat = zero
axmat(1,1) = aux_fld(1,k)-one
axmat(2,2) = aux_fld(2,k)-one

!temp2 = MATMUL(Psi_2s,iO)
call dgemm('N','N',2,Npart,Npart,one,Psi_2s,2,iO, &
        &  Npart,zero,temp2,2)
!temp4 = MATMUL(temp2,PsiT_2s) ! This is the Green-Function
call dgemm('N','N',2,2,Npart,one,temp2,2,PsiT_2s,Npart, &
        &  zero,temp4,2)

! invert the 2x2 matrix axmat
at = axmat(1,1)
bt = axmat(1,2)
ct = axmat(2,1)
dt = axmat(2,2)

det = at*dt-bt*ct
axmat(1,1)=dt/det
axmat(1,2)=-bt/det
axmat(2,1)=-ct/det
axmat(2,2)=at/det

temp4 = temp4 + axmat

at = temp4(1,1)
bt = temp4(1,2)
ct = temp4(2,1)
dt = temp4(2,2)

dets = at*dt-bt*ct

o_ratio_temp(k) = dets*det
o_ratio_real(k) = max(o_ratio_temp(k),zero)
end do

! calculate the normalization of the pdf

sum_o = o_ratio_real(1)+o_ratio_real(2)

if (sum_o.le.zero) w(kwlk)=zero

if (w(kwlk).gt.zero) then

w(kwlk)=w(kwlk)*sum_o*half !The 1/2 acounts for the pdf before import samp

! sampling the auxiliary fields
if (o_ratio_real(1)/sum_o.gt.rannyu()) then
xspin=1
else
xspin=2
end if ! x_spin

! store auxiliary field for back propagation

afs(ibp,ksite,kwlk) = xspin

! propagate the walker

do k=1,Npart
Psi(ksite,k)=Psi(ksite,k)*aux_fld(1,xspin)
Psi(ksite+Nsites,k)=Psi(ksite+Nsites,k)*aux_fld(2,xspin)
end do

! update the overlap
!call dgemm('C','N',Npart,Npart,2*Nsites,one,Psi_T,2*Nsites,Psi, &
!        &  2*Nsites,zero,iO,Npart)
!call matinv(iO,Npart,det)

!id = zero
!id(1,1) = one
!id(2,2) = one
axmat = zero
axmat(1,1) = aux_fld(1,xspin)-one
axmat(2,2) = aux_fld(2,xspin)-one

!temp2 = MATMUL(Psi_2s,iO)
call dgemm('N','N',2,Npart,Npart,one,Psi_2s,2,iO, &
       &   Npart,zero,temp2,2)
!temp3 = MATMUL(iO,PsiT_2s)
call dgemm('N','N',Npart,2,Npart,one,iO,Npart,PsiT_2s, &
       &   Npart,zero,temp3,Npart)
!temp4 = MATMUL(temp2,PsiT_2s) ! This is the Green-Function
call dgemm('N','N',2,2,Npart,one,temp2,2,PsiT_2s, &
       &   Npart,zero,temp4,2)

! invert the 2x2 matrix axmat
at = axmat(1,1)
bt = axmat(1,2)
ct = axmat(2,1)
dt = axmat(2,2)

det = at*dt-bt*ct
axmat(1,1)=dt/det
axmat(1,2)=-bt/det
axmat(2,1)=-ct/det
axmat(2,2)=at/det

temp4 = temp4 + axmat
! invert temp4
at = temp4(1,1)
bt = temp4(1,2)
ct = temp4(2,1)
dt = temp4(2,2)

det = at*dt-bt*ct
temp4(1,1)=dt/det
temp4(1,2)=-bt/det
temp4(2,1)=-ct/det
temp4(2,2)=at/det


!temp5 = MATMUL(temp3,temp4)
call dgemm('N','N',Npart,2,2,one,temp3,Npart,temp4, &
       &   2,zero,temp5,Npart)
!temp6 = MATMUL(temp5,temp2)
call dgemm('N','N',Npart,Npart,2,one,temp5,Npart,temp2, &
       &   2,zero,temp6,Npart)

iO = iO - temp6

O_T(kwlk)=O_T(kwlk)*o_ratio_temp(xspin)

end if !w(kwlk)



deallocate(temp1)
deallocate(Psi_2s)
deallocate(PsiT_2s)
deallocate(temp2)
deallocate(temp3)
deallocate(temp5)
deallocate(temp6)


end 

subroutine back_prop(jwlk)
Implicit none
Real(8), Allocatable:: temp(:,:)
Integer:: ii,jj,kk,xspin,jwlk

eta = Psi_T

do ii = 2, N_bp
! k propagate spin-up sector
!eta_up=matmul(rho_0,eta_up)
Allocate(temp(1:2*Nsites,1:Npart))
temp = eta
call dgemm('C','N',2*Nsites,Npart,2*Nsites,one,rho_0,2*Nsites,temp, &
        &  2*Nsites,zero,eta,2*Nsites)

        do jj = 1,Nsites
                xspin = afs(N_bp+1-ii,jj,jwlk)
                do kk=1,Npart
                        ! v propagate spin-up sector
                        eta(jj,kk)=eta(jj,kk)*aux_fld(1,xspin)
                        ! v propagate spin-dwn sector
                        eta(jj+Nsites,kk)=eta(jj+Nsites,kk)*aux_fld(2,xspin)
                end do
        end do

! k propagate spin-up sector
!eta_up=matmul(rho_0,eta_up)
temp = eta
call dgemm('C','N',2*Nsites,Npart,2*Nsites,one,rho_0,2*Nsites,temp, &
        &  2*Nsites,zero,eta,2*Nsites)


deallocate(temp)
if (mod(ii-1,Nspo).eq.0) call stblz_bp

end do

end

subroutine measure_mix
Implicit none
Real(8), Dimension(:,:), Allocatable:: temp, G, HG
Real(8), Dimension(:), Allocatable:: sumHG
Real(8):: E_pot, E_kin, n_int
Integer:: k,l



Allocate(temp(1:2*Nsites,1:Npart))
Allocate(G(1:2*Nsites,1:2*Nsites))
Allocate(HG(1:2*Nsites,1:2*Nsites))
Allocate(sumHG(1:2*Nsites))

call dgemm('N','N',2*Nsites,Npart,Npart,one,Psi,2*Nsites,iO, &
        &  Npart,zero,temp,2*Nsites)

call dgemm('N','C',2*Nsites,2*Nsites,Npart,one,temp,2*Nsites,Psi_T, &
        &  2*Nsites,zero,G,2*Nsites)


! potential energy
n_int=zero
do k=1,Nsites
   n_int = n_int + G(k,k)*G(k+Nsites,k+Nsites)
   n_int = n_int - G(k,k+Nsites)*G(k+Nsites,k)
end do
E_pot = U*n_int
! kinetic energy

do k=1,2*Nsites
do l=1,2*Nsites
HG(k,l) = H_0(l,k)*(G(k,l))
end do
end do
sumHG=zero
do k=1,2*Nsites
do l=1,2*Nsites
sumHG(k)=sumHG(k)+HG(k,l)
end do
end do
E_kin=0
do k=1,2*Nsites
E_kin=E_kin+sumHG(k)
end do

!total energy

En=E_kin+E_pot

deallocate(temp)
deallocate(G)
deallocate(HG)
deallocate(sumHG)

end
subroutine measure_bp(lwlk)
Implicit none
Real(8), Dimension(:,:), Allocatable:: temp, G, HG, &
                                    & iep, Psi_old
Real(8), Dimension(:), Allocatable:: sumHG
Real(8) :: c, sidex, sidey
Real(8):: E_pot, E_kin, n_int, det, rl(6,2), & 
        & rl4(12,2), rl5(6,2), rkl(2), a1(2), a2(2)
Integer, Allocatable:: xkl(:), ykl(:)
Integer:: kk,ll,lwlk,isite(2),irkl(2), ii, jj, ic
Integer:: alpha, beta, gamm, delta, eps, xi
Integer:: ia, ib, jg, jd, ke, kx, Ns
Complex(8):: sisjsk
Real(8):: r, cr(ncmx), c2
Integer:: iy, j1, j2, j3, ir

Allocate(Psi_old(1:2*Nsites,1:Npart))
Allocate(iep(1:Npart,1:Npart))
Allocate(xkl(1:Nsites))
Allocate(ykl(1:Nsites))
Allocate(temp(1:2*Nsites,1:Npart))
Allocate(G(1:2*Nsites,1:2*Nsites))
Allocate(HG(1:2*Nsites,1:2*Nsites))
Allocate(sumHG(1:2*Nsites))

Psi_old(:,:) = parents(:,:,lwlk)

call dgemm('C','N',Npart,Npart,2*Nsites,one,eta,2*Nsites,Psi_old, &
        &  2*Nsites,zero,iep,Npart)

call matinv(iep,Npart,det)

call dgemm('N','N',2*Nsites,Npart,Npart,one,Psi_old,2*Nsites,iep, &
        &  Npart,zero,temp,2*Nsites)

call dgemm('N','C',2*Nsites,2*Nsites,Npart,one,temp,2*Nsites,eta, &
        &  2*Nsites,zero,G,2*Nsites)

! potential energy
n_int=zero
do kk=1,Nsites
   n_int = n_int + G(kk,kk)*G(kk+Nsites,kk+Nsites)
   n_int = n_int - G(kk,kk+Nsites)*G(kk+Nsites,kk)
end do
E_pot = U*n_int
! kinetic energy

do kk=1,2*Nsites
do ll=1,2*Nsites
HG(kk,ll) = H_0(ll,kk)*(G(kk,ll))
end do
end do
sumHG=zero
do kk=1,2*Nsites
do ll=1,2*Nsites
sumHG(kk)=sumHG(kk)+HG(kk,ll)
end do
end do
E_kin=0
do kk=1,2*Nsites
E_kin=E_kin+sumHG(kk)
end do

Enp=E_pot

Enk=E_kin

 !correlations between pairs

ssij = zero
nnij = zero
nij = zero

sidex = float(Lx)
sidey = float(Ly)*dy

        do ll = 1, Nsites
        do kk = 1, Nsites
! we are computing the modulus squared of nq, no closest image convention is required
                rkl(1) = sites(ll,1) - sites(kk,1)
                rkl(2) = sites(ll,2) - sites(kk,2)
                irkl(1) = nint(rkl(1)/dx) + 2*Lx + 1
                irkl(2) = nint(rkl(2)/dy) + Ly + 1

                if (ll.eq.kk) then

                c = G(ll,ll)+G(ll+Nsites,ll+Nsites)
                ssij(irkl(1),irkl(2)) = ssij(irkl(1),irkl(2)) + c  &
              &         +two*G(ll,ll+Nsites)*G(ll+Nsites,ll)-two*G(ll,ll)*G(ll+Nsites,ll+Nsites)

                nnij(irkl(1),irkl(2)) = nnij(irkl(1),irkl(2)) + c - c*c &
              &         -two*G(ll,ll+Nsites)*G(ll+Nsites,ll)+two*G(ll,ll)*G(ll+Nsites,ll+Nsites)
                nij(irkl(1),irkl(2)) = nij(irkl(1),irkl(2)) + c
                else
                c = G(ll,ll)*G(kk,kk)-G(ll,kk)*G(kk,ll) &
              &   + G(ll+Nsites,ll+Nsites)*G(kk+Nsites,kk+Nsites) &
              &   - G(ll+Nsites,kk+Nsites)*G(kk+Nsites,ll+Nsites)

                ssij(irkl(1),irkl(2)) = ssij(irkl(1),irkl(2)) +  c - G(ll,ll)*G(kk+Nsites,kk+Nsites) &
              &          +  G(ll,kk+Nsites)*G(kk+Nsites,ll) &
              &          -  G(ll+Nsites,ll+Nsites)*G(kk,kk) &
              &          +  G(ll+Nsites,kk)*G(kk,ll+Nsites)
                nnij(irkl(1),irkl(2)) =  nnij(irkl(1),irkl(2)) - G(ll,kk)*G(kk,ll) &
              &          - G(ll+Nsites,kk+Nsites)*G(kk+Nsites,ll+Nsites) &
              &          - G(ll,kk+Nsites)*G(kk+Nsites,ll) &
              &          - G(ll+Nsites,kk)*G(kk,ll+Nsites)
               nij(irkl(1),irkl(2)) = nij(irkl(1),irkl(2)) + G(ll,kk) + &
              &                       G(ll+Nsites,kk+Nsites)
                end if

        end do
        end do

        ! chiral order parameter

! sum over upsidedown triangles

chi = zero
kk = 0
Ns = Nsites
do jj = 1, Ly-1
        do ii = (jj-1)*Lx + 1  , jj*Lx -1
                do alpha = 1,2
                do beta = 1,2
                do gamm = 1,2
                do delta = 1,2
                do eps = 1,2
                do xi = 1,2
                ia = (alpha-1)*Ns + ii+Lx
                ib = (beta-1)*Ns + ii+Lx
                jg = (gamm-1)*Ns + ii
                jd = (delta-1)*Ns + ii
                ke = (eps-1)*Ns + ii+1
                kx = (xi-1)*Ns + ii+1

                sisjsk = sgmx(alpha,beta)*(sgmy(gamm,delta)*sgmz(eps,xi)-sgmz(gamm,delta)*sgmy(eps,xi)) &
              &        + sgmy(alpha,beta)*(sgmz(gamm,delta)*sgmx(eps,xi)-sgmx(gamm,delta)*sgmz(eps,xi)) &
              &        + sgmz(alpha,beta)*(sgmx(gamm,delta)*sgmy(eps,xi)-sgmy(gamm,delta)*sgmx(eps,xi))

                c = G(ia,ib)*(G(jg,jd)*G(ke,kx)-G(jg,kx)*G(ke,jd)) &
              &   + G(ia,jd)*(G(ke,ib)*G(jg,kx)-G(jg,ib)*G(ke,kx)) &
              &   + G(ia,kx)*(G(jg,ib)*G(ke,jd)-G(ke,ib)*G(jg,jd))

                chi = chi + sisjsk*c
                end do
                end do
                end do
                end do
                end do
                end do
               kk = kk+1
              end do
              end do

! sum over up triangles

do jj = 2 , Ly
        do ii = (jj-1)*Lx + 1  , jj*Lx -1
                do alpha = 1,2
                do beta = 1,2
                do gamm = 1,2
                do delta = 1,2
                do eps = 1,2
                do xi = 1,2
                ia = (alpha-1)*Ns + ii+1-Lx
                ib = (beta-1)*Ns + ii+1-Lx
                jg = (gamm-1)*Ns + ii+1
                jd = (delta-1)*Ns + ii+1
                ke = (eps-1)*Ns + ii
                kx = (xi-1)*Ns + ii

                sisjsk = sgmx(alpha,beta)*(sgmy(gamm,delta)*sgmz(eps,xi)-sgmz(gamm,delta)*sgmy(eps,xi)) &
              &        + sgmy(alpha,beta)*(sgmz(gamm,delta)*sgmx(eps,xi)-sgmx(gamm,delta)*sgmz(eps,xi)) &
              &        + sgmz(alpha,beta)*(sgmx(gamm,delta)*sgmy(eps,xi)-sgmy(gamm,delta)*sgmx(eps,xi))

                c = G(ia,ib)*(G(jg,jd)*G(ke,kx)-G(jg,kx)*G(ke,jd)) &
              &   + G(ia,jd)*(G(ke,ib)*G(jg,kx)-G(jg,ib)*G(ke,kx)) &
              &   + G(ia,kx)*(G(jg,ib)*G(ke,jd)-G(ke,ib)*G(jg,jd))

                chi = chi + sisjsk*c
                end do
                end do
                end do
                end do
                end do
                end do
               kk = kk+1
              end do
              end do

                chi = chi/float(kk)

spc = zero
tpc = zero
cr = zero
do kk = 2 , Ly-1
        do ii = (kk-1)*Lx + 1  , kk*Lx -1
                do ll = 2, Ly-1
                        do jj = (ll-1)*Lx+1, ll*Lx-1
                        
                        if (mod(kk,2).eq.0) then
                        iy = ii-Lx
                        else
                        iy = ii-Lx+1
                        end if
                        if (mod(ll,2).eq.0) then
                        j1 = jj-Lx
                        j3 = jj+Lx
                        else
                        j1 = jj-Lx+1
                        j3 = jj+Lx+1
                        end if
                        j2 = jj+1

                        rkl(1) = sites(ii,1) - sites(jj,1)
                        rkl(2) = sites(ii,2) - sites(jj,2)
                        
                        !PBC
                        rkl(1) = rkl(1) - sidex*nint(rkl(1)/sidex)
                        rkl(2) = rkl(2) - sidey*nint(rkl(2)/sidey)

                        r = sqrt(rkl(1)**2+rkl(2)**2)
                        ir = nint(r/dr) + 1
                        cr(ir) = cr(ir) + one

                        c = G(iy+Ns,j1+Ns)*G(ii,jj) - G(iy+Ns,jj)*G(ii,j1+Ns) &
                       &  + G(iy,j1)*G(ii+Ns,jj+Ns) - G(iy,jj+Ns)*G(ii+Ns,j1)
                        c2= G(iy+Ns,j1)*G(ii,jj+Ns) - G(iy+Ns,jj+Ns)*G(ii,j1) &
                       &  + G(iy,j1+Ns)*G(ii+Ns,jj) - G(iy,jj)*G(ii+Ns,j1+Ns)

                        spc(ir) = spc(ir) + abs(c-c2)
                        tpc(ir) = tpc(ir) + abs(c+c2)

                        c = G(iy+Ns,j2+Ns)*G(ii,jj) - G(iy+Ns,jj)*G(ii,j2+Ns) &
                       &  + G(iy,j2)*G(ii+Ns,jj+Ns) - G(iy,jj+Ns)*G(ii+Ns,j2)
                        c2= G(iy+Ns,j2)*G(ii,jj+Ns) - G(iy+Ns,jj+Ns)*G(ii,j2) &
                       &  + G(iy,j2+Ns)*G(ii+Ns,jj) - G(iy,jj)*G(ii+Ns,j2+Ns)

                        spc(ir) = spc(ir) + abs(c-c2)
                        tpc(ir) = tpc(ir) + abs(c+c2)

                        c = G(iy+Ns,j3+Ns)*G(ii,jj) - G(iy+Ns,jj)*G(ii,j3+Ns) &
                       &  + G(iy,j3)*G(ii+Ns,jj+Ns) - G(iy,jj+Ns)*G(ii+Ns,j3)
                        c2= G(iy+Ns,j3)*G(ii,jj+Ns) - G(iy+Ns,jj+Ns)*G(ii,j3) &
                       &  + G(iy,j3+Ns)*G(ii+Ns,jj) - G(iy,jj)*G(ii+Ns,j3+Ns)

                        spc(ir) = spc(ir) + abs(c-c2)
                        tpc(ir) = tpc(ir) + abs(c+c2)
 
 
                        end do
                end do
        end do
end do

do ii=1, ncmx
        if (cr(ii).gt.zero) then
                spc(ii) = spc(ii)/(six*cr(ii))
                tpc(ii) = tpc(ii)/(six*cr(ii))
        end if
end do

deallocate(temp)
deallocate(G)
deallocate(HG)
deallocate(sumHG)
deallocate(iep)
deallocate(xkl)
deallocate(ykl)

end

subroutine pop_cntrl
Implicit none
Real(8), Allocatable:: Psi_n(:,:,:), o_new(:), parents_n(:,:,:)
Real(8):: sum_w, d
Integer, Allocatable:: afs_n(:,:,:)
Integer:: l,k,n,nw

Allocate(Psi_n(1:2*Nsites,1:Npart,1:N_wlk))
Allocate(parents_n(1:2*Nsites,1:Npart,1:N_wlk))
Allocate(afs_n(1:N_bp,1:Nsites,1:N_wlk))
Allocate(o_new(1:N_wlk))

sum_w=zero
do k=1,N_wlk
sum_w=sum_w+w(k)
end do

d = float(N_wlk)/sum_w

sum_w=-rannyu()
nw=0

do k=1, N_wlk
sum_w=sum_w+w(k)*d
n=ceiling(sum_w)
do l=nw+1,n
Psi_n(:,:,l)=Psi_wlk(:,:,k)
! update the back-propagation storage
parents_n(:,:,l)=parents(:,:,k)
afs_n(:,:,l) = afs(:,:,k)
o_new(l)=O_T(k)
end do
nw=n
end do

Psi_wlk=Psi_n
parents = parents_n
afs = afs_n
O_T=o_new
w=one

deallocate(Psi_n)
deallocate(parents_n)
deallocate(afs_n)
deallocate(o_new)


end

subroutine stblz
Implicit none
Real(8), Allocatable:: H(:,:),QR(:,:),t(:),w(:),Q(:,:),R(:,:)
Integer:: k,l,lwlk,lwork,info
Real(8):: dr


Allocate(H(1:2*Nsites,1:Npart))
Allocate(QR(1:2*Nsites,1:Npart))
Allocate(Q(1:2*Nsites,1:2*Nsites))
Allocate(R(1:Npart,1:Npart))
Allocate(t(1:Npart))

! finding defining lworku and lword
Allocate(w(1:Npart))
QR(:,:) = Psi_wlk(:,:,1)
call dgeqrf(2*Nsites,Npart,QR,2*Nsites,t,w,-1,info)
lwork = int(w(1))
deallocate(w)

Allocate(w(1:lwork))

!>>>>

do lwlk=1,N_wlk

QR(:,:) = Psi_wlk(:,:,lwlk)
call dgeqrf(2*Nsites,Npart,QR,2*Nsites,t,w,lwork,info)

H=zero
R=zero

do k=1,Npart-1
do l=k+1,Npart
R(k,l)=QR(k,l)
end do
end do
do k=1,Npart
R(k,k)=QR(k,k)
end do

do k=1,Npart
do l=1,2*Nsites

if (l.lt.k) then
H(l,k)=zero
else if (l.eq.k) then
H(l,k)=one
else
H(l,k)=QR(l,k)
end if
end do
end do

call dorg2r(2*Nsites,Npart,Npart,H,2*Nsites,t,w,info)

Psi_wlk(:,:,lwlk)=H(:,:)

call matinv(R,Npart,dr)
O_T(lwlk)=O_T(lwlk)/(dr)

end do!lwlk

deallocate(H)
deallocate(QR)
deallocate(Q)
deallocate(R)
deallocate(t)
deallocate(w)
end

subroutine stblz_bp
Implicit none
Real(8), Allocatable:: H(:,:),QR(:,:),t(:),w(:),Q(:,:),R(:,:)
Integer:: k,l,lwlk,lwork,info
Real(8):: dr


Allocate(H(1:2*Nsites,1:Npart))
Allocate(QR(1:2*Nsites,1:Npart))
Allocate(Q(1:2*Nsites,1:2*Nsites))
Allocate(R(1:Npart,1:Npart))
Allocate(t(1:Npart))

! finding defining lworku and lword
Allocate(w(1:Npart))
QR(:,:) = eta(:,:)
call dgeqrf(2*Nsites,Npart,QR,2*Nsites,t,w,-1,info)
lwork = int(w(1))
deallocate(w)

Allocate(w(1:lwork))

!>>>>

QR = eta

call dgeqrf(2*Nsites,Npart,QR,2*Nsites,t,w,lwork,info)

H=zero
R=zero

do k=1,Npart-1
do l=k+1,Npart
R(k,l)=QR(k,l)
end do
end do
do k=1,Npart
R(k,k)=QR(k,k)
end do

do k=1,Npart
do l=1,2*Nsites

if (l.lt.k) then
H(l,k)=zero
else if (l.eq.k) then
H(l,k)=one
else
H(l,k)=QR(l,k)
end if
end do
end do

call dorg2r(2*Nsites,Npart,Npart,H,2*Nsites,t,w,info)

eta=H


deallocate(H)
deallocate(QR)
deallocate(Q)
deallocate(R)
deallocate(t)
deallocate(w)

end


end module
