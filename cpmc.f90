program cpmc
use SYSTEM
use ONE_BODY
use INITIALIZATION
use WALK
IMPLICIT NONE
Real(8) :: Eloc, Eploc, Ekloc, sidex, sidey, chiloc
Real(8) :: wloc, a1(2), a2(2), rkl(2), r
Real(8), Allocatable:: sscorav(:,:), nncorav(:,:), ncorav(:,:)
Real(8), Dimension(ncmx):: spcloc, tpcloc
Integer :: ib,js, i, j, irkl(2)
      

      open(100,file='cpmc.in')

      read(100,*) U
      read(100,*) Tx
      read(100,*) Ty
      read(100,*) Tz
      read(100,*) Kx
      read(100,*) Ky
      read(100,*) Kz
      read(100,*) Lx
      read(100,*) Ly
      read(100,*) Lz
      read(100,*) N_u
      read(100,*) N_d
      read(100,*) tau
      read(100,*) irn(1)
      read(100,*) irn(2)
      read(100,*) irn(3)
      read(100,*) irn(4)
      read(100,*) N_wlk
      read(100,*) Nbeq
      read(100,*) Nbac
      read(100,*) Nspb
      read(100,*) Nspm
      read(100,*) Nsppc
      read(100,*) Nspo
      read(100,*) N_bp

      call setrn(irn)
      call init
      open(200,file='cpmc.out')
      open(300,file='ss.out')
      open(400,file='nn.out')


      write(200,*) 'CMPC run'
      write(200,*) 'Equilibration blocks'
      write(200,*) 'E', 'Err', 'Epot', 'Err', 'wsum'

      Allocate(sscorav(1:4*Lx+1,1:2*Ly+1))
      Allocate(nncorav(1:4*Lx+1,1:2*Ly+1))
      Allocate(sscor(1:4*Lx+1,1:2*Ly+1))
      Allocate(nncor(1:4*Lx+1,1:2*Ly+1))
      Allocate(ncorav(1:4*Lx+1,1:2*Ly+1))
      Allocate(ncor(1:4*Lx+1,1:2*Ly+1))




      ! initializing correlations

!lattice vectors

      dx = 0.5_8
      dy = dx*sqrt(3._8)


a1(1) = 2*dx
a1(2) = zero
a2(1) = dx
a2(2) = dy


      Eacum=zero
      Eacum2=zero
      Epacum=zero
      Epacum2=zero
      Ekacum=zero
      Ekacum2=zero
      chiacum=zero
      chiacum2=zero
      spcacum = zero
      spc2 = zero
      tpcacum = zero
      tpc2 = zero

      sscorav = zero
      nncorav = zero
      ncorav = zero
      do ib = 1,Nbeq
      Esum=zero
      Epsum=zero
      Eksum=zero
      chisum=zero
      spcsum = zero
      tpcsum = zero
      ibp = 1 ! initiate back-propagation storage
      parents = Psi_wlk
      do js = 1,Nspb

      if (mod(js,Nspm).eq.0) then
      imeasure_mix=1
      else
      imeasure_mix=0
      end if
      if (mod(js,N_bp).eq.0) then
      imeasure_bp=1
      else
      imeasure_bp=0
      end if
      if (mod(js,Nsppc).eq.0) then
      ipop_cntrl=1
      else
      ipop_cntrl=0
      end if
      if (mod(js,Nspo).eq.0) then
      istblz=1
      else
      istblz=0
      end if
      call stpwlk(Eloc,Eploc,Ekloc,chiloc,spcloc,tpcloc,wloc)
      if (imeasure_mix.eq.1) then
        Esum=Esum+Eloc
        fac_norm=(Eloc-half*U*(N_u+N_d))*tau
      end if
      if (imeasure_bp.eq.1) then
        Epsum = Epsum + Eploc
        Eksum = Eksum + Ekloc
        chisum = chisum + chiloc
        spcsum = spcsum + spcloc
        tpcsum = tpcsum + tpcloc
        ibp = 1
        parents = Psi_wlk
      else
        ibp = ibp + 1
      end if

      end do
      
      Esum=Esum/float(Nspb/Nspm)
      Eacum=Eacum+Esum
      Eacum2=Eacum2+Esum**2
      Eav=Eacum/float(ib)

      Epsum=Epsum/float(Nspb/N_bp)
      Epacum=Epacum+Epsum
      Epacum2=Epacum2+Epsum**2
      Epav=Epacum/float(ib)

      Eksum=Eksum/float(Nspb/N_bp)
      Ekacum=Ekacum+Eksum
      Ekacum2=Ekacum2+Eksum**2
      Ekav=Ekacum/float(ib)

      chisum=chisum/float(Nspb/N_bp)
      chiacum=chiacum+chisum
      chiacum2=chiacum2+chisum**2
      chiav=chiacum/float(ib)

      spcsum=spcsum/float(Nspb/N_bp)
      spcacum=spcacum+spcsum
      spc2(:)=spc2(:)+spcsum(:)**2
      spairc=spcacum/float(ib)

      tpcsum=tpcsum/float(Nspb/N_bp)
      tpcacum=tpcacum+tpcsum
      tpc2(:)=tpc2(:)+tpcsum(:)**2
      tpairc=tpcacum/float(ib)

      if (ib.eq.1) then
      Eerr=zero
      Eperr=zero
      Ekerr=zero      
      else
      Eav2=Eacum2/float(ib)
      Eerr = sqrt(abs(Eav2-Eav**2)/(ib-1))
 
      Epav2=Epacum2/float(ib)
      Eperr = sqrt(abs(Epav2-Epav**2)/(ib-1))
      
      Ekav2=Ekacum2/float(ib)
      Ekerr = sqrt(abs(Ekav2-Ekav**2)/(ib-1))

      chiav2=chiacum2/float(ib)
      chierr = sqrt(abs(chiav2-chiav**2)/(ib-1))

      spcav2 = spc2/float(ib)
      spcerr(:) = sqrt(abs(spcav2(:)-spairc(:)**2)/(ib-1))

      tpcav2 = tpc2/float(ib)
      tpcerr(:) = sqrt(abs(tpcav2(:)-tpairc(:)**2)/(ib-1))
 
      end if
      write(200,'(7f10.5)') Eav, Eerr, Epav, Eperr, Ekav, Ekerr, wloc/1000
      end do   

      write(200,*) 'Accumulation blocks'
      write(200,*) 'E', 'Err', 'Ep', 'Err', 'wsum'

      Eacum=zero
      Eacum2=zero
      Epacum=zero
      Epacum2=zero
      Ekacum=zero
      Ekacum2=zero
      chiacum=zero
      chiacum2=zero
      spcacum = zero
      spc2 = zero
      tpcacum = zero
      tpc2 = zero


      ssacum = zero
      nnacum = zero
      ssacum2 = zero
      nnacum2 = zero
      sserr = zero
      nnerr = zero
      sscorav = zero
      nncorav = zero
      ncorav = zero


      do ib = 1,Nbac

      Esum=zero
      Epsum=zero
      Eksum=zero
      sssum = zero
      nnsum = zero
      chisum = zero
      spcsum = zero
      tpcsum = zero
      ! initiate back-propagation storage
      ibp = 1
      parents = Psi_wlk

      do js = 1,Nspb

      if (mod(js,Nspm).eq.0) then
      imeasure_mix=1
      else
      imeasure_mix=0
      end if
      if (mod(js,N_bp).eq.0) then
      imeasure_bp=1
      else
      imeasure_bp=0
      end if
      if (mod(js,Nsppc).eq.0) then
      ipop_cntrl=1
      else
      ipop_cntrl=0
      end if
      if (mod(js,Nspo).eq.0) then
      istblz=1
      else
      istblz=0
      end if

      call stpwlk(Eloc,Eploc,Ekloc,chiloc,spcloc,tpcloc,wloc)
      if (imeasure_mix.eq.1) then
        Esum=Esum+Eloc
        fac_norm=(Eloc-half*U*(N_u+N_d))*tau
      end if
      if (imeasure_bp.eq.1) then
        Epsum = Epsum + Eploc
        Eksum = Eksum + Ekloc
        chisum = chisum + chiloc
        spcsum = spcsum + spcloc
        tpcsum = tpcsum + tpcloc
        sssum = sssum + ssloc
        nnsum = nnsum + nnloc
        sscorav = sscorav+sscor
        nncorav = nncorav+nncor
        ncorav = ncorav + ncor
        ibp = 1
        parents = Psi_wlk
      else
        ibp = ibp + 1
      end if



      end do
      
      !if (mod(ib,5).eq.0) call snapshot(ib/5)
      !call snapshot(ib)

      Esum=Esum/float(Nspb/Nspm)
      Eacum=Eacum+Esum
      Eacum2=Eacum2+Esum**2
      Eav=Eacum/float(ib)

      Epsum=Epsum/float(Nspb/N_bp)
      Epacum=Epacum+Epsum
      Epacum2=Epacum2+Epsum**2
      Epav=Epacum/float(ib)

      Eksum=Eksum/float(Nspb/N_bp)
      Ekacum=Ekacum+Eksum
      Ekacum2=Ekacum2+Eksum**2
      Ekav=Ekacum/float(ib)

      chisum=chisum/float(Nspb/N_bp)
      chiacum=chiacum+chisum
      chiacum2=chiacum2+chisum**2
      chiav=chiacum/float(ib)

      spcsum=spcsum/float(Nspb/N_bp)
      spcacum=spcacum+spcsum
      spc2(:)=spc2(:)+spcsum(:)**2
      spairc=spcacum/float(ib)

      tpcsum=tpcsum/float(Nspb/N_bp)
      tpcacum=tpcacum+tpcsum
      tpc2(:)=tpc2(:)+tpcsum(:)**2
      tpairc=tpcacum/float(ib)

      sssum=sssum/float(6*Nsites*Nspb/N_bp)
      nnsum=nnsum/float(6*Nsites*Nspb/N_bp)
      ssacum=ssacum+sssum
      nnacum=nnacum+nnsum
      do i=1,5
      ssacum2(i)=ssacum2(i)+sssum(i)**2
      nnacum2(i)=nnacum2(i)+nnsum(i)**2
      end do
      ssav=ssacum/float(ib)
      nnav=nnacum/float(ib)

      if (ib.eq.1) then
      Eerr=zero
      Eperr=zero
      Ekerr=zero
      sserr=zero
      else
      Eav2=Eacum2/float(ib)
      Eerr = sqrt(abs(Eav2-Eav**2)/(ib-1))
      Epav2=Epacum2/float(ib)
      Eperr = sqrt(abs(Epav2-Epav**2)/(ib-1))
      Ekav2=Ekacum2/float(ib)
      Ekerr = sqrt(abs(Ekav2-Ekav**2)/(ib-1))
      ssav2=ssacum2/float(ib)
      nnav2=nnacum2/float(ib)
      spcav2 = spc2/float(ib)
      spcerr(:) = sqrt(abs(spcav2(:)-spairc(:)**2)/(ib-1))

      tpcav2 = tpc2/float(ib)
      tpcerr(:) = sqrt(abs(tpcav2(:)-tpairc(:)**2)/(ib-1))

      do i=1,5
      sserr(i) = sqrt(abs(ssav2(i)-ssav(i)**2)/(ib-1))
      nnerr(i) = sqrt(abs(nnav2(i)-nnav(i)**2)/(ib-1))
      end do
      chiav2=chiacum2/float(ib)
      chierr = sqrt(abs(chiav2-chiav**2)/(ib-1))
      end if
      write(200,'(7f10.5)') Eav, Eerr, Epav, Eperr, chiav, chierr, wloc/1000
      end do   

      write(200,*) 'FINAL RESULTS'
      write(200,*) 'Energy =', Eav, '+-', Eerr
      write(200,*) 'Doublon dens =', Epav/(U*Nsites), '+-', Eperr/(U*Nsites)
      do i=1,5
      write(200,*) 'SzSz corr =', ssav(i), '+-', sserr(i), i
      write(200,*) 'nn corr=', nnav(i), '+-', nnerr(i), i
      end do
      write(200,*) 'Kin energy =', Ekav, '+-', Ekerr
      close(200)


      sscorav = sscorav/float(Nbac*Nspb/N_bp)
      nncorav = nncorav/float(Nbac*Nspb/N_bp)
      ncorav = ncorav/float(Nbac*Nspb/N_bp)

      open(333,file='nk.out')

      do j = 1, 2*Ly+1

        rkl(2) = (j - Ly - 1)*dy
        do i = 1, 4*Lx+1
                ! no closest image here
                rkl(1) = (i - 2*Lx - 1)*dx

                write(300,*) rkl(1), rkl(2), sscorav(i,j)
                write(400,*) rkl(1), rkl(2), nncorav(i,j)
                write(333,*) rkl(1), rkl(2), ncorav(i,j)
        end do
      end do

      open(3, file='spair.out')
      open(4, file='tpair.out')

      do i=1,ncmx
        if (spairc(i) .ne. zero .OR. tpairc(i) .ne. zero) then
        r = float(i-1)*dr
        write(3,*) r, spairc(i), spcerr(i)
        write(4,*) r, tpairc(i), tpcerr(i)
        end if
      end do

      end program
