Module SYSTEM
Implicit none
Integer, Parameter:: ncmx = 500
Real(8), Parameter:: dr = 0.1
Real(8), Parameter:: zero=0._8, one=1._8, two=2._8, three=3._8, four=4._8, &
                   & five=5._8, six=6._8, seven=7._8, eight=8._8, nine=9._8, &
                   & ten=10._8, half=one/two
Real(8), Parameter:: pi=4._8*atan(one)
Real(8) :: U, Tx, Ty, Tz, Kx, Ky, Kz, tau, fac_norm, hgamma, aux_fld(2,2)
Real(8), Dimension(:,:), Allocatable:: H_0, rho_0, Psi_T
Real(8), Dimension(:,:,:), Allocatable:: Psi_wlk
Real(8), Dimension(:), Allocatable:: O_T
Real(8), Dimension(:), Allocatable:: w
Real(8) :: E_T, Esum, Epsum, Eksum, chisum
Real(8) :: Eacum, Eacum2, Eav, Eav2, Eerr
Real(8) :: Epacum, Epacum2, Epav, Epav2, Eperr
Real(8) :: Ekacum, Ekacum2, Ekav, Ekav2, Ekerr
Real(8) :: chiacum, chiacum2, chiav, chiav2, chierr
Real(8), Dimension (5) :: nnsum, sssum, nnacum, ssacum
Real(8), Dimension (5) :: nnav, nnav2, nnerr, nnacum2
Real(8), Dimension (5) :: ssav, ssav2, sserr, ssacum2
Real(8), Dimension (5) :: ssloc, nnloc
Real(8), Allocatable :: sscor(:,:) , nncor(:,:), ncor(:,:)
Complex(8), Dimension(2,2):: sgmx, sgmy, sgmz
Integer :: Lx, Ly, Lz, N_u, N_d, N_wlk, irn(4)
Integer :: Nbeq,Nbac,Nspb,Nspm,Nsppc,Nspo
Integer :: iuhf
Integer, Dimension(:), Allocatable:: sup_old, sdwn_old
Real(8) :: p_old
Real(8), Allocatable, Dimension(:):: pup_old, pdwn_old, oru_old, ord_old
! back-propagation variables
Integer, Allocatable, Dimension(:,:,:) :: afs ! store auxiliary fields
Real(8), Allocatable, Dimension(:,:,:):: parents  ! store population of parents
Integer:: N_bp ! number of back-propagation projections
Real(8), Dimension(:,:), Allocatable:: sites
Real(8), Dimension(ncmx):: spairc, spcsum, spcacum, spc2, spcerr, spcav2
Real(8), Dimension(ncmx):: tpairc, tpcsum, tpcacum, tpc2, tpcerr, tpcav2


end Module
