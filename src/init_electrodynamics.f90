! Copyright 2021, the GITM Development Team (see srcDoc/dev_team.md for members)
! Full license can be found in LICENSE
subroutine init_electrodynamics

  use ModGITM
  use ModInputs
  use ModConstants
  use ModElectrodynamics
  use ModLinearSolver
  use ModMPI
  use ModTime
  use ModMagTrace

  implicit none

  integer, external :: jday

  integer :: i,j,k,bs, iError, iDir, iLon, iLat, iAlt, ip, im, iOff

  integer :: iEquator
  
  real :: GeoLat, GeoLon, GeoAlt, xAlt, len, ped, hal
  real :: sp_d1d1_d, sp_d2d2_d, sp_d1d2_d, sh
  real :: xmag, ymag, zmag, bmag, signz, magpot, lShell
  real :: mlatMC, mltMC, jul, shl, spl, length, kdpm_s, kdlm_s, je1_s, je2_s
  real :: kpm_s, klm_s, xstretch, ystretch   
  real :: sinIm, spp, sll, shh, scc, sccline, sppline, sllline, shhline, be3

  real :: q2, dju,dl, cD, ReferenceAlt

  real :: magloctime_local(0:nLons+1,0:nLats+1)

  real, dimension(-1:nLons+2, -1:nLats+2, -1:nAlts+2) :: &
       e_density, Vi, Ve, MeVen, MeVei, MiVin, VeOe, ViOi, &
       JuDotB, sigmap_d2d2_d, sigmap_d1d1_d, sigmap_d1d2_d, sigmah, &
       ue1, ue2, kmp, kml     !, je1, je2, ed1, ed2

  real, dimension(-1:nLons+2, -1:nLats+2, -1:nAlts+2, 3) :: &
       Gradient_GC

  real, dimension(-1:nLons+2, -1:nLats+2, -1:nAlts+2, nBlocksMax) :: tmp3D
  real, dimension(-1:nLons+2, -1:nLats+2, nBlocksMax)             :: tmp2D
  
  real :: aLat, aLon, gLat, gLon, Date, sLat, sLon, gLatMC, gLonMC

  real :: residual, oldresidual, a, tmp

  logical :: IsDone, IsFirstTime = .true., DoTestMe, Debug=.False.

  integer :: iLm, iLp, jLat, iI, MaxIteration, nIteration, iLonNoon
  integer :: nX
  integer :: iStart, iEnd, iAve

  logical :: UseNewTrace = .false.
  
  external :: matvec_gitm2

  if (Debug) &
       write(*,*)'DBG: entered init_electrodynamics: ',&
       DipoleStrength, UseDynamo, Is1D
  
  if (DipoleStrength == 0) return

  ! Used to have a return here for UseDynamo=F, but sort of need this for
  ! coupling if we are in the SWMF.

  if (Is1D) return

  call report("init_electrodynamics",1)
  call start_timing("init_electrodyn")

  iEquator = DynamoHighLatBoundary/MagLatRes + 1
  iAve = (DynamoLonAverage/2) / MagLonRes
     
  ! JMB IsFirstTime
  if (IsFirstTime) then
     IsFirstTime = .false.

     nMagLats = (2*DynamoHighLatBoundary)/MagLatRes+1
     nMagLons = 360.0 / MagLonRes

     allocate(DivJuAltMC(nMagLons+1,nMagLats), &
          SigmaHallMC(nMagLons+1,nMagLats), &
          SigmaPedersenMC(nMagLons+1,nMagLats), &
          SigmaLLMC(nMagLons+1,nMagLats), &
          SigmaPPMC(nMagLons+1,nMagLats), &
          AverageMC(nMagLons+1,nMagLats), &
          SigmaHHMC(nMagLons+1,nMagLats), &
          SigmaCCMC(nMagLons+1,nMagLats), &
          SigmaLPMC(nMagLons+1,nMagLats), &
          SigmaPLMC(nMagLons+1,nMagLats), &
          KpmMC(nMagLons+1,nMagLats), KlmMC(nMagLons+1,nMagLats), &
          KDpmMC(nMagLons+1,nMagLats), &
          KDlmMC(nMagLons+1,nMagLats), &
          MagBufferMC(nMagLons+1,nMagLats), &
          LengthMC(nMagLons+1,nMagLats), &
          GeoLatMC(nMagLons+1,nMagLats), &
          GeoLonMC(nMagLons+1,nMagLats), &
          MagLocTimeMC(nMagLons+1,nMagLats), &
          MagLonMC(nMagLons+1,nMagLats), MagLatMC(nMagLons+1,nMagLats), &
          solver_a_mc(nMagLons+1,nMagLats), &
          solver_b_mc(nMagLons+1,nMagLats), &
          solver_c_mc(nMagLons+1,nMagLats), &
          solver_d_mc(nMagLons+1,nMagLats), &
          solver_e_mc(nMagLons+1,nMagLats), &
          solver_s_mc(nMagLons+1,nMagLats), &
          deltalmc(nMagLons+1,nMagLats), deltapmc(nMagLons+1,nMagLats),  &
          dSigmaLLdlMC(nMagLons+1,nMagLats),  &
          dSigmaLPdlMC(nMagLons+1,nMagLats),  &
          dSigmaPLdpMC(nMagLons+1,nMagLats),  &
          dSigmaLLdpMC(nMagLons+1,nMagLats), SigmaCowlingMC(nMagLons+1,nMagLats), &
          dSigmaCowlingdpMC(nMagLons+1,nMagLats), dKDlmdpMC(nMagLons+1,nMagLats), &
          dSigmaPPdpMC(nMagLons+1,nMagLats), &
          dKDpmdpMC(nMagLons+1,nMagLats), dKlmdlMC(nMagLons+1,nMagLats), &
          dKDlmdlMC(nMagLons+1,nMagLats),  dKpmdpMC(nMagLons+1,nMagLats), &
          DynamoPotentialMC(nMagLons+1,nMagLats), &
          Ed1new(nMagLons+1,nMagLats), Ed2new(nMagLons+1,nMagLats), &
          OldPotMC(nMagLons+1,nMagLats), &
          stat = iError)

     allocate(SmallMagLocTimeMC(nMagLons+1,2), &
          SmallMagLatMC(nMagLons+1,2), &
          SmallPotentialMC(nMagLons+1,2), &
          stat = iError)

     DivJuAltMC = 0.0
     SigmaHallMC = 0.0
     SigmaPedersenMC = 0.0
     SigmaLLMC = 0.0
     SigmaPPMC = 0.0
     SigmaHHMC = 0.0
     SigmaCCMC = 0.0
     SigmaLPMC = 0.0
     SigmaPLMC = 0.0
     KDpmMC = 0.0
     KDlmMC = 0.0
     KpmMC = 0.0
     KlmMC = 0.0
     MagBufferMC = 0.0
     LengthMC = 0.0
     GeoLatMC = 0.0
     GeoLonMC = 0.0
     MagLocTimeMC = 0.0
     MagLonMC = 0.0
     MagLatMC = 0.0
     solver_a_mc = 0.0
     solver_b_mc = 0.0
     solver_c_mc = 0.0
     solver_d_mc = 0.0
     solver_e_mc = 0.0
     solver_s_mc = 0.0
     deltalmc = 0.0
     deltapmc = 0.0
     dSigmaLLdlMC = 0.0
     dSigmaLPdlMC = 0.0
     dSigmaPLdpMC = 0.0
     dSigmaPPdpMC = 0.0
     dKDpmdpMC = 0.0
     dKDlmdlMC = 0.0
     dKpmdpMC = 0.0
     dKlmdlMC = 0.0
     SigmaCowlingMC = 0.0 
     dSigmaCowlingdpMC = 0.0
     dSigmaLLdpMC = 0.0
     dKDlmdpMC = 0.0
     Ed1new = 0.0
     Ed2new = 0.0
     DynamoPotentialMC = 0.0
     OldPotMC = 0.0

     if (iError /= 0) then
        call CON_stop("Error allocating array DivJuAltMC")
     endif

     date = iStartTime(1) + float(iJulianDay)/float(jday(iStartTime(1),12,31))

     iStart = float(iProc)/nProcs * (nMagLons+1) + 1
     iEnd   = float(iProc+1)/nProcs * (nMagLons+1)

     GeoLatMC = -1.0e32
     GeoLonMC = -1.0e32

     MagLatMC = -1.0e32
     MagLonMC = -1.0e32

     do i=iStart,iEnd
        if (iDebugLevel > 1) &
             write(*,*) "==> Calculating Apex->Geo", i, iStart, iEnd
        do j=1,nMagLats

           MagLatMC(i,j)     = float(j-1) * MagLatRes &
                - DynamoHighLatBoundary

           MagLonMC(i,j)     = 360.0 * float(i-1) / float(nMagLons)

           if (UseApex) then 

              aLat = MagLatMC(i,j)
              aLon = MagLonMC(i,j)

              if (i > 1) then
                 sLat = GeoLatMC(i-1,j)
                 sLon = GeoLonMC(i-1,j)
              elseif (j > 1) then
                 sLat = GeoLatMC(i,j-1)
                 sLon = GeoLonMC(i,j-1)
              else
                 sLat = -100.0
              endif

              ReferenceAlt = Altitude_GB(1,1,0,1)/1000.0
              AltMinIono = ReferenceAlt

              call apex_to_geo(date, aLat, aLon, ReferenceAlt, &
                   gLat, gLon, sLat, sLon)

              GeoLatMC(i,j) = gLat*pi/180.0
              GeoLonMC(i,j) = gLon*pi/180.0

           else

              aLat = MagLatMC(i,j)
              aLon = MagLonMC(i,j)

              ReferenceAlt = Altitude_GB(1,1,0,1)/1000.0
              AltMinIono = ReferenceAlt

              call dipole_to_geo(aLat, aLon, ReferenceAlt, gLat, gLon)
              GeoLatMC(i,j) = gLat*pi/180.0
              GeoLonMC(i,j) = gLon*pi/180.0

           endif

        enddo
     enddo

     if (iDebugLevel > 3) &
          write(*,*) "====> Starting Messagepass for apex_to_geo"

     bs = nMagLats * (nMagLons+1)

     MagBufferMC = GeoLatMC
     call MPI_AllREDUCE(MagBufferMC, GeoLatMC,  &
          bs, MPI_REAL, MPI_MAX, iCommGITM, iError)
     MagBufferMC = GeoLonMC
     call MPI_AllREDUCE(MagBufferMC, GeoLonMC,  &
          bs, MPI_REAL, MPI_MAX, iCommGITM, iError)
     MagBufferMC = MagLatMC
     call MPI_AllREDUCE(MagBufferMC, MagLatMC,  &
          bs, MPI_REAL, MPI_MAX, iCommGITM, iError)
     MagBufferMC = MagLonMC
     call MPI_AllREDUCE(MagBufferMC, MagLonMC,  &
          bs, MPI_REAL, MPI_MAX, iCommGITM, iError)

     do i=1,nMagLons+1
        do j=2,nMagLats-1
           deltalmc(i,j) = MagLatMC(i,j+1) - MagLatMC(i,j-1) 
        enddo
        deltalmc(i,1) = 2*(MagLatMC(i,2) - MagLatMC(i,1))
        deltalmc(i,nMagLats) = 2*(MagLatMC(i,nMagLats) - MagLatMC(i,nMagLats-1))
     enddo

     do j=1,nMagLats
        do i=2,nMagLons
           deltapmc(i,j) = MagLonMC(i+1,j) - MagLonMC(i-1,j) 
        enddo
        deltapmc(1,j) = 2*(MagLonMC(2,j) - MagLonMC(1,j))
        deltapmc(nMagLons+1,j) = deltapmc(1,j)
     enddo

     deltapmc = deltapmc * Pi / 180.0 / 2.0
     deltalmc = deltalmc * Pi / 180.0 / 2.0

     if (UseNewTrace) then 
        LengthFieldLine   = 0.0
        call MMT_Init
     endif
     
  endif

  if(UseApex .and. IsEarth) then
     do i=1,nMagLons+1
        do j=1,nMagLats

           call magloctm(MagLonMC(i,j),SubsolarLatitude,   &
                SubsolarLongitude,  &
                MagneticPoleColat, &
                MagneticPoleLon,mltMC)

           MagLocTimeMC(i,j) = mod(mltMC+24.0,24.0)

        enddo
     enddo
  end if

end subroutine init_electrodynamics

subroutine init_electrodynamics_1d(iBlock)

  use ModGITM
  use ModInputs
  use ModConstants
  use ModElectrodynamics
  use ModLinearSolver

  implicit none
  integer, intent(in ) :: iBlock

  integer :: k

  real :: q2

  real, dimension(-1:nLons+2, -1:nLats+2, -1:nAlts+2) :: &
       e_density, Vi, Ve, MeVen, MeVei, MiVin, VeOe, ViOi

  call report("init_electrodynamics_1d",1)
  call start_timing("init_electrodyn_1d")

  q2 = Element_Charge * Element_Charge

  e_density = IDensityS(:,:,:,ie_,iBlock)  

  Vi = Collisions(:,:,:,iVIN_)
  Ve = Collisions(:,:,:,iVEN_) ! + Collisions(:,:,:,iVEI_)

  MeVen = Mass_Electron * Collisions(:,:,:,iVEN_)
  MeVei = Mass_Electron * Collisions(:,:,:,iVEI_)
  MiVin = MeanIonMass * Collisions(:,:,:,iVIN_)

  VeOe = Ve**2 + e_gyro**2
  ViOi = Vi**2 + i_gyro**2

  Sigma_0 = q2 * E_Density / (1.0/MeVen + 1.0/MiVin)

  Sigma_Pedersen = ((1.0/MeVen) * (Ve*Ve/VeOe) + &
       (1.0/MiVin) * (Vi*Vi/ViOi)) * E_Density * q2

  Sigma_Hall = ((1.0/MeVen) * (Ve*e_gyro/VeOe) - &
       (1.0/MiVin) * (Vi*i_gyro/ViOi)) * E_Density * q2

  PedersenConductance(:,:,iBlock) = 0.0
  HallConductance(:,:,iBlock)     = 0.0

  do k=1,nAlts
     PedersenConductance(:,:,iBlock) = PedersenConductance(:,:,iBlock) + &
          Sigma_Pedersen(:,:,k)*dAlt_GB(:,:,k,iBlock)
     HallConductance(:,:,iBlock)     = HallConductance(:,:,iBlock)     + &
          Sigma_Hall(:,:,k)    *dAlt_GB(:,:,k,iBlock)
  enddo

  call end_timing("init_electrodyn_1d")

end subroutine init_electrodynamics_1d
