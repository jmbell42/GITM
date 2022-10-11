! Copyright 2021, the GITM Development Team (see srcDoc/dev_team.md for members)
! Full license can be found in LICENSE

subroutine initialize_gitm(TimeIn)

  use ModGITM
  use ModInputs
  use ModConstants
  use ModPlanet
  use ModRates
  use ModSphereInterface
  use ModTime
  use ModEUV
  use ModIndicesInterfaces
  implicit none

  type (UAM_ITER) :: r_iter

  real(Real8_), intent(in) :: TimeIn
  integer :: iLat, iAlt, iBlock, iSpecies, iLon, iError,nAltsMean,iilon,iilat

  real :: TempAve
  real :: TempDiff
  real :: InvScaleHeightS(-1:nLons+2,-1:nLats+2)

  real :: LogRho(-1:nLons+2,-1:nLats+2), NewSumRho(-1:nLons+2,-1:nLats+2)
  real :: GradAlt_CD(nLons, nLats, nAlts, 3)

  real :: TempUnit_const, t, h,meanaltzero,altdiff, dAlt

  logical :: IsThere, IsOk, IsDone, IsFirstTime = .true.

  real :: DistM, DistP, Ratio2, InvDenom
  !----------------------------------------------------------------------------

!! Sorry but this is a double-negative
!! This checks to see if the planet is not, not Titan.
!! That is to say, it checks to see if this is actually Titan.

  if (   .not. (index(cPlanet,"Titan") == 0)  ) then 
     call init_radcooling
  endif

  if (.not.IsFirstTime) return

  IsFirstTime = .false.

  call start_timing("initialize")
  call report("initialize",1)

  ! ------------------------------------------------------------------------
  ! initialize stuff
  ! ------------------------------------------------------------------------

  if (TimeIn /= CurrentTime) then
     CurrentTime = TimeIn
     call time_real_to_int(CurrentTime, iTimeArray)
     call fix_vernal_time
  endif

  call get_f107(CurrentTime, f107, iError)
  call get_f107a(CurrentTime, f107a, iError)

  call init_grid

  if (DoRestart) then
     call read_inputs(trim(restartInDir)//"/header.rst")
     call set_inputs
     call read_restart(restartInDir)
     call init_msis
!     if (UsePerturbation) call user_create_perturbation
  endif

  call set_RrTempInd

  inquire(file='GITM.STOP',EXIST=IsThere)
  if (IsThere .and. iProc == 0) then
     open(iOutputUnit_, file = 'GITM.STOP', status = 'OLD')
     close(iOutputUnit_, status = 'DELETE')
  endif

  !\
  ! Initialize the EUV wave spectrum
  !/

  call init_euv

  Gravity_GB = 0.0

  if (.not. DoRestart) then

     if (UseStretchedAltitude) then
        call init_altitude
     else
        if (UseTopography) then
           if (AltMin > 1.0) then 
              write(*,*) 'When using topography, the minimum altitude'
              write(*,*) 'must be zero.  Stopping...'
              call stop_gitm('Incorrect minimum altitude')
           endif
           altzero = 0.0
           call init_topography

           meanAltZero = sum(altzero)/(nlons*nlats*nblocks)
           altdiff = AltMinUniform - meanAltZero
           naltsmean = floor(nalts - nalts*((AltMax-AltMinUniform)/(AltMax)))

           do iBlock = 1, nBlocks
              do iLon = -1,nLons +2
                 do iLat = -1, nLats +2
                     iiLon = min(max(iLon,1),nLons)
                     iilat = min(max(iLat,1),nLats)

                    do iAlt = -1,nAltsmean
                       
                       dAlt = (AltMinUniform-AltZero(iiLon,iiLat,iBlock))/nAltsmean
                       Altitude_GB(iLon,iLat,iAlt,iBlock) = &
                            AltZero(iiLon,iiLat,iBlock) + ialt*dAlt
                    enddo
!                  if (ilon ==1) then
!                     write(93,*) ilat,altzero(iilon,iilat,iblock)
!                     write(93,*)ilat, iblock,latitude(ilat,iblock)*180/pi,altitude_GB(1,ilat,0,iblock)
!                  endif
                                                         
                 enddo
              enddo

              dAlt = (AltMax-AltMinUniform)/(nAlts-naltsmean)

              do iAlt = 1, (nAlts+2)-naltsmean
                 Altitude_GB(:,:,iAlt+naltsmean,1:nBlocks) = &
                      AltMinUniform + iAlt * dAlt
              enddo
           enddo

        else
           ! Uniform grid
           do iAlt=-1,nAlts+2
              Altitude_GB(:,:,iAlt,1:nBlocks) = &
                   AltMin + (iAlt-0.5)*(AltMax-AltMin)/nAlts
           enddo
        endif
        
     end if
  endif

  ! Calculate vertical cell sizes
  do iAlt = 0,nAlts+1
     ! Cell interface is taken to be half way between cell centers
     ! so the cell size is half of the cell center distance
     ! between cells i-1 and i+1: 
     dAlt_GB(:,:,iAlt,1:nBlocks) = 0.5* &
          ( Altitude_GB(:,:,iAlt+1,1:nBlocks) &
          - Altitude_GB(:,:,iAlt-1,1:nBlocks))
  enddo
  dAlt_GB(:,:,-1,1:nBlocks)      = dAlt_GB(:,:,0,1:nBlocks)
  dAlt_GB(:,:,nAlts+2,1:nBlocks) = dAlt_GB(:,:,nAlts+1,1:nBlocks)

  RadialDistance_GB(:,:,:,1:nBlocks)    = rBody + Altitude_GB(:,:,:,1:nBlocks)
  InvRadialDistance_GB(:,:,:,1:nBlocks) = 1/RadialDistance_GB(:,:,:,1:nBlocks)

  Gravity_GB(:,:,:,1:nBlocks) = -Gravitational_Constant &
       *(rBody*InvRadialDistance_GB(:,:,:,1:nBlocks)) ** 2

  if (UseStretchedAltitude) then
     do iAlt=1,nAlts
        Gravity_GB(:,:,iAlt,:) = -Gravitational_Constant &
             *(rBody/RadialDistance_GB(:,:,iAlt,:))**2 
     enddo
  endif

  if (iDebugLevel > 2) then
     do iAlt=-1,nAlts+2
        write(*,*) "===>Altitude : ", &
             iAlt, Altitude_GB(1,1,iAlt,1), RadialDistance_GB(1,1,iAlt,1), &
             Gravity_GB(1,1,iAlt,1)
     end do
  endif

  if (Is1D) then
     Latitude(0,1)  = Latitude(1,1) - 1.0 * pi/180.0
     Latitude(-1,1) = Latitude(0,1) - 1.0 * pi/180.0
     Latitude(2,1)  = Latitude(1,1) + 1.0 * pi/180.0
     Latitude(3,1)  = Latitude(2,1) + 1.0 * pi/180.0

     Longitude(0,1)  = Longitude(1,1) - 1.0 * pi/180.0
     Longitude(-1,1) = Longitude(0,1) - 1.0 * pi/180.0
     Longitude(2,1)  = Longitude(1,1) + 1.0 * pi/180.0
     Longitude(3,1)  = Longitude(2,1) + 1.0 * pi/180.0
  endif

  ! Precalculate the limited tangent and unlimited cosine of the latitude
  TanLatitude(:,1:nBlocks) = min(abs(tan(Latitude(:,1:nBlocks))),100.0) * &
       sign(1.0,Latitude(:,1:nBlocks))

  CosLatitude(:,1:nBlocks) = cos(Latitude(:,1:nBlocks))

  ! This is done so we don't get a /0 below.

  dLatDist_GB = 1.0
  dLatDist_FB = 1.0
  dLonDist_GB = 1.0
  dLonDist_FB = 1.0

  ! Precalculate the physical size of cells in the Lat and Lon directions
  do iBlock = 1, nBlocks
     do iAlt = -1, nAlts+2
        do iLat = 0, nLats+1
           do iLon = 0, nLons+1
              ! This is the cell size assuming that cell interface is half way
              dLatDist_GB(iLon, iLat, iAlt, iBlock) = 0.5 * &
                (Latitude(iLat+1, iBlock) - Latitude(iLat-1, iBlock)) * &
                RadialDistance_GB(iLon, iLat, iAlt, iBlock)

              ! This is the distance between neighboring cells
              ! Note that face(i) is between cells i and i-1 (like in BATSRUS)
              dLatDist_FB(iLon, iLat, iAlt, iBlock) = &
                   (Latitude(iLat, iBlock) - Latitude(iLat-1, iBlock)) * &
                   0.5*(RadialDistance_GB(iLon, iLat  , iAlt, iBlock) &
                   +    RadialDistance_GB(iLon, iLat-1, iAlt, iBlock))
              
              ! This is the cell size assuming that cell interface is half way
              dLonDist_GB(iLon, iLat, iAlt, iBlock) = 0.5 * &
                   (Longitude(iLon+1,iBlock) - Longitude(iLon-1,iBlock)) * &
                   RadialDistance_GB(iLon, iLat, iAlt, iBlock)* &
                   max(abs(CosLatitude(iLat,iBlock)),0.01)

              ! This is the distance between neighboring cells
              dLonDist_FB(iLon, iLat, iAlt, iBlock) = &
                   (Longitude(iLon,iBlock) - Longitude(iLon-1,iBlock)) * &
                   0.5*(RadialDistance_GB(iLon,   iLat, iAlt, iBlock) &
                   +    RadialDistance_GB(iLon-1, iLat, iAlt, iBlock)) &
                   *max(abs(CosLatitude(iLat,iBlock)),0.01)

              CellVolume(iLon,iLat,iAlt,iBlock) = &
                   dLonDist_GB(iLon, iLat, iAlt, iBlock) * &
                   dLatDist_GB(iLon, iLat, iAlt, iBlock) * &
                   dAlt_GB(iLon,iLat,iAlt,iBlock) 

           enddo

           ! Fill in longitude ghost cells
           dLonDist_FB(-1, iLat, iAlt, iBlock) = &
                dLonDist_FB(0, iLat, iAlt, iBlock)
           dLonDist_FB(nLons+2, iLat, iAlt, iBlock) = &
                dLonDist_FB(nLons+1, iLat, iAlt, iBlock)
           dLatDist_FB(-1, iLat, iAlt, iBlock) = &
                dLatDist_FB(0, iLat, iAlt, iBlock)
           dLatDist_FB(nLons+2, iLat, iAlt, iBlock) = &
                dLatDist_FB(nLons+1, iLat, iAlt, iBlock)
        enddo
        ! Fill in latitude ghost cells
        dLonDist_FB(:, -1, iAlt, iBlock) = &
             dLonDist_FB(:, 0, iAlt, iBlock)
        dLonDist_FB(:, nLats+2, iAlt, iBlock) = &
             dLonDist_FB(:, nLats+1, iAlt, iBlock)
        dLatDist_FB(:, -1, iAlt, iBlock) = &
             dLatDist_FB(:, 0, iAlt, iBlock)
        dLatDist_FB(:, nLats+2, iAlt, iBlock) = &
             dLatDist_FB(:, nLats+1, iAlt, iBlock)
     enddo
  enddo

  InvDLatDist_GB = 1.0/dLatDist_GB
  InvDLatDist_FB = 1.0/dLatDist_FB
  InvDLonDist_GB = 1.0/dLonDist_GB
  InvDLonDist_FB = 1.0/dLonDist_FB

  ! Precalculate the coefficients for the gradient calculation
  do iBlock = 1, nBlocks
     do iLon = 1, nLons
        DistM    = Longitude(iLon,iBlock) - Longitude(iLon-1,iBlock)
        DistP    = Longitude(iLon+1,iBlock) - Longitude(iLon,iBlock)
        Ratio2   = (DistM / DistP)**2
        InvDenom = 1.0/(Ratio2*DistP + DistM)

        GradLonP_CB(iLon, iBlock) =  InvDenom*Ratio2
        GradLon0_CB(iLon, iBlock) =  InvDenom*(1-Ratio2)
        GradLonM_CB(iLon, iBlock) = -InvDenom
     enddo

     do iLat = 1, nLats
        DistM    = Latitude(iLat,iBlock) - Latitude(iLat-1,iBlock)
        DistP    = Latitude(iLat+1,iBlock) - Latitude(iLat,iBlock)
        Ratio2   = (DistM / DistP)**2
        InvDenom = 1.0/(Ratio2*DistP + DistM)

        GradLatP_CB(iLat, iBlock) =  InvDenom*Ratio2
        GradLat0_CB(iLat, iBlock) =  InvDenom*(1-Ratio2)
        GradLatM_CB(iLat, iBlock) = -InvDenom
     enddo
  end do

  if(UseTopography)then
     !!! What about the maxi tricks ??? Why is that there???
     do iBlock = 1, nBlocks
        call UAM_gradient(Altitude_GB, GradAlt_CD, iBlock)
        dAltDLon_CB(:,:,:,iBlock) = GradAlt_CD(:,:,:,iEast_)
        dAltDLat_CB(:,:,:,iBlock) = GradAlt_CD(:,:,:,iNorth_)
     end do
  else
     dAltDLon_CB = 0.
     dAltDLat_CB = 0.
  end if

  call init_heating_efficiency

!! Some Titan-Specific Startup Routines here (Regardless of Restart or Not)

  if (   .not. (index(cPlanet,"Titan") == 0)  ) then 
     call init_magheat
     call init_isochem
     call init_aerosol
  endif

  if (   .not. (index(cPlanet,"Mars") == 0)  ) then 
     call init_isochem
  endif

  if (.not. DoRestart) then

     Potential = 0.0
     Velocity = 0.0
     IVelocity = 0.0
     VerticalVelocity = 0.0

     if (UseMsis) then
        call init_msis
     else

        TempUnit_const = 1. * Mass(1) / Boltzmanns_Constant
        TempAve  = (TempMax+TempMin)/2/TempUnit_const
        TempDiff = (TempMax-TempMin)/2/TempUnit_const

        do iBlock = 1, nBlocks

           do iAlt=-1,nAlts+2
              call get_temperature(0.0, 0.0, Altitude_GB(:,:,iAlt,iBlock), t, h)
              Temperature(:,:,iAlt,iBlock)  = t/TempUnit_const
              eTemperature(:,:,iAlt,iBlock) = t
              iTemperature(:,:,iAlt,iBlock) = t
           enddo
           
           do iAlt=-1,nAlts+2

              Rho(:,:,iAlt,iBlock) = 0.0

              NewSumRho = 0.0

              do iSpecies = 1, nSpecies

                 InvScaleHeightS = -Gravity_GB(:,:,iAlt,iBlock) * &
                      Mass(iSpecies) / &
                      (Temperature(:,:,iAlt,iBlock)*TempUnit_const* &
                      Boltzmanns_Constant)

                 if(iAlt < 2)then
                    LogNS(:,:,iAlt,iSpecies,iBlock) = &
                         - (Altitude_GB(:,:,iAlt,iBlock)-AltMin) &
                         * InvScaleHeightS + LogNS0(iSpecies)
                 elseif(iAlt > 1 .and. iAlt < nAlts+2)then
                    LogNS(:,:,iAlt,iSpecies,iBlock) = &
                         LogNS(:,:,iAlt-1,iSpecies,iBlock) &
                         -(Altitude_GB(:,:,iAlt  ,iBlock) &
                         - Altitude_GB(:,:,iAlt-1,iBlock))*InvScaleHeightS &
                         -(Temperature(:,:,iAlt+1,iBlock)   &
                         - Temperature(:,:,iAlt-1,iBlock))  &
                         /(2.0*Temperature(:,:,iAlt,iBlock))
                 else
                    LogNS(:,:,iAlt,iSpecies,iBlock) = &
                         LogNS(:,:,iAlt-1,iSpecies,iBlock) &
                         -(Altitude_GB(:,:,iAlt  ,iBlock) &
                         - Altitude_GB(:,:,iAlt-1,iBlock))*InvScaleHeightS &
                         -(Temperature(:,:,iAlt  ,iBlock)  &
                         - Temperature(:,:,iAlt-1,iBlock)) &
                         /(Temperature(:,:,iAlt,iBlock))
                 endif
              
                 NewSumRho      = NewSumRho + &
                      Mass(iSpecies)*exp(LogNS(:,:,iAlt,iSpecies,iBlock))
              
              enddo
              
              do iSpecies=1,nSpecies

                 NDensityS(:,:,iAlt,iSpecies,iBlock) = &
                      exp(LogNS(:,:,iAlt,iSpecies,iBlock))

                 NDensity(:,:,iAlt,iBlock) = NDensity(:,:,iAlt,iBlock) + &
                      NDensityS(:,:,iAlt,iSpecies,iBlock)

                 Rho(:,:,iAlt,iBlock) = Rho(:,:,iAlt,iBlock) + &
                      Mass(iSpecies) * NDensityS(:,:,iAlt,iSpecies,iBlock)

              enddo

           enddo
           
        enddo

     endif

     if (UseIRI .and. IsEarth) then
        call init_iri
     else
        if (IsEarth) then
           do iBlock = 1, nBlocks
              IDensityS(:,:,:,:,iBlock)    = 1.00e8
              IDensityS(:,:,:,ie_,iBlock)  = 1.00e8*(nIons-1)
              eTemperature(:,:,:,iBlock) = 500.0
              iTemperature(:,:,:,iBlock) = 500.0
           enddo
        endif
     endif

  endif

  if (UseWACCMTides) then
     call read_waccm_tides
     call update_waccm_tides
  endif

  if (UseGSWMTides) then
     call read_tides
     call update_tides
  endif

  call init_b0
  if (IsEarth) call init_energy_deposition

  if (UseApex .and. IsEarth) then
     call report("subsolr",2)
     call SUBSOLR(iTimeArray(1),iJulianDay,iTimeArray(4),&
          iTimeArray(5),iTimeArray(6),SubsolarLatitude, &
          SubsolarLongitude)
  endif

  if (.not.Is1D) call exchange_messages_sphere

  call calc_pressure

  ! The iLon and iLat are dummy variables...
  call UA_calc_electrodynamics(iLon,iLat)
  do iBlock = 1, nBlocks
     call init_mesh_lats(iBlock)
     call init_mesh_lons(iBlock)
     call init_mesh_alts(iBlock)
  enddo 

  do iBlock = 1, nBlocks
    call calc_eddy_diffusion_coefficient(iBlock)
    call calc_rates(iBlock)
    call calc_viscosity(iBlock)
  enddo

  call end_timing("initialize")

end subroutine initialize_gitm



subroutine init_mesh_lats(iBlock)

   use ModGITM
   use ModMesh
   use ModSizeGITM, only: nLats, nLons, nAlts, nBlocks
   implicit none
   integer, intent(in) :: iBlock

   integer :: iLat, iAlt, iLon
 
   real :: CL_P120, CL_P121, CL_P122
   real :: CR_M120, CR_M121, CR_M122
   real :: X_P12, X_P32, X_P52
   real :: X_M12, X_M32, X_M52
   real :: IS0, IS1, IS2
   real :: IS0Z, IS1Z, IS2Z

   real :: LocaldLat(-2:nLats+3)  ! Need this for extrapolation


   do iAlt = 1, nAlts
     do iLon = 1, nLons
        LocaldLat(-1:nLats+2) = dLatDist_FB(iLon,-1:nLats+2,iAlt,iBlock)
        LocaldLat(-2)         = dLatDist_FB(iLon,-1,iAlt,iBlock)
        LocaldLat(nLats+3)         = dLatDist_FB(iLon,nLats+2,iAlt,iBlock)

        !                                                   X+3/2 = X+1/2 
        !                                       X+1/2 = 0.5dLonFB(i+1)
        !                                0.0 ---|
        !   o-----X-----o-----X-----o-----X-----o-----X-----o-----X-----o-----X-----o
        !   |    i-2    |    i-1    |     i     |    i+1    |    i+2    |    i+3    |
        ! i-5/2   |   i-3/2   |   i-1/2   |   i+1/2   |   i+3/2   |   i+5/2   |
        !dFB(i-2) |  dFB(i-1) |   dFB(i)  |  dFB(i+1) |  dFB(i+2) |  dFB(i+3) |
        !                     |--     2 x dGB(i)    --|
        do iLat=0,nLats+1  
           X_P12 = 0.5*LocaldLat(iLat+1)
           X_P32 =     LocaldLat(iLat+1) + 0.5*LocaldLat(iLat+2)
           if(iLat .eq. nLats + 1) then
              X_P52 = 1.0*LocaldLat(iLat+1) + LocaldLat(iLat+2) + &
                      0.5*LocaldLat(iLat+2)
           else
              X_P52 = 1.0*LocaldLat(iLat+1) + LocaldLat(iLat+2) + &
                      0.5*LocaldLat(iLat+3)
           endif

           X_M12 = -0.5*LocaldLat(iLat)
           X_M32 = -1.0*LocaldLat(iLat) - 0.5*LocaldLat(iLat-1)
           X_M52 = -1.0*LocaldLat(iLat  ) - LocaldLat(iLat-1) - &
                    0.5*LocaldLat(iLat-2)

           !------------------
           ! Uniform grid (X+1/2) =  0.5*h
           ! Uniform grid (X+3/2) =  1.5*h
           ! Uniform grid (X+5/2) =  2.5*h
           !------------------
           ! Uniform grid (X-1/2) = -0.5*h
           ! Uniform grid (X-3/2) = -1.5*h
           ! Uniform grid (X-5/2) = -2.5*h
           !------------------
      LatMesh_ULP120(iLon,iLat,iAlt,1,iBlock) = 1.0
      LatMesh_ULP120(iLon,iLat,iAlt,2,iBlock) = &
        ( (X_P32 - X_P12)/(X_P52 - X_M12))*((X_P52 - X_P12)/(X_P32 - X_M12))
      LatMesh_ULP120(iLon,iLat,iAlt,3,iBlock) = &
        ( (X_P32 - X_P12)/(X_P52 - X_M12))*((X_P12 - X_M12)/(X_P52 - X_P12))

      LatMesh_ULP121(iLon,iLat,iAlt,1,iBlock) = 1.0
      LatMesh_ULP121(iLon,iLat,iAlt,2,iBlock) = &
        ( (X_P12 - X_M12)/(X_P32 - X_M32))*((X_P12 - X_M32)/(X_P32 - X_M12))
      LatMesh_ULP121(iLon,iLat,iAlt,3,iBlock) = &
        ( (X_P12 - X_M12)/(X_P32 - X_M32))*((X_P32 - X_P12)/(X_P12 - X_M32))
!          
      LatMesh_ULP122(iLon,iLat,iAlt,1,iBlock) = 1.0
      LatMesh_ULP122(iLon,iLat,iAlt,2,iBlock) = &
        ( (X_P12 - X_M12)/(X_M12 - X_M52))*((X_P12 - X_M32)/(X_P12 - X_M52))
      LatMesh_ULP122(iLon,iLat,iAlt,3,iBlock) = &
        (1.0 +  (X_P12 - X_M12)/(X_P12 - X_M32) +  (X_P12 - X_M12)/(X_P12 - X_M52))


      LatMesh_URM120(iLon,iLat,iAlt,1,iBlock) = 1.0
      LatMesh_URM120(iLon,iLat,iAlt,2,iBlock) = &
        (1.0 +  (X_P12 - X_M12)/(X_P32 - X_M12) +  (X_P12 - X_M12)/(X_P52 - X_M12))
      LatMesh_URM120(iLon,iLat,iAlt,3,iBlock) = &
        ( (X_P12 - X_M12)/(X_P52 - X_M12))*((X_P32 - X_M12)/(X_P52 - X_P12))

      LatMesh_URM121(iLon,iLat,iAlt,1,iBlock) = 1.0
      LatMesh_URM121(iLon,iLat,iAlt,2,iBlock) = &
        ( (X_P12 - X_M12)/(X_P32 - X_M32))*((X_P32 - X_M12)/(X_P12 - X_M32))
      LatMesh_URM121(iLon,iLat,iAlt,3,iBlock) = &
        ( (X_P12 - X_M12)/(X_P32 - X_M32))*((X_M12 - X_M32)/(X_P32 - X_M12))

      LatMesh_URM122(iLon,iLat,iAlt,1,iBlock) = 1.0
      LatMesh_URM122(iLon,iLat,iAlt,2,iBlock) = &
        ( (X_M12 - X_M32)/(X_P12 - X_M52))*((X_M12 - X_M52)/(X_P12 - X_M32))
      LatMesh_URM122(iLon,iLat,iAlt,3,iBlock) = &
        ( (X_M12 - X_M32)/(X_P12 - X_M52))*((X_P12 - X_M12)/(X_M12 - X_M52))

      LatMesh_CLP120(iLon,iLat,iAlt,iBlock) = &
          (X_P12 - X_M52)/(X_P52 - X_M52)*&
          (X_P12 - X_M32)/(X_P52 - X_M32)

      LatMesh_CLP121(iLon,iLat,iAlt,iBlock) = &
          (X_P12 - X_M52)/(X_P52 - X_M52)*&
          (X_P52 - X_P12)/(X_P52 - X_M32)*&
        ( (X_P52 - X_M32)/(X_P32 - X_M52) + 1.0)

      LatMesh_CLP122(iLon,iLat,iAlt,iBlock) = &
          (X_P32 - X_P12)/(X_P52 - X_M52)*&
          (X_P52 - X_P12)/(X_P32 - X_M52)

  !------
      LatMesh_CRM120(iLon,iLat,iAlt,iBlock) = &
          (X_M12 - X_M32)/(X_P52 - X_M52)*&
          (X_M12 - X_M52)/(X_P52 - X_M32)

      LatMesh_CRM121(iLon,iLat,iAlt,iBlock) = &
          (X_M12 - X_M52)/(X_P52 - X_M52)*&
          (X_P52 - X_M12)/(X_P52 - X_M32)*&
        ( (X_P52 - X_M32)/(X_P32 - X_M52) + 1.0)

      LatMesh_CRM122(iLon,iLat,iAlt,iBlock) = &
          (X_P32 - X_M12)/(X_P52 - X_M52)*&
          (X_P52 - X_M12)/(X_P32 - X_M52)

!  real :: LatMesh_IS0(1:nLons,0:nLats+1,1:nAlts,1:3,nBlocksMax)=0.0
!  real :: LatMesh_IS1(1:nLons,0:nLats+1,1:nAlts,1:3,nBlocksMax)=0.0
!  real :: LatMesh_IS2(1:nLons,0:nLats+1,1:nAlts,1:3,nBlocksMax)=0.0

      LatMesh_IS0(iLon,iLat,iAlt,1,iBlock) = &
         4.0*( ( (X_P12 - X_M12)/(X_P52 - X_M12))**2.0)*&
               (( 1.0/(X_P52 - X_P12))**2.0)* &
               (10.0*(X_P12 - X_M12)**2.0 + (X_P32 - X_M12)*(X_P32 - X_P12)) 
      LatMesh_IS0(iLon,iLat,iAlt,2,iBlock) = &
         4.0*( ( (X_P12 - X_M12)/(X_P52 - X_M12))**2.0)*&
               ( 1.0/ ( (X_P52 - X_P12)*(X_P32 - X_M12)))*&
               (20.0*(X_P12 - X_M12)**2.0 + 2.0*(X_P32 - X_M12)*(X_P32 - X_P12) + &
                    (X_P52 - X_M12)*(2.0*X_P32 - X_P12 - X_M12) ) 
      LatMesh_IS0(iLon,iLat,iAlt,3,iBlock) = &
         4.0*( ( (X_P12 - X_M12)/(X_P52 - X_M12))**2.0)*&
               ( (1.0/(X_P32 - X_M12))**2.0)* &
               (10.0*(X_P12 - X_M12)**2.0 + (X_P52 + X_P32 - 2.0*X_M12)* &
                    (X_P52 + X_P32 - X_P12 - X_M12) ) 

      LatMesh_IS1(iLon,iLat,iAlt,1,iBlock) = &
         4.0*( ( (X_P12 - X_M12)/(X_P32 - X_M32))**2.0)*&
               ( (1.0/(X_P12 - X_M32))**2.0)* &
               (10.0*(X_P12 - X_M12)**2.0 + (X_P32 - X_M12)*(X_P32 - X_P12)) 
      LatMesh_IS1(iLon,iLat,iAlt,2,iBlock) = &
         4.0*( ( (X_P12 - X_M12)/(X_P32 - X_M32))**2.0)*&
               ( 1.0/( (X_P32 - X_M12)*(X_P12 - X_M32)))*&
               (20.0*(X_P12 - X_M12)**2.0 - (X_P32 - X_P12)*(X_M12 - X_M32) - &
                    (X_P32 - X_M12)*(X_P12 - X_M32) ) 
      LatMesh_IS1(iLon,iLat,iAlt,3,iBlock) = &
         4.0*( ( (X_P12 - X_M12)/(X_P32 - X_M32))**2.0)*&
               ( (1.0/(X_P32 - X_M12))**2.0)* &
               (10.0*(X_P12 - X_M12)**2.0 + (X_M12 - X_M32)*(X_P12 - X_M32) ) 


      LatMesh_IS2(iLon,iLat,iAlt,1,iBlock) = &
          4.0*( ( (X_P12 - X_M12)/(X_P12 - X_M52))**2.0)*&
               ( (1.0/(X_M12 - X_M52))**2.0)* &
               (10.0*(X_P12 - X_M12)**2.0 + (X_P12 - X_M32)*(X_M12 - X_M32)) 
      LatMesh_IS2(iLon,iLat,iAlt,2,iBlock) = &
          4.0*( ( (X_P12 - X_M12)/(X_P12 - X_M52))**2.0)*&
               (1.0/((X_P12 - X_M32)*(X_M12 - X_M52)))*&
               (20.0*(X_P12 - X_M12)**2.0 + 2.0*(X_P12 - X_M32)*(X_M12 - X_M32) + &
                    (X_P12 - X_M52)*(X_P12 + X_M12 - 2.0*X_M32) ) 
      LatMesh_IS2(iLon,iLat,iAlt,3,iBlock) = &
          4.0*( ( (X_P12 - X_M12)/(X_P12 - X_M52))**2.0)*&
               ( (1.0/(X_P12 - X_M32))**2.0)* &
               (10.0*(X_P12 - X_M12)**2.0 + (2.0*X_P12 - X_M52 - X_M32)* &
                    (X_P12 + X_M12 - X_M32 - X_M52) ) 

        enddo !iLat=0,nLats+1  
     enddo !iLon = 1, nLons
   enddo !iAlt = 1, nAlts


end subroutine !init_mesh_lats

subroutine init_mesh_lons(iBlock)

   use ModGITM
   use ModMesh
   use ModSizeGITM, only: nLats, nLons, nAlts, nBlocks
   implicit none
   integer, intent(in) :: iBlock

   integer :: iLat, iAlt, iLon
 
   real :: CL_P120, CL_P121, CL_P122
   real :: CR_M120, CR_M121, CR_M122
   real :: X_P12, X_P32, X_P52
   real :: X_M12, X_M32, X_M52
   real :: IS0, IS1, IS2
   real :: IS0Z, IS1Z, IS2Z

   real :: LocaldLon(-2:nLons+3)  ! Need this for extrapolation


   do iAlt = 1, nAlts
     do iLat = 1, nLats
        !LocaldLon(-2:nLons+3) = dLonDist_FB(-2:nLons+3,iLat,iAlt,iBlock)
        LocaldLon(-1:nLons+2) = dLonDist_FB(-1:nLons+2,iLat,iAlt,iBlock)
        LocaldLon(-2) = dLonDist_FB(-1,iLat,iAlt,iBlock)
        LocaldLon(nLons+3) = dLonDist_FB(nLons+2,iLat,iAlt,iBlock)

        !                                                   X+3/2 = X+1/2 
        !                                       X+1/2 = 0.5dLonFB(i+1)
        !                                0.0 ---|
        !   o-----X-----o-----X-----o-----X-----o-----X-----o-----X-----o-----X-----o
        !   |    i-2    |    i-1    |     i     |    i+1    |    i+2    |    i+3    |
        ! i-5/2   |   i-3/2   |   i-1/2   |   i+1/2   |   i+3/2   |   i+5/2   |
        !dFB(i-2) |  dFB(i-1) |   dFB(i)  |  dFB(i+1) |  dFB(i+2) |  dFB(i+3) |
        !                     |--     2 x dGB(i)    --|
        do iLon=0,nLons+1  
           X_P12 = 0.5*LocaldLon(iLon+1)
           X_P32 =     LocaldLon(iLon+1) + 0.5*LocaldLon(iLon+2)
           if(iLon .eq. nLons + 1) then
              X_P52 = 1.0*LocaldLon(iLon+1) + LocaldLon(iLon+2) + &
                      0.5*LocaldLon(iLon+2)
           else
              X_P52 = 1.0*LocaldLon(iLon+1) + LocaldLon(iLon+2) + &
                      0.5*LocaldLon(iLon+3)
           endif

           X_M12 = -0.5*LocaldLon(iLon)
           X_M32 = -1.0*LocaldLon(iLon) - 0.5*LocaldLon(iLon-1)
           X_M52 = -1.0*LocaldLon(iLon  ) - LocaldLon(iLon-1) - &
                    0.5*LocaldLon(iLon-2)

           !------------------
           ! Uniform grid (X+1/2) =  0.5*h
           ! Uniform grid (X+3/2) =  1.5*h
           ! Uniform grid (X+5/2) =  2.5*h
           !------------------
           ! Uniform grid (X-1/2) = -0.5*h
           ! Uniform grid (X-3/2) = -1.5*h
           ! Uniform grid (X-5/2) = -2.5*h
           !------------------
      LonMesh_ULP120(iLon,iLat,iAlt,1,iBlock) = 1.0
      LonMesh_ULP120(iLon,iLat,iAlt,2,iBlock) = &
        ( (X_P32 - X_P12)/(X_P52 - X_M12))*((X_P52 - X_P12)/(X_P32 - X_M12))
      LonMesh_ULP120(iLon,iLat,iAlt,3,iBlock) = &
        ( (X_P32 - X_P12)/(X_P52 - X_M12))*((X_P12 - X_M12)/(X_P52 - X_P12))

      LonMesh_ULP121(iLon,iLat,iAlt,1,iBlock) = 1.0
      LonMesh_ULP121(iLon,iLat,iAlt,2,iBlock) = &
        ( (X_P12 - X_M12)/(X_P32 - X_M32))*((X_P12 - X_M32)/(X_P32 - X_M12))
      LonMesh_ULP121(iLon,iLat,iAlt,3,iBlock) = &
        ( (X_P12 - X_M12)/(X_P32 - X_M32))*((X_P32 - X_P12)/(X_P12 - X_M32))
!          
      LonMesh_ULP122(iLon,iLat,iAlt,1,iBlock) = 1.0
      LonMesh_ULP122(iLon,iLat,iAlt,2,iBlock) = &
        ( (X_P12 - X_M12)/(X_M12 - X_M52))*((X_P12 - X_M32)/(X_P12 - X_M52))
      LonMesh_ULP122(iLon,iLat,iAlt,3,iBlock) = &
        (1.0 +  (X_P12 - X_M12)/(X_P12 - X_M32) +  (X_P12 - X_M12)/(X_P12 - X_M52))


      LonMesh_URM120(iLon,iLat,iAlt,1,iBlock) = 1.0
      LonMesh_URM120(iLon,iLat,iAlt,2,iBlock) = &
        (1.0 +  (X_P12 - X_M12)/(X_P32 - X_M12) +  (X_P12 - X_M12)/(X_P52 - X_M12))
      LonMesh_URM120(iLon,iLat,iAlt,3,iBlock) = &
        ( (X_P12 - X_M12)/(X_P52 - X_M12))*((X_P32 - X_M12)/(X_P52 - X_P12))

      LonMesh_URM121(iLon,iLat,iAlt,1,iBlock) = 1.0
      LonMesh_URM121(iLon,iLat,iAlt,2,iBlock) = &
        ( (X_P12 - X_M12)/(X_P32 - X_M32))*((X_P32 - X_M12)/(X_P12 - X_M32))
      LonMesh_URM121(iLon,iLat,iAlt,3,iBlock) = &
        ( (X_P12 - X_M12)/(X_P32 - X_M32))*((X_M12 - X_M32)/(X_P32 - X_M12))

      LonMesh_URM122(iLon,iLat,iAlt,1,iBlock) = 1.0
      LonMesh_URM122(iLon,iLat,iAlt,2,iBlock) = &
        ( (X_M12 - X_M32)/(X_P12 - X_M52))*((X_M12 - X_M52)/(X_P12 - X_M32))
      LonMesh_URM122(iLon,iLat,iAlt,3,iBlock) = &
        ( (X_M12 - X_M32)/(X_P12 - X_M52))*((X_P12 - X_M12)/(X_M12 - X_M52))

      LonMesh_CLP120(iLon,iLat,iAlt,iBlock) = &
          (X_P12 - X_M52)/(X_P52 - X_M52)*&
          (X_P12 - X_M32)/(X_P52 - X_M32)

      LonMesh_CLP121(iLon,iLat,iAlt,iBlock) = &
          (X_P12 - X_M52)/(X_P52 - X_M52)*&
          (X_P52 - X_P12)/(X_P52 - X_M32)*&
        ( (X_P52 - X_M32)/(X_P32 - X_M52) + 1.0)

      LonMesh_CLP122(iLon,iLat,iAlt,iBlock) = &
          (X_P32 - X_P12)/(X_P52 - X_M52)*&
          (X_P52 - X_P12)/(X_P32 - X_M52)

  !------
      LonMesh_CRM120(iLon,iLat,iAlt,iBlock) = &
          (X_M12 - X_M32)/(X_P52 - X_M52)*&
          (X_M12 - X_M52)/(X_P52 - X_M32)

      LonMesh_CRM121(iLon,iLat,iAlt,iBlock) = &
          (X_M12 - X_M52)/(X_P52 - X_M52)*&
          (X_P52 - X_M12)/(X_P52 - X_M32)*&
        ( (X_P52 - X_M32)/(X_P32 - X_M52) + 1.0)

      LonMesh_CRM122(iLon,iLat,iAlt,iBlock) = &
          (X_P32 - X_M12)/(X_P52 - X_M52)*&
          (X_P52 - X_M12)/(X_P32 - X_M52)

!  real :: LonMesh_IS0(1:nLons,0:nLats+1,1:nAlts,1:3,nBlocksMax)=0.0
!  real :: LonMesh_IS1(1:nLons,0:nLats+1,1:nAlts,1:3,nBlocksMax)=0.0
!  real :: LonMesh_IS2(1:nLons,0:nLats+1,1:nAlts,1:3,nBlocksMax)=0.0

      LonMesh_IS0(iLon,iLat,iAlt,1,iBlock) = &
         4.0*( ( (X_P12 - X_M12)/(X_P52 - X_M12))**2.0)*&
               (( 1.0/(X_P52 - X_P12))**2.0)* &
               (10.0*(X_P12 - X_M12)**2.0 + (X_P32 - X_M12)*(X_P32 - X_P12)) 
      LonMesh_IS0(iLon,iLat,iAlt,2,iBlock) = &
         4.0*( ( (X_P12 - X_M12)/(X_P52 - X_M12))**2.0)*&
               ( 1.0/ ( (X_P52 - X_P12)*(X_P32 - X_M12)))*&
               (20.0*(X_P12 - X_M12)**2.0 + 2.0*(X_P32 - X_M12)*(X_P32 - X_P12) + &
                    (X_P52 - X_M12)*(2.0*X_P32 - X_P12 - X_M12) ) 
      LonMesh_IS0(iLon,iLat,iAlt,3,iBlock) = &
         4.0*( ( (X_P12 - X_M12)/(X_P52 - X_M12))**2.0)*&
               ( (1.0/(X_P32 - X_M12))**2.0)* &
               (10.0*(X_P12 - X_M12)**2.0 + (X_P52 + X_P32 - 2.0*X_M12)* &
                    (X_P52 + X_P32 - X_P12 - X_M12) ) 

      LonMesh_IS1(iLon,iLat,iAlt,1,iBlock) = &
         4.0*( ( (X_P12 - X_M12)/(X_P32 - X_M32))**2.0)*&
               ( (1.0/(X_P12 - X_M32))**2.0)* &
               (10.0*(X_P12 - X_M12)**2.0 + (X_P32 - X_M12)*(X_P32 - X_P12)) 
      LonMesh_IS1(iLon,iLat,iAlt,2,iBlock) = &
         4.0*( ( (X_P12 - X_M12)/(X_P32 - X_M32))**2.0)*&
               ( 1.0/( (X_P32 - X_M12)*(X_P12 - X_M32)))*&
               (20.0*(X_P12 - X_M12)**2.0 - (X_P32 - X_P12)*(X_M12 - X_M32) - &
                    (X_P32 - X_M12)*(X_P12 - X_M32) ) 
      LonMesh_IS1(iLon,iLat,iAlt,3,iBlock) = &
         4.0*( ( (X_P12 - X_M12)/(X_P32 - X_M32))**2.0)*&
               ( (1.0/(X_P32 - X_M12))**2.0)* &
               (10.0*(X_P12 - X_M12)**2.0 + (X_M12 - X_M32)*(X_P12 - X_M32) ) 


      LonMesh_IS2(iLon,iLat,iAlt,1,iBlock) = &
          4.0*( ( (X_P12 - X_M12)/(X_P12 - X_M52))**2.0)*&
               ( (1.0/(X_M12 - X_M52))**2.0)* &
               (10.0*(X_P12 - X_M12)**2.0 + (X_P12 - X_M32)*(X_M12 - X_M32)) 
      LonMesh_IS2(iLon,iLat,iAlt,2,iBlock) = &
          4.0*( ( (X_P12 - X_M12)/(X_P12 - X_M52))**2.0)*&
               (1.0/((X_P12 - X_M32)*(X_M12 - X_M52)))*&
               (20.0*(X_P12 - X_M12)**2.0 + 2.0*(X_P12 - X_M32)*(X_M12 - X_M32) + &
                    (X_P12 - X_M52)*(X_P12 + X_M12 - 2.0*X_M32) ) 
      LonMesh_IS2(iLon,iLat,iAlt,3,iBlock) = &
          4.0*( ( (X_P12 - X_M12)/(X_P12 - X_M52))**2.0)*&
               ( (1.0/(X_P12 - X_M32))**2.0)* &
               (10.0*(X_P12 - X_M12)**2.0 + (2.0*X_P12 - X_M52 - X_M32)* &
                    (X_P12 + X_M12 - X_M32 - X_M52) ) 

        enddo !iLat=0,nLats+1  
     enddo !iLon = 1, nLons
   enddo !iAlt = 1, nAlts


end subroutine !init_mesh_lats

subroutine init_mesh_alts(iBlock)

   use ModGITM
   use ModMesh
   use ModSizeGITM, only: nLats, nLons, nAlts, nBlocks
   implicit none
   integer, intent(in) :: iBlock

   integer :: iLat, iAlt, iLon
   integer :: iiAlt
 
   real :: CL_P120, CL_P121, CL_P122
   real :: CR_M120, CR_M121, CR_M122
   real :: X_P12, X_P32, X_P52
   real :: X_M12, X_M32, X_M52
   real :: IS0, IS1, IS2
   real :: IS0Z, IS1Z, IS2Z
   !----------------------------------------
   real :: LocaldAlt(-2:nAlts+3)  ! Need this for extrapolation
   real :: dZ
   real :: Alt_M52, Alt_P52
   real :: Alt_M32, Alt_P32
   !----------------------------------------
   real :: h1, h2, h3, h4
   real :: MeshH1, MeshH2, MeshH3, MeshH4
   real :: MeshCoef0, MeshCoef1, &
           MeshCoef2, MeshCoef3, &
           MeshCoef4
   real :: hm1, hm2, hm3, hm4
   real :: MeshHm1, MeshHm2, MeshHm3, MeshHm4


   do iLon = 1, nLons
     do iLat = 1, nLats
        !                                                   X+3/2 = X+1/2 
        !                                       X+1/2 = 0.5dLonFB(i+1)
        !                                0.0 ---|
        !                               Z(i)    
        !   o-----X-----o-----X-----o-----X-----o-----X-----o-----X-----o-----X-----o
        !   |    i-2    |    i-1    |     i     |    i+1    |    i+2    |    i+3    |
        ! i-5/2   |   i-3/2   |   i-1/2   |   i+1/2   |   i+3/2   |   i+5/2   |
        !dFB(i-2) |  dFB(i-1) |   dFB(i)  |  dFB(i+1) |  dFB(i+2) |  dFB(i+3) |
        !                     |--     2 x dGB(i)    --|
        do iAlt=0,nAlts+1  
           X_M12 = 0.5*(Altitude_GB(iLon,iLat,iAlt  ,iBlock) + &
                        Altitude_GB(iLon,iLat,iAlt-1,iBlock))

           dZ = Altitude_GB(iLon,iLat,0,iBlock) - Altitude_GB(iLon,iLat,-1,iBlock)
           if(iAlt .eq. 0) then
  !                                    |<------ dZ  ------>|
  !           |---- (-2)    ---|----  -1     ----|----     0   
              !   ( iAlt-2)    |     iAlt-1      |        iAlt  
  !           M52              M32               M12          
              X_M32 = Altitude_GB(iLon,iLat,-1,iBlock) - 0.5*dZ
              X_M52 = Altitude_GB(iLon,iLat,-1,iBlock) - 1.5*dZ 
           elseif(iAlt .eq. 1) then
  !                                    |<------ dZ  ------>|
  !           |----  -1     ---|----   0     ----|----     1   
              !     iAlt-2     |     iAlt-1      |        iAlt  
  !           M52              M32               M12          
              X_M32 = 0.5*(Altitude_GB(iLon,iLat,iAlt-1 ,iBlock) + &
                           Altitude_GB(iLon,iLat,iAlt-2 ,iBlock))
              X_M52 = Altitude_GB(iLon,iLat,-1,iBlock) - 0.5*dZ 
           else
              X_M32 = 0.5*(Altitude_GB(iLon,iLat,iAlt-1 ,iBlock) + &
                           Altitude_GB(iLon,iLat,iAlt-2 ,iBlock))
              X_M52 = 0.5*(Altitude_GB(iLon,iLat,iAlt-2 ,iBlock) + &
                           Altitude_GB(iLon,iLat,iAlt-3 ,iBlock))
           endif

           X_P12 = 0.5*(Altitude_GB(iLon,iLat,iAlt+1,iBlock) + &
                        Altitude_GB(iLon,iLat,iAlt  ,iBlock))
           dZ = Altitude_GB(iLon,iLat,nAlts+2,iBlock) - Altitude_GB(iLon,iLat,nAlts+1,iBlock)
           if(iAlt .eq. nAlts+1) then
  !                                   |<------dZ------>|<------------
  !           |---- nAlts+1 ---|---- nAlts+2 ----|---- (nAlts+3 ---- |
              !     iAlt       |     iAlt+1      |      iAlt+2      |
              !               P12               P32                  P52
              X_P32 = Altitude_GB(iLon,iLat,nAlts+2,iBlock) + 0.5*dZ
              X_P52 = Altitude_GB(iLon,iLat,nAlts+2,iBlock) + 0.5*dZ + 1.0*dZ
           elseif(iAlt .eq. nAlts) then
  !                                   |<------dZ------>|<------------
  !           |---- nAlts ---|---- nAlts+1 ----|---- nAlts + 2 ---- |
              !     iAlt     |      iAlt + 1   |      iAlt + 2      |
              !             P12               P32                  P52
              X_P32 = 0.5*(Altitude_GB(iLon,iLat,iAlt+1 ,iBlock) + &
                           Altitude_GB(iLon,iLat,iAlt+2 ,iBlock))
              X_P52 = Altitude_GB(iLon,iLat,nAlts+2,iBlock) + 0.5*dZ 
           else
              X_P32 = 0.5*(Altitude_GB(iLon,iLat,iAlt+1 ,iBlock) + &
                           Altitude_GB(iLon,iLat,iAlt+2 ,iBlock))
              X_P52 = 0.5*(Altitude_GB(iLon,iLat,iAlt+2 ,iBlock) + &
                           Altitude_GB(iLon,iLat,iAlt+3 ,iBlock))
           endif

           !------------------
           ! Uniform grid (X+1/2) =  0.5*h
           ! Uniform grid (X+3/2) =  1.5*h
           ! Uniform grid (X+5/2) =  2.5*h
           !------------------
           ! Uniform grid (X-1/2) = -0.5*h
           ! Uniform grid (X-3/2) = -1.5*h
           ! Uniform grid (X-5/2) = -2.5*h
           !------------------
      AltMesh_ULP120(iLon,iLat,iAlt,1,iBlock) = 1.0
      AltMesh_ULP120(iLon,iLat,iAlt,2,iBlock) = &
        ( (X_P32 - X_P12)/(X_P52 - X_M12))*((X_P52 - X_P12)/(X_P32 - X_M12))
      AltMesh_ULP120(iLon,iLat,iAlt,3,iBlock) = &
        ( (X_P32 - X_P12)/(X_P52 - X_M12))*((X_P12 - X_M12)/(X_P52 - X_P12))

      AltMesh_ULP121(iLon,iLat,iAlt,1,iBlock) = 1.0
      AltMesh_ULP121(iLon,iLat,iAlt,2,iBlock) = &
        ( (X_P12 - X_M12)/(X_P32 - X_M32))*((X_P12 - X_M32)/(X_P32 - X_M12))
      AltMesh_ULP121(iLon,iLat,iAlt,3,iBlock) = &
        ( (X_P12 - X_M12)/(X_P32 - X_M32))*((X_P32 - X_P12)/(X_P12 - X_M32))
!          
      AltMesh_ULP122(iLon,iLat,iAlt,1,iBlock) = 1.0
      AltMesh_ULP122(iLon,iLat,iAlt,2,iBlock) = &
        ( (X_P12 - X_M12)/(X_M12 - X_M52))*((X_P12 - X_M32)/(X_P12 - X_M52))
      AltMesh_ULP122(iLon,iLat,iAlt,3,iBlock) = &
        (1.0 +  (X_P12 - X_M12)/(X_P12 - X_M32) +  (X_P12 - X_M12)/(X_P12 - X_M52))


      AltMesh_URM120(iLon,iLat,iAlt,1,iBlock) = 1.0
      AltMesh_URM120(iLon,iLat,iAlt,2,iBlock) = &
        (1.0 +  (X_P12 - X_M12)/(X_P32 - X_M12) +  (X_P12 - X_M12)/(X_P52 - X_M12))
      AltMesh_URM120(iLon,iLat,iAlt,3,iBlock) = &
        ( (X_P12 - X_M12)/(X_P52 - X_M12))*((X_P32 - X_M12)/(X_P52 - X_P12))

      AltMesh_URM121(iLon,iLat,iAlt,1,iBlock) = 1.0
      AltMesh_URM121(iLon,iLat,iAlt,2,iBlock) = &
        ( (X_P12 - X_M12)/(X_P32 - X_M32))*((X_P32 - X_M12)/(X_P12 - X_M32))
      AltMesh_URM121(iLon,iLat,iAlt,3,iBlock) = &
        ( (X_P12 - X_M12)/(X_P32 - X_M32))*((X_M12 - X_M32)/(X_P32 - X_M12))

      AltMesh_URM122(iLon,iLat,iAlt,1,iBlock) = 1.0
      AltMesh_URM122(iLon,iLat,iAlt,2,iBlock) = &
        ( (X_M12 - X_M32)/(X_P12 - X_M52))*((X_M12 - X_M52)/(X_P12 - X_M32))
      AltMesh_URM122(iLon,iLat,iAlt,3,iBlock) = &
        ( (X_M12 - X_M32)/(X_P12 - X_M52))*((X_P12 - X_M12)/(X_M12 - X_M52))

      AltMesh_CLP120(iLon,iLat,iAlt,iBlock) = &
          (X_P12 - X_M52)/(X_P52 - X_M52)*&
          (X_P12 - X_M32)/(X_P52 - X_M32)

      AltMesh_CLP121(iLon,iLat,iAlt,iBlock) = &
          (X_P12 - X_M52)/(X_P52 - X_M52)*&
          (X_P52 - X_P12)/(X_P52 - X_M32)*&
        ( (X_P52 - X_M32)/(X_P32 - X_M52) + 1.0)

      AltMesh_CLP122(iLon,iLat,iAlt,iBlock) = &
          (X_P32 - X_P12)/(X_P52 - X_M52)*&
          (X_P52 - X_P12)/(X_P32 - X_M52)

  !------
      AltMesh_CRM120(iLon,iLat,iAlt,iBlock) = &
          (X_M12 - X_M32)/(X_P52 - X_M52)*&
          (X_M12 - X_M52)/(X_P52 - X_M32)

      AltMesh_CRM121(iLon,iLat,iAlt,iBlock) = &
          (X_M12 - X_M52)/(X_P52 - X_M52)*&
          (X_P52 - X_M12)/(X_P52 - X_M32)*&
        ( (X_P52 - X_M32)/(X_P32 - X_M52) + 1.0)

      AltMesh_CRM122(iLon,iLat,iAlt,iBlock) = &
          (X_P32 - X_M12)/(X_P52 - X_M52)*&
          (X_P52 - X_M12)/(X_P32 - X_M52)

!  real :: LonMesh_IS0(1:nLons,0:nLats+1,1:nAlts,1:3,nBlocksMax)=0.0
!  real :: LonMesh_IS1(1:nLons,0:nLats+1,1:nAlts,1:3,nBlocksMax)=0.0
!  real :: LonMesh_IS2(1:nLons,0:nLats+1,1:nAlts,1:3,nBlocksMax)=0.0

      AltMesh_IS0(iLon,iLat,iAlt,1,iBlock) = &
         4.0*( ( (X_P12 - X_M12)/(X_P52 - X_M12))**2.0)*&
               (( 1.0/(X_P52 - X_P12))**2.0)* &
               (10.0*(X_P12 - X_M12)**2.0 + (X_P32 - X_M12)*(X_P32 - X_P12)) 
      AltMesh_IS0(iLon,iLat,iAlt,2,iBlock) = &
         4.0*( ( (X_P12 - X_M12)/(X_P52 - X_M12))**2.0)*&
               ( 1.0/ ( (X_P52 - X_P12)*(X_P32 - X_M12)))*&
               (20.0*(X_P12 - X_M12)**2.0 + 2.0*(X_P32 - X_M12)*(X_P32 - X_P12) + &
                    (X_P52 - X_M12)*(2.0*X_P32 - X_P12 - X_M12) ) 
      AltMesh_IS0(iLon,iLat,iAlt,3,iBlock) = &
         4.0*( ( (X_P12 - X_M12)/(X_P52 - X_M12))**2.0)*&
               ( (1.0/(X_P32 - X_M12))**2.0)* &
               (10.0*(X_P12 - X_M12)**2.0 + (X_P52 + X_P32 - 2.0*X_M12)* &
                    (X_P52 + X_P32 - X_P12 - X_M12) ) 

      AltMesh_IS1(iLon,iLat,iAlt,1,iBlock) = &
         4.0*( ( (X_P12 - X_M12)/(X_P32 - X_M32))**2.0)*&
               ( (1.0/(X_P12 - X_M32))**2.0)* &
               (10.0*(X_P12 - X_M12)**2.0 + (X_P32 - X_M12)*(X_P32 - X_P12)) 
      AltMesh_IS1(iLon,iLat,iAlt,2,iBlock) = &
         4.0*( ( (X_P12 - X_M12)/(X_P32 - X_M32))**2.0)*&
               ( 1.0/( (X_P32 - X_M12)*(X_P12 - X_M32)))*&
               (20.0*(X_P12 - X_M12)**2.0 - (X_P32 - X_P12)*(X_M12 - X_M32) - &
                    (X_P32 - X_M12)*(X_P12 - X_M32) ) 
      AltMesh_IS1(iLon,iLat,iAlt,3,iBlock) = &
         4.0*( ( (X_P12 - X_M12)/(X_P32 - X_M32))**2.0)*&
               ( (1.0/(X_P32 - X_M12))**2.0)* &
               (10.0*(X_P12 - X_M12)**2.0 + (X_M12 - X_M32)*(X_P12 - X_M32) ) 


      AltMesh_IS2(iLon,iLat,iAlt,1,iBlock) = &
          4.0*( ( (X_P12 - X_M12)/(X_P12 - X_M52))**2.0)*&
               ( (1.0/(X_M12 - X_M52))**2.0)* &
               (10.0*(X_P12 - X_M12)**2.0 + (X_P12 - X_M32)*(X_M12 - X_M32)) 
      AltMesh_IS2(iLon,iLat,iAlt,2,iBlock) = &
          4.0*( ( (X_P12 - X_M12)/(X_P12 - X_M52))**2.0)*&
               (1.0/((X_P12 - X_M32)*(X_M12 - X_M52)))*&
               (20.0*(X_P12 - X_M12)**2.0 + 2.0*(X_P12 - X_M32)*(X_M12 - X_M32) + &
                    (X_P12 - X_M52)*(X_P12 + X_M12 - 2.0*X_M32) ) 
      AltMesh_IS2(iLon,iLat,iAlt,3,iBlock) = &
          4.0*( ( (X_P12 - X_M12)/(X_P12 - X_M52))**2.0)*&
               ( (1.0/(X_P12 - X_M32))**2.0)* &
               (10.0*(X_P12 - X_M12)**2.0 + (2.0*X_P12 - X_M52 - X_M32)* &
                    (X_P12 + X_M12 - X_M32 - X_M52) ) 

        enddo !iLat=0,nLats+1  
     enddo !iLon = 1, nLons
   enddo !iAlt = 1, nAlts

  ! Calculate the Central Difference and One-Sided Gradient Coefficients

  ! Use these for Central difference mesh coefficients
  !real ::     AltCDMeshCoef(1:nLons,1:nLats,1:nAlts,1:5,nBlocksMax)=0.0
  !real :: AltOSDownMeshCoef(1:nLons,1:nLats,1:nAlts,1:5,nBlocksMax)=0.0
  !real ::   AltOSUpMeshCoef(1:nLons,1:nLats,1:nAlts,1:5,nBlocksMax)=0.0
   do iLon = 1, nLons
     do iLat = 1, nLats
       do iAlt = 1, nAlts

         h1 = Altitude_GB(iLon,iLat,iAlt-1,iBlock) - Altitude_GB(iLon,iLat,iAlt-2,iBlock)
         h2 = Altitude_GB(iLon,iLat,iAlt  ,iBlock) - Altitude_GB(iLon,iLat,iAlt-1,iBlock)
         h3 = Altitude_GB(iLon,iLat,iAlt+1,iBlock) - Altitude_GB(iLon,iLat,iAlt  ,iBlock)
         h4 = Altitude_GB(iLon,iLat,iAlt+2,iBlock) - Altitude_GB(iLon,iLat,iAlt+1,iBlock)

         MeshH2 = h2 + h1
         MeshH3 = h3 + h2 + h1
         MeshH4 = h4 + h3 + h2 + h1
!
         MeshCoef0 = (h2*h3*(h3+h4))/(h1*MeshH2*MeshH3*MeshH4)
         MeshCoef1 = -1.0*(MeshH2*h3*(h3 + h4))/(h1*h2*(h2+h3)*(h2+h3+h4))
         MeshCoef3 = MeshH2*h2*(h4 + h3)/(MeshH3*(h2+h3)*h3*h4) 
         MeshCoef4 = -1.0*MeshH2*h2*h3/(MeshH4*(h2+h3+h4)*(h3+h4)*h4)
 
         MeshCoef2 = (h2*h3*(h3+h4) + &
                      MeshH2*h3*(h3+h4) - &
                      MeshH2*h2*(h3+h4) - &
                      MeshH2*h2*h3)/&
                      (MeshH2*h2*h3*(h3+h4))

         AltCDMeshCoef(iLon,iLat,iAlt,1,iBlock)=MeshCoef0
         AltCDMeshCoef(iLon,iLat,iAlt,2,iBlock)=MeshCoef1
         AltCDMeshCoef(iLon,iLat,iAlt,3,iBlock)=MeshCoef2
         AltCDMeshCoef(iLon,iLat,iAlt,4,iBlock)=MeshCoef3
         AltCDMeshCoef(iLon,iLat,iAlt,5,iBlock)=MeshCoef4



       enddo !iAlt = 1, nAlts
     enddo !iLat = 1, nLats
   enddo !iLon = 1, nLons

  ! One-Sided Derivative Coeficients

   do iLon = 1, nLons
     do iLat = 1, nLats

       ! One Sided Gradient Mesh Coefficients at the bottom of the model
       do iAlt = -2, 0
        iiAlt = iAlt + 3
        h1 = Altitude_GB(iLon,iLat,iAlt+2,iBlock) - Altitude_GB(iLon,iLat,iAlt+1,iBlock)
        h2 = Altitude_GB(iLon,iLat,iAlt+3,iBlock) - Altitude_GB(iLon,iLat,iAlt+2,iBlock)
        h3 = Altitude_GB(iLon,iLat,iAlt+4,iBlock) - Altitude_GB(iLon,iLat,iAlt+3,iBlock)
        h4 = Altitude_GB(iLon,iLat,iAlt+5,iBlock) - Altitude_GB(iLon,iLat,iAlt+4,iBlock)
    !h1 = dAlt_F(iAlt+2) ! dAlt_F(1) = Alt(1) - Alt(0);  h1 in notes  
    !h2 = dAlt_F(iAlt+3) ! dAlt_F(2) = Alt(2) - Alt(1);  h2 in notes
    !h3 = dAlt_F(iAlt+4) ! dAlt_F(3) = Alt(3) - Alt(2);  h3 in notes
    !h4 = dAlt_F(iAlt+5) ! dAlt_F(4) = Alt(4) - Alt(3);  h4 in notes

    !! Mesh Coefficients are summations over the individual mesh scales
    !MeshH1 = h1                 
    !MeshH2 = h1 + h2            
    !MeshH3 = h1 + h2 + h3
    !MeshH4 = h1 + h2 + h3 + h4

    !!! 3rd Order Mesh Coef
    !MeshCoef0 = -1.0*( MeshH2*MeshH3*MeshH4 + MeshH1*MeshH3*MeshH4 + &
    !                   MeshH1*MeshH2*MeshH4 + MeshH1*MeshH2*MeshH3)/&
    !                 (MeshH1*MeshH2*MeshH3*MeshH4) 
    !MeshCoef1 =  1.0*( MeshH2*MeshH3*MeshH4)/(h1*h2*(h2 + h3)*(h2 + h3 + h4))
    !MeshCoef2 = -1.0*( MeshH1*MeshH3*MeshH4)/(MeshH2*h2*h3*(h3+h4))
    !MeshCoef3 =  1.0*( MeshH1*MeshH2*MeshH4)/(MeshH3*(h3+h2)*h3*h4)
    !MeshCoef4 = -1.0*( MeshH1*MeshH2*MeshH3)/(MeshH4*(h2+h3+h4)*(h3+h4)*h4)
!
!
        MeshH1 = h1                 
        MeshH2 = h1 + h2            
        MeshH3 = h1 + h2 + h3
        MeshH4 = h1 + h2 + h3 + h4

        OS_LBMeshCoef(iLon,iLat,iiAlt,1,iBlock) = &
                          -1.0*( MeshH2*MeshH3*MeshH4 + MeshH1*MeshH3*MeshH4 + &
                          MeshH1*MeshH2*MeshH4 + MeshH1*MeshH2*MeshH3)/&
                         (MeshH1*MeshH2*MeshH3*MeshH4) 
        OS_LBMeshCoef(iLon,iLat,iiAlt,2,iBlock) = &
                     1.0*( MeshH2*MeshH3*MeshH4)/(h1*h2*(h2 + h3)*(h2 + h3 + h4))
        OS_LBMeshCoef(iLon,iLat,iiAlt,3,iBlock) = &
                    -1.0*( MeshH1*MeshH3*MeshH4)/(MeshH2*h2*h3*(h3+h4))
        OS_LBMeshCoef(iLon,iLat,iiAlt,4,iBlock) = &
                     1.0*( MeshH1*MeshH2*MeshH4)/(MeshH3*(h3+h2)*h3*h4)
        OS_LBMeshCoef(iLon,iLat,iiAlt,5,iBlock) = &
                    -1.0*( MeshH1*MeshH2*MeshH3)/(MeshH4*(h2+h3+h4)*(h3+h4)*h4)
       enddo 
       do iAlt = nAlts+1, nAlts+3
          iiAlt = iAlt - nAlts
          hm1 = Altitude_GB(iLon,iLat,iAlt-1,iBlock) - Altitude_GB(iLon,iLat,iAlt-2,iBlock)
          hm2 = Altitude_GB(iLon,iLat,iAlt-2,iBlock) - Altitude_GB(iLon,iLat,iAlt-3,iBlock)
          hm3 = Altitude_GB(iLon,iLat,iAlt-3,iBlock) - Altitude_GB(iLon,iLat,iAlt-4,iBlock)
          hm4 = Altitude_GB(iLon,iLat,iAlt-4,iBlock) - Altitude_GB(iLon,iLat,iAlt-5,iBlock)

! 
          MeshHm1 = hm1                 
          MeshHm2 = hm1 + hm2            
          MeshHm3 = hm1 + hm2 + hm3
          MeshHm4 = hm1 + hm2 + hm3 + hm4
       !!! 3rd Order Mesh Coef
        OS_UBMeshCoef(iLon,iLat,iiAlt,1,iBlock) = &
              1.0*( MeshHm2*MeshHm3*MeshHm4 + MeshHm1*MeshHm3*MeshHm4 + &
                    MeshHm1*MeshHm2*MeshHm4 + MeshHm1*MeshHm2*MeshHm3)/&
                   (MeshHm1*MeshHm2*MeshHm3*MeshHm4) 
!!
        OS_UBMeshCoef(iLon,iLat,iiAlt,2,iBlock) = &
                -1.0*( MeshHm2*MeshHm3*MeshHm4)/&
                     (hm1*hm2*(hm2 + hm3)*(hm2 + hm3 + hm4))
        OS_UBMeshCoef(iLon,iLat,iiAlt,3,iBlock) = &
                 1.0*( MeshHm1*MeshHm3*MeshHm4)/(MeshHm2*hm2*hm3*(hm3+hm4))
        OS_UBMeshCoef(iLon,iLat,iiAlt,4,iBlock) = &
                -1.0*( MeshHm1*MeshHm2*MeshHm4)/(MeshHm3*(hm3+hm2)*hm3*hm4)
        OS_UBMeshCoef(iLon,iLat,iiAlt,5,iBlock) = &
                 1.0*( MeshHm1*MeshHm2*MeshHm3)/&
                     (MeshHm4*(hm2+hm3+hm4)*(hm3+hm4)*hm4)
       enddo 

     enddo !iLat = 1, nLats
   enddo !iLon = 1, nLons

end subroutine !init_mesh_alts

