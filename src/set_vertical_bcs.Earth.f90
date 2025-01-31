! Copyright 2021, the GITM Development Team (see srcDoc/dev_team.md for members)
! Full license can be found in LICENSE

!  ------------------------------------------------------------
!  set boundary conditions in the vertical direction
!  ------------------------------------------------------------ 

subroutine set_vertical_bcs(LogRho,LogNS,Vel_GD,Temp, LogINS, iVel, VertVel)

  use ModSizeGitm, only: nAlts
  use ModPlanet, only: nSpecies, nIonsAdvect, Mass, nIons, IsEarth,&
                       iN2_,iNO_,iO_3P_, iN_4S_
  use ModGITM, only: iEast_, iNorth_, iUp_
  use ModInputs
  use ModConstants
  use ModTime, only: UTime, iJulianDay,currenttime
  use ModVertical, only: &
       Lat, Lon, Gravity_G, &
       Altitude_G, dAlt_F, &
       MeanMajorMass_1d, &
       iLon1D, iLat1D, iBlock1D, &
       Centrifugal, InvRadialDistance_C, Coriolis, &
       MLatVertical, SZAVertical
  use ModIndicesInterfaces, only: get_HPI
  use ModTides, only: TidesNorth, TidesEast, TidesTemp

  use EUA_ModMsis90, ONLY: meter6

  implicit none

  real, intent(inout) :: &
       LogRho(-1:nAlts+2), &
       LogNS(-1:nAlts+2,nSpecies), &
       LogINS(-1:nAlts+2,nIons), &
       Vel_GD(-1:nAlts+2,3), &
       IVel(-1:nAlts+2,3), &
       Temp(-1:nAlts+2), &
       VertVel(-1:nAlts+2,nSpecies)

  integer :: iSpecies, iAlt, iDir
  real    :: InvScaleHeightS, InvScaleHgt, Alt, Lst, Ap = 4.0, dn, dt
  logical :: IsFirstTime = .true., UseMsisBCs = .false.
  real    :: HP, v(2)
  integer :: ierror
  real    :: temptemp
  real    :: logNS_Species(nSpecies)
  real :: n0, n1, n2, n3, n4, n5
  real :: tec
  
  integer, dimension(25) :: sw

  ! Useful for Upper Boundary Conditions--Hydrostatic
  real :: InvAtmScaleHeight
  real :: NS(-1:nAlts+2,1:nSpecies), NT(-1:nAlts+2)
  real :: SumRho
  real :: EffectiveGravity(-1:nAlts+2)
  real :: MeanGravity, MeanMass, MeanTemp
  logical :: IsHydrostatic(1:nSpecies)

  logical :: UsePlasmasphereBC
  
  ! Gradient Terms
  real :: dLogNS, dTemp, dVertVel
  real :: dLogINS

  ! Use for 4-th Order Forward Differences
  ! Need a 5-point Stencil
  real :: h1, h2, h3, h4
  real :: MeshH1, MeshH2, MeshH3, MeshH4
  real :: MeshCoef0, MeshCoef1, &
          MeshCoef2, MeshCoef3, &
          MeshCoef4

  ! Use for 4-th Order Backward Differences
  ! Need a 5-point Stencil
  real :: hm1, hm2, hm3, hm4
  real :: MeshHm1, MeshHm2, MeshHm3, MeshHm4
  real :: MeshCoefm0, MeshCoefm1, MeshCoefm2, MeshCoefm3, MeshCoefm4

  ! Temporary arrays of the lower boundary condition to pass to WP (WP-GITM)
  real :: LogNS_LBC(2,nSpecies), VelGD_LBC(2,3)

  !------------------------------------------------------------------------
  !------------------------------------------------------------------------
  ! Bottom
  !------------------------------------------------------------------------
  !------------------------------------------------------------------------

  NS(-1:nAlts+2,1:nSpecies) = exp(LogNS(-1:nAlts+2,1:nSpecies))

  do iAlt = -1, nAlts + 2
     EffectiveGravity(iAlt) = &
        Gravity_G(iAlt) + &
        Centrifugal / InvRadialDistance_C(iAlt) + & 
        (Vel_GD(iAlt,iNorth_)**2 + Vel_GD(iAlt,iEast_)**2) &
        * InvRadialDistance_C(iAlt) + & 
        Coriolis * Vel_GD(iAlt,iEast_)
  enddo 

  !------------------------------------------------------------
  ! This is so we can have a fixed lower BC (MSIS specified):
  !------------------------------------------------------------

  if (IsEarth) UseMsisBCs = UseMsis

  if (IsFirstTime .and. UseMsisBCs) then
     call meter6(.true.)
     sw = 1
     IsFirstTime = .true.
  endif

  if (UseMsisBCs) then

     ! We don't often have Ap, but have hemispheric power. So, use that:
     call get_HPI(CurrentTime, HP, iError)  
     if (iError > 0) hp = 40.0
     Ap = min(200.,max(-40.72 + 1.3 * HP, 10.))

     do iAlt = -1, 0
        Alt = Altitude_G(iAlt)/1000.0 * &
             (1.0 - &
             MsisOblateFactor/2.0 + &
             MsisOblateFactor * cos(Lat*3.1415/180.0))
        Lst = mod(UTime/3600.0+Lon/15.0,24.0)

        call msis_bcs(iJulianDay,UTime,Alt,Lat,Lon,Lst, &
             F107A,F107,AP,LogNS_Species, temptemp, &
             LogRho(iAlt),v)

        LogNS(iAlt,:) = logNS_Species
        
        if (.not. DuringPerturb) temp(iAlt) = temptemp

        vel_gd(iAlt,iEast_) = v(iEast_)
        vel_gd(iAlt,iNorth_) = v(iNorth_)

     enddo
  else
  endif ! Use MSIS BCs
!
  if (UseGSWMTides) then
     Vel_GD(-1:0,iEast_)  = TidesEast(iLon1D,iLat1D,1:2,iBlock1D)
     Vel_GD(-1:0,iNorth_) = TidesNorth(iLon1D,iLat1D,1:2,iBlock1D)
     Temp(-1:0)           = TidesTemp(iLon1D,iLat1D,1:2,iBlock1D) + Temp(-1:0)
  endif
  if (UseWACCMTides) then
     Vel_GD(-1:0,iEast_)  = TidesEast(iLon1D,iLat1D,1:2,iBlock1D)
     Vel_GD(-1:0,iNorth_) = TidesNorth(iLon1D,iLat1D,1:2,iBlock1D)
     Temp(-1:0)           = TidesTemp(iLon1D,iLat1D,1:2,iBlock1D)
  endif

  ! For WP-GITM: add neutral atmospheric perturbations caused by
  ! tsunami or earthquake
  if (UseBcPerturbation) then
     ! some compilers do not work with passing non-continuous sub-arrays,
     ! extract sub-arrays to temporary arrays to be safe.
     LogNS_LBC = LogNS(-1:0, :)
     VelGD_LBC = Vel_GD(-1:0, iEast_ : iUp_)
     call user_bc_perturbation(LogRho(-1:0), &
                               LogNS_LBC, &
                               VelGD_LBC, &
                               Temp(-1:0))
     LogNS(-1:0, :) = LogNS_LBC
     Vel_GD(-1:0, iEast_ : iUp_) = VelGD_LBC
  endif

  ! If we use MSIS or Perturbation Settings, then just allow the data to specify the fields
  !if (UseBcPerturbation) then
  if(UseMSISBCs .or. UseGSWMTides .or. &
     UseWACCMTides .or. UseBcPerturbation) then
    ! Don't update the lower boundary cells of the Densities
    ! Let them be specified by the data

    do iSpecies = 1, nSpecies
       if (IsPhotoChemical(iSpecies)) then
          do iAlt = 0, -1, -1
             h1 = dAlt_F(iAlt+2) ! dAlt_F(1) = Alt(1) - Alt(0);  h1 in notes  
             h2 = dAlt_F(iAlt+3) ! dAlt_F(2) = Alt(2) - Alt(1);  h2 in notes
             h3 = dAlt_F(iAlt+4) ! dAlt_F(3) = Alt(3) - Alt(2);  h3 in notes
             h4 = dAlt_F(iAlt+5) ! dAlt_F(4) = Alt(4) - Alt(3);  h4 in notes

             ! Mesh Coefficients are summations over the individual mesh scales
             MeshH1 = h1                 
             MeshH2 = h1 + h2            
             MeshH3 = h1 + h2 + h3
             MeshH4 = h1 + h2 + h3 + h4

             !!! 3rd Order Mesh Coef
             MeshCoef0 = -1.0*( MeshH2*MeshH3*MeshH4 + MeshH1*MeshH3*MeshH4 + &
                                MeshH1*MeshH2*MeshH4 + MeshH1*MeshH2*MeshH3)/&
                              (MeshH1*MeshH2*MeshH3*MeshH4) 
             MeshCoef1 =  1.0*( MeshH2*MeshH3*MeshH4)/&
                               (h1*h2*(h2 + h3)*(h2 + h3 + h4))
             MeshCoef2 = -1.0*( MeshH1*MeshH3*MeshH4)/(MeshH2*h2*h3*(h3+h4))
             MeshCoef3 =  1.0*( MeshH1*MeshH2*MeshH4)/(MeshH3*(h3+h2)*h3*h4)
             MeshCoef4 = -1.0*( MeshH1*MeshH2*MeshH3)/&
                               (MeshH4*(h2+h3+h4)*(h3+h4)*h4)

             dLogNS = LogNS(iAlt+1,iSpecies)*MeshCoef0 + & 
                      LogNS(iAlt+2,iSpecies)*MeshCoef1 + &
                      LogNS(iAlt+3,iSpecies)*MeshCoef2 + &
                      LogNS(iAlt+4,iSpecies)*MeshCoef3 + &
                      LogNS(iAlt+5,iSpecies)*MeshCoef4

             LogNS(iAlt, iSpecies) = &
                  LogNS(iAlt+1,iSpecies) - dAlt_F(iAlt+1) * dLogNS 

          ! Only allow the NO density to decrease by roughly 1/3 scale
          ! height max
          !if (iSpecies == iNO_) then
          !   if (LogNS(iAlt,iSpecies) < LogNS(iAlt+1,iSpecies) - 0.333) then
          !      LogNS(iAlt,iSpecies) = LogNS(iAlt+1,iSpecies) - 0.333
          !   endif
          !endif
          enddo !iAlt = 0,-1,-1
       endif !(.not. IsPhotoChemical(iSpecies)) then
    enddo !iSpecies = 1, nSpecies

  else ! Use MSIS/WACCM/GSWM

    !----------- CASE: No Specified LBC --------------------
    ! No Data Specification (keep 0 cell fixed)
    ! Extend the densities downward
    ! Neutral Densities
    do iSpecies = 1, nSpecies

       if (IsPhotoChemical(iSpecies)) then
           do iAlt = 0, -1, -1
              h1 = dAlt_F(iAlt+2) ! dAlt_F(1) = Alt(1) - Alt(0);  h1 in notes  
              h2 = dAlt_F(iAlt+3) ! dAlt_F(2) = Alt(2) - Alt(1);  h2 in notes
              h3 = dAlt_F(iAlt+4) ! dAlt_F(3) = Alt(3) - Alt(2);  h3 in notes
              h4 = dAlt_F(iAlt+5) ! dAlt_F(4) = Alt(4) - Alt(3);  h4 in notes

              ! Mesh Coefficients are summations over the individual mesh scales
              MeshH1 = h1                 
              MeshH2 = h1 + h2            
              MeshH3 = h1 + h2 + h3
              MeshH4 = h1 + h2 + h3 + h4

              !!! 3rd Order Mesh Coef
              MeshCoef0 = -1.0*( MeshH2*MeshH3*MeshH4 + MeshH1*MeshH3*MeshH4 + &
                                 MeshH1*MeshH2*MeshH4 + MeshH1*MeshH2*MeshH3)/&
                              (MeshH1*MeshH2*MeshH3*MeshH4) 
              MeshCoef1 =  1.0*( MeshH2*MeshH3*MeshH4)/&
                                (h1*h2*(h2 + h3)*(h2 + h3 + h4))
              MeshCoef2 = -1.0*( MeshH1*MeshH3*MeshH4)/(MeshH2*h2*h3*(h3+h4))
              MeshCoef3 =  1.0*( MeshH1*MeshH2*MeshH4)/(MeshH3*(h3+h2)*h3*h4)
              MeshCoef4 = -1.0*( MeshH1*MeshH2*MeshH3)/&
                                (MeshH4*(h2+h3+h4)*(h3+h4)*h4)

              dLogNS = LogNS(iAlt+1,iSpecies)*MeshCoef0 + & 
                       LogNS(iAlt+2,iSpecies)*MeshCoef1 + &
                       LogNS(iAlt+3,iSpecies)*MeshCoef2 + &
                       LogNS(iAlt+4,iSpecies)*MeshCoef3 + &
                       LogNS(iAlt+5,iSpecies)*MeshCoef4

              LogNS(iAlt, iSpecies) = &
                   LogNS(iAlt+1,iSpecies) - dAlt_F(iAlt+1) * dLogNS 
           enddo !iAlt = 0,-1,-1
       else  ! Is PhotoChemical Check
          iAlt = -1
          h1 = dAlt_F(iAlt+2) ! dAlt_F(1) = Alt(1) - Alt(0);  h1 in notes  
          h2 = dAlt_F(iAlt+3) ! dAlt_F(2) = Alt(2) - Alt(1);  h2 in notes
          h3 = dAlt_F(iAlt+4) ! dAlt_F(3) = Alt(3) - Alt(2);  h3 in notes
          h4 = dAlt_F(iAlt+5) ! dAlt_F(4) = Alt(4) - Alt(3);  h4 in notes

          ! Mesh Coefficients are summations over the individual mesh scales
          MeshH1 = h1                 
          MeshH2 = h1 + h2            
          MeshH3 = h1 + h2 + h3
          MeshH4 = h1 + h2 + h3 + h4

          !!! 3rd Order Mesh Coef
          MeshCoef0 = -1.0*( MeshH2*MeshH3*MeshH4 + MeshH1*MeshH3*MeshH4 + &
                             MeshH1*MeshH2*MeshH4 + MeshH1*MeshH2*MeshH3)/&
                          (MeshH1*MeshH2*MeshH3*MeshH4) 
          MeshCoef1 =  1.0*( MeshH2*MeshH3*MeshH4)/&
                            (h1*h2*(h2 + h3)*(h2 + h3 + h4))
          MeshCoef2 = -1.0*( MeshH1*MeshH3*MeshH4)/(MeshH2*h2*h3*(h3+h4))
          MeshCoef3 =  1.0*( MeshH1*MeshH2*MeshH4)/(MeshH3*(h3+h2)*h3*h4)
          MeshCoef4 = -1.0*( MeshH1*MeshH2*MeshH3)/&
                            (MeshH4*(h2+h3+h4)*(h3+h4)*h4)

          dLogNS = LogNS(iAlt+1,iSpecies)*MeshCoef0 + & 
                   LogNS(iAlt+2,iSpecies)*MeshCoef1 + &
                   LogNS(iAlt+3,iSpecies)*MeshCoef2 + &
                   LogNS(iAlt+4,iSpecies)*MeshCoef3 + &
                   LogNS(iAlt+5,iSpecies)*MeshCoef4

          ! Only extend down to -1 Cells
          LogNS(iAlt, iSpecies) = &
               LogNS(iAlt+1,iSpecies) - dAlt_F(iAlt+1) * dLogNS 
       endif ! PhotoChemical Check

   enddo  ! iSpecies loop

   ! If there is no LBC Specification
   ! Set Vel_GD(East,North) = 0.0  at the LBC
    do iAlt = 0, -1, -1
       Vel_GD(iAlt,iEast_ ) = 0.0
       Vel_GD(iAlt,iNorth_) = 0.0
    enddo !iAlt = 0, -1, -1

  endif ! Use MSIS/WACCM/GSWM

  ! JMB: CONTINUE LBC SETTINGS
  do iAlt = 0, -1, -1
     h1 = dAlt_F(iAlt+2) ! dAlt_F(1) = Alt(1) - Alt(0);  h1 in notes  
     h2 = dAlt_F(iAlt+3) ! dAlt_F(2) = Alt(2) - Alt(1);  h2 in notes
     h3 = dAlt_F(iAlt+4) ! dAlt_F(3) = Alt(3) - Alt(2);  h3 in notes
     h4 = dAlt_F(iAlt+5) ! dAlt_F(4) = Alt(4) - Alt(3);  h4 in notes

     ! Mesh Coefficients are summations over the individual mesh scales
     MeshH1 = h1                 
     MeshH2 = h1 + h2            
     MeshH3 = h1 + h2 + h3
     MeshH4 = h1 + h2 + h3 + h4

     !!! 3rd Order Mesh Coef
     MeshCoef0 = -1.0*( MeshH2*MeshH3*MeshH4 + MeshH1*MeshH3*MeshH4 + &
                        MeshH1*MeshH2*MeshH4 + MeshH1*MeshH2*MeshH3)/&
                      (MeshH1*MeshH2*MeshH3*MeshH4) 
     MeshCoef1 =  1.0*( MeshH2*MeshH3*MeshH4)/&
                       (h1*h2*(h2 + h3)*(h2 + h3 + h4))
     MeshCoef2 = -1.0*( MeshH1*MeshH3*MeshH4)/(MeshH2*h2*h3*(h3+h4))
     MeshCoef3 =  1.0*( MeshH1*MeshH2*MeshH4)/(MeshH3*(h3+h2)*h3*h4)
     MeshCoef4 = -1.0*( MeshH1*MeshH2*MeshH3)/&
                       (MeshH4*(h2+h3+h4)*(h3+h4)*h4)

    ! -------------------------------------------------
    ! Ions
    ! -------------------------------------------------
    ! Ions Float at the LBC
    do iSpecies = 1, nIons-1
       ! Simply assume no gradient, since this is a region in which
       ! chemitry is dominant.  It shouldn't matter, really:
       LogINS(iAlt,iSpecies) = LogINS(iAlt+1,iSpecies)
    enddo ! Ions

    ! Electon Density:
    if (UseImprovedIonAdvection) then
       LogINS(iAlt,nIons) = sum(LogINS(iAlt,1:nIons-1))
    else
       LogINS(iAlt,nIons) = alog(sum(exp(LogINS(iAlt,1:nIons-1))))
    endif
    ! Set the Ion Bulk Winds
    do iDir = 1, 3
       !dVertVel = IVel(iAlt+1,iDir)*MeshCoef0 + & 
       !           IVel(iAlt+2,iDir)*MeshCoef1 + &
       !           IVel(iAlt+3,iDir)*MeshCoef2 + &
       !           IVel(iAlt+4,iDir)*MeshCoef3 + &
       !           IVel(iAlt+5,iDir)*MeshCoef4
!
!       IVel(iAlt  ,iDir) = IVel(iAlt+1,iDir) - dAlt_F(iAlt+1)*dVertVel 
       ! Fix the winds at the LBC
       IVel(iAlt  ,iDir) = 0.0
    enddo  ! iDir

    ! -------------------------------------------------
    ! Return to the neutrals, to deal with (vertical) winds now...
    ! -------------------------------------------------
    
    do iSpecies = 1, nSpecies
         dVertVel = VertVel(iAlt+1, iSpecies) * MeshCoef0 + & 
                    VertVel(iAlt+2, iSpecies) * MeshCoef1 + &
                    VertVel(iAlt+3, iSpecies) * MeshCoef2 + &
                    VertVel(iAlt+4, iSpecies) * MeshCoef3 + &
                    VertVel(iAlt+5, iSpecies) * MeshCoef4

         VertVel(iAlt, iSpecies) = VertVel(iAlt+1,iSpecies) - &
                                   dAlt_F(iAlt+1) * dVertVel 
    enddo ! nSpecies

    ! Update NS & LogRho to use in the vertical wind calculation:
    NS(iAlt,1:nSpecies) = exp(LogNS(iAlt,1:nSpecies))
    SumRho = 0.0     
    do iSpecies = 1, nSpecies
       SumRho  = SumRho  + Mass(iSpecies) * NS(iAlt, iSpecies)
    enddo
    LogRho(iAlt) = alog(SumRho)

    ! Calculate the bulk vertical winds:
    Vel_GD(iAlt,iUp_) = 0.0
    do iSpecies = 1, nSpecies
       Vel_GD(iAlt,iUp_) = Vel_GD(iAlt,iUp_) + &
            NS(iAlt,iSpecies) * Mass(iSpecies) * VertVel(iAlt,iSpecies) / &
            exp(LogRho(iAlt))
    enddo 

  enddo ! End Outer IAlt Loop (0, -1, -1) 

  ! Special Case for PhotoChemical Species
  ! Allow O to flow upward
  do iSpecies = 1, nSpecies
     if ( (IsPhotoChemical(iSpecies)) .and. &
          (VertVel(1, iSpecies) .gt. 0.0) .and. &
          (iSpecies .ne. iO_3P_) ) then 
        ! Do NOT allow photochemical species to flow upward!
        ! The problem becomes unconstrained
        VertVel( 0, iSpecies) = -1.0 * VertVel(1, iSpecies)
        VertVel(-1, iSpecies) = -1.0 * VertVel(1, iSpecies)
     endif  
  enddo 

  !------------------------------------------------------------------------
  !------------------------------------------------------------------------
  ! TOP
  !------------------------------------------------------------------------
  !------------------------------------------------------------------------

  ! Slip flow at the top
  ! Assume zero gradients in the velocities & temps
  !IVel(nAlts+1:nAlts+2, iEast_)  = IVel(nAlts, iEast_)
  !IVel(nAlts+1:nAlts+2, iNorth_) = IVel(nAlts, iNorth_)


  !! Vertical wind upper BCs:
  !do iSpecies = 1, nSpecies
  !   VertVel(nAlts+1, iSpecies) =  VertVel(nAlts, iSpecies)
  !   VertVel(nAlts+2, iSpecies) =  VertVel(nAlts, iSpecies)
  !   ! Can't have photochemical species flowing downward
  !   if (IsPhotoChemical(iSpecies) .and. &
  !        (VertVel(nAlts,iSpecies) .lt. 0.0) ) then
  !      VertVel(nAlts+1, iSpecies) = -VertVel(nAlts, iSpecies)
  !      VertVel(nAlts+2, iSpecies) = -VertVel(nAlts-1, iSpecies)
  !   endif
  !enddo 

  !if (Vel_GD(nAlts,iUp_) < 0.0) then
  !   Vel_GD(nAlts+1, iUp_) = -Vel_GD(nAlts, iUp_)
  !   Vel_GD(nAlts+2, iUp_) = -Vel_GD(nAlts, iUp_)
  !else
  !   Vel_GD(nAlts+1, iUp_) =  Vel_GD(nAlts, iUp_)
  !   Vel_GD(nAlts+2, iUp_) =  Vel_GD(nAlts, iUp_)
  !endif

  ! Vertical Ion Drifts:
  !if (IVel(nAlts,iUp_) < 0.0) then
  !   IVel(nAlts+1,iUp_) = -IVel(nAlts,iUp_)
  !   IVel(nAlts+2,iUp_) = -IVel(nAlts,iUp_)
  !else
  !   IVel(nAlts+1,iUp_) = IVel(nAlts,iUp_)
  !   IVel(nAlts+2,iUp_) = IVel(nAlts,iUp_)
  !endif

  ! Constant temperature (zero gradient)

  do iAlt = nAlts+1, nAlts+2
     hm1 = dAlt_F(iAlt  )
     hm2 = dAlt_F(iAlt-1)
     MeshCoefm0 = (2.0*hm1 + hm2)/(hm1*(hm1+hm2))
     MeshCoefm1 = -1.0*(hm1 + hm2)/(hm1*hm2)
     MeshCoefm2 = hm1/(hm2*(hm1+hm2))

    ! dT/dr = 0
     Temp(iAlt) = &
        -1.0*(MeshCoefm1*Temp(iAlt-1) + & 
              MeshCoefm2*Temp(iAlt-2))/ & 
              MeshCoefm0
     ! dVs/dr = 0.0
     do iSpecies =1, nSpecies
        VertVel(iAlt,iSpecies) = &
           -1.0*(MeshCoefm1*VertVel(iAlt-1,iSpecies) + & 
                 MeshCoefm2*VertVel(iAlt-2,iSpecies))/ & 
                 MeshCoefm0
     enddo !iSpecies =1, nSpecies

     do iDir = iEast_, iNorth_
        Vel_GD(iAlt,iDir) = &
           -1.0*(MeshCoefm1*Vel_GD(iAlt-1,iDir) + & 
                 MeshCoefm2*Vel_GD(iAlt-2,iDir))/ & 
                 MeshCoefm0
        IVel(iAlt,iDir) = &
           -1.0*(MeshCoefm1*IVel(iAlt-1,iDir) + & 
                 MeshCoefm2*IVel(iAlt-2,iDir))/ & 
                 MeshCoefm0
     enddo !iDir = 1, 3
  enddo !iAlt = nAlts+1, nAlts+2

  !! Vertical wind upper BCs:
  do iSpecies = 1, nSpecies
     ! Can't have photochemical species flowing downward
     if (IsPhotoChemical(iSpecies) .and. &
          (VertVel(nAlts,iSpecies) .lt. 0.0) ) then
        VertVel(nAlts+1, iSpecies) = 0.0
        VertVel(nAlts+2, iSpecies) = 0.0
     endif
  enddo 

  ! Limit Ion Downwelling
  if (IVel(nAlts,iUp_) < 0.0) then
     IVel(nAlts+1,iUp_) = 0.0
     IVel(nAlts+2,iUp_) = 0.0
  endif


  ! Hydrostatic pressure for the neutrals
  do iSpecies=1,nSpecies
     do iAlt = nAlts+1, nAlts+2
        MeanGravity = -0.5*(EffectiveGravity(iAlt) + EffectiveGravity(iAlt-1))
        MeanTemp =  0.5*( Temp(iAlt-1) + Temp(iAlt) )
        InvScaleHeightS =  MeanGravity * Mass(iSpecies) / &
                           (MeanTemp * Boltzmanns_Constant)

        NS(iAlt, iSpecies) = &
             NS(iAlt-1, iSpecies) * (Temp(iAlt-1) / Temp(iAlt)) * &
             exp(-InvScaleHeightS * dAlt_F(iAlt)) 

        LogNS(iAlt,iSpecies) = alog(NS(iAlt,iSpecies))
     enddo
  enddo

  UsePlasmasphereBC = .false.
  if (UseNighttimeIonBCs) then
     if ( SZAVertical > Pi/2 .and. &
          (( MLatVertical > 30.0 .and. MLatVertical <70.0) .or. &
          ( MLatVertical > -60.0 .and. MLatVertical < -30.0 ))) then
        UsePlasmasphereBC = .true.
     endif
  endif

  if (UseImprovedIonAdvection) then
     tec = sum(dAlt_f(1:nAlts) * LogINS(1:nAlts,1))/1e16
  else
     tec = sum(dAlt_f(1:nAlts) * exp(LogINS(1:nAlts,1)))/1e16
  endif

  do iAlt = nAlts + 1, nAlts + 2

     hm1 = dAlt_F(iAlt-1) ! 
     hm2 = dAlt_F(iAlt-2) ! 
     hm3 = dAlt_F(iAlt-3) ! 
     hm4 = dAlt_F(iAlt-4) ! 

     ! Mesh Coefficients are summations over the individual mesh scales
     MeshHm1 = hm1                 
     MeshHm2 = hm1 + hm2            
     MeshHm3 = hm1 + hm2 + hm3
     MeshHm4 = hm1 + hm2 + hm3 + hm4

     ! 3rd Order Mesh Coef
     MeshCoefm0 =  1.0*( MeshHm2*MeshHm3*MeshHm4 + MeshHm1*MeshHm3*MeshHm4 + &
          MeshHm1*MeshHm2*MeshHm4 + MeshHm1*MeshHm2*MeshHm3)/&
          (MeshHm1*MeshHm2*MeshHm3*MeshHm4) 
     MeshCoefm1 = -1.0*( MeshHm2*MeshHm3*MeshHm4)/&
          (hm1*hm2*(hm2 + hm3)*(hm2 + hm3 + hm4))
     MeshCoefm2 =  1.0*( MeshHm1*MeshHm3*MeshHm4)/(MeshHm2*hm2*hm3*(hm3+hm4))
     MeshCoefm3 = -1.0*( MeshHm1*MeshHm2*MeshHm4)/(MeshHm3*(hm3+hm2)*hm3*hm4)
     MeshCoefm4 =  1.0*( MeshHm1*MeshHm2*MeshHm3)/&
          (MeshHm4*(hm2+hm3+hm4)*(hm3+hm4)*hm4)
           
     do iSpecies=1,nIons-1

        if (UseImprovedIonAdvection) then

           call check_for_negative_densities

           n0 = alog(LogINS(iAlt  ,iSpecies))
           n1 = alog(LogINS(iAlt-1,iSpecies))
           n2 = alog(LogINS(iAlt-2,iSpecies))
           n3 = alog(LogINS(iAlt-3,iSpecies))
           n4 = alog(LogINS(iAlt-4,iSpecies))
           n5 = alog(LogINS(iAlt-5,iSpecies))
           if (DoCheckForNans) then
              if (isnan(n0)) write(*,*) 'n0 :',iAlt, LogINS(iAlt  ,iSpecies)
              if (isnan(n1)) write(*,*) 'n1 :',iAlt-1, LogINS(iAlt-1  ,iSpecies)
              if (isnan(n2)) write(*,*) 'n2 :',iAlt-2, LogINS(iAlt-2  ,iSpecies)
              if (isnan(n3)) write(*,*) 'n3 :',iAlt-3, LogINS(iAlt-3  ,iSpecies)
              if (isnan(n4)) write(*,*) 'n4 :',iAlt-4, LogINS(iAlt-4  ,iSpecies)
              if (isnan(n5)) write(*,*) 'n5 :',iAlt-5, LogINS(iAlt-5  ,iSpecies)
           endif
        else
           n0 = LogINS(iAlt  ,iSpecies)
           n1 = LogINS(iAlt-1,iSpecies)
           n2 = LogINS(iAlt-2,iSpecies)
           n3 = LogINS(iAlt-3,iSpecies)
           n4 = LogINS(iAlt-4,iSpecies)
           n5 = LogINS(iAlt-5,iSpecies)
        endif

        dn = MeshCoefm0*n1 + &  
             MeshCoefm1*n2 + &  
             MeshCoefm2*n3 + &  
             MeshCoefm3*n4 + &  
             MeshCoefm4*n5      

        if (tec < MinTEC .and. UsePlasmasphereBC) then
           ! do nothing
        else

           ! Limit the slope of the gradient to be negative, since the
           ! ion density should be decreasing at the top of the model.
           
           if (dn > -0.001*n0/dAlt_F(nAlts)) dn = -0.001*n0/dAlt_F(nAlts)
           n0 = n1 + dn*dAlt_F(iAlt)

           if (UseImprovedIonAdvection) then
              LogINS(iAlt,iSpecies) = exp(n0)
           else
              LogINS(iAlt,iSpecies) = n0
           endif

        endif
              
        if (DoCheckForNans) then
           if (isnan(LogINS(iAlt,1))) &
                write(*,*) 'svbc ',iAlt,LogINS(iAlt,1), n0, dn, &
                n1, n2, n3, n4, n5, tec, MinTEC, UsePlasmasphereBC, &
                dAlt_F(iAlt)
        endif

     enddo

     ! Electon Density:
     if (UseImprovedIonAdvection) then
        LogINS(iAlt,nIons) = sum(LogINS(iAlt,1:nIons-1))
     else
        LogINS(iAlt,nIons) = alog(sum(exp(LogINS(iAlt,1:nIons-1))))
     endif
           
  enddo

  ! JMB: Final Neutral Wind Specifications
  do iAlt = nAlts+1, nAlts+2
    ! Update NS & LogRho to use in the vertical wind calculation:
    NS(iAlt,1:nSpecies) = exp(LogNS(iAlt,1:nSpecies))
    SumRho = 0.0     
    do iSpecies = 1, nSpecies
       SumRho  = SumRho  + Mass(iSpecies) * NS(iAlt, iSpecies)
    enddo
    LogRho(iAlt) = alog(SumRho)

    ! Calculate the bulk vertical winds:
    Vel_GD(iAlt,iUp_) = 0.0
    do iSpecies = 1, nSpecies
       Vel_GD(iAlt,iUp_) = Vel_GD(iAlt,iUp_) + &
            NS(iAlt,iSpecies) * Mass(iSpecies) * VertVel(iAlt,iSpecies) / &
            exp(LogRho(iAlt))
    enddo 
  enddo !iAlt = nAlts+1, nAlts+2

contains

  subroutine check_for_negative_densities

    integer :: iAltSub
    logical :: IsFound

    if (minval(LogINS(iAlt-5:iAlt, iSpecies)) < 0.0) then

       if (DoCheckForNaNs) write(*,*) "Correcting negative ion density in set_vertical_bcs : ", iAlt, iSpecies

       if (LogINS(iAlt-6, iSpecies) < 0.0) then
          write(*,*) 'Negative Ion density too close to the upper boundary: ', iSpecies
          call stop_gitm('Stopping in set_vertical_bcs')
       endif

       do iAltSub = iAlt-5, iAlt-1

          if (LogINS(iAltSub,iSpecies) < 0.0) then

             if (LogINS(iAltSub-1,iSpecies) > 0.0 .and. &
                  LogINS(iAltSub+1,iSpecies) > 0.0) then
                ! Average two points around it:
                LogINS(iAltSub,iSpecies) = &
                     (LogINS(iAltSub-1,iSpecies) + LogINS(iAltSub+1,iSpecies))/2.0
             else
                ! Or take previous point and decrease it by 5 percent:
                LogINS(iAltSub,iSpecies) = LogINS(iAltSub-1,iSpecies) * 0.95
             endif

          endif

       enddo

       ! Top point:
       if (LogINS(iAlt,iSpecies) < 0.0) then
          ! Take previous point and decrease it by 5 percent:
          LogINS(iAlt,iSpecies) = LogINS(iAlt-1,iSpecies) * 0.95
       endif

    endif


  end subroutine check_for_negative_densities


end subroutine set_vertical_bcs

