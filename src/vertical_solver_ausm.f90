! Copyright 2021, the GITM Development Team (see srcDoc/dev_team.md for members)
! Full license can be found in LICENSE

!\
! ------------------------------------------------------------
! advance
! ------------------------------------------------------------
!/

subroutine check_ion_densities(iDen)

  use ModSizeGitm
  use ModPlanet, only: nIonsAdvect, nIons
  real, intent(inout) :: iDen(-1:nAlts+2,nIons)

  do iIon = 1, nIonsAdvect
     do iAlt = 1, nAlts+2
        if (iDen(iAlt, iIon) < 0.0) then
           iDen(iAlt,iIon) = max(iDen(iAlt-1,iIon)*0.99,10.0)
        endif
     enddo
  enddo

end subroutine check_ion_densities

subroutine advance_vertical_1d_ausm

  use ModVertical
  use ModPlanet, only: Mass, cSpecies
  use ModGITM, ONLY : Dt, iCommGITM, iProc, iUp_
  use ModInputs, only: UseBarriers, iDebugLevel, IsPhotoChemical
  implicit none
  !-----------------------------------------------------------
  integer :: iError, iAlt, iSpecies, iDir
  !!!!! Variables for the Runga-Kutta 4th Order Time-stepping
  real :: OrigLogNS(-1:nAlts+2,1:nSpecies)
  real :: OrigLogINS(-1:nAlts+2,1:nIons)
  real :: OrigLogRho(-1:nAlts+2)
  real :: OrigVel_GD(-1:nAlts+2,1:3)
  real :: OrigTemp(-1:nAlts+2)
  real :: OrigVS(-1:nAlts+2,1:nSpecies)

  real :: UpdatedLogNS(-1:nAlts+2,1:nSpecies)
  real :: UpdatedLogINS(-1:nAlts+2,1:nIons)
  real :: UpdatedLogRho(-1:nAlts+2)
  real :: UpdatedVel_GD(-1:nAlts+2,1:3)
  real :: UpdatedTemp(-1:nAlts+2)
  real :: UpdatedVS(-1:nAlts+2,1:nSpecies)

  real :: FinalLogNS(-1:nAlts+2,1:nSpecies)
  real :: FinalLogINS(-1:nAlts+2,1:nIons)
  real :: FinalLogRho(-1:nAlts+2)
  real :: FinalVel_GD(-1:nAlts+2,1:3)
  real :: FinalTemp(-1:nAlts+2)
  real :: FinalVS(-1:nAlts+2,1:nSpecies)

!!! RStage-4 Coefficients
  real :: Stage1LogNS(-1:nAlts+2,1:nSpecies)
  real :: Stage1LogINS(-1:nAlts+2,1:nIons)
  real :: Stage1LogRho(-1:nAlts+2)
  real :: Stage1Vel_GD(-1:nAlts+2,1:3)
  real :: Stage1Temp(-1:nAlts+2)
  real :: Stage1VS(-1:nAlts+2,1:nSpecies)

  real :: Stage2LogNS(-1:nAlts+2,1:nSpecies)
  real :: Stage2LogINS(-1:nAlts+2,1:nIons)
  real :: Stage2LogRho(-1:nAlts+2)
  real :: Stage2Vel_GD(-1:nAlts+2,1:3)
  real :: Stage2Temp(-1:nAlts+2)
  real :: Stage2VS(-1:nAlts+2,1:nSpecies)

  real :: Stage3LogNS(-1:nAlts+2,1:nSpecies)
  real :: Stage3LogINS(-1:nAlts+2,1:nIons)
  real :: Stage3LogRho(-1:nAlts+2)
  real :: Stage3Vel_GD(-1:nAlts+2,1:3)
  real :: Stage3Temp(-1:nAlts+2)
  real :: Stage3VS(-1:nAlts+2,1:nSpecies)

  real :: Stage4LogNS(-1:nAlts+2,1:nSpecies)
  real :: Stage4LogINS(-1:nAlts+2,1:nIons)
  real :: Stage4LogRho(-1:nAlts+2)
  real :: Stage4Vel_GD(-1:nAlts+2,1:3)
  real :: Stage4Temp(-1:nAlts+2)
  real :: Stage4VS(-1:nAlts+2,1:nSpecies)

  real :: DtIn

  if (UseBarriers) call MPI_BARRIER(iCommGITM,iError)
  if (iDebugLevel > 6) write(*,*) "=======> vertical bcs 1", iproc

!!! =================
!!! General RK4 Update:
!!! Y(n+1) = Y(n) + Dt/6*(k1 + 2k2 + 2k3 + k4)
!!! Time(n+1) = Time(n) + Dt
!!! 
!!! k1 = f(tn,yn)
!!! k2 = f(tn + Dt/2, Yn + Dt/2*k1)
!!! k3 = f(tn + Dt/2, Yn + Dt/2*k2)
!!! k4 = f(tn + Dt, Yn + Dt*k3)
!!! =================

  ! Step 1, Fill in Ghost Cells
  call set_vertical_bcs(LogRho,LogNS,Vel_GD,Temp,LogINS,IVel,VertVel)
  ! Store our original time step from GITM (CFL-limited).

!!! Set the Original State -> Orig = Y(n)
   OrigLogNS(-1:nAlts+2,1:nSpecies)    =   LogNS(-1:nAlts+2,1:nSpecies)
  OrigLogINS(-1:nAlts+2,1:nIons) =  LogINS(-1:nAlts+2,1:nIons)
  OrigLogRho(-1:nAlts+2)               =  LogRho(-1:nAlts+2)
  OrigVel_GD(-1:nAlts+2,1:3)           =  Vel_GD(-1:nAlts+2,1:3)
    OrigTemp(-1:nAlts+2)               =    Temp(-1:nAlts+2)
      OrigVS(-1:nAlts+2,1:nSpecies)    = VertVel(-1:nAlts+2,1:nSpecies)

   Stage1LogNS(-1:nAlts+2,1:nSpecies)    =   LogNS(-1:nAlts+2,1:nSpecies)
  Stage1LogINS(-1:nAlts+2,1:nIons) =  LogINS(-1:nAlts+2,1:nIons)
  Stage1LogRho(-1:nAlts+2)               =  LogRho(-1:nAlts+2)
  Stage1Vel_GD(-1:nAlts+2,1:3)           =  Vel_GD(-1:nAlts+2,1:3)
    Stage1Temp(-1:nAlts+2)               =    Temp(-1:nAlts+2)
      Stage1VS(-1:nAlts+2,1:nSpecies)    = VertVel(-1:nAlts+2,1:nSpecies)

   Stage2LogNS(-1:nAlts+2,1:nSpecies)    =   LogNS(-1:nAlts+2,1:nSpecies)
  Stage2LogINS(-1:nAlts+2,1:nIons) =  LogINS(-1:nAlts+2,1:nIons)
  Stage2LogRho(-1:nAlts+2)               =  LogRho(-1:nAlts+2)
  Stage2Vel_GD(-1:nAlts+2,1:3)           =  Vel_GD(-1:nAlts+2,1:3)
    Stage2Temp(-1:nAlts+2)               =    Temp(-1:nAlts+2)
      Stage2VS(-1:nAlts+2,1:nSpecies)    = VertVel(-1:nAlts+2,1:nSpecies)

   Stage3LogNS(-1:nAlts+2,1:nSpecies)    =   LogNS(-1:nAlts+2,1:nSpecies)
  Stage3LogINS(-1:nAlts+2,1:nIons) =  LogINS(-1:nAlts+2,1:nIons)
  Stage3LogRho(-1:nAlts+2)               =  LogRho(-1:nAlts+2)
  Stage3Vel_GD(-1:nAlts+2,1:3)           =  Vel_GD(-1:nAlts+2,1:3)
    Stage3Temp(-1:nAlts+2)               =    Temp(-1:nAlts+2)
      Stage3VS(-1:nAlts+2,1:nSpecies)    = VertVel(-1:nAlts+2,1:nSpecies)

   Stage4LogNS(-1:nAlts+2,1:nSpecies)    =   LogNS(-1:nAlts+2,1:nSpecies)
  Stage4LogINS(-1:nAlts+2,1:nIons) =  LogINS(-1:nAlts+2,1:nIons)
  Stage4LogRho(-1:nAlts+2)               =  LogRho(-1:nAlts+2)
  Stage4Vel_GD(-1:nAlts+2,1:3)           =  Vel_GD(-1:nAlts+2,1:3)
    Stage4Temp(-1:nAlts+2)               =    Temp(-1:nAlts+2)
      Stage4VS(-1:nAlts+2,1:nSpecies)    = VertVel(-1:nAlts+2,1:nSpecies)

  NewLogNS = LogNS
  NewLogINS = LogINS
  NewLogRho = LogRho
  NewVel_GD = Vel_GD
  NewTemp = Temp
  NewVertVel = VertVel

  DtIn = 0.5*Dt  !!! Store this so that it doesn't change
  call advance_vertical_1stage_ausm(DtIn, &
       LogRho, LogNS, Vel_GD, Temp, NewLogRho, NewLogNS, NewVel_GD, NewTemp, &
       LogINS, NewLogINS, IVel, VertVel, NewVertVel)

!!! note that Stage 1 -> updated by a 1/2 step
!!! (NewLogNS - LogNS) = f(tn + Dt/2, yn + dt/2)
       Stage2LogNS(1:nAlts,1:nSpecies)      = &
       Stage1LogNS(1:nAlts,1:nSpecies)      + &
     (NewLogNS(1:nAlts,1:nSpecies) - LogNS(1:nAlts,1:nSpecies))

      Stage2LogINS(1:nAlts,1:nIons)  = &
      Stage1LogINS(1:nAlts,1:nIons)  + &
    (NewLogINS(1:nAlts,1:nIons) - LogINS(1:nAlts,1:nIons))

      Stage2LogRho(1:nAlts)                = &
      Stage1LogRho(1:nAlts)                + &
    (NewLogRho(1:nAlts) - LogRho(1:nAlts))

      Stage2Vel_GD(1:nAlts,1:3)            = & 
      Stage1Vel_GD(1:nAlts,1:3)            + & 
    (NewVel_GD(1:nAlts,1:3) - Vel_GD(1:nAlts,1:3))

        Stage2Temp(1:nAlts)                  = &
        Stage1Temp(1:nAlts)                  + &
     (NewTemp(1:nAlts) -  Temp(1:nAlts))

          Stage2VS(1:nAlts,1:nSpecies)      = &
          Stage1VS(1:nAlts,1:nSpecies)      + &
   (NewVertVel(1:nAlts,1:nSpecies) - VertVel(1:nAlts,1:nSpecies))

  !!! Now Calculate the Next Update Stage
  !!! We need Y(Updated) = Y(n) + 0.5*Stage2
   UpdatedVel_GD(-1:nAlts+2,1:3) = &
        Stage2Vel_GD(-1:nAlts+2,1:3) 

    UpdatedLogNS(-1:nAlts+2,1:nSpecies) = &
         Stage2LogNS(-1:nAlts+2,1:nSpecies)

   UpdatedLogINS(-1:nAlts+2,1:nIons) = &
         Stage2LogINS(-1:nAlts+2,1:nIons) 

   UpdatedLogRho(-1:nAlts+2) = &
         Stage2LogRho(-1:nAlts+2) 

     UpdatedTemp(-1:nAlts+2) = &
          Stage2Temp(-1:nAlts+2) 

       UpdatedVS(-1:nAlts+2,1:nSpecies) = &
            Stage2VS(-1:nAlts+2,1:nSpecies) 

  DtIn = 0.5*Dt
!!! UpdateStage 1 Upper Boundary
  call set_vertical_bcs(UpdatedLogRho, UpdatedLogNS, UpdatedVel_GD, &
                        UpdatedTemp, UpdatedLogINS, IVel, UpdatedVS)

   LogNS  = UpdatedLogNS
   LogINS = UpdatedLogINS
   LogRho = UpdatedLogRho
   Vel_GD = UpdatedVel_GD
     Temp = UpdatedTemp
  VertVel = UpdatedVS

  NewLogNS  = LogNS
  NewLogINS = LogINS
  NewLogRho = LogRho
  NewVel_GD = Vel_GD
     NewTemp = Temp
  NewVertVel = VertVel

!!!!! Calculate K2

  call advance_vertical_1stage_ausm(DtIn, &
       LogRho, LogNS, Vel_GD, Temp, NewLogRho, NewLogNS, NewVel_GD, NewTemp, &
       LogINS, NewLogINS, IVel, VertVel, NewVertVel)

       Stage3LogNS(1:nAlts,1:nSpecies)      = &
       Stage1LogNS(1:nAlts,1:nSpecies)      + &
     (NewLogNS(1:nAlts,1:nSpecies) - LogNS(1:nAlts,1:nSpecies))

      Stage3LogINS(1:nAlts,1:nIons)  = &
      Stage1LogINS(1:nAlts,1:nIons)  + &
    (NewLogINS(1:nAlts,1:nIons) - LogINS(1:nAlts,1:nIons))

      Stage3LogRho(1:nAlts)                = &
      Stage1LogRho(1:nAlts)                + &
    (NewLogRho(1:nAlts) - LogRho(1:nAlts))

      Stage3Vel_GD(1:nAlts,1:3)            = & 
      Stage1Vel_GD(1:nAlts,1:3)            + & 
    (NewVel_GD(1:nAlts,1:3) - Vel_GD(1:nAlts,1:3))

        Stage3Temp(1:nAlts)                  = &
        Stage1Temp(1:nAlts)                  + &
      (NewTemp(1:nAlts) -  Temp(1:nAlts))

          Stage3VS(1:nAlts,1:nSpecies)      = &
          Stage1VS(1:nAlts,1:nSpecies)      + &
   (NewVertVel(1:nAlts,1:nSpecies) - VertVel(1:nAlts,1:nSpecies))

  !!! Now Calculate the Next Update Stage
  !!! We need Y(Updated) = Y(n) + 0.5*Stage3
   UpdatedVel_GD(-1:nAlts+2,1:3) = &
        Stage3Vel_GD(-1:nAlts+2,1:3) 

    UpdatedLogNS(-1:nAlts+2,1:nSpecies) = &
         Stage3LogNS(-1:nAlts+2,1:nSpecies)

   UpdatedLogINS(-1:nAlts+2,1:nIons) = &
         Stage3LogINS(-1:nAlts+2,1:nIons) 

   UpdatedLogRho(-1:nAlts+2) = &
         Stage3LogRho(-1:nAlts+2) 

     UpdatedTemp(-1:nAlts+2) = &
          Stage3Temp(-1:nAlts+2) 

       UpdatedVS(-1:nAlts+2,1:nSpecies) = &
            Stage3VS(-1:nAlts+2,1:nSpecies) 

  DtIn = Dt
!! Update Boundary Conditions

  call set_vertical_bcs(UpdatedLogRho, UpdatedLogNS, UpdatedVel_GD, &
                          UpdatedTemp, UpdatedLogINS, IVel, UpdatedVS)

!
  LogNS  = UpdatedLogNS
  LogINS = UpdatedLogINS
  LogRho = UpdatedLogRho
  Vel_GD = UpdatedVel_GD
  Temp = UpdatedTemp
  VertVel = UpdatedVS

  NewLogNS = LogNS
  NewLogINS = LogINS
  NewLogRho = LogRho
  NewVel_GD = Vel_GD
  NewTemp = Temp
  NewVertVel = VertVel
!
!
!!!!!! Calculate K3

  call advance_vertical_1stage_ausm(DtIn, &
       LogRho, LogNS, Vel_GD, Temp, NewLogRho, NewLogNS, NewVel_GD, NewTemp, &
       LogINS, NewLogINS, IVel, VertVel, NewVertVel)
!!!! K3 Coefficients for RK-4

       Stage4LogNS(1:nAlts,1:nSpecies)      = &
       Stage1LogNS(1:nAlts,1:nSpecies)      + &
     (NewLogNS(1:nAlts,1:nSpecies) - LogNS(1:nAlts,1:nSpecies))

      Stage4LogINS(1:nAlts,1:nIons)  = &
      Stage1LogINS(1:nAlts,1:nIons)  + &
    (NewLogINS(1:nAlts,1:nIons) - LogINS(1:nAlts,1:nIons))

      Stage4LogRho(1:nAlts)                = &
      Stage1LogRho(1:nAlts)                + &
    (NewLogRho(1:nAlts) - LogRho(1:nAlts))

      Stage4Vel_GD(1:nAlts,1:3)            = & 
      Stage1Vel_GD(1:nAlts,1:3)            + & 
    (NewVel_GD(1:nAlts,1:3) - Vel_GD(1:nAlts,1:3))

        Stage4Temp(1:nAlts)                  = &
        Stage1Temp(1:nAlts)                  + &
      (NewTemp(1:nAlts) -  Temp(1:nAlts))

          Stage4VS(1:nAlts,1:nSpecies)      = &
          Stage1VS(1:nAlts,1:nSpecies)      + &
   (NewVertVel(1:nAlts,1:nSpecies) - VertVel(1:nAlts,1:nSpecies))

  !!! Now Calculate the Next Update Stage
  !!! We need Y(Updated) = Y(n) + 0.5*Stage4
   UpdatedVel_GD(-1:nAlts+2,1:3) = &
        Stage4Vel_GD(-1:nAlts+2,1:3) 

    UpdatedLogNS(-1:nAlts+2,1:nSpecies) = &
         Stage4LogNS(-1:nAlts+2,1:nSpecies)

   UpdatedLogINS(-1:nAlts+2,1:nIons) = &
         Stage4LogINS(-1:nAlts+2,1:nIons) 

   UpdatedLogRho(-1:nAlts+2) = &
         Stage4LogRho(-1:nAlts+2) 

     UpdatedTemp(-1:nAlts+2) = &
          Stage4Temp(-1:nAlts+2) 

       UpdatedVS(-1:nAlts+2,1:nSpecies) = &
          Stage4VS(-1:nAlts+2,1:nSpecies) 


!!!! Update Boundary Conditions
  call set_vertical_bcs(UpdatedLogRho, UpdatedLogNS, UpdatedVel_GD, &
                          UpdatedTemp, UpdatedLogINS, IVel, UpdatedVS)

  LogNS  = UpdatedLogNS
  LogINS = UpdatedLogINS
  LogRho = UpdatedLogRho
  Vel_GD = UpdatedVel_GD
  Temp = UpdatedTemp
  VertVel = UpdatedVS

  NewLogNS = LogNS
  NewLogINS = LogINS
  NewLogRho = LogRho
  NewVel_GD = Vel_GD
  NewTemp = Temp
  NewVertVel = VertVel
  
!! Calculate K4 (Final Coefficient)

  DtIn = 0.5*Dt
  call advance_vertical_1stage_ausm(DtIn, &
       LogRho, LogNS, Vel_GD, Temp, NewLogRho, NewLogNS, NewVel_GD, NewTemp, &
       LogINS, NewLogINS, IVel, VertVel, NewVertVel)

  ! This section ensures that our lower boundary conditions are maintained
  ! and not overwritten.
  FinalLogNS = Stage1LogNS
  FinalLogINS = Stage1LogINS
  FinalLogRho = Stage1LogRho
  FinalVel_GD = Stage1Vel_GD
  FinalTemp   = Stage1Temp
  FinalVS     = Stage1VS

!!! Set the Updated State:  Stage 2
  FinalLogNS(1:nAlts,:) = &
     (1.0/3.0)*(-1.0*Stage1LogNS(1:nAlts,:) + Stage2LogNS(1:nAlts,:)  +  &
                 2.0*Stage3LogNS(1:nAlts,:) + Stage4LogNS(1:nAlts,:)  +  &
                       (NewLogNS(1:nAlts,:) - LogNS(1:nAlts,:)) )

  FinalLogINS(1:nAlts,:) = &
     (1.0/3.0)*(-1.0*Stage1LogINS(1:nAlts,:) + Stage2LogINS(1:nAlts,:)  +  &
                 2.0*Stage3LogINS(1:nAlts,:) + Stage4LogINS(1:nAlts,:)  +  &
                       (NewLogINS(1:nAlts,:) - LogINS(1:nAlts,:)) )

  FinalLogRho(1:nAlts) = &
     (1.0/3.0)*(-1.0*Stage1LogRho(1:nAlts) + Stage2LogRho(1:nAlts)  +  &
                 2.0*Stage3LogRho(1:nAlts) + Stage4LogRho(1:nAlts)  +  &
                       (NewLogRho(1:nAlts) - LogRho(1:nAlts)) )

  FinalVel_GD(1:nAlts,:) = &
     (1.0/3.0)*(-1.0*Stage1Vel_GD(1:nAlts,:) + Stage2Vel_GD(1:nAlts,:)  +  &
                 2.0*Stage3Vel_GD(1:nAlts,:) + Stage4Vel_GD(1:nAlts,:)  +  &
                       (NewVel_GD(1:nAlts,:) - Vel_GD(1:nAlts,:)) )

  FinalTemp(1:nAlts) = &
     (1.0/3.0)*(-1.0*Stage1Temp(1:nAlts) + Stage2Temp(1:nAlts)  +  &
                 2.0*Stage3Temp(1:nAlts) + Stage4Temp(1:nAlts)  +  &
                       (NewTemp(1:nAlts) - Temp(1:nAlts)) )

  FinalVS(1:nAlts,:) = &
     (1.0/3.0)*(-1.0*Stage1VS(1:nAlts,:) + Stage2VS(1:nAlts,:)  +  &
                 2.0*Stage3VS(1:nAlts,:) + Stage4VS(1:nAlts,:)  +  &
                  (NewVertVel(1:nAlts,:) - VertVel(1:nAlts,:)) )

  call set_vertical_bcs(FinalLogRho, FinalLogNS, FinalVel_GD, &
                          FinalTemp, FinalLogINS, IVel, FinalVS)

   LogNS = FinalLogNS
  LogINS = FinalLogINS
  LogRho = FinalLogRho
  Vel_GD = FinalVel_GD
    Temp = FinalTemp
 VertVel = FinalVS

  if (UseBarriers) call MPI_BARRIER(iCommGITM,iError)
  if (iDebugLevel > 7) &
       write(*,*) "========> Done with advance_vertical_1d", iproc

end subroutine advance_vertical_1d_ausm

!=============================================================================
subroutine advance_vertical_1stage_ausm( DtIn, &
     LogRho, LogNS, Vel_GD, Temp, NewLogRho, NewLogNS, NewVel_GD, NewTemp, &
     LogINS, NewLogINS, IVel, VertVel, NewVertVel)

  ! With fluxes and sources based on LogRho..Temp, update NewLogRho..NewTemp

  use ModGITM, only: &
       Dt, iEast_, iNorth_, iUp_, ThermalDiffCoefS
  use ModPlanet
  use ModSizeGitm
  use ModVertical, only : &
       Heating, EddyCoef_1D, Coriolis, &
       MeanMajorMass_1d, Gamma_1d, InvRadialDistance_C, &
       !ChemSources_1d, &
       Centrifugal, &
       CD_MeshCoefs, &
       Gravity_G, Altitude_G,Cv_1D, dAlt_F
  use ModTime
  use ModInputs
  use ModConstants
  use ModSources, only : EddyCondAdia
  implicit none

  real, intent(in) :: DtIn
  real, intent(in) :: LogRho(-1:nAlts+2)
  real, intent(in) :: LogNS(-1:nAlts+2,nSpecies)
  real, intent(in) :: LogINS(-1:nAlts+2,nIons)
  real, intent(in) :: Vel_GD(-1:nAlts+2,3)
  real, intent(in) :: IVel(-1:nAlts+2,3)
  real, intent(in) :: Temp(-1:nAlts+2)
  real, intent(in) :: VertVel(-1:nAlts+2,nSpecies)

  real, intent(inout) :: NewLogRho(-1:nAlts+2)
  real, intent(inout) :: NewLogNS(-1:nAlts+2,nSpecies)
  real, intent(inout) :: NewLogINS(-1:nAlts+2,nIons)
  real, intent(inout) :: NewVel_GD(-1:nAlts+2,3)
  real :: NewVel2_G(-1:nAlts+2)
  real, intent(inout) :: NewTemp(-1:nAlts+2)
  real, intent(out) :: NewVertVel(-1:nAlts+2,nSpecies)
  real :: NS(-1:nAlts+2,nSpecies), Pressure1D(-1:nAlts+2)
  real :: Rho(-1:nAlts+2)

  real :: LogNum(-1:nAlts+2)

  real, dimension(1:nAlts)    :: GradLogRho, DivVel, GradTemp, GradTempKoM, &
       DiffLogRho, DiffTemp, GradTmp, DiffTmp, DiffLogNum, GradLogNum, &
       DiviVel
  real, dimension(1:nAlts,3) :: GradVel_CD, DiffVel_CD
  real, dimension(1:nAlts,3) :: GradiVel_CD, DiffiVel_CD

  real, dimension(1:nAlts,nSpecies)    :: GradLogNS, DiffLogNS, &
       GradVertVel, DiffVertVel, DivVertVel
  real, dimension(1:nAlts,nIons) :: GradLogINS, DiffLogINS
  real :: NewSumRho, NewLogSumRho, rat, ed

  integer :: iAlt, iSpecies, jSpecies, iDim

  real, dimension(-1:nAlts+2)    :: NT
  real, dimension(-1:nAlts+2)    :: Press, LogPress
  real, dimension(1:nAlts)    :: DiffLogPress, GradLogPress
  real, dimension(1:nAlts,nSpecies)    :: EddyDiffusionVel

  real :: nVel(1:nAlts,1:nSpecies)
  integer :: nFilter, iFilter
  real :: LowFilter

! Parameters Used for the Sponge
! This Sponge is useful to dampen out spurious modes
! oscillating between the bottom and top of the model.

  integer :: nAltsSponge = 12
  real :: kSP, NuSP, AmpSP

  ! JMB:  Adding Eddy Diffusion Variables here
  ! Note:  These are used in the calc_neutral_friction
  !--------------------------------------------------------------------------
  !! Eddy Diffusion Variables
  real, dimension(1:nAlts,nSpecies)    :: GradLogConS
  real, dimension(-1:nAlts+2,nSpecies)    :: ConS, LogConS
  real, dimension(1:nAlts,nSpecies)    :: EddyCoefRatio_1d
  !--------------------------------------------------------------------------
  ! 4th Order Gradients on a Non-Uniform Mesh (5-point Stencil)
  ! Used for calculating the d(ln[Chi])/dr -> Log of the concentration gradient
  !--------------------------------------------------------------------------
  real :: h1, h2, h3, h4
  real :: MeshH1, MeshH2, MeshH3, MeshH4
  real :: MeshCoef0, MeshCoef1, &
          MeshCoef2, MeshCoef3, &
          MeshCoef4
  ! ----------------------------------------------------
  ! JMB:  AUSM Variables
  ! ----------------------------------------------------
  real ::    RhoS(-1:nAlts+2,1:nSpecies),&
          NewRhoS(-1:nAlts+2,1:nSpecies),&
      AUSMRhoSFluxes(1:nAlts,1:nSpecies)

  real :: InvScaleHeight, MeanGravity, MeanTemp, MeanMass

  real ::   HydroNS(-1:nAlts+2,1:nSpecies),&
          HydroRhoS(-1:nAlts+2,1:nSpecies), &
            HydroNT(-1:nAlts+2),&
     HydroPressureS(-1:nAlts+2,1:nSpecies),&
      DeviationRhoS(-1:nAlts+2,1:nSpecies),&
 DeviationRhoSRatio(-1:nAlts+2,1:nSpecies),&
      HydroPressure(-1:nAlts+2), &
           HydroRho(-1:nAlts+2)

  ! Momentum Fluxes and Variables
  real ::    MomentumS(-1:nAlts+2,1:nSpecies),&
          NewMomentumS(-1:nAlts+2,1:nSpecies),&
              Momentum(-1:nAlts+2,1:3),&        ! Bulk Momentum
           NewMomentum(-1:nAlts+2,1:3),&        ! Bulk Momentum
       AUSMMomentumSFluxes(1:nAlts,1:nSpecies), &
       AUSMMomentumFluxes(1:nAlts,3)

  real :: PressureS(-1:nAlts+2,1:nSpecies), &
          NewNS(-1:nAlts+2,1:nSpecies), &
          NewNT(-1:nAlts+2)

  real :: NS_small(nAlts,nSpecies)
  
  real :: TotalEnergy(-1:nAlts+2),&
       NewTotalEnergy(-1:nAlts+2),&
       AUSMTotalEnergyFluxes(1:nAlts), &
       NewPress(-1:nAlts+2), NewRho(-1:nAlts+2)

  real :: RadialDistance_C(-1:nAlts+2)
  real :: EffectiveGravity(-1:nAlts+2)
  real :: EffectiveGravity_Energy(-1:nAlts+2)

  ! JMB:  Use these as Limiters on Winds for an initial startup
  real :: TimeFactor, Vel0, DeltaV, VelocityCap
  real, dimension(-1:nAlts+2,nSpecies) :: DijS, LambdaS
  real :: InvDij, TempDij, denscale, kTOverM
  real :: InvScaleHeightAtm, ExpArg

  !--------------------------------------------------------------------------
  !--------------------------------------------------------------------------
  Vel0 = 1.0 ! initial velocity in m/s
  TimeFactor = exp(-tSimulation/1800.0)
  DeltaV = (MaximumVerticalVelocity - Vel0)*(1.0 - TimeFactor)
  VelocityCap = Vel0 + DeltaV*(1.0 - TimeFactor)

  do iAlt = -1, nAlts + 2
     EffectiveGravity(iAlt) = &
        Gravity_G(iAlt) + &
        Centrifugal / InvRadialDistance_C(iAlt) + & 
        (Vel_GD(iAlt,iNorth_)**2 + Vel_GD(iAlt,iEast_)**2) &
        * InvRadialDistance_C(iAlt) + & 
        Coriolis * Vel_GD(iAlt,iEast_)

     ! Use this version in energy conservation.
     ! Curvature and Coriolis forces disappear when calculating energy
     ! conservation.
     EffectiveGravity_Energy(iAlt) = &
        Gravity_G(iAlt) + &
        Centrifugal / InvRadialDistance_C(iAlt) 
  enddo 

  NS = exp(LogNS)

  Rho = exp(LogRho)
  LogNum = alog(sum(NS,dim=2))
  nFilter = 10
  do iAlt = -1, nAlts + 2
     RadialDistance_C(iAlt) = 1.0/InvRadialDistance_C(iAlt)
  enddo 
  
  NT(-1:nAlts+2) = 0.0
  do iSpecies = 1, nSpecies
     NT(-1:nAlts+2) = NT(-1:nAlts+2) + &
                      NS(-1:nAlts+2,iSpecies)
  enddo 

  do iAlt = -1, nAlts + 2
    Press(iAlt) = NT(iAlt)*Boltzmanns_Constant*Temp(iAlt)
    LogPress(iAlt) = alog(Press(iAlt))
  enddo

  do iAlt = -1, nAlts+2
     do iSpecies = 1, nSpecies
        InvDij = 0.0
        kTOverM = Boltzmanns_Constant * Temp(iAlt) / Mass(iSpecies)
        denscale = 1.0/NT(iAlt) 
        do jSpecies = 1, nSpecies
           if (jSpecies == iSpecies) cycle
              TempDij = (1.0e-04)*&              ! Scales the Dij from cm^2/s -> m^2/s
                (   Diff0(iSpecies,jSpecies)*( Temp(iAlt)**DiffExp(iSpecies,jSpecies) )   ) / &
                (    NT(iAlt)*(1.0e-06) )     ! Converts to #/cm^-3
           InvDij = InvDij + &
                denscale*NS(iAlt, jSpecies)/ &
                ( TempDij )
        enddo  ! End DO over jSpecies
        DijS(iAlt,iSpecies) = 1.0/InvDij
        LambdaS(iAlt,iSpecies) = EddyCoef_1d(iAlt)/DijS(iAlt,iSpecies)
     enddo  !End DO Over iSpecies
  enddo !iAlt = 1, nAlts


  call calc_rusanov_alts_ausm(LogPress ,GradLogPress,  DiffLogPress)
  call calc_rusanov_alts_ausm(LogRho ,GradLogRho,  DiffLogRho)
  call calc_rusanov_alts_ausm(LogNum ,GradLogNum,  DiffLogNum)
  call calc_rusanov_alts_ausm(Temp   ,GradTemp,    DiffTemp)
  do iDim = 1, 3
     call calc_rusanov_alts_ausm(Vel_GD(:,iDim), &
          GradVel_CD(:,iDim),DiffVel_CD(:,iDim))
     ! call calc_rusanov_alts_ausm(iVel(:,iDim), &
     call calc_rusanov_alts_rusanov(iVel(:,iDim), &
          GradiVel_CD(:,iDim),DiffiVel_CD(:,iDim))
  enddo

  ! Add geometrical correction to gradient and obtain divergence
  DivVel = GradVel_CD(:,iUp_) + 2*Vel_GD(1:nAlts,iUp_)*InvRadialDistance_C(1:nAlts)
  !DiviVel = DivIonVelCoef * (GradiVel_CD(:,iUp_) + 2*iVel(1:nAlts,iUp_)*InvRadialDistance_C(1:nAlts))

  do iSpecies=1,nSpecies

     call calc_rusanov_alts_ausm(LogNS(:,iSpecies),GradTmp, DiffTmp)
     GradLogNS(:,iSpecies) = GradTmp
     DiffLogNS(:,iSpecies) = DiffTmp

     call calc_rusanov_alts_ausm(VertVel(:,iSpecies),GradTmp, DiffTmp)
     GradVertVel(:,iSpecies) = GradTmp
     DiffVertVel(:,iSpecies) = DiffTmp
     DivVertVel(:,iSpecies) = GradVertVel(:,iSpecies) + &
          2*VertVel(1:nAlts,iSpecies)*InvRadialDistance_C(1:nAlts)

  enddo

  do iSpecies=1,nIons-1
     ! call calc_rusanov_alts_ausm(LogINS(:,iSpecies), GradTmp, DiffTmp)
     call calc_rusanov_alts_rusanov(LogINS(:,iSpecies), GradTmp, DiffTmp)
     GradLogINS(:,iSpecies) = GradTmp
     DiffLogINS(:,iSpecies) = DiffTmp
  enddo

  ! JMB 2022 Update: Use CD Coefficients
  do iSpecies = 1, nSpecies
    LogConS(-1:nAlts+2,iSpecies) = &
         alog(Mass(iSpecies)*NS(-1:nAlts+2,iSpecies)/Rho(-1:nAlts+2))
  enddo 

  do iAlt = 1, nAlts
   do iSpecies = 1, nSpecies
         GradLogConS(iAlt,iSpecies) =  &
            CD_MeshCoefs(iAlt,1)*LogConS(iAlt-2,iSpecies)&
         +  CD_MeshCoefs(iAlt,2)*LogConS(iAlt-1,iSpecies)&
         +  CD_MeshCoefs(iAlt,3)*LogConS(iAlt  ,iSpecies)&
         +  CD_MeshCoefs(iAlt,4)*LogConS(iAlt+1,iSpecies)&
         +  CD_MeshCoefs(iAlt,5)*LogConS(iAlt+2,iSpecies)

      enddo  ! iSpecies Loop
  enddo ! iAlt Loop

  !!!! JMB AUSM: BEGIN THE HYDROSTATIC BACKGROUND
  !!!! We define a hydrostatic background to subtract off
  !!!! this removes errors introduced in the Grad(P) - rho*g
  !!!! that reduces the accuracy of non-hydrostatic calculations.
  HydroNS(-1:nAlts+2,1:nSpecies) = NS(-1:nAlts+2,1:nSpecies)
  HydroNT(-1:nAlts+2)            = NT(-1:nAlts+2)
  do iSpecies =  1, nSpecies
     do iAlt = 1, nAlts + 2
     MeanMass = Mass(iSpecies)
     MeanGravity = 0.5*( EffectiveGravity(iAlt) + EffectiveGravity(iAlt-1) )
     MeanTemp = 0.5*( Temp(iAlt) + Temp(iAlt-1) )
     InvScaleHeight = -1.0* MeanMass*MeanGravity/&
                      (Boltzmanns_Constant*MeanTemp)

     MeanMass = 0.5*(MeanMajorMass_1d(iAlt-1) + MeanMajorMass_1d(iAlt))
     InvScaleHeightAtm = -1.0* MeanMass*MeanGravity/&
                          (Boltzmanns_Constant*MeanTemp)
     ExpArg = ( LambdaS(iAlt,iSpecies)*InvScaleHeightAtm + &
                 InvScaleHeight)/&
              ( 1.0 + LambdaS(iAlt,iSpecies))

     HydroNS(iAlt,iSpecies) = &
           HydroNS(iAlt-1,iSpecies)*(Temp(iAlt-1)/Temp(iAlt))*&
           exp (-1.0*dAlt_F(iAlt)*ExpArg)

     enddo 
  enddo 
!--------------
  do iAlt = -1, nAlts + 2
     do iSpecies =  1, nSpecies
       HydroPressureS(iAlt,iSpecies) = &
           HydroNS(iAlt,iSpecies)*Boltzmanns_Constant*Temp(iAlt)
       HydroRhoS(iAlt,iSpecies) = Mass(iSpecies)*HydroNS(iAlt,iSpecies)
     enddo 
  enddo 
!--------------
  do iAlt = 1, nAlts + 2
    HydroNT(iAlt) = 0.0
    HydroRho(iAlt) = 0.0
    do iSpecies = 1, nSpecies
       HydroNT(iAlt) = HydroNT(iAlt) + &
        HydroNS(iAlt,iSpecies)

       HydroRho(iAlt) = HydroRho(iAlt) + &
        HydroRhoS(iAlt,iSpecies)
    enddo 
  enddo 
  do iSpecies = 1, nSpecies
     RhoS(-1:nAlts+2,iSpecies) =  &
        Mass(iSpecies)*NS(-1:nAlts+2,iSpecies)
     NewRhoS(-1:nAlts+2,iSpecies) = RhoS(-1:nAlts+2,iSpecies)
     MomentumS(-1:nAlts+2,iSpecies) =  &
        Mass(iSpecies)*NS(-1:nAlts+2,iSpecies)*&
         VertVel(-1:nAlts+2,iSpecies)
     PressureS(-1:nAlts+2,iSpecies) =  &
        NS(-1:nAlts+2,iSpecies)*Temp(-1:nAlts+2)*&
        Boltzmanns_Constant
  enddo 

  do iAlt = -1, nAlts + 2
   TotalEnergy(iAlt) = &
        Press(iAlt)/(Gamma_1d(iAlt) - 1.0) + &
        0.5*Rho(iAlt)*(Vel_GD(iAlt,iUp_)**2.0 + &
                       Vel_GD(iAlt,iNorth_)**2.0 + &
                       Vel_GD(iAlt,iEast_)**2.0) 
  enddo 
  do iDim = 1, 3
     Momentum(-1:nAlts+2,iDim) = Rho(-1:nAlts+2)*Vel_GD(-1:nAlts+2,iDim)
  enddo 

  DeviationRhoSRatio(-1:nAlts+2,1:nSpecies) = & 
            abs(RhoS(-1:nAlts+2,1:nSpecies) - &
           HydroRhoS(-1:nAlts+2,1:nSpecies))/&
                RhoS(-1:nAlts+2,1:nSpecies)

  DeviationRhoS(-1:nAlts+2,1:nSpecies) = &
           RhoS(-1:nAlts+2, 1:nSpecies) - &
      HydroRhoS(-1:nAlts+2,1:nSpecies)

  NewRho = Rho
  NewPress = Press
  NewTotalEnergy = TotalEnergy
  ! Call the AUSM Solvers
  call calc_all_fluxes_hydro(DtIn, RhoS, PressureS, HydroPressureS, HydroRhoS, &
        HydroPressure,  HydroRho, AUSMRhoSFluxes,AUSMMomentumSFluxes, &
        AUSMTotalEnergyFluxes, AUSMMomentumFluxes)

  AmpSP = (1.0/(10.0*DtIn))
  kSP = nAltsSponge + 1

  do iAlt = 1,nAlts
     NewLogRho(iAlt) = NewLogRho(iAlt) - DtIn * &
          (DivVel(iAlt) + Vel_GD(iAlt,iUp_) * GradLogRho(iAlt) ) &
          + DtIn * DiffLogRho(iAlt)

     do iSpecies=1,nSpecies
        NewRhoS(iAlt,iSpecies) = RhoS(iAlt,iSpecies) &
               - DtIn*(AUSMRhoSFluxes(iAlt,iSpecies))
        NewLogNS(iAlt,iSpecies) = alog( NewRhoS(iAlt,iSpecies)/Mass(iSpecies) )
     enddo

     do iSpecies=1,nIonsAdvect
!        if (UseImprovedIonAdvection) then
!           ! This works for non-log ion densities:
!           NewLogINS(iAlt,iSpecies) = NewLogINS(iAlt,iSpecies) - DtIn * &
!                (DiviVel(iAlt) * LogINS(iAlt,iSpecies) + &
!                IVel(iAlt,iUp_) * GradLogINS(iAlt,iSpecies) ) &
!                + DtIn * DiffLogINS(iAlt,iSpecies)
!        else
           NewLogINS(iAlt,iSpecies) = NewLogINS(iAlt,iSpecies) - DtIn * &
                (IVel(iAlt,iUp_) * GradLogINS(iAlt,iSpecies) ) &
                + DtIn * DiffLogINS(iAlt,iSpecies)
!        endif
     enddo

  enddo !iAlt = 1,nAlts


  NewNS  = 0.0
  NewNT  = 0.0
  NewRho = 0.0

  do iAlt = -1, nAlts+2
     do iSpecies = 1, nSpecies
         NewNS(iAlt,iSpecies) = NewRhoS(iAlt,iSpecies)/Mass(iSpecies)

         NewRho(iAlt) = NewRho(iAlt) + &
                  NewRhoS(iAlt,iSpecies)

         NewNT(iAlt) = NewNT(iAlt) + &
               NewNS(iAlt,iSpecies)
     enddo 
  enddo 


  if (iAlt >= (nAlts - nAltsSponge)) then
     NuSP = AmpSP*(1.0 - cos( pi*(kSP - (nAlts - iAlt))/kSP) )
  else
     NuSP = 0.0
  endif

  if (UseDamping) then
     VertTau(iAlt) = &
          15 - (1 - exp(-1.0*altitude_G(ialt)/1000.0/40.0))*5.0
  endif

  do iAlt = 1,nAlts
     do iSpecies=1,nSpecies

        NewMomentumS(iAlt,iSpecies) = MomentumS(iAlt,iSpecies) - &
              DtIn*(AUSMMomentumSFluxes(iAlt,iSpecies)) + &
              DtIn*DeviationRhoS(iAlt,iSpecies)*EffectiveGravity(iAlt) 

!        NewMomentumS(iAlt,iSpecies) = NewMomentumS(iAlt,iSpecies) + &
!              DtIn*ChemSources_1d(iAlt,iSpecies)*&
!              VertVel(iAlt,iSpecies)*Mass(iSpecies)
!------------------
        NewVertVel(iAlt,iSpecies) = &
           NewMomentumS(iAlt,iSpecies)/NewRhoS(iAlt,iSpecies) 
        ! Correction for making the background atmosphere use the full
        ! diffusive equilibrium formulation
        MeanGravity = 0.5*(EffectiveGravity(iAlt-1) + &
                           EffectiveGravity(iAlt  ) )
        MeanMass = 0.5*(MeanMajorMass_1d(iAlt-1) + &
                        MeanMajorMass_1d(iAlt  ) )
      
        NewVertVel(iAlt,iSpecies) = NewVertVel(iAlt,iSpecies) + &
             DtIn*(HydroRhoS(iAlt,iSpecies)/RhoS(iAlt,iSpecies))*MeanGravity *&
             (LambdaS(iAlt,iSpecies)/(1.0 + LambdaS(iAlt,iSpecies)))*&
               ( 1.0 - MeanMass/Mass(iSpecies)) 

        ! Thermal Diffusion Effects (For Light Species H2, H, and He) 
        ! ThermalDiffCoefS is set in calc_rates
        ! Note:  ThermalDiffCoefS is < 0.0 for light species
        ! This forces light species into hot zones and heavy species into cold zones
        NewVertVel(iAlt,iSpecies) = NewVertVel(iAlt,iSpecies) - &
          DtIn*(ThermalDiffCoefS(iSpecies)*Boltzmanns_Constant*&
             GradTemp(iAlt))/Mass(iSpecies)
     enddo ! iSpecies

  enddo ! iAlt

  if (UseNeutralFriction) then
     nVel(1:nAlts,1:nSpecies) = NewVertVel(1:nAlts,1:nSpecies)
     NS_small = NewNS(1:nAlts,1:nSpecies)

     call calc_neutral_friction(DtIn,nVel(1:nAlts,1:nSpecies), &
                         EddyCoef_1d(1:nAlts), &
                               NewNT(1:nAlts), &
                               NS_small, &
                         GradLogConS(1:nAlts,1:nSpecies), &
                                Temp(1:nAlts))
     NewVertVel(1:nAlts,1:nSpecies) = nVel(1:nAlts,1:nSpecies)
  endif 


  do iAlt = -1, nAlts+2
     do iSpecies=1,nSpecies
        NewVertVel(iAlt, iSpecies) = max(-MaximumVerticalVelocity, &
             NewVertVel(iAlt, iSpecies))
        NewVertVel(iAlt, iSpecies) = min( MaximumVerticalVelocity, &
             NewVertVel(iAlt, iSpecies))
     enddo
  enddo

  NewVel_GD(-1:nAlts+2,iUp_) = 0.0
  do iAlt = -1, nAlts+2
     do iSpecies=1,nSpecies
!        NewVertVel(iAlt, iSpecies) = max(-MaximumVerticalVelocity, &
!             NewVertVel(iAlt, iSpecies))
!        NewVertVel(iAlt, iSpecies) = min( MaximumVerticalVelocity, &
!             NewVertVel(iAlt, iSpecies))
!
!        NewVertVel(iAlt, iSpecies) = max(-VelocityCap, &
!             NewVertVel(iAlt, iSpecies))
!        NewVertVel(iAlt, iSpecies) = min( VelocityCap, &
!             NewVertVel(iAlt, iSpecies))
!
        NewVel_GD(iAlt,iUp_) = NewVel_GD(iAlt,iUp_) + &
             NewVertVel(iAlt, iSpecies)*NewRhoS(iAlt,iSpecies)/NewRho(iAlt)
     enddo
  enddo


  do iAlt = 1, nAlts

      NewMomentum(iAlt,iEast_) = Momentum(iAlt,iEast_) - &
           DtIn*(AUSMMomentumFluxes(iAlt,iEast_)) 
 
      NewVel_GD(iAlt,iEast_) = NewMomentum(iAlt,iEast_)/NewRho(iAlt)
      NewVel_GD(iAlt,iEast_) = NewVel_GD(iAlt,iEast_) &
          + 0.35*exp(-1.0*(Altitude_G(iAlt) - Altitude_G(0))/&
            (45.0e+03)) * DtIn*DiffVel_CD(iAlt,iEast_)
 
      NewMomentum(iAlt,iNorth_) = Momentum(iAlt,iNorth_) - &
           DtIn*(AUSMMomentumFluxes(iAlt,iNorth_)) 
 
      NewVel_GD(iAlt,iNorth_) = NewMomentum(iAlt,iNorth_)/NewRho(iAlt)
      NewVel_GD(iAlt,iNorth_) = NewVel_GD(iAlt,iNorth_) &
          + 0.35*exp(-1.0*(Altitude_G(iAlt) - Altitude_G(0))/&
            (45.0e+03)) * DtIn*DiffVel_CD(iAlt,iNorth_)

      ! dVphi/dt = - (V grad V)_phi
!       NewVel_GD(iAlt,iEast_) = Vel_GD(iAlt,iEast_) - DtIn * &
!            Vel_GD(iAlt,iUp_)*GradVel_CD(iAlt,iEast_) &
!            + DtIn * DiffVel_CD(iAlt,iEast_)
! 
!       ! dVtheta/dt = - (V grad V)_theta
!       NewVel_GD(iAlt,iNorth_) = Vel_GD(iAlt,iNorth_) - DtIn * &
!            Vel_GD(iAlt,iUp_)*GradVel_CD(iAlt,iNorth_) &
!            + DtIn * DiffVel_CD(iAlt,iNorth_)
 
     ! dT/dt = -(V.grad T + (gamma - 1) T div V +  &
     !        (gamma - 1) * g  * grad (KeH^2  * rho) /rho 
     ! AUSM Method
     NewTotalEnergy(iAlt)   = TotalEnergy(iAlt) - &
         DtIn*AUSMTotalEnergyFluxes(iAlt) + &
         DtIn*Rho(iAlt)*Vel_GD(iAlt,iUp_)*EffectiveGravity_Energy(iAlt) 
 
     NewPress(iAlt) = &
        (NewTotalEnergy(iAlt) - &
         0.5*NewRho(iAlt)*(NewVel_GD(iAlt,iUp_)**2.0 + &
                           NewVel_GD(iAlt,iNorth_)**2.0 + &
                           NewVel_GD(iAlt,iEast_ )**2.0 ) )*&
                           (Gamma_1d(iAlt) - 1.0)
 
     NewTemp(iAlt) = NewPress(iAlt)/(Boltzmanns_Constant*NewNT(iAlt))

     NewTemp(iAlt) = NewTemp(iAlt) + &
       0.25*exp(-1.0*(Altitude_G(iAlt) - Altitude_G(0))/&
       (10.0e+03)) * DtIn*DiffTemp(iAlt)

  enddo ! iAlt

  do iAlt = 1, nAlts
     NewSumRho    = sum( Mass(1:nSpecies)*exp(NewLogNS(iAlt,1:nSpecies)) )
     NewLogRho(iAlt) = alog(NewSumRho)
  enddo

end subroutine advance_vertical_1stage_ausm

!\
! ------------------------------------------------------------
! calc_rusanov
! ------------------------------------------------------------
!/

subroutine calc_rusanov_alts_ausm(Var, GradVar, DiffVar)

  use ModSizeGitm
  use ModVertical, only : dAlt_C, cMax
  implicit none

  real, intent(in) :: Var(-1:nAlts+2)
  real, intent(out):: GradVar(1:nAlts), DiffVar(1:nAlts)

  real, dimension(1:nAlts+1) :: VarLeft, VarRight, DiffFlux
  !------------------------------------------------------------

  call calc_facevalues_alts_ausm(Var, VarLeft, VarRight)

  ! Gradient based on averaged Left/Right values
  GradVar = 0.5 * &
       (VarLeft(2:nAlts+1)+VarRight(2:nAlts+1) - &
       VarLeft(1:nAlts)-VarRight(1:nAlts))/dAlt_C(1:nAlts)

  ! Rusanov/Lax-Friedrichs diffusive term
  DiffFlux = 0.5 * max(cMax(0:nAlts),cMax(1:nAlts+1)) * (VarRight - VarLeft)

  DiffVar = (DiffFlux(2:nAlts+1) - DiffFlux(1:nAlts))/dAlt_C(1:nAlts)

end subroutine calc_rusanov_alts_ausm

!\
! ------------------------------------------------------------
! calc_facevalues_alts_ausm
! ------------------------------------------------------------
!/

subroutine calc_facevalues_alts_ausm(Var, VarLeft, VarRight)

  use ModVertical, only: dAlt_F, InvDAlt_F
  use ModSizeGITM, only: nAlts
  use ModLimiterGitm

  implicit none
  
  real, intent(in) :: Var(-1:nAlts+2)
  real, intent(out):: VarLeft(1:nAlts+1), VarRight(1:nAlts+1)

  real :: dVarUp, dVarDown, dVarLimited(0:nAlts+1)

  real, parameter :: Factor1=0.6250000 ! 15/24
  real, parameter :: Factor2=0.0416667 !  1/24
  real :: h

  integer :: i

  do i=1,nAlts

     ! 4th order scheme for calculating face values

     h  = InvDAlt_F(i+1)*2.0
     dVarUp   = h*(Factor1*(Var(i+1)-Var(i)  ) - Factor2*(Var(i+2)-Var(i-1)))
     h  = InvDAlt_F(i)*2.0
     dVarDown = h*(Factor1*(Var(i)  -Var(i-1)) - Factor2*(Var(i+1)-Var(i-2)))

!     ! This is Gabor's scheme
!     dVarUp            = (Var(i+1) - Var(i))   * InvDAlt_F(i+1)
!     dVarDown          = (Var(i)   - Var(i-1)) * InvDAlt_F(i)

     dVarLimited(i) = Limiter_mc(dVarUp, dVarDown)

  end do

  i = 0
  dVarUp            = (Var(i+1) - Var(i))   * InvDAlt_F(i+1)
  dVarDown          = (Var(i)   - Var(i-1)) * InvDAlt_F(i)
  dVarLimited(i) = Limiter_mc(dVarUp, dVarDown)

  i = nAlts+1
  dVarUp            = (Var(i+1) - Var(i))   * InvDAlt_F(i+1)
  dVarDown          = (Var(i)   - Var(i-1)) * InvDAlt_F(i)
  dVarLimited(i) = Limiter_mc(dVarUp, dVarDown)

  do i=1,nAlts+1
     VarLeft(i)  = Var(i-1) + 0.5*dVarLimited(i-1) * dAlt_F(i)
     VarRight(i) = Var(i)   - 0.5*dVarLimited(i)   * dAlt_F(i)
  end do

end subroutine calc_facevalues_alts_ausm


subroutine calc_all_fluxes_hydro(DtIn, RhoS, PressureS, HydroPressureS, HydroRhoS, &
            HydroPressure, HydroRho, DivRhoSFlux, DivMomentumSFlux, &
            DivEnergyFlux, DivMomentumFlux)

  use ModSizeGitm
  use ModVertical, only : dAlt_C, cMax, VertVel, Gamma_1D, &
                          LogRho, Vel_GD, MeanMajorMass_1d, &
                          CellVol1D, Area_P12, Area_M12, &
                          Temp, Altitude_G, dAlt_F


  use ModPlanet, only : nSpecies, nIonsAdvect, Mass, RBody, &
                        iN2_, cSpecies
  use ModGITM, only : Dt,iUp_, iEast_, iNorth_
  use ModConstants, only : Boltzmanns_Constant

  implicit none

! Passed from Vertical_Solver In
  real, intent(in) :: DtIn
  real, intent(in) ::           RhoS(-1:nAlts+2, 1:nSpecies)
  real, intent(in) ::      PressureS(-1:nAlts+2,1:nSpecies)
  real, intent(in) ::      HydroRhoS(-1:nAlts+2, 1:nSpecies)
  real, intent(in) :: HydroPressureS(-1:nAlts+2,1:nSpecies)
  real, intent(in) ::       HydroRho(-1:nAlts+2) 
  real, intent(in) ::  HydroPressure(-1:nAlts+2)
  real, intent(out)::      DivRhoSFlux(1:nAlts,1:nSpecies)
  real, intent(out):: DivMomentumSFlux(1:nAlts,1:nSpecies)
  real, intent(out)::    DivEnergyFlux(1:nAlts)
  real, intent(out)::  DivMomentumFlux(1:nAlts,1:3)

  real, dimension(-1:nAlts+2, 1:nSpecies) :: LogRhoS
  real, dimension(-1:nAlts+2, 1:nSpecies) :: LogPS

  real, dimension(-1:nAlts+2,1:nIonsAdvect) :: LogRhoI
  real, dimension(1:nAlts,3) :: IVelLeft_M12, IVelRight_M12
  real, dimension(1:nAlts,3) :: IVelLeft_P12, IVelRight_P12
  real, dimension(1:nAlts,1:nIonsAdvect) :: RhoIFlux_M12, RhoIFlux_P12

!!! Hydrostatic Variables
  real, dimension(-1:nAlts+2,1:nSpecies) :: LogHydroPressureS
  real, dimension(-1:nAlts+2,1:nSpecies) :: LogHydroRhoS

  real :: SubCs
  integer :: iSpecies, iAlt, iDim
  !------------------------------------------------------------

  ! ==================== AUSM Flux Variables
  real, dimension( 1:nAlts) :: BulkIVel_P12, BulkIVel_M12    
  real :: Kp(1:nSpecies), Ku(1:nSpecies)

  real, dimension( 1:nAlts) :: AreaFunction_P12, AreaFunction_M12    
  real, dimension( 1:nAlts) :: LocalCellVolume

  real:: MInf, LiouBeta, LiouAlpha


  ! New Streamlined Variables
  real, dimension(0:nAlts,1:nSpecies) ::         RhoSLeft,         RhoSRight
  real, dimension(0:nAlts,1:nSpecies) ::      LogRhoSLeft,      LogRhoSRight
  real, dimension(0:nAlts,1:nSpecies) ::    HydroRhoSLeft,    HydroRhoSRight
  real, dimension(0:nAlts,1:nSpecies) :: LogHydroRhoSLeft, LogHydroRhoSRight
  real, dimension(0:nAlts,1:nSpecies) ::     MeanRhoS
  ! New Pressure Variables
  real, dimension(0:nAlts,1:nSpecies) :: PressureSLeft, PressureSRight
  real, dimension(0:nAlts,1:nSpecies) :: HydroPressureSLeft, HydroPressureSRight
  real, dimension(0:nAlts,1:nSpecies) :: MeanPressureS
  real, dimension(0:nAlts,1:nSpecies) :: MeanHydroPressureS
  ! New Species Winds 
  real, dimension(0:nAlts,1:nSpecies) ::     VertVelSLeft,     VertVelSRight
  ! Bulk Fields
  real, dimension(0:nAlts)            :: RhoLeft, RhoRight
  real, dimension(0:nAlts)            :: TemperatureLeft, TemperatureRight
  real, dimension(0:nAlts)            :: PressureLeft, PressureRight
  real, dimension(0:nAlts,1:3)        :: VelGDLeft, VelGDRight
  real, dimension(0:nAlts)            :: GammaLeft, GammaRight
  ! Derived Fields
  real, dimension(0:nAlts)            :: EnergyLeft, EnergyRight
  real, dimension(0:nAlts)            :: EnthalpyLeft, EnthalpyRight
  real, dimension(0:nAlts,1:nSpecies) :: EnthalpySLeft, EnthalpySRight

  ! Sound Speed Fields
  real, dimension(0:nAlts)            :: SoundSpeedLeft, SoundSpeedRight
  real, dimension(0:nAlts)            :: InterfaceSoundSpeed
  real, dimension(0:nAlts,1:nSpecies) :: SoundSpeedSLeft, SoundSpeedSRight
  real, dimension(0:nAlts,1:nSpecies) :: InterfaceSoundSpeedS

  ! Mach Numbers 
  real, dimension(0:nAlts,1:nSpecies) :: MachSLeft, MachSRight
  real, dimension(0:nAlts,1:nSpecies) :: MachSBar2
  real, dimension(0:nAlts,1:nSpecies) :: MachSZero2, MachSZero
  ! AUSM+-up Parameters
  real, dimension(0:nAlts,1:nSpecies) :: KpSParam, KuSParam
  
  ! Mach Functions 
  ! M(1)
  real, dimension(0:nAlts,1:nSpecies) ::  MachSFn1PLeft,  MachSFn1NLeft
  real, dimension(0:nAlts,1:nSpecies) :: MachSFn1PRight, MachSFn1NRight
  ! M(2)
  real, dimension(0:nAlts,1:nSpecies) ::  MachSFn2PLeft,  MachSFn2NLeft
  real, dimension(0:nAlts,1:nSpecies) :: MachSFn2PRight, MachSFn2NRight
  ! M(4)
  real, dimension(0:nAlts,1:nSpecies) ::  MachSFn4PLeft,  MachSFn4NLeft
  real, dimension(0:nAlts,1:nSpecies) :: MachSFn4PRight, MachSFn4NRight

  ! MPress
  real, dimension(0:nAlts,1:nSpecies) ::  AUSMMachPressureParamS
  ! Numerical Velocity and Pressure
  real, dimension(0:nAlts,1:nSpecies) ::  NumericalVertVelS
  real, dimension(0:nAlts,1:nSpecies) ::  NumericalPressureS
  real, dimension(0:nAlts           ) ::  BulkNumericalVertVel
  ! Interface Flux
  real, dimension(0:nAlts,1:nSpecies) ::  RhoSFlux
  real, dimension(0:nAlts,1:nSpecies) ::  MomentumSFlux
  real, dimension(0:nAlts,1:3       ) ::  MomentumFlux
  real, dimension(0:nAlts           ) ::  EnergyFlux
  ! Thornber Correction
  real :: ZVar, VL, VR, CBar, ML, MR
  real :: MachScaling

  real, dimension(0:nAlts,1:nSpecies) :: PsiP_Left, PsiN_Right
  ! Chen [2022], g(M) function
  real, dimension(0:nAlts,1:nSpecies) :: gMfuncS
  real, dimension(0:nAlts,1:nSpecies) :: kL, kR
  integer :: jAlt
  real :: DeltaPU, DeltaPU_P1, DeltaPU_M1
  real :: DeltaMU, DeltaMU_P1, DeltaMU_M1

  real :: NonHydrostaticPressureSLeft, NonHydrostaticPressureSRight

  ! New Streamlined Variables
  real, dimension(0:nAlts,1:nSpecies) ::         RhoSLeftWENO,         RhoSRightWENO
  real, dimension(0:nAlts,1:nSpecies) ::      LogRhoSLeftWENO,      LogRhoSRightWENO
  real, dimension(0:nAlts,1:nSpecies) ::    HydroRhoSLeftWENO,    HydroRhoSRightWENO
  real, dimension(0:nAlts,1:nSpecies) :: LogHydroRhoSLeftWENO, LogHydroRhoSRightWENO

  ! New Pressure Variables
  real, dimension(0:nAlts,1:nSpecies) :: PressureSLeftWENO, PressureSRightWENO
  real, dimension(0:nAlts,1:nSpecies) :: HydroPressureSLeftWENO, HydroPressureSRightWENO
  ! New Species Winds 
  real, dimension(0:nAlts,1:nSpecies) ::     VertVelSLeftWENO,     VertVelSRightWENO
  ! Bulk Fields
  real, dimension(0:nAlts)            :: RhoLeftWENO, RhoRightWENO
  real, dimension(0:nAlts)            :: TemperatureLeftWENO, TemperatureRightWENO
  real, dimension(0:nAlts)            :: PressureLeftWENO, PressureRightWENO
  real, dimension(0:nAlts,1:3)        :: VelGDLeftWENO, VelGDRightWENO
  real, dimension(0:nAlts)            :: GammaLeftWENO, GammaRightWENO
  ! Derived Fields
  real, dimension(0:nAlts)            :: EnergyLeftWENO, EnergyRightWENO
  real, dimension(0:nAlts)            :: EnthalpyLeftWENO, EnthalpyRightWENO
  real, dimension(0:nAlts,1:nSpecies) :: EnthalpySLeftWENO, EnthalpySRightWENO

  ! Sound Speed Fields
  real, dimension(0:nAlts)            :: SoundSpeedLeftWENO, SoundSpeedRightWENO
  real, dimension(0:nAlts,1:nSpecies) :: SoundSpeedSLeftWENO, SoundSpeedSRightWENO

  ! Mach Numbers 
  real, dimension(0:nAlts,1:nSpecies) :: MachSLeftWENO, MachSRightWENO
  real, dimension(0:nAlts,1:nSpecies) :: MachSZero2WENO, MachSZeroWENO
  ! AUSM+-up Parameters
  
  ! Mach Functions 
  ! M(1)
  real, dimension(0:nAlts,1:nSpecies) ::  MachSFn1PLeftWENO,  MachSFn1NLeftWENO
  real, dimension(0:nAlts,1:nSpecies) :: MachSFn1PRightWENO, MachSFn1NRightWENO
  ! M(2)
  real, dimension(0:nAlts,1:nSpecies) ::  MachSFn2PLeftWENO,  MachSFn2NLeftWENO
  real, dimension(0:nAlts,1:nSpecies) :: MachSFn2PRightWENO, MachSFn2NRightWENO
  ! M(4)
  real, dimension(0:nAlts,1:nSpecies) ::  MachSFn4PLeftWENO,  MachSFn4NLeftWENO
  real, dimension(0:nAlts,1:nSpecies) :: MachSFn4PRightWENO, MachSFn4NRightWENO


  real, dimension(0:nAlts,1:nSpecies) :: PsiP_LeftWENO, PsiN_RightWENO
  real :: NonHydrostaticPressureSLeftWENO, NonHydrostaticPressureSRightWENO

  ! Numerical Velocity and Pressure
  real, dimension(0:nAlts,1:nSpecies) ::  NumericalVertVelSWENO
  real, dimension(0:nAlts,1:nSpecies) ::  NumericalPressureSWENO
  real, dimension(0:nAlts           ) ::  BulkNumericalVertVelWENO
  real, dimension(0:nAlts,1:nSpecies) ::  AUSMMachPressureParamSWENO
  real, dimension(0:nAlts,1:nSpecies) :: InterfaceSoundSpeedSWENO
  real, dimension(0:nAlts)            :: InterfaceSoundSpeedWENO
  real, dimension(0:nAlts,1:nSpecies) ::     MeanRhoSWENO
  real, dimension(0:nAlts,1:nSpecies) :: MeanPressureSWENO
  real, dimension(0:nAlts,1:nSpecies) :: MeanHydroPressureSWENO
  real, dimension(0:nAlts,1:nSpecies) :: MachSBar2WENO
  real, dimension(0:nAlts,1:nSpecies) :: KpSParamWENO, KuSParamWENO
  ! Interface Flux
  real, dimension(0:nAlts,1:nSpecies) ::  RhoSFluxWENO
  real, dimension(0:nAlts,1:nSpecies) ::  MomentumSFluxWENO
  real, dimension(0:nAlts,1:3       ) ::  MomentumFluxWENO
  real, dimension(0:nAlts           ) ::  EnergyFluxWENO
  ! Thornber Correction
  real :: ZVarWENO, VLWENO, VRWENO, CBarWENO, MLWENO, MRWENO
  real :: MachScalingWENO
  ! Chen [2022], g(M) function
  real, dimension(0:nAlts,1:nSpecies) :: gMfuncSWENO
  real, dimension(0:nAlts,1:nSpecies) :: kLWENO, kRWENO

  real, dimension(0:nAlts,1:nSpecies) ::  LogNumericalPressureS
  real, dimension(0:nAlts,1:nSpecies) ::  LogNumericalHydroPressureS
  real, dimension(0:nAlts,1:nSpecies) ::     NumericalHydroPressureS
  real, dimension(1:nAlts,1:nSpecies) ::  CellCenteredLogNumericalPressureS
  real, dimension(1:nAlts,1:nSpecies) ::  CellCenteredLogNumericalHydroPressureS
  real, dimension(1:nAlts,1:nSpecies) ::  CellCenteredNumericalPressureS
  real, dimension(1:nAlts,1:nSpecies) ::  CellCenteredNumericalHydroPressureS

  real, dimension(0:nAlts,1:nSpecies) ::  LogNumericalPressureSWENO
  real, dimension(0:nAlts,1:nSpecies) ::  LogNumericalHydroPressureSWENO
  real, dimension(0:nAlts,1:nSpecies) ::     NumericalHydroPressureSWENO
  real, dimension(1:nAlts,1:nSpecies) ::  CellCenteredLogNumericalPressureSWENO
  real, dimension(1:nAlts,1:nSpecies) ::  CellCenteredLogNumericalHydroPressureSWENO
  real, dimension(1:nAlts,1:nSpecies) ::  CellCenteredNumericalPressureSWENO
  real, dimension(1:nAlts,1:nSpecies) ::  CellCenteredNumericalHydroPressureSWENO



  ! BEGIN SETTING UP VARIABLES
  ! Basic Geometry of the Grid
  !   |------------------|--------- X ---------|------------------|
  !   |-----(i-1)--------|--------- i ---------|------(i+1)-------|
  !                   (i-1/2)               (i+1/2)
  !                    L | R                 L | R
  !            UL_M12(i) | UR_M12(i)  UL_P12(i) | UR_P12(i)
  !                     i ranges from 1 -> nAlts
  !
  ! Current Variable Setup focuses only on the +1/2 Values from 0 to nAlts
  !         |           UL(i-1)| UR(i-1)       UL(i)| UR(i)
  !                           i ranges from 0 -> nAlts
  ! U(i)  |---0---|---- 1 ------|  ...(i-1) -- | -- (i) -- | ... (nAlts) -- | -- (nAlts+1) 
  ! UI(i) |--   UI(0) -----   UI(1) ...      UI(i-1)         ...        UI(nAlts) 
  ! L/R   | [UL(0)|UR(0)] [UL(1)|UR(1)] [UL(i-1)|UR(i-1)]    ...  [UL(nAlts)|UR(nAlts)]


  MInf = 1.0e-19
  LiouBeta  = 1.0/8.0
  LiouAlpha = 3.0/16.0 ! just set FA = 1.0

       LogRhoS(-1:nAlts+2,1:nSpecies) =      alog(RhoS(-1:nAlts+2,1:nSpecies))
  LogHydroRhoS(-1:nAlts+2,1:nSpecies) = alog(HydroRhoS(-1:nAlts+2,1:nSpecies))
  ! BEGINWENO
  ! Establish Cell Interface values 
  ! Get the RhoS facevaluess and Bulk Rho facevalues
   RhoLeftWENO(0:nAlts) = 0.0
  RhoRightWENO(0:nAlts) = 0.0
  do iSpecies = 1, nSpecies
     call calc_weno_facevalues(LogRhoS(-1:nAlts+2,iSpecies), &
                      LogRhoSLeftWENO( 0:nAlts  ,iSpecies), &
                     LogRhoSRightWENO( 0:nAlts  ,iSpecies) )
      RhoSLeftWENO(:,iSpecies) = exp( LogRhoSLeftWENO(:,iSpecies))
     RhoSRightWENO(:,iSpecies) = exp(LogRhoSRightWENO(:,iSpecies))
     ! Bulk Rho
      RhoLeftWENO(:) =  RhoLeftWENO(:) +    RhoSLeftWENO(:,iSpecies)
     RhoRightWENO(:) = RhoRightWENO(:) +   RhoSRightWENO(:,iSpecies)
  enddo 
  ! End RhoS, Rho Facevalues

  ! Begin HydroRhoS, HydroRho Facevalues:
  ! Note: Hydro denotes hydrostatic background
  do iSpecies = 1, nSpecies
     call calc_weno_facevalues(LogHydroRhoS(-1:nAlts+2,iSpecies), &
                      LogHydroRhoSLeftWENO( 0:nAlts  ,iSpecies), &
                     LogHydroRhoSRightWENO( 0:nAlts  ,iSpecies) )
      HydroRhoSLeftWENO(:,iSpecies) = exp( LogHydroRhoSLeftWENO(:,iSpecies))
     HydroRhoSRightWENO(:,iSpecies) = exp(LogHydroRhoSRightWENO(:,iSpecies))
  enddo 

  ! Calculate the LeftWENO and RightWENO Faces of the Temperatures 
  call calc_weno_facevalues(Temp(-1:nAlts+2), &
                TemperatureLeftWENO( 0:nAlts  ), &
               TemperatureRightWENO( 0:nAlts  ) )

  PressureRightWENO = 0.0
   PressureLeftWENO = 0.0
  do iSpecies = 1, nSpecies
     PressureSLeftWENO(0:nAlts,iSpecies) = &
         (RhoSLeftWENO(0:nAlts,iSpecies)/Mass(iSpecies))*&
         Boltzmanns_Constant*TemperatureLeftWENO(0:nAlts)
     PressureLeftWENO = PressureLeftWENO + PressureSLeftWENO(:,iSpecies)

     HydroPressureSLeftWENO(0:nAlts,iSpecies) = &
         (HydroRhoSLeftWENO(0:nAlts,iSpecies)/Mass(iSpecies))*&
        Boltzmanns_Constant*TemperatureLeftWENO(0:nAlts)

     PressureSRightWENO(0:nAlts,iSpecies) = &
        (RhoSRightWENO(0:nAlts,iSpecies)/Mass(iSpecies))*&
         Boltzmanns_Constant*TemperatureRightWENO(0:nAlts)
     PressureRightWENO = PressureRightWENO + PressureSRightWENO(:,iSpecies)

    HydroPressureSRightWENO(0:nAlts,iSpecies) = &
        (HydroRhoSRightWENO(0:nAlts,iSpecies)/Mass(iSpecies))*&
         Boltzmanns_Constant*TemperatureRightWENO(0:nAlts)
  enddo 

  ! Begin Vertical, Horizontal Species and Bulk Winds
  do iSpecies = 1, nSpecies
    call calc_weno_facevalues(VertVel(:,iSpecies), &
                        VertVelSLeftWENO(:,iSpecies), &
                       VertVelSRightWENO(:,iSpecies))
  enddo 

  do iDim = 1, 3
     call calc_weno_facevalues(Vel_GD(:,iDim), &
                           VelGDLeftWENO(:,iDim), &
                          VelGDRightWENO(:,iDim))
  enddo 

   VelGDLeftWENO(:,iUp_) = 0.0
  VelGDRightWENO(:,iUp_) = 0.0
  do iSpecies = 1, nSpecies
      VelGDLeftWENO(:,iUp_) = VelGDLeftWENO(:,iUp_) + &
           RhoSLeftWENO(:,iSpecies)*VertVelSLeftWENO(:,iSpecies)/&
            RhoLeftWENO(:)

      VelGDRightWENO(:,iUp_) = VelGDRightWENO(:,iUp_) + &
           RhoSRightWENO(:,iSpecies)*VertVelSRightWENO(:,iSpecies)/&
            RhoRightWENO(:)
  enddo !iSpecies = 1, nSpecies

  call calc_weno_facevalues(Gamma_1d(:), GammaLeftWENO(:), GammaRightWENO(:))

  do iAlt = 0, nAlts
     EnergyLeftWENO(iAlt) = &
         ( 1.0/(GammaLeftWENO(iAlt) - 1.0))*PressureLeftWENO(iAlt) + &
           0.5*RhoLeftWENO(iAlt)*&
          (VelGDLeftWENO(iAlt,iUp_)**2.0 + VelGDLeftWENO(iAlt,iEast_)**2.0 + &
           VelGDLeftWENO(iAlt,iNorth_)**2.0)

     EnergyRightWENO(iAlt) = &
         ( 1.0/(GammaRightWENO(iAlt) - 1.0))*PressureRightWENO(iAlt) + &
           0.5*RhoRightWENO(iAlt)*&
          (VelGDRightWENO(iAlt,iUp_)**2.0 + VelGDRightWENO(iAlt,iEast_)**2.0 + &
           VelGDRightWENO(iAlt,iNorth_)**2.0)
  enddo 

  do iAlt = 0, nAlts 
     EnthalpyLeftWENO(iAlt) = &
       0.5*(VelGDLeftWENO(iAlt,iUp_)**2.0 + &
            VelGDLeftWENO(iAlt,iEast_)**2.0 + &
            VelGDLeftWENO(iAlt,iNorth_)**2.0) + &
          (GammaLeftWENO(iAlt)/(GammaLeftWENO(iAlt) - 1.0))*&
          PressureLeftWENO(iAlt)/RhoLeftWENO(iAlt)

     EnthalpyRightWENO(iAlt) = &
       0.5*(VelGDRightWENO(iAlt,iUp_)**2.0 + &
            VelGDRightWENO(iAlt,iEast_)**2.0 + &
            VelGDRightWENO(iAlt,iNorth_)**2.0) + &
          (GammaRightWENO(iAlt)/(GammaRightWENO(iAlt) - 1.0))*&
          PressureRightWENO(iAlt)/RhoRightWENO(iAlt)

     do iSpecies = 1, nSpecies

        EnthalpySLeftWENO(iAlt,iSpecies) = &
          0.5*(VertVelSLeftWENO(iAlt,iSpecies)**2.0 + &
               VelGDLeftWENO(iAlt,iEast_)**2.0 + &
               VelGDLeftWENO(iAlt,iNorth_)**2.0) + &
             (GammaLeftWENO(iAlt)/(GammaLeftWENO(iAlt) - 1.0))*&
             PressureSLeftWENO(iAlt,iSpecies)/RhoSLeftWENO(iAlt,iSpecies)

        EnthalpySRightWENO(iAlt,iSpecies) = &
          0.5*(VertVelSRightWENO(iAlt,iSpecies)**2.0 + &
               VelGDRightWENO(iAlt,iEast_)**2.0 + &
               VelGDRightWENO(iAlt,iNorth_)**2.0) + &
             (GammaRightWENO(iAlt)/(GammaRightWENO(iAlt) - 1.0))*&
             PressureSRightWENO(iAlt,iSpecies)/RhoSRightWENO(iAlt,iSpecies)

     enddo !iSpecies = 1, nSpecies
  enddo 

  do iAlt = 0, nAlts
     SubCs = &
        sqrt(2.0*( (GammaLeftWENO(iAlt) - 1.0 )/(GammaLeftWENO(iAlt) + 1.0)) *&
                 EnthalpyLeftWENO(iAlt) )
     SoundSpeedLeftWENO(iAlt) = (SubCs**2.0)/max(SubCs, VelGDLeftWENO(iAlt,iUp_))

     ! Note the -1.0 factor in the RightWENO face
     SubCs = &
        sqrt(2.0*( (GammaRightWENO(iAlt) - 1.0 )/(GammaRightWENO(iAlt) + 1.0)) *&
                 EnthalpyRightWENO(iAlt) )
     SoundSpeedRightWENO(iAlt) = (SubCs**2.0)/max(SubCs, -1.0*VelGDRightWENO(iAlt,iUp_))

     InterfaceSoundSpeedWENO(iAlt) = min(SoundSpeedLeftWENO(iAlt), SoundSpeedRightWENO(iAlt))

     do iSpecies = 1, nSpecies
        SubCs = &
           sqrt(2.0*( (GammaLeftWENO(iAlt) - 1.0 )/(GammaLeftWENO(iAlt) + 1.0)) *&
                    EnthalpySLeftWENO(iAlt,iSpecies) )
        SoundSpeedSLeftWENO(iAlt,iSpecies) = &
           (SubCs**2.0)/max(SubCs, VertVelSLeftWENO(iAlt,iSpecies))

        SubCs = &
           sqrt(2.0*( (GammaRightWENO(iAlt) - 1.0 )/(GammaRightWENO(iAlt) + 1.0)) *&
                    EnthalpySRightWENO(iAlt,iSpecies) )
        ! Note the -1.0 factor in the RightWENO face
        SoundSpeedSRightWENO(iAlt,iSpecies) = &
           (SubCs**2.0)/max(SubCs, -1.0*VertVelSRightWENO(iAlt,iSpecies))
        InterfaceSoundSpeedSWENO(iAlt,iSpecies) = &
           min(SoundSpeedSLeftWENO(iAlt,iSpecies), SoundSpeedSRightWENO(iAlt,iSpecies))
     enddo !iSpecies = 1, nSpecies
  enddo !iAlt = 0, nAlts

  do iSpecies = 1, nSpecies
     do iAlt = 0, nAlts
        ! Implement Thornber Entropy Correction at low Mach Numbers
        ! see. Thornber et al. [2008], J. Comp. Phys, 227, 4873-4894
        ! Implementation here taken from Houim et al. [2011], J Comp.Phys, 230, 8527-8553
        CBarWENO = InterfaceSoundSpeedSWENO(iAlt,iSpecies)
        VLWENO   =         VertVelSLeftWENO(iAlt,iSpecies)
        VRWENO   =        VertVelSRightWENO(iAlt,iSpecies)
        MLWENO   = sqrt(VLWENO*VLWENO)/CBarWENO
        MRWENO   = sqrt(VRWENO*VRWENO)/CBarWENO
        ZVarWENO = min(1.0, max(MLWENO,MRWENO))
         VertVelSLeftWENO(iAlt,iSpecies) = 0.5*(VLWENO + VRWENO) + 0.5*ZVarWENO*(VLWENO - VRWENO)
        VertVelSRightWENO(iAlt,iSpecies) = 0.5*(VLWENO + VRWENO) + 0.5*ZVarWENO*(VRWENO - VLWENO)
     enddo !iAlt = 1, nAlts
  enddo !iSpecies = 1, nSpecies

  MeanPressureSWENO = 0.5*(PressureSLeftWENO + PressureSRightWENO)
  MeanHydroPressureSWENO = 0.5*(HydroPressureSLeftWENO + HydroPressureSRightWENO)
  MeanRhoSWENO = 0.5*(RhoSLeftWENO + RhoSRightWENO)

  do iAlt = 0, nAlts 
    do iSpecies = 1, nSpecies
       MachSLeftWENO(iAlt,iSpecies) = &
         VertVelSLeftWENO(iAlt,iSpecies)/InterfaceSoundSpeedSWENO(iAlt,iSpecies)
       MachSRightWENO(iAlt,iSpecies) = &
         VertVelSRightWENO(iAlt,iSpecies)/InterfaceSoundSpeedSWENO(iAlt,iSpecies)
    enddo !iSpecies = 1, nSpecies
  enddo !iAlt = 1, nAlts 

  do iAlt = 0, nAlts
     do iSpecies = 1, nSpecies
        MachSBar2WENO(iAlt,iSpecies) = &
           0.5*( MachSLeftWENO(iAlt,iSpecies)**2.0 + &
                MachSRightWENO(iAlt,iSpecies)**2.0 )
        MachSZero2WENO(iAlt,iSpecies) = &
           min(1.0, max(MachSBar2WENO(iAlt,iSpecies), MInf))
        MachSZeroWENO(iAlt,iSpecies) = sqrt(MachSZero2WENO(iAlt,iSpecies))
     enddo !iSpecies = 1, nSpecies
  enddo !iAlt = 1, nAlts

  ! Create Mach-Scaling Pressure-Correction Functions
  ! Establish Mach Number Scaling Functions g(M)
  ! We specify a g(M) such that it goes away as |M|->1
  ! g(M) ~ M^2 as M -> 0.0
  Kp(1:nSpecies) = 0.10             !! Ullrich et al. [2011]
  Ku(1:nSpecies) = 10.00             !! Ullrich et al. [2011]

  do iAlt = 0, nAlts
     do iSpecies = 1, nSpecies
        ! Note that MZero is scaled to be within [0,1]
        ! Specify pressure correction fluxes (LiouKpS)
        ! Note: MZero -> [0,1]
        ! Note: M2Zero  = (MZero)^2.0-> [0,1]
        ! As M -> 0.0, MachScaling ~ M^2 - M^4 ~ M^2 -> Varies like a pressure at Low Mach
        ! As M >> 1.0, MZero = 1.0, MachScaling ~ 1 - 1 -> 0.0
        !    This shuts our correction down at high mach number
        ! Mathematically, we want AUSM+-up Numerical Flux(Mp):
        !        Mp ~ d(P)*/(rho*cs*(fa*cs + zeta*dz/dt)) ~ dP/rho*cs^2

        MachScalingWENO = MachSZeroWENO(iAlt,iSpecies)**2.0
        KpSParamWENO(iAlt,iSpecies) = MachScalingWENO*Kp(iSpecies)

        ! Note that MZero ranges from [0, 1.0] (cap at 1.0)
        !MachScalingWENO = 1.0 - MachSZeroWENO(iAlt,iSpecies)
        !KuSParamWENO(iAlt,iSpecies) = MachScaling*Ku(iSpecies)
        !MachScalingWENO = MachSZeroWENO(iAlt,iSpecies)*(1.0 - MachSZeroWENO(iAlt,iSpecies))
        MachScalingWENO = (1.0 - MachSZeroWENO(iAlt,iSpecies))
        KuSParamWENO(iAlt,iSpecies) = MachScaling*Ku(iSpecies)
     enddo 
  enddo 

  ! BEGIN Mach Functions 
  ! M(1)
  do iAlt = 0, nAlts
     do iSpecies = 1, nSpecies
        MachSFn1PLeftWENO(iAlt,iSpecies) = &
            0.5*(MachSLeftWENO(iAlt,iSpecies) + abs(MachSLeftWENO(iAlt,iSpecies)) )
        MachSFn1PRightWENO(iAlt,iSpecies) = &
            0.5*(MachSRightWENO(iAlt,iSpecies) + abs(MachSRightWENO(iAlt,iSpecies)) )

        MachSFn1NLeftWENO(iAlt,iSpecies) = &
            0.5*(MachSLeftWENO(iAlt,iSpecies) - abs(MachSLeftWENO(iAlt,iSpecies)) )
        MachSFn1NRightWENO(iAlt,iSpecies) = &
            0.5*(MachSRightWENO(iAlt,iSpecies) - abs(MachSRightWENO(iAlt,iSpecies)) )

     enddo !iSpecies = 1, nSpecies
  enddo !iAlt = 1, nAlts
  ! M(2)
  do iAlt = 0, nAlts
     do iSpecies = 1, nSpecies
         MachSFn2PLeftWENO(iAlt,iSpecies) =  0.25*( MachSLeftWENO(iAlt,iSpecies) + 1.0)**2.0
        MachSFn2PRightWENO(iAlt,iSpecies) =  0.25*(MachSRightWENO(iAlt,iSpecies) + 1.0)**2.0

         MachSFn2NLeftWENO(iAlt,iSpecies) = -0.25*( MachSLeftWENO(iAlt,iSpecies) - 1.0)**2.0
        MachSFn2NRightWENO(iAlt,iSpecies) = -0.25*(MachSRightWENO(iAlt,iSpecies) - 1.0)**2.0
     enddo 
  enddo 
  ! M(4)
  do iAlt = 0, nAlts
    do iSpecies = 1, nSpecies

      if ( abs(MachSLeftWENO(iAlt,iSpecies)) .ge. 1.0) then 
         MachSFn4PLeftWENO(iAlt,iSpecies) = MachSFn1PLeftWENO(iAlt,iSpecies)

         ! These are the Slau2 versions (alpha = 0.0)
         PsiP_LeftWENO(iAlt,iSpecies) = &
            0.5*(1.0 + abs(MachSLeftWENO(iAlt,iSpecies))/MachSLeftWENO(iAlt,iSpecies) )
      else
         MachSFn4PLeftWENO(iAlt,iSpecies) = MachSFn2PLeftWENO(iAlt,iSpecies)*&
                   (1.0 - 16.0*LiouBeta*MachSFn2NLeftWENO(iAlt,iSpecies))

         ! These are the Slau2 versions (alpha = 0.0)
         PsiP_LeftWENO(iAlt,iSpecies) = &
            0.25*( (MachSLeftWENO(iAlt,iSpecies) + 1.0)**2.0)*(2.0 - MachSLeftWENO(iAlt,iSpecies)) + &
                  LiouAlpha*MachSLeftWENO(iAlt,iSpecies)*(MachSLeftWENO(iAlt,iSpecies)**2.0 - 1.0)**2.0 
      endif 
!
      if ( abs(MachSRightWENO(iAlt,iSpecies)) .ge. 1.0) then 
         MachSFn4NRightWENO(iAlt,iSpecies) = MachSFn1NRightWENO(iAlt,iSpecies)
         ! These are the Slau2 versions (alpha = 0.0)
         PsiN_RightWENO(iAlt,iSpecies) = &
            0.5*(1.0 + abs(MachSRightWENO(iAlt,iSpecies))/MachSRightWENO(iAlt,iSpecies))
      else
         MachSFn4NRightWENO(iAlt,iSpecies) = MachSFn2NRightWENO(iAlt,iSpecies)*&
                   (1.0 + 16.0*LiouBeta*MachSFn2PRightWENO(iAlt,iSpecies))
         ! These are the Slau2 versions (alpha = 0.0)
         !PsiN_RightWENO(iAlt,iSpecies) = &
         !   0.25*( (MachSRightWENO(iAlt,iSpecies) - 1.0)**2.0)*(2.0 + MachSRightWENO(iAlt,iSpecies)) 
         PsiN_RightWENO(iAlt,iSpecies) = &
            0.25*( (MachSRightWENO(iAlt,iSpecies) - 1.0)**2.0)*(2.0 + MachSRightWENO(iAlt,iSpecies)) - &
                  LiouAlpha*MachSRightWENO(iAlt,iSpecies)*(MachSRightWENO(iAlt,iSpecies)**2.0 - 1.0)**2.0 
      endif 
    enddo 
  enddo 

  !!! Begin 2nd Order Mach Number Polynomials
  !real, dimension(0:nAlts,1:nSpecies) ::  MachPressureSParamAUSM
  do iAlt = 0, nAlts 
   do iSpecies = 1, nSpecies 
      ! Note that MPress varies as Kp*(g(M))* dPs/(rhos*as^2)
      ! g(M) is a "shut off" for M > 1
      ! in all cases rhos*as^2 ~ rhos*(Gamma*Ps/rhos) ~ Gamma*Ps 
      ! now:  In the vertical direction
      !
      ! dPs/Ps ~ (1/Ps)*(dP/dr)*dr ~ -rhos*(grav)/(ns*kb*T)*dr  
      ! Thus as M -> 0 and P approaches hydrostatic
      ! Therefore MPress ~ Kp*ms*g*dr/(kb*T) (no units)
      ! But we want this to vary as M^2 as a pressure flux as M->0
      ! So we want APC Scalar as f(M) ~ Ms^2
      ! Mp_APC ~ f(M)*Kp*dPs/(rhos*cs^2) -> Ms^2*Kp*(dr/Hs)
!!!!! Hydrostatic Species Use This one

      AUSMMachPressureParamSWENO(iAlt,iSpecies) = &
          KpSParamWENO(iAlt,iSpecies)*max( (1.0 - MachSBar2WENO(iAlt,iSpecies)), 0.0)*&
         ( (PressureSRightWENO(iAlt, iSpecies) - HydroPressureSRightWENO(iAlt,iSpecies) ) - &
           ( PressureSLeftWENO(iAlt, iSpecies) -  HydroPressureSLeftWENO(iAlt,iSpecies) ) )/&
         ( MeanRhoSWENO(iAlt,iSpecies)*InterfaceSoundSpeedSWENO(iAlt,iSpecies)**2.0)

   enddo !iSpecies = 1, nSpecies 
  enddo !iAlt = 0, nAlts 


  do iAlt = 0, nAlts 
     do iSpecies = 1, nSpecies
        NumericalVertVelSWENO(iAlt,iSpecies) = &
          InterfaceSoundSpeedSWENO(iAlt,iSpecies)*&
          ( MachSFn4PLeftWENO(iAlt,iSpecies) + &  
           MachSFn4NRightWENO(iAlt,iSpecies) - &  
          AUSMMachPressureParamSWENO(iAlt,iSpecies))
     enddo !iSpecies = 1, nSpecies
  enddo !iAlt = 0, nAlts 

!  do iAlt = 0, nAlts
!     do iSpecies = 1, nSpecies
!         NonHydrostaticPressureSLeftWENO  =  PressureSLeftWENO(iAlt,iSpecies) -  HydroPressureSLeftWENO(iAlt,iSpecies)
!        NonHydrostaticPressureSRightWENO  = PressureSRightWENO(iAlt,iSpecies) - HydroPressureSRightWENO(iAlt,iSpecies)
!
!        NumericalPressureSWENO(iAlt,iSpecies) = &
!             PsiP_LeftWENO(iAlt,iSpecies)*NonHydrostaticPressureSLeftWENO +  &
!            PsiN_RightWENO(iAlt,iSpecies)*NonHydrostaticPressureSRightWENO - &
!            KuSParamWENO(iAlt,iSpecies)*PsiP_LeftWENO(iAlt,iSpecies)*PsiN_RightWENO(iAlt,iSpecies)*&
!            (RhoSLeftWENO(iAlt,iSpecies) + RhoSRightWENO(iAlt,iSpecies))*&
!             InterfaceSoundSpeedSWENO(iAlt,iSpecies)*&
!            (VertVelSLeftWENO(iAlt,iSpecies) - VertVelSRightWENO(iAlt,iSpecies))
!     enddo !iSpecies = 1, nSpecies
!  enddo !iAlt = 1, nAlts

  do iAlt = 0, nAlts
     do iSpecies = 1, nSpecies
        NumericalPressureSWENO(iAlt,iSpecies) = &
             PsiP_LeftWENO(iAlt,iSpecies)*PressureSLeftWENO(iAlt,iSpecies) +  &
            PsiN_RightWENO(iAlt,iSpecies)*PressureSRightWENO(iAlt,iSpecies) - &
            KuSParamWENO(iAlt,iSpecies)*PsiP_LeftWENO(iAlt,iSpecies)*PsiN_RightWENO(iAlt,iSpecies)*&
            (RhoSLeftWENO(iAlt,iSpecies) + RhoSRightWENO(iAlt,iSpecies))*&
             InterfaceSoundSpeedSWENO(iAlt,iSpecies)*&
            (VertVelSLeftWENO(iAlt,iSpecies) - VertVelSRightWENO(iAlt,iSpecies))

        LogNumericalPressureSWENO(iAlt,iSpecies) = dlog(NumericalPressureSWENO(iAlt,iSpecies))

        NumericalHydroPressureSWENO(iAlt,iSpecies) = &
             PsiP_LeftWENO(iAlt,iSpecies)*HydroPressureSLeftWENO(iAlt,iSpecies) +  &
            PsiN_RightWENO(iAlt,iSpecies)*HydroPressureSRightWENO(iAlt,iSpecies) 

        LogNumericalHydroPressureSWENO(iAlt,iSpecies) = dlog(NumericalHydroPressureSWENO(iAlt,iSpecies))

     enddo !iSpecies = 1, nSpecies
  enddo !iAlt = 1, nAlts

  do iAlt = 1, nAlts
    do iSpecies = 1, nSpecies
       CellCenteredLogNumericalPressureSWENO(iAlt,iSpecies) = &
           0.5*( LogNumericalPressureSWENO(iAlt  ,iSpecies) + & 
                 LogNumericalPressureSWENO(iAlt-1,iSpecies) )
       CellCenteredLogNumericalHydroPressureSWENO(iAlt,iSpecies) = &
           0.5*( LogNumericalHydroPressureSWENO(iAlt  ,iSpecies) + & 
                 LogNumericalHydroPressureSWENO(iAlt-1,iSpecies) )

       CellCenteredNumericalPressureSWENO(iAlt,iSpecies) = &
           exp(CellCenteredLogNumericalPressureSWENO(iAlt,iSpecies))
       CellCenteredNumericalHydroPressureSWENO(iAlt,iSpecies) = &
           exp(CellCenteredLogNumericalHydroPressureSWENO(iAlt,iSpecies))
    enddo !iSpecies = 1, nSpecies
  enddo 
  !ENDWENO STUFF
  !ENDWENO STUFF

  !BEGINWENOFLUX 
  do iSpecies = 1, nSpecies
     do iAlt = 0, nAlts 

        ! Note if V > 0.0, then |V| =  V:  (1/2)*(V + |V|) = V, (1/2)*(V - |V|) = 0.0
        ! Note if V < 0.0, then |V| = -V:  (1/2)*(V + |V|) = V, (1/2)*(V - |V|) = 0.0
        RhoSFluxWENO(iAlt,iSpecies) = &
           0.5* RhoSLeftWENO(iAlt,iSpecies)*&
                (     NumericalVertVelSWENO(iAlt,iSpecies) +  & 
                  dabs(NumericalVertVelSWENO(iAlt,iSpecies)) ) + &
           0.5*RhoSRightWENO(iAlt,iSpecies)*&
                (     NumericalVertVelSWENO(iAlt,iSpecies) -  & 
                  dabs(NumericalVertVelSWENO(iAlt,iSpecies)) ) 

        MomentumSFluxWENO(iAlt,iSpecies) = &
           0.5* RhoSLeftWENO(iAlt,iSpecies)*VertVelSLeftWENO(iAlt,iSpecies)*&
                (     NumericalVertVelSWENO(iAlt,iSpecies) +    &
                  dabs(NumericalVertVelSWENO(iAlt,iSpecies)) ) + &
           0.5*RhoSRightWENO(iAlt,iSpecies)*VertVelSRightWENO(iAlt,iSpecies)*&
                (     NumericalVertVelSWENO(iAlt,iSpecies) -    &
                  dabs(NumericalVertVelSWENO(iAlt,iSpecies)) ) 
     enddo 
  enddo 

  BulkNumericalVertVelWENO(0:nAlts) = 0.0
  do iSpecies = 1, nSpecies
     BulkNumericalVertVelWENO = &
     BulkNumericalVertVelWENO + &
           ( RhoSLeftWENO(:,iSpecies) + RhoSRightWENO(:,iSpecies) )*&
             NumericalVertVelSWENO(:,iSpecies)/&
             (RhoLeftWENO + RhoRightWENO )
  enddo 

  do iAlt = 0, nAlts 
     EnergyFluxWENO(iAlt) = &
        0.5*( EnergyLeftWENO(iAlt) + PressureLeftWENO(iAlt) )* &
            ( BulkNumericalVertVelWENO(iAlt) +  &
          abs(BulkNumericalVertVelWENO(iAlt)) ) + & 
        0.5*( EnergyRightWENO(iAlt) + PressureRightWENO(iAlt) )* &
            ( BulkNumericalVertVelWENO(iAlt) -  &
          abs(BulkNumericalVertVelWENO(iAlt)) ) 

     MomentumFluxWENO(iAlt,1:3) = &
        0.5*( RhoLeftWENO(iAlt)*VelGDLeftWENO(iAlt,1:3))*&
            ( BulkNumericalVertVelWENO(iAlt) +  &
          abs(BulkNumericalVertVelWENO(iAlt)) ) + & 
        0.5*( RhoRightWENO(iAlt)*VelGDRightWENO(iAlt,1:3))*&
            ( BulkNumericalVertVelWENO(iAlt) -  &
          abs(BulkNumericalVertVelWENO(iAlt)) ) 
  enddo 
  ! ENDWENOFLUX

  ! BEGINTVD
  ! Establish Cell Interface values 
  ! Get the RhoS facevaluess and Bulk Rho facevalues
   RhoLeft(0:nAlts) = 0.0
  RhoRight(0:nAlts) = 0.0
  do iSpecies = 1, nSpecies
     call calc_tvd_facevalues(LogRhoS(-1:nAlts+2,iSpecies), &
                          LogRhoSLeft( 0:nAlts  ,iSpecies), &
                         LogRhoSRight( 0:nAlts  ,iSpecies) )
      RhoSLeft(:,iSpecies) = exp( LogRhoSLeft(:,iSpecies))
     RhoSRight(:,iSpecies) = exp(LogRhoSRight(:,iSpecies))
     ! Bulk Rho
      RhoLeft(:) =  RhoLeft(:) +    RhoSLeft(:,iSpecies)
     RhoRight(:) = RhoRight(:) +   RhoSRight(:,iSpecies)
  enddo 
  ! End RhoS, Rho Facevalues

  ! Begin HydroRhoS, HydroRho Facevalues:
  ! Note: Hydro denotes hydrostatic background
  do iSpecies = 1, nSpecies
     call calc_tvd_facevalues(LogHydroRhoS(-1:nAlts+2,iSpecies), &
                          LogHydroRhoSLeft( 0:nAlts  ,iSpecies), &
                         LogHydroRhoSRight( 0:nAlts  ,iSpecies) )
      HydroRhoSLeft(:,iSpecies) = exp( LogHydroRhoSLeft(:,iSpecies))
     HydroRhoSRight(:,iSpecies) = exp(LogHydroRhoSRight(:,iSpecies))
  enddo 
  ! End HydroRhoS, HydroRho Facevalues

  ! Calculate the Left and Right Faces of the Temperatures 
  call calc_tvd_facevalues(Temp(-1:nAlts+2), &
                TemperatureLeft( 0:nAlts  ), &
               TemperatureRight( 0:nAlts  ) )

  PressureRight = 0.0
   PressureLeft = 0.0
  do iSpecies = 1, nSpecies
     PressureSLeft(0:nAlts,iSpecies) = &
         (RhoSLeft(0:nAlts,iSpecies)/Mass(iSpecies))*&
         Boltzmanns_Constant*TemperatureLeft(0:nAlts)
     PressureLeft = PressureLeft + PressureSLeft(:,iSpecies)

     HydroPressureSLeft(0:nAlts,iSpecies) = &
         (HydroRhoSLeft(0:nAlts,iSpecies)/Mass(iSpecies))*&
        Boltzmanns_Constant*TemperatureLeft(0:nAlts)

     PressureSRight(0:nAlts,iSpecies) = &
        (RhoSRight(0:nAlts,iSpecies)/Mass(iSpecies))*&
         Boltzmanns_Constant*TemperatureRight(0:nAlts)
     PressureRight = PressureRight + PressureSRight(:,iSpecies)

    HydroPressureSRight(0:nAlts,iSpecies) = &
        (HydroRhoSRight(0:nAlts,iSpecies)/Mass(iSpecies))*&
         Boltzmanns_Constant*TemperatureRight(0:nAlts)
  enddo 

  ! Begin Vertical, Horizontal Species and Bulk Winds
  do iSpecies = 1, nSpecies
    call calc_tvd_facevalues(VertVel(:,iSpecies), &
                        VertVelSLeft(:,iSpecies), &
                       VertVelSRight(:,iSpecies))
  enddo 

  do iDim = 1, 3
     call calc_tvd_facevalues(Vel_GD(:,iDim), &
                           VelGDLeft(:,iDim), &
                          VelGDRight(:,iDim))
  enddo 

   VelGDLeft(:,iUp_) = 0.0
  VelGDRight(:,iUp_) = 0.0
  do iSpecies = 1, nSpecies
      VelGDLeft(:,iUp_) = VelGDLeft(:,iUp_) + &
           RhoSLeft(:,iSpecies)*VertVelSLeft(:,iSpecies)/&
            RhoLeft(:)

      VelGDRight(:,iUp_) = VelGDRight(:,iUp_) + &
           RhoSRight(:,iSpecies)*VertVelSRight(:,iSpecies)/&
            RhoRight(:)
  enddo !iSpecies = 1, nSpecies

  call calc_tvd_facevalues(Gamma_1d(:), GammaLeft(:), GammaRight(:))

  do iAlt = 0, nAlts
     EnergyLeft(iAlt) = &
         ( 1.0/(GammaLeft(iAlt) - 1.0))*PressureLeft(iAlt) + &
           0.5*RhoLeft(iAlt)*&
          (VelGDLeft(iAlt,iUp_)**2.0 + VelGDLeft(iAlt,iEast_)**2.0 + &
           VelGDLeft(iAlt,iNorth_)**2.0)

     EnergyRight(iAlt) = &
         ( 1.0/(GammaRight(iAlt) - 1.0))*PressureRight(iAlt) + &
           0.5*RhoRight(iAlt)*&
          (VelGDRight(iAlt,iUp_)**2.0 + VelGDRight(iAlt,iEast_)**2.0 + &
           VelGDRight(iAlt,iNorth_)**2.0)
  enddo 

  do iAlt = 0, nAlts 
     EnthalpyLeft(iAlt) = &
       0.5*(VelGDLeft(iAlt,iUp_)**2.0 + &
            VelGDLeft(iAlt,iEast_)**2.0 + &
            VelGDLeft(iAlt,iNorth_)**2.0) + &
          (GammaLeft(iAlt)/(GammaLeft(iAlt) - 1.0))*&
          PressureLeft(iAlt)/RhoLeft(iAlt)

     EnthalpyRight(iAlt) = &
       0.5*(VelGDRight(iAlt,iUp_)**2.0 + &
            VelGDRight(iAlt,iEast_)**2.0 + &
            VelGDRight(iAlt,iNorth_)**2.0) + &
          (GammaRight(iAlt)/(GammaRight(iAlt) - 1.0))*&
          PressureRight(iAlt)/RhoRight(iAlt)

     do iSpecies = 1, nSpecies

        EnthalpySLeft(iAlt,iSpecies) = &
          0.5*(VertVelSLeft(iAlt,iSpecies)**2.0 + &
               VelGDLeft(iAlt,iEast_)**2.0 + &
               VelGDLeft(iAlt,iNorth_)**2.0) + &
             (GammaLeft(iAlt)/(GammaLeft(iAlt) - 1.0))*&
             PressureSLeft(iAlt,iSpecies)/RhoSLeft(iAlt,iSpecies)

        EnthalpySRight(iAlt,iSpecies) = &
          0.5*(VertVelSRight(iAlt,iSpecies)**2.0 + &
               VelGDRight(iAlt,iEast_)**2.0 + &
               VelGDRight(iAlt,iNorth_)**2.0) + &
             (GammaRight(iAlt)/(GammaRight(iAlt) - 1.0))*&
             PressureSRight(iAlt,iSpecies)/RhoSRight(iAlt,iSpecies)

     enddo !iSpecies = 1, nSpecies
  enddo 

  do iAlt = 0, nAlts
     SubCs = &
        sqrt(2.0*( (GammaLeft(iAlt) - 1.0 )/(GammaLeft(iAlt) + 1.0)) *&
                 EnthalpyLeft(iAlt) )
     SoundSpeedLeft(iAlt) = (SubCs**2.0)/max(SubCs, VelGDLeft(iAlt,iUp_))

     ! Note the -1.0 factor in the Right face
     SubCs = &
        sqrt(2.0*( (GammaRight(iAlt) - 1.0 )/(GammaRight(iAlt) + 1.0)) *&
                 EnthalpyRight(iAlt) )
     SoundSpeedRight(iAlt) = (SubCs**2.0)/max(SubCs, -1.0*VelGDRight(iAlt,iUp_))

     InterfaceSoundSpeed(iAlt) = min(SoundSpeedLeft(iAlt), SoundSpeedRight(iAlt))
     do iSpecies = 1, nSpecies
        SubCs = &
           sqrt(2.0*( (GammaLeft(iAlt) - 1.0 )/(GammaLeft(iAlt) + 1.0)) *&
                    EnthalpySLeft(iAlt,iSpecies) )
        SoundSpeedSLeft(iAlt,iSpecies) = &
           (SubCs**2.0)/max(SubCs, VertVelSLeft(iAlt,iSpecies))

        SubCs = &
           sqrt(2.0*( (GammaRight(iAlt) - 1.0 )/(GammaRight(iAlt) + 1.0)) *&
                    EnthalpySRight(iAlt,iSpecies) )
        ! Note the -1.0 factor in the Right face
        SoundSpeedSRight(iAlt,iSpecies) = &
           (SubCs**2.0)/max(SubCs, -1.0*VertVelSRight(iAlt,iSpecies))
        InterfaceSoundSpeedS(iAlt,iSpecies) = &
           min(SoundSpeedSLeft(iAlt,iSpecies), SoundSpeedSRight(iAlt,iSpecies))
     enddo !iSpecies = 1, nSpecies
  enddo !iAlt = 0, nAlts

  do iSpecies = 1, nSpecies
     do iAlt = 0, nAlts
        ! Implement Thornber Entropy Correction at low Mach Numbers
        ! see. Thornber et al. [2008], J. Comp. Phys, 227, 4873-4894
        ! Implementation here taken from Houim et al. [2011], J Comp.Phys, 230, 8527-8553
        CBar = InterfaceSoundSpeedS(iAlt,iSpecies)
        VL   =         VertVelSLeft(iAlt,iSpecies)
        VR   =        VertVelSRight(iAlt,iSpecies)
        ML   = sqrt(VL*VL)/CBar
        MR   = sqrt(VR*VR)/CBar
        ZVar = min(1.0, max(ML,MR))
         VertVelSLeft(iAlt,iSpecies) = 0.5*(VL + VR) + 0.5*ZVar*(VL - VR)
        VertVelSRight(iAlt,iSpecies) = 0.5*(VL + VR) + 0.5*ZVar*(VR - VL)
     enddo !iAlt = 1, nAlts
  enddo !iSpecies = 1, nSpecies

  MeanPressureS = 0.5*(PressureSLeft + PressureSRight)
  MeanHydroPressureS = 0.5*(HydroPressureSLeft + HydroPressureSRight)
  MeanRhoS = 0.5*(RhoSLeft + RhoSRight)

  do iAlt = 0, nAlts 
    do iSpecies = 1, nSpecies
       MachSLeft(iAlt,iSpecies) = &
         VertVelSLeft(iAlt,iSpecies)/InterfaceSoundSpeedS(iAlt,iSpecies)
       MachSRight(iAlt,iSpecies) = &
         VertVelSRight(iAlt,iSpecies)/InterfaceSoundSpeedS(iAlt,iSpecies)
    enddo !iSpecies = 1, nSpecies
  enddo !iAlt = 1, nAlts 

  do iAlt = 0, nAlts
     do iSpecies = 1, nSpecies
        MachSBar2(iAlt,iSpecies) = &
           0.5*( MachSLeft(iAlt,iSpecies)**2.0 + &
                MachSRight(iAlt,iSpecies)**2.0 )
        MachSZero2(iAlt,iSpecies) = &
           min(1.0, max(MachSBar2(iAlt,iSpecies), MInf))
        MachSZero(iAlt,iSpecies) = sqrt(MachSZero2(iAlt,iSpecies))
     enddo !iSpecies = 1, nSpecies
  enddo !iAlt = 1, nAlts

  ! Create Mach-Scaling Pressure-Correction Functions
  ! Establish Mach Number Scaling Functions g(M)
  ! We specify a g(M) such that it goes away as |M|->1
  ! g(M) ~ M^2 as M -> 0.0
  Kp(1:nSpecies) = 0.10             !! Ullrich et al. [2011]
  Ku(1:nSpecies) = 0.25             !! Ullrich et al. [2011]

  do iAlt = 0, nAlts
     do iSpecies = 1, nSpecies
        ! Note that MZero is scaled to be within [0,1]
        ! Specify pressure correction fluxes (LiouKpS)
        ! Note: MZero -> [0,1]
        ! Note: M2Zero  = (MZero)^2.0-> [0,1]
        ! As M -> 0.0, MachScaling ~ M^2 - M^4 ~ M^2 -> Varies like a pressure at Low Mach
        ! As M >> 1.0, MZero = 1.0, MachScaling ~ 1 - 1 -> 0.0
        !    This shuts our correction down at high mach number
        ! Mathematically, we want AUSM+-up Numerical Flux(Mp):
        !        Mp ~ d(P)*/(rho*cs*(fa*cs + zeta*dz/dt)) ~ dP/rho*cs^2

        MachScaling = MachSZero(iAlt,iSpecies)**2.0
        KpSParam(iAlt,iSpecies) = MachScaling*Kp(iSpecies)

        ! Note that MZero ranges from [0, 1.0] (cap at 1.0)
        MachScaling = MachSZero(iAlt,iSpecies)*(1.0 - MachSZero(iAlt,iSpecies))
        KuSParam(iAlt,iSpecies) = MachScaling*Ku(iSpecies)
     enddo 
  enddo 

  ! BEGIN Mach Functions 
  ! M(1)
  do iAlt = 0, nAlts
     do iSpecies = 1, nSpecies
        MachSFn1PLeft(iAlt,iSpecies) = &
            0.5*(MachSLeft(iAlt,iSpecies) + abs(MachSLeft(iAlt,iSpecies)) )
        MachSFn1PRight(iAlt,iSpecies) = &
            0.5*(MachSRight(iAlt,iSpecies) + abs(MachSRight(iAlt,iSpecies)) )

        MachSFn1NLeft(iAlt,iSpecies) = &
            0.5*(MachSLeft(iAlt,iSpecies) - abs(MachSLeft(iAlt,iSpecies)) )
        MachSFn1NRight(iAlt,iSpecies) = &
            0.5*(MachSRight(iAlt,iSpecies) - abs(MachSRight(iAlt,iSpecies)) )

     enddo !iSpecies = 1, nSpecies
  enddo !iAlt = 1, nAlts
  ! M(2)
  do iAlt = 0, nAlts
     do iSpecies = 1, nSpecies
         MachSFn2PLeft(iAlt,iSpecies) =  0.25*( MachSLeft(iAlt,iSpecies) + 1.0)**2.0
        MachSFn2PRight(iAlt,iSpecies) =  0.25*(MachSRight(iAlt,iSpecies) + 1.0)**2.0

         MachSFn2NLeft(iAlt,iSpecies) = -0.25*( MachSLeft(iAlt,iSpecies) - 1.0)**2.0
        MachSFn2NRight(iAlt,iSpecies) = -0.25*(MachSRight(iAlt,iSpecies) - 1.0)**2.0
     enddo 
  enddo 
  ! M(4)
  do iAlt = 0, nAlts
    do iSpecies = 1, nSpecies

      if ( abs(MachSLeft(iAlt,iSpecies)) .ge. 1.0) then 
         MachSFn4PLeft(iAlt,iSpecies) = MachSFn1PLeft(iAlt,iSpecies)

         ! These are the Slau2 versions (alpha = 0.0)
         PsiP_Left(iAlt,iSpecies) = &
            0.5*(1.0 + abs(MachSLeft(iAlt,iSpecies))/MachSLeft(iAlt,iSpecies) )
      else
         MachSFn4PLeft(iAlt,iSpecies) = MachSFn2PLeft(iAlt,iSpecies)*&
                   (1.0 - 16.0*LiouBeta*MachSFn2NLeft(iAlt,iSpecies))

         ! These are the Slau2 versions (alpha = 0.0)
         PsiP_Left(iAlt,iSpecies) = &
            0.25*( (MachSLeft(iAlt,iSpecies) + 1.0)**2.0)*(2.0 - MachSLeft(iAlt,iSpecies)) + &
                  LiouAlpha*MachSLeft(iAlt,iSpecies)*(MachSLeft(iAlt,iSpecies)**2.0 - 1.0)**2.0 
      endif 
!
      if ( abs(MachSRight(iAlt,iSpecies)) .ge. 1.0) then 
         MachSFn4NRight(iAlt,iSpecies) = MachSFn1NRight(iAlt,iSpecies)
         ! These are the Slau2 versions (alpha = 0.0)
         PsiN_Right(iAlt,iSpecies) = &
            0.5*(1.0 + abs(MachSRight(iAlt,iSpecies))/MachSRight(iAlt,iSpecies))
      else
         MachSFn4NRight(iAlt,iSpecies) = MachSFn2NRight(iAlt,iSpecies)*&
                   (1.0 + 16.0*LiouBeta*MachSFn2PRight(iAlt,iSpecies))
         ! These are the Slau2 versions (alpha = 0.0)
         !PsiN_Right(iAlt,iSpecies) = &
         !   0.25*( (MachSRight(iAlt,iSpecies) - 1.0)**2.0)*(2.0 + MachSRight(iAlt,iSpecies)) 
         PsiN_Right(iAlt,iSpecies) = &
            0.25*( (MachSRight(iAlt,iSpecies) - 1.0)**2.0)*(2.0 + MachSRight(iAlt,iSpecies)) - &
                  LiouAlpha*MachSRight(iAlt,iSpecies)*(MachSRight(iAlt,iSpecies)**2.0 - 1.0)**2.0 
      endif 
    enddo 
  enddo 

  !!! Begin 2nd Order Mach Number Polynomials
  !real, dimension(0:nAlts,1:nSpecies) ::  MachPressureSParamAUSM
  do iAlt = 0, nAlts 
   do iSpecies = 1, nSpecies 
      ! Note that MPress varies as Kp*(g(M))* dPs/(rhos*as^2)
      ! g(M) is a "shut off" for M > 1
      ! in all cases rhos*as^2 ~ rhos*(Gamma*Ps/rhos) ~ Gamma*Ps 
      ! now:  In the vertical direction
      !
      ! dPs/Ps ~ (1/Ps)*(dP/dr)*dr ~ -rhos*(grav)/(ns*kb*T)*dr  
      ! Thus as M -> 0 and P approaches hydrostatic
      ! Therefore MPress ~ Kp*ms*g*dr/(kb*T) (no units)
      ! But we want this to vary as M^2 as a pressure flux as M->0
      ! So we want APC Scalar as f(M) ~ Ms^2
      ! Mp_APC ~ f(M)*Kp*dPs/(rhos*cs^2) -> Ms^2*Kp*(dr/Hs)
!!!!! Hydrostatic Species Use This one

      AUSMMachPressureParamS(iAlt,iSpecies) = &
          KpSParam(iAlt,iSpecies)*max( (1.0 - MachSBar2(iAlt,iSpecies)), 0.0)*&
         ( (PressureSRight(iAlt, iSpecies) - HydroPressureSRight(iAlt,iSpecies) ) - &
           ( PressureSLeft(iAlt, iSpecies) -  HydroPressureSLeft(iAlt,iSpecies) ) )/&
         ( MeanRhoS(iAlt,iSpecies)*InterfaceSoundSpeedS(iAlt,iSpecies)**2.0)

   enddo !iSpecies = 1, nSpecies 
  enddo !iAlt = 0, nAlts 

  do iAlt = 0, nAlts 
     do iSpecies = 1, nSpecies
        NumericalVertVelS(iAlt,iSpecies) = &
          InterfaceSoundSpeedS(iAlt,iSpecies)*&
          ( MachSFn4PLeft(iAlt,iSpecies) + &  
           MachSFn4NRight(iAlt,iSpecies) - &  
          AUSMMachPressureParamS(iAlt,iSpecies))
     enddo !iSpecies = 1, nSpecies
  enddo !iAlt = 0, nAlts 

  do iAlt = 0, nAlts
     do iSpecies = 1, nSpecies
        NumericalPressureS(iAlt,iSpecies) = &
             PsiP_Left(iAlt,iSpecies)*PressureSLeft(iAlt,iSpecies) +  &
            PsiN_Right(iAlt,iSpecies)*PressureSRight(iAlt,iSpecies) - &
            KuSParam(iAlt,iSpecies)*PsiP_Left(iAlt,iSpecies)*PsiN_Right(iAlt,iSpecies)*&
            (RhoSLeft(iAlt,iSpecies) + RhoSRight(iAlt,iSpecies))*&
             InterfaceSoundSpeedS(iAlt,iSpecies)*&
            (VertVelSLeft(iAlt,iSpecies) - VertVelSRight(iAlt,iSpecies))

        LogNumericalPressureS(iAlt,iSpecies) = dlog(NumericalPressureS(iAlt,iSpecies))

        NumericalHydroPressureS(iAlt,iSpecies) = &
             PsiP_Left(iAlt,iSpecies)*HydroPressureSLeft(iAlt,iSpecies) +  &
            PsiN_Right(iAlt,iSpecies)*HydroPressureSRight(iAlt,iSpecies) 

        LogNumericalHydroPressureS(iAlt,iSpecies) = dlog(NumericalHydroPressureS(iAlt,iSpecies))

     enddo !iSpecies = 1, nSpecies
  enddo !iAlt = 1, nAlts

  do iAlt = 1, nAlts
    do iSpecies = 1, nSpecies
       CellCenteredLogNumericalPressureS(iAlt,iSpecies) = &
           0.5*( LogNumericalPressureS(iAlt  ,iSpecies) + & 
                 LogNumericalPressureS(iAlt-1,iSpecies) )
       CellCenteredLogNumericalHydroPressureS(iAlt,iSpecies) = &
           0.5*( LogNumericalHydroPressureS(iAlt  ,iSpecies) + & 
                 LogNumericalHydroPressureS(iAlt-1,iSpecies) )

       CellCenteredNumericalPressureS(iAlt,iSpecies) = &
           exp(CellCenteredLogNumericalPressureS(iAlt,iSpecies))
       CellCenteredNumericalHydroPressureS(iAlt,iSpecies) = &
           exp(CellCenteredLogNumericalHydroPressureS(iAlt,iSpecies))
    enddo !iSpecies = 1, nSpecies
  enddo 
! ENDTVD STUFF
!
!  do iAlt = 0, nAlts
!     do iSpecies = 1, nSpecies
!! AUSM+-up Version
!         NonHydrostaticPressureSLeft  =  PressureSLeft(iAlt,iSpecies) -  HydroPressureSLeft(iAlt,iSpecies)
!        NonHydrostaticPressureSRight  = PressureSRight(iAlt,iSpecies) - HydroPressureSRight(iAlt,iSpecies)
!
!        NumericalPressureS(iAlt,iSpecies) = &
!             PsiP_Left(iAlt,iSpecies)*NonHydrostaticPressureSLeft +  &
!            PsiN_Right(iAlt,iSpecies)*NonHydrostaticPressureSRight - &
!            KuSParam(iAlt,iSpecies)*PsiP_Left(iAlt,iSpecies)*PsiN_Right(iAlt,iSpecies)*&
!            (RhoSLeft(iAlt,iSpecies) + RhoSRight(iAlt,iSpecies))*&
!             InterfaceSoundSpeedS(iAlt,iSpecies)*&
!            (VertVelSLeft(iAlt,iSpecies) - VertVelSRight(iAlt,iSpecies))
!!
!     enddo !iSpecies = 1, nSpecies
!  enddo !iAlt = 1, nAlts
! ENDTVD STUFF

  !BEGINTVDFLUX 
  do iSpecies = 1, nSpecies
     do iAlt = 0, nAlts 

        ! Note if V > 0.0, then |V| =  V:  (1/2)*(V + |V|) = V, (1/2)*(V - |V|) = 0.0
        ! Note if V < 0.0, then |V| = -V:  (1/2)*(V + |V|) = V, (1/2)*(V - |V|) = 0.0
        RhoSFlux(iAlt,iSpecies) = &
           0.5* RhoSLeft(iAlt,iSpecies)*&
                (     NumericalVertVelS(iAlt,iSpecies) +  & 
                  dabs(NumericalVertVelS(iAlt,iSpecies)) ) + &
           0.5*RhoSRight(iAlt,iSpecies)*&
                (     NumericalVertVelS(iAlt,iSpecies) -  & 
                  dabs(NumericalVertVelS(iAlt,iSpecies)) ) 

        MomentumSFlux(iAlt,iSpecies) = &
           0.5* RhoSLeft(iAlt,iSpecies)*VertVelSLeft(iAlt,iSpecies)*&
                (     NumericalVertVelS(iAlt,iSpecies) +    &
                  dabs(NumericalVertVelS(iAlt,iSpecies)) ) + &
           0.5*RhoSRight(iAlt,iSpecies)*VertVelSRight(iAlt,iSpecies)*&
                (     NumericalVertVelS(iAlt,iSpecies) -    &
                  dabs(NumericalVertVelS(iAlt,iSpecies)) ) 
     enddo 
  enddo 

  BulkNumericalVertVel(0:nAlts) = 0.0
  do iSpecies = 1, nSpecies
     BulkNumericalVertVel = &
     BulkNumericalVertVel + &
           ( RhoSLeft(:,iSpecies) + RhoSRight(:,iSpecies) )*&
             NumericalVertVelS(:,iSpecies)/&
             (RhoLeft + RhoRight )
  enddo 

  do iAlt = 0, nAlts 
     EnergyFlux(iAlt) = &
        0.5*( EnergyLeft(iAlt) + PressureLeft(iAlt) )* &
            ( BulkNumericalVertVel(iAlt) +  &
          abs(BulkNumericalVertVel(iAlt)) ) + & 
        0.5*( EnergyRight(iAlt) + PressureRight(iAlt) )* &
            ( BulkNumericalVertVel(iAlt) -  &
          abs(BulkNumericalVertVel(iAlt)) ) 

     MomentumFlux(iAlt,1:3) = &
        0.5*( RhoLeft(iAlt)*VelGDLeft(iAlt,1:3))*&
            ( BulkNumericalVertVel(iAlt) +  &
          abs(BulkNumericalVertVel(iAlt)) ) + & 
        0.5*( RhoRight(iAlt)*VelGDRight(iAlt,1:3))*&
            ( BulkNumericalVertVel(iAlt) -  &
          abs(BulkNumericalVertVel(iAlt)) ) 

  enddo 
  !ENDTVDFLUX


  !BEGINDIVFLUX 
  ! Define Cell Volumes and radial Areas
  do iAlt = 1, nAlts
     AreaFunction_P12(iAlt) = Area_P12(iAlt)
     AreaFunction_M12(iAlt) = Area_M12(iAlt)
     LocalCellVolume(iAlt) = CellVol1D(iAlt)
  enddo 

  ! Divergence of Fluxes
  do iAlt = 1, nAlts
     do iSpecies = 1, nSpecies
        DivRhoSFlux(iAlt,iSpecies) = &
          ( AreaFunction_P12(iAlt)*RhoSFluxWENO(iAlt  ,iSpecies) - &
            AreaFunction_M12(iAlt)*RhoSFluxWENO(iAlt-1,iSpecies) )/&
            LocalCellVolume(iAlt)

        DivMomentumSFlux(iAlt,iSpecies) = &
          ( AreaFunction_P12(iAlt)*MomentumSFluxWENO(iAlt  ,iSpecies) - &
            AreaFunction_M12(iAlt)*MomentumSFluxWENO(iAlt-1,iSpecies) )/&
            LocalCellVolume(iAlt)

        ! Add PRessure Gradient
        ! Note that we want Del(Pressure') = Del(Pressure) - Del(PressureH)
        !  --> Del(P ) = dP/dr = P*d[Ln(P)]/dr
        !  --> Del(Ph) = dPh/dr = Ph*d[Ln(Ph)]/dr

        ! Step1 Add on the Del(P)
        DivMomentumSFlux(iAlt,iSpecies) = &
        DivMomentumSFlux(iAlt,iSpecies) + &
           CellCenteredNumericalPressureSWENO(iAlt,iSpecies)*&
           (LogNumericalPressureSWENO(iAlt  ,iSpecies) - &
            LogNumericalPressureSWENO(iAlt-1,iSpecies))/dAlt_C(iAlt)

        ! Subtract Away Hydrostatic Part
        DivMomentumSFlux(iAlt,iSpecies) = &
        DivMomentumSFlux(iAlt,iSpecies) - &
           CellCenteredNumericalHydroPressureSWENO(iAlt,iSpecies)*&
           (LogNumericalHydroPressureSWENO(iAlt  ,iSpecies) - &
            LogNumericalHydroPressureSWENO(iAlt-1,iSpecies))/dAlt_C(iAlt)
     enddo !iAlt


     DivMomentumFlux(iAlt,1:3) = &
          ( AreaFunction_P12(iAlt)*MomentumFluxWENO(iAlt  ,1:3) - &
            AreaFunction_M12(iAlt)*MomentumFluxWENO(iAlt-1,1:3) )/&
            LocalCellVolume(iAlt)
     DivEnergyFlux(iAlt) = &
          ( AreaFunction_P12(iAlt)*EnergyFluxWENO(iAlt  ) - &
            AreaFunction_M12(iAlt)*EnergyFluxWENO(iAlt-1) )/&
            LocalCellVolume(iAlt)

  enddo !iSpecies = 1, nSpecies

end subroutine calc_all_fluxes_hydro


 subroutine calc_tvd_facevalues(Var, VarLeft, VarRight)
   use ModVertical, only: dAlt_F, InvDAlt_F
   use ModSizeGITM, only: nAlts
   use ModLimiterGitm
 
   implicit none
   
   real, intent(in) :: Var(-1:nAlts+2)
   real, intent(out):: VarLeft(0:nAlts), VarRight(0:nAlts)
 
   real :: dVarUp, dVarDown, dVarLimited(0:nAlts+1)
 
   integer :: i
   ! Currently we have Var(i) corresponds to the interface between (i) and (i+1)
   !   |         |                  |                          |                  |
   !   |  UL(0)->|<-UR(0); UL(1) -> |<- UR(1) ....   UL(N-1)-> |<-UR(N-1); UL(N)->|<-UR(N) 
   !   |         |                  |              i=nAlts-1   |                  |
   !      i = 0        i = 1

   do i=0,nAlts+1
        dVarUp            = (Var(i+1) - Var(i  )) * InvDAlt_F(i+1)
        dVarDown          = (Var(i  ) - Var(i-1)) * InvDAlt_F(i  )
        dVarLimited(i) = Limiter_mc(dVarUp, dVarDown)
   end do

   do i = 0, nAlts
       VarLeft(i) = Var(i  ) + 0.5*dVarLimited(i  )*dAlt_F(i+1)
      VarRight(i) = Var(i+1) - 0.5*dVarLimited(i+1)*dAlt_F(i+1)
   enddo !i = 0, nAlts

 end subroutine calc_tvd_facevalues

 subroutine calc_weno_facevalues(Var, VarLeft, VarRight)
 
   use ModVertical, only: &
       dAlt_F, dAlt_C, InvDAlt_F, Altitude_G, &
       Mesh_ULP120, Mesh_ULP121, Mesh_ULP122, &
       Mesh_CLP120, Mesh_CLP121, Mesh_CLP122, &
       Mesh_URM120, Mesh_URM121, Mesh_URM122, &
       Mesh_CRM120, Mesh_CRM121, Mesh_CRM122, &
       UB_MeshCoefs, LB_MeshCoefs, &
       Mesh_IS0, Mesh_IS1, Mesh_IS2
     
   use ModSizeGITM, only: nAlts
   use ModLimiterGitm
 
   implicit none
   
   real, intent(in) :: Var(-1:nAlts+2)
   real, intent(out):: VarLeft(0:nAlts), VarRight(0:nAlts)
 
   real :: UL_P120, UL_P121, UL_P122
   real :: UR_M120, UR_M121, UR_M122
   real :: CL_P120, CL_P121, CL_P122
   real :: CR_M120, CR_M121, CR_M122
   real :: X_P12, X_P32, X_P52
   real :: X_M12, X_M32, X_M52
   real :: IS0, IS1, IS2
   real :: IS0Z, IS1Z, IS2Z
   real :: TRis
   real :: Alpha_P0, Alpha_P1, Alpha_P2
   real :: Alpha_M0, Alpha_M1, Alpha_M2
   real :: W_P0, W_P1, W_P2
   real :: W_M0, W_M1, W_M2
   integer :: i, k,iAlt
   real :: RISk
   real :: WENOEpsilon
   !------------ USE THESE FOR THE FAST VERSION --------------
   real :: C_P120, C_P121, C_P122
   real :: C_M120, C_M121, C_M122
   !-------
   real :: DenominatorL, DenominatorR

!!! Use for 4-th Order Backward Differences
!!! Need a 5-point Stencil
  real :: dVar
  real :: LocalVar(-2:nAlts+3)  ! Need this for extrapolation
  real :: LocalVarLeft(0:nAlts), LocalVarRight(0:nAlts)
!  real :: LocalVarLeft_M12(0:nAlts+1), LocalVarRight_M12(0:nAlts+1)
!  real :: LocalVarLeft_P12(0:nAlts+1), LocalVarRight_P12(0:nAlts+1)
  ! WENOZ  - Factors
  real :: Tau5Z
  LocalVar(-1:nAlts+2) = Var(-1:nAlts+2)
  ! =\ 
  ! ==\
  ! Extend Local Vars upward and downward
  iAlt = -2
  dVar  = LB_MeshCoefs(1,1)*Var(iAlt+1) + &  
          LB_MeshCoefs(1,2)*Var(iAlt+2) + &  
          LB_MeshCoefs(1,3)*Var(iAlt+3) + &  
          LB_MeshCoefs(1,4)*Var(iAlt+4) + &  
          LB_MeshCoefs(1,5)*Var(iAlt+5)      
  LocalVar(iAlt) = Var(iAlt+1) - dAlt_F(iAlt+1)*dVar 

  iAlt = nAlts + 3
  dVar  = UB_MeshCoefs(3,1)*Var(iAlt-1) + &  
          UB_MeshCoefs(3,2)*Var(iAlt-2) + &  
          UB_MeshCoefs(3,3)*Var(iAlt-3) + &  
          UB_MeshCoefs(3,4)*Var(iAlt-4) + &  
          UB_MeshCoefs(3,5)*Var(iAlt-5)      
  LocalVar(iAlt) = LocalVar(iAlt-1) + dVar*dAlt_F(iAlt-1)

  ! Use WENO Reconstruction
  !WENOEpsilon = 1.0e-6
  WENOEpsilon = 1.0e-2
  do iAlt=1,nAlts
      UL_P120 = &
        Mesh_ULP120(iAlt,1)*   & 
                LocalVar(iAlt+1)  + &
        Mesh_ULP120(iAlt,2)*   &
                (LocalVar(iAlt) - LocalVar(iAlt+1)) - &
        Mesh_ULP120(iAlt,3)*   &
                (LocalVar(iAlt+2) - LocalVar(iAlt+1))
!          
      UL_P121 = &
        Mesh_ULP121(iAlt,1)*   & 
                   LocalVar(iAlt  ) + &
        Mesh_ULP121(iAlt,2)*   &
               (LocalVar(iAlt+1) - LocalVar(iAlt  )) - &
        Mesh_ULP121(iAlt,3)*   &
                (LocalVar(iAlt-1) - LocalVar(iAlt  )) 

      UL_P122 = &
        Mesh_ULP122(iAlt,1)*   & 
                   LocalVar(iAlt-1) + &
        Mesh_ULP122(iAlt,2)*   &
               (LocalVar(iAlt-2) - LocalVar(iAlt-1)) + &
        Mesh_ULP122(iAlt,3)*   &
               (LocalVar(iAlt  ) - LocalVar(iAlt-1)) 
       
      UR_M120 = &
        Mesh_URM120(iAlt,1)*   & 
                   LocalVar(iAlt+1) + &
        Mesh_URM120(iAlt,2)*   &
                (LocalVar(iAlt  ) - LocalVar(iAlt+1)) + & 
        Mesh_URM120(iAlt,3)*   &
                (LocalVar(iAlt+2) - LocalVar(iAlt+1)) 

      UR_M121 = &
        Mesh_URM121(iAlt,1)*   & 
                   LocalVar(iAlt  ) + &
        Mesh_URM121(iAlt,2)*   &
                  (LocalVar(iAlt-1) - LocalVar(iAlt  )) - &
        Mesh_URM121(iAlt,3)*   &
                  (LocalVar(iAlt+1) - LocalVar(iAlt  )) 
       
      UR_M122 = &
        Mesh_URM122(iAlt,1)*   & 
                   LocalVar(iAlt-1) + &
        Mesh_URM122(iAlt,2)*   &
                  (LocalVar(iAlt) - LocalVar(iAlt-1)) - &
        Mesh_URM122(iAlt,3)*   &
                  (LocalVar(iAlt-2) - LocalVar(iAlt-1)) 

      CL_P120 = Mesh_CLP120(iAlt)
      CL_P121 = Mesh_CLP121(iAlt)
      CL_P122 = Mesh_CLP122(iAlt)

      CR_M120 = Mesh_CRM120(iAlt)
      CR_M121 = Mesh_CRM121(iAlt)
      CR_M122 = Mesh_CRM122(iAlt)

      IS0 = Mesh_IS0(iAlt,1)*&
            (LocalVar(iAlt+2) - LocalVar(iAlt+1))**2.0 + &
            Mesh_IS0(iAlt,2)*&
            ((LocalVar(iAlt+2) - LocalVar(iAlt+1))*(LocalVar(iAlt) - LocalVar(iAlt+1))) + &
            Mesh_IS0(iAlt,3)*&
            ((LocalVar(iAlt  ) - LocalVar(iAlt+1))**2.0)
            
      IS1 = Mesh_IS1(iAlt,1)*&
            (LocalVar(iAlt-1) - LocalVar(iAlt  ))**2.0 + &
            Mesh_IS1(iAlt,2)*&
            ((LocalVar(iAlt+1) - LocalVar(iAlt  ))*(LocalVar(iAlt-1) - LocalVar(iAlt))) + &
            Mesh_IS1(iAlt,3)*&
            ((LocalVar(iAlt+1) - LocalVar(iAlt  ))**2.0)

      IS2 = Mesh_IS2(iAlt,1)*&
            (LocalVar(iAlt-2) - LocalVar(iAlt-1))**2.0 + &
            Mesh_IS2(iAlt,2)*&
            ((LocalVar(iAlt  ) - LocalVar(iAlt-1))*(LocalVar(iAlt-2) - LocalVar(iAlt-1))) + &
            Mesh_IS2(iAlt,3)*&
            ((LocalVar(iAlt  ) - LocalVar(iAlt-1))**2.0)

      Tau5Z = abs(IS0 - IS2)

      Alpha_P0 = CL_P120*(1.0 + (Tau5Z/(IS0  + WENOEpsilon))**2.0)
      Alpha_P1 = CL_P121*(1.0 + (Tau5Z/(IS1  + WENOEpsilon))**2.0)
      Alpha_P2 = CL_P122*(1.0 + (Tau5Z/(IS2  + WENOEpsilon))**2.0)

      Alpha_M0 = CR_M120*(1.0 + (Tau5Z/(IS0  + WENOEpsilon))**2.0)
      Alpha_M1 = CR_M121*(1.0 + (Tau5Z/(IS1  + WENOEpsilon))**2.0)
      Alpha_M2 = CR_M122*(1.0 + (Tau5Z/(IS2  + WENOEpsilon))**2.0)

      DenominatorL = Alpha_P0 + Alpha_P1 + Alpha_P2
      DenominatorR = Alpha_M0 + Alpha_M1 + Alpha_M2

      W_P0 = Alpha_P0/DenominatorL
      W_P1 = Alpha_P1/DenominatorL
      W_P2 = Alpha_P2/DenominatorL

      W_M0 = Alpha_M0/DenominatorR
      W_M1 = Alpha_M1/DenominatorR
      W_M2 = Alpha_M2/DenominatorR

      ! For a given "i" we calculate the inner faces of the cell
      ! That is we interpolate to the right-face of the -1/2 edge
      ! That is we interpolate to the left-face  of the +1/2 edge
      ! Let's define
      !     U(i)          = LocalVar(i)
      !     UL_P12/M12(i) = LocalVarLeft_P12(i)/LocalVarLeft_M12(i)
      !     UR_P12/M12(i) = LocalVarRight_P12(i)/LocalVarRight_M12(i)
      !      (i-1)        |               (i)                |       (i+1)
      !                   | UR(i)  <----- U(i) ----->  UL(i) |
      !  Old Notation
      !                   | UR_M12(i)  <--U(i) --> UL_P12(i) |
      !   Old Scheme (i) ranged from i = 1, nAlts
      !
      !  New Notation
      !                   | UR(i-1)    <--U(i) --> UL(i)     |
      !   New Scheme (i) ranges from i = 0, nAlts
      !   In New Scheme (Special Cases):  i = 0:  

!       LocalVarLeft_P12(iAlt  ) = W_P0*UL_P120 + W_P1*UL_P121 + W_P2*UL_P122
!      LocalVarRight_M12(iAlt  ) = W_M0*UR_M120 + W_M1*UR_M121 + W_M2*UR_M122

       LocalVarLeft(iAlt  ) = W_P0*UL_P120 + W_P1*UL_P121 + W_P2*UL_P122
      LocalVarRight(iAlt-1) = W_M0*UR_M120 + W_M1*UR_M121 + W_M2*UR_M122
   enddo !i=1,nAlts

   ! Special Cases
   iAlt = 0
   UL_P120 = &
        Mesh_ULP120(iAlt,1)*   & 
                LocalVar(iAlt+1)  + &
        Mesh_ULP120(iAlt,2)*   &
                (LocalVar(iAlt) - LocalVar(iAlt+1)) - &
        Mesh_ULP120(iAlt,3)*   &
                (LocalVar(iAlt+2) - LocalVar(iAlt+1))
!          
   UL_P121 = &
        Mesh_ULP121(iAlt,1)*   & 
                   LocalVar(iAlt  ) + &
        Mesh_ULP121(iAlt,2)*   &
               (LocalVar(iAlt+1) - LocalVar(iAlt  )) - &
        Mesh_ULP121(iAlt,3)*   &
                (LocalVar(iAlt-1) - LocalVar(iAlt  )) 

   UL_P122 = &
        Mesh_ULP122(iAlt,1)*   &
                   LocalVar(iAlt-1) + &
        Mesh_ULP122(iAlt,2)*   &
               (LocalVar(iAlt-2) - LocalVar(iAlt-1)) + &
        Mesh_ULP122(iAlt,3)*   &
               (LocalVar(iAlt  ) - LocalVar(iAlt-1)) 
       
   CL_P120 = Mesh_CLP120(iAlt)
   CL_P121 = Mesh_CLP121(iAlt)
   CL_P122 = Mesh_CLP122(iAlt)

   IS0 = Mesh_IS0(iAlt,1)*&
         (LocalVar(iAlt+2) - LocalVar(iAlt+1))**2.0 + &
         Mesh_IS0(iAlt,2)*&
         ((LocalVar(iAlt+2) - LocalVar(iAlt+1))*(LocalVar(iAlt) - LocalVar(iAlt+1))) + &
         Mesh_IS0(iAlt,3)*&
         ((LocalVar(iAlt  ) - LocalVar(iAlt+1))**2.0)
            
   IS1 = Mesh_IS1(iAlt,1)*&
         (LocalVar(iAlt-1) - LocalVar(iAlt  ))**2.0 + &
         Mesh_IS1(iAlt,2)*&
         ((LocalVar(iAlt+1) - LocalVar(iAlt  ))*(LocalVar(iAlt-1) - LocalVar(iAlt))) + &
         Mesh_IS1(iAlt,3)*&
         ((LocalVar(iAlt+1) - LocalVar(iAlt  ))**2.0)

   IS2 = Mesh_IS2(iAlt,1)*&
         (LocalVar(iAlt-2) - LocalVar(iAlt-1))**2.0 + &
         Mesh_IS2(iAlt,2)*&
         ((LocalVar(iAlt  ) - LocalVar(iAlt-1))*(LocalVar(iAlt-2) - LocalVar(iAlt-1))) + &
         Mesh_IS2(iAlt,3)*&
         ((LocalVar(iAlt  ) - LocalVar(iAlt-1))**2.0)

   Tau5Z = abs(IS0 - IS2)

   Alpha_P0 = CL_P120*(1.0 + (Tau5Z/(IS0  + WENOEpsilon))**2.0)
   Alpha_P1 = CL_P121*(1.0 + (Tau5Z/(IS1  + WENOEpsilon))**2.0)
   Alpha_P2 = CL_P122*(1.0 + (Tau5Z/(IS2  + WENOEpsilon))**2.0)

   DenominatorL = Alpha_P0 + Alpha_P1 + Alpha_P2

   W_P0 = Alpha_P0/DenominatorL
   W_P1 = Alpha_P1/DenominatorL
   W_P2 = Alpha_P2/DenominatorL

   LocalVarLeft(iAlt  ) = W_P0*UL_P120 + W_P1*UL_P121 + W_P2*UL_P122


   iAlt = nAlts+1
   UR_M120 = &
     Mesh_URM120(iAlt,1)*   & 
                LocalVar(iAlt+1) + &
     Mesh_URM120(iAlt,2)*   &
             (LocalVar(iAlt  ) - LocalVar(iAlt+1)) + & 
     Mesh_URM120(iAlt,3)*   &
             (LocalVar(iAlt+2) - LocalVar(iAlt+1)) 

   UR_M121 = &
     Mesh_URM121(iAlt,1)*   & 
                LocalVar(iAlt  ) + &
     Mesh_URM121(iAlt,2)*   &
               (LocalVar(iAlt-1) - LocalVar(iAlt  )) - &
     Mesh_URM121(iAlt,3)*   &
               (LocalVar(iAlt+1) - LocalVar(iAlt  )) 
       
   UR_M122 = &
     Mesh_URM122(iAlt,1)*   & 
                LocalVar(iAlt-1) + &
     Mesh_URM122(iAlt,2)*   &
               (LocalVar(iAlt) - LocalVar(iAlt-1)) - &
     Mesh_URM122(iAlt,3)*   &
               (LocalVar(iAlt-2) - LocalVar(iAlt-1)) 

   CR_M120 = Mesh_CRM120(iAlt)
   CR_M121 = Mesh_CRM121(iAlt)
   CR_M122 = Mesh_CRM122(iAlt)

   IS0 = Mesh_IS0(iAlt,1)*&
         (LocalVar(iAlt+2) - LocalVar(iAlt+1))**2.0 + &
         Mesh_IS0(iAlt,2)*&
         ((LocalVar(iAlt+2) - LocalVar(iAlt+1))*(LocalVar(iAlt) - LocalVar(iAlt+1))) + &
         Mesh_IS0(iAlt,3)*&
         ((LocalVar(iAlt  ) - LocalVar(iAlt+1))**2.0)
            
   IS1 = Mesh_IS1(iAlt,1)*&
         (LocalVar(iAlt-1) - LocalVar(iAlt  ))**2.0 + &
         Mesh_IS1(iAlt,2)*&
         ((LocalVar(iAlt+1) - LocalVar(iAlt  ))*(LocalVar(iAlt-1) - LocalVar(iAlt))) + &
         Mesh_IS1(iAlt,3)*&
         ((LocalVar(iAlt+1) - LocalVar(iAlt  ))**2.0)

   IS2 = Mesh_IS2(iAlt,1)*&
         (LocalVar(iAlt-2) - LocalVar(iAlt-1))**2.0 + &
         Mesh_IS2(iAlt,2)*&
         ((LocalVar(iAlt  ) - LocalVar(iAlt-1))*(LocalVar(iAlt-2) - LocalVar(iAlt-1))) + &
         Mesh_IS2(iAlt,3)*&
         ((LocalVar(iAlt  ) - LocalVar(iAlt-1))**2.0)

   Tau5Z = abs(IS0 - IS2)

   Alpha_M0 = CR_M120*(1.0 + (Tau5Z/(IS0  + WENOEpsilon))**2.0)
   Alpha_M1 = CR_M121*(1.0 + (Tau5Z/(IS1  + WENOEpsilon))**2.0)
   Alpha_M2 = CR_M122*(1.0 + (Tau5Z/(IS2  + WENOEpsilon))**2.0)

   DenominatorR = Alpha_M0 + Alpha_M1 + Alpha_M2

   W_M0 = Alpha_M0/DenominatorR
   W_M1 = Alpha_M1/DenominatorR
   W_M2 = Alpha_M2/DenominatorR

   ! This isi LocalVarRight_P12(nAlts)
   LocalVarRight(iAlt-1) = W_M0*UR_M120 + W_M1*UR_M121 + W_M2*UR_M122

   do iAlt = 0, nAlts
       VarLeft(iAlt) = LocalVarLeft(iAlt)
      VarRight(iAlt) = LocalVarRight(iAlt)
   enddo 

 end subroutine !calc_weno_facevalues(Var, VarLeft, VarRight)


