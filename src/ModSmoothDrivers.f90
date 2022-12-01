!  Copyright (C) 2002 Regents of the University of Michigan, portions used with permission 
!  For more information, see http://csem.engin.umich.edu/tools/swmf

! JMB: 11/2022:
! Use these variables to smooth out drivers
module ModSmoothDrivers

  use ModGITM, only: nBlocksmax,nLons,nLats,nAlts,nIons,nSpecies,nSpeciesTotal
  use ModKind, only: Real8_

  implicit none

  ! Electric Potential Variables
  real ::  PreviousPotential(-1:nLons+2,-1:nLats+2,-1:nAlts+2,1:nBlocksMax)=0.0
  real :: PreviousPotentialY(-1:nLons+2,-1:nLats+2,-1:nAlts+2,1:nBlocksMax)=0.0
  real ::      NextPotential(-1:nLons+2,-1:nLats+2,-1:nAlts+2,1:nBlocksMax)=0.0
  real ::     NextPotentialY(-1:nLons+2,-1:nLats+2,-1:nAlts+2,1:nBlocksMax)=0.0
  real(Real8_)          :: LastPotentialTime
  real ::     tSimLastPotential, tSimNextPotential
  real ::     tSimCurrentPotential

  ! Auroral Variables: Original Variables in ModGITM.f90
  real :: PreviousElectronAverageEnergy(-1:nLons+2,-1:nLats+2,1:nBlocksMax)=0.0
  real ::     NextElectronAverageEnergy(-1:nLons+2,-1:nLats+2,1:nBlocksMax)=0.0
  real ::    PreviousElectronEnergyFlux(-1:nLons+2,-1:nLats+2,1:nBlocksMax)=0.0
  real ::        NextElectronEnergyFlux(-1:nLons+2,-1:nLats+2,1:nBlocksMax)=0.0
  ! These are Optional: Set in UAM.in
  real ::      PreviousIonAverageEnergy(-1:nLons+2,-1:nLats+2,1:nBlocksMax)=0.0
  real ::          NextIonAverageEnergy(-1:nLons+2,-1:nLats+2,1:nBlocksMax)=0.0
  real ::         PreviousIonEnergyFlux(-1:nLons+2,-1:nLats+2,1:nBlocksMax)=0.0
  real ::             NextIonEnergyFlux(-1:nLons+2,-1:nLats+2,1:nBlocksMax)=0.0
  ! Auroral Variables: Different Auroral Models
  real :: PreviousElectronEnergyFluxMono(-1:nLons+2,-1:nLats+2,1:nBlocksMax)=0.0
  real :: PreviousElectronEnergyFluxWave(-1:nLons+2,-1:nLats+2,1:nBlocksMax)=0.0
  real ::     NextElectronEnergyFluxMono(-1:nLons+2,-1:nLats+2,1:nBlocksMax)=0.0
  real ::     NextElectronEnergyFluxWave(-1:nLons+2,-1:nLats+2,1:nBlocksMax)=0.0
  real :: PreviousElectronNumberFluxMono(-1:nLons+2,-1:nLats+2,1:nBlocksMax)=0.0
  real :: PreviousElectronNumberFluxWave(-1:nLons+2,-1:nLats+2,1:nBlocksMax)=0.0
  real ::     NextElectronNumberFluxMono(-1:nLons+2,-1:nLats+2,1:nBlocksMax)=0.0
  real ::     NextElectronNumberFluxWave(-1:nLons+2,-1:nLats+2,1:nBlocksMax)=0.0
  real(Real8_)          :: LastAuroraTime
  real ::     tSimLastAurora, tSimNextAurora
  real ::     tSimCurrentAurora

  ! Electric Potential Variables
  real :: PreviousAuroralIonRateS(1:nLons,1:nLats,1:nAlts,1:nSpecies,1:nBlocksMax)=0.0
  real :: PreviousAuroralHeatingRate(1:nLons,1:nLats,1:nAlts,1:nBlocksMax)=0.0
  real ::        NextAuroralIonRateS(1:nLons,1:nLats,1:nAlts,1:nSpecies,1:nBlocksMax)=0.0
  real ::     NextAuroralHeatingRate(1:nLons,1:nLats,1:nAlts,1:nBlocksMax)=0.0
  real(Real8_)          :: LastAuroraIonHeatTime
  real ::     tSimLastAuroraIonHeat, tSimNextAuroraIonHeat
  real ::     tSimCurrentAuroraIonHeat
  real ::        FrozenTempUnit(1:nLons,1:nLats,1:nAlts,1:nBlocksMax)=0.0
  real ::              FrozenCp(1:nLons,1:nLats,1:nAlts,1:nBlocksMax)=0.0
  real ::             FrozenRho(1:nLons,1:nLats,1:nAlts,1:nBlocksMax)=0.0

  !EUV Variables
  real ::       PreviousEuvHeatingRate(1:nLons,1:nLats,1:nAlts,1:nBlocksMax)=0.0
  real :: PreviousPhotoElectronHeatingRate(1:nLons,1:nLats,1:nAlts,1:nBlocksMax)=0.0
  real ::          PreviouseEuvHeatingRate(1:nLons,1:nLats,1:nAlts,1:nBlocksMax)=0.0
    !allocate(EuvHeating(nLons, nLats, nAlts,nBlocks))
    !allocate(eEuvHeating(nLons, nLats, nAlts,nBlocks))
    !allocate(PhotoElectronHeating(nLons, nLats, nAlts,nBlocks))
  real ::  PreviousEuvIonRateS(1:nLons,1:nLats,1:nAlts,1:nIons,1:nBlocksMax)=0.0
  real :: PreviousEuvDissRateS(1:nLons,1:nLats,1:nAlts,1:nSpeciesTotal,1:nBlocksMax)=0.0
    !allocate(EuvIonRateS(nLons, nLats, nAlts, nIons,nBlocks))
    !allocate(EuvDissRateS(nLons, nLats, nAlts, nSpeciesTotal,nBlocks))

  real ::       NextEuvHeatingRate(1:nLons,1:nLats,1:nAlts,1:nBlocksMax)=0.0
  real :: NextPhotoElectronHeatingRate(1:nLons,1:nLats,1:nAlts,1:nBlocksMax)=0.0
  real ::          NexteEuvHeatingRate(1:nLons,1:nLats,1:nAlts,1:nBlocksMax)=0.0
  real ::  NextEuvIonRateS(1:nLons,1:nLats,1:nAlts,1:nIons,1:nBlocksMax)=0.0
  real :: NextEuvDissRateS(1:nLons,1:nLats,1:nAlts,1:nSpeciesTotal,1:nBlocksMax)=0.0
  real(Real8_)          :: LastEuvTime
  real ::     tSimLastEuv, tSimNextEuv
  real ::     tSimCurrentEuv


end module ModSmoothDrivers
