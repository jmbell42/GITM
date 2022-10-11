!  Copyright (C) 2002 Regents of the University of Michigan, portions used with permission 
!  For more information, see http://csem.engin.umich.edu/tools/swmf
!-----------------------------------------------------------------------------
! $Id: Earth.f90,v 1.21 2017/08/05 01:33:09 ridley Exp $
!
! Author: Aaron Ridley, UMichigan
!
! Modified: AGB Oct 2013 - Corrected spelling of photoelectron heating
!                          efficiency variable
!-----------------------------------------------------------------------------

subroutine fill_photo

  use ModPlanet
  use ModEUV

  implicit none

  integer :: iSpecies, iWave

  PhotoAbs = 0.0
  PhotoIon = 0.0
  PhotoDis = 0.0
  PhotoElecIon = 0.0
  PhotoElecDiss = 0.0

  night_photoion = 0.0
  night_photoabs = 0.0

  photoabs(:, iO_3P_) = PhotoAbs_O
  photoabs(:, iO2_) = PhotoAbs_O2
  photoabs(:, iN2_) = PhotoAbs_N2

  if (nSpecies > 3) then
     iSpecies = iN_4S_
     photoabs(:,min(iSpecies,nSpecies)) = PhotoIon_N
  endif

  ! JMB:  06/25/2016
  if (nSpecies > 5) then
     iSpecies = iHe_
     photoabs(:,min(iSpecies,nSpecies)) = PhotoAbs_He
  endif

  ! This may need to be as defined below....
  photoion(:,iN2P_) = PhotoIon_N2
  photoion(:,iO2P_) = PhotoIon_O2
  photoion(:,iNP_) = PhotoIon_N
  photoion(:,iO_4SP_) = PhotoIon_OPlus4S
  photoion(:,iO_2DP_) = PhotoIon_OPlus2D
  photoion(:,iO_2PP_) = PhotoIon_OPlus2P
  photoion(:,iHeP_) = PhotoAbs_He

  PhotoIonFrom(iN2P_) = iN2_
  PhotoIonFrom(iO2P_) = iO2_
  PhotoIonFrom(iNP_) = iN_4S_
  PhotoIonFrom(iO_4SP_) = iO_3P_
  PhotoIonFrom(iO_2DP_) = iO_3P_ 
  PhotoIonFrom(iO_2PP_) = iO_3P_ 
  PhotoIonFrom(iHeP_) = iHe_

  ! Photoelectrons:
  ! N2:
  ! PE Ratio:  N2 + e- -> N2+
  PhotoElecIon(:,iN2P_) = PhotoElec_N2_N2Plus

  ! O2:
  ! PE Ratio:  O2 + e- -> O2+
  PhotoElecIon(:,iO2P_) = PhotoElec_O2_O2Plus

  ! O:
  ! O + e* -> O(4S)+ + e
  PhotoElecIon(:,iO_4SP_) = PhotoElec_O_O4SPlus
  ! O + e* -> O(2D)+ + e
  PhotoElecIon(:,iO_2DP_) = PhotoElec_O_O2DPlus
  ! O + e* -> O(2P)+ + e
  PhotoElecIon(:,iO_2PP_) = PhotoElec_O_O2PPlus

  ! Dissociation:
  PhotoElecDiss(:,iN2_) = PhotoElec_N2_N4S
  PhotoElecDiss(:,iO2_) = PhotoElec_O2_O3P
  
  ! PE Ratio:  N2 + e- -> N(2D) + N(4S)
  pelecratio_N2(:,1) = PhotoElec_N2_N4S
  ! PE Ratio:  N2 + e- -> N+ + N(4S)
  pelecratio_N2(:,2) = PhotoElec_N2_NPlus
  ! PE Ratio:  N2 + e- -> N2+
  pelecratio_N2(:,3) = PhotoElec_N2_N2Plus

  ! PE Ratio:  O2 + e- -> O(4S) + O(3P)
  pelecratio_O2(:,1) = PhotoElec_O2_O3P
  ! PE Ratio:  O2 + e- -> O2+
  pelecratio_O2(:,2) = PhotoElec_O2_O2Plus
  ! PE Ratio:  O2 + e- -> O(4S)+ + O(3P)
  pelecratio_O2(:,3) = PhotoElec_O2_OPlus
  
  do iWave = 1, Num_WaveLengths_High
     if (waves(iWave) >= 1250.0 .and. wavel(iWave) <= 1750.0) then
        PhotoDis(iWave, iO2_) = &
             photoabs(iWave,iO2_) - PhotoIon(iWave, iO2P_)
     endif
     if (waves(iWave) >= 800.0 .and. wavel(iWave) <= 1250.0) then
        PhotoDis(iWave, iN2_) = &
             photoabs(iWave,iN2_) - PhotoIon(iWave, iN2P_)
     endif
  enddo

  ! Night time ionization:

  night_photoion(:,iN2P_) = Night_PhotoIon_N2 *1.e-4
  night_photoion(:,iO2P_) = Night_PhotoIon_O2 *1.e-4
  night_photoion(:,iNOP_) = Night_PhotoIon_NO *1.e-4
  night_photoion(:,iO_4SP_) = Night_PhotoIon_OPlus4S *1.e-4
  night_photoion(:,iO_2DP_) = Night_PhotoIon_OPlus2D *1.e-4
  night_photoion(:,iO_2PP_) = Night_PhotoIon_OPlus2P *1.e-4
  night_photoion(:,iNP_) = Night_PhotoIon_N *1.e-4

  night_photoabs(:,iN2_) = Night_PhotoAbs_N2 *1.e-4
  night_photoabs(:,iO2_) = Night_PhotoAbs_O2 *1.e-4
  night_photoabs(:,iNO_) = Night_PhotoAbs_NO *1.e-4
  night_photoabs(:,iO_3P_) = Night_PhotoAbs_O *1.e-4
  night_photoabs(:,iN_4S_) = Night_PhotoAbs_N *1.e-4
  
end subroutine fill_photo

subroutine calc_planet_sources(iBlock)

  use ModInputs
  use ModSources
  use ModEUV
  use ModGITM
  use ModTime
  
  implicit none

  integer, intent(in) :: iBlock

  integer :: iAlt, iError, iDir, iLat, iLon

  real :: tmp2(nLons, nLats, nAlts)
  real :: tmp3(nLons, nLats, nAlts)
  real :: Omega(nLons, nLats, nAlts)

  LowAtmosRadRate = 0.0

  !\
  ! Cooling ----------------------------------------------------------
  !/

  if (UseBarriers) call MPI_BARRIER(iCommGITM,iError)
  if (iDebugLevel > 4) write(*,*) "=====> NO cooling", iproc, UseNOCooling

  call calc_co2(iBlock)

  if (UseCO2Cooling) then

     call calc_co2_cooling(iBlock)
     ! fomichev version uses cool-to-space only
     !call calc_co2_cooling_fomichev(iBlock)

  else
     CO2Cooling = 0.0
  endif

  if (UseNOCooling) then

     call calc_no_cooling(iBlock)

  else

     NOCooling = 0.0

  endif

  if (UseBarriers) call MPI_BARRIER(iCommGITM,iError)
  if (iDebugLevel > 4) write(*,*) "=====> UseOCooling", iproc, UseOCooling

  if (UseOCooling) then 

     call calc_o3p_cooling(iBlock)

  else

     OCooling = 0.0

  endif

  RadCooling(1:nLons,1:nLats,1:nAlts,iBlock) = &
       OCooling + NOCooling + CO2Cooling

  PhotoElectronHeating(:,:,:,iBlock) = 0.0
  PhotoElectronHeating(:,:,:,iBlock) = &
       PhotoElectronHeatingEfficiency * &
       35.0*1.602e-19*&
       ( &
       EuvIonRateS(:,:,:,iO2P_,iBlock) + &
       EuvIonRateS(:,:,:,iN2P_,iBlock) + & !* &
       EuvIonRateS(:,:,:,iO_4SP_,iBlock) + &
       EuvIonRateS(:,:,:,iO_2DP_,iBlock) + &
       EuvIonRateS(:,:,:,iO_2PP_,iBlock))

  PhotoElectronHeating2d = 0.0
  do iAlt = 1, nAlts
     PhotoElectronHeating2d(1:nLons,1:nLats) = &
          PhotoElectronHeating2d(1:nLons,1:nLats) + &
          PhotoElectronHeating(:,:,iAlt,iBlock)  * &
          dAlt_GB(1:nLons,1:nLats,iAlt,iBlock)
  enddo
  
  PhotoElectronHeating(:,:,:,iBlock) = &
       PhotoElectronHeating(:,:,:,iBlock) / &
       Rho(1:nLons,1:nLats,1:nAlts,iBlock) / &
       cp(1:nLons,1:nLats,1:nAlts,iBlock) / &
       TempUnit(1:nLons,1:nLats,1:nAlts)

!--------------------------------------------------------------------
! GLOW
!--------------------------------------------------------------------

if (UseGlow) then
     if (dt < 10000.) then
        if  (floor((tSimulation-dt)/DtGlow) /= &
             floor(tsimulation/DtGlow)) then   

           call start_timing("glow")
           isInitialGlow = .True.

           if (iDebugLevel > 4) write(*,*) "=====> going into get_glow", iproc

           do iLat = 1, nLats
              do iLon = 1, nLons

                 call get_glow(iLon,iLat,iBlock)
                 
              enddo
           enddo

           call end_timing("glow")

        endif
     endif
     PhotoElectronDensity(:,:,:,:,iBlock) = PhotoElectronRate(:,:,:,:,iBlock) * dt
  endif


end subroutine calc_planet_sources

!---------------------------------------------------------------------
! Initialize Heating Rates
!---------------------------------------------------------------------

subroutine init_heating_efficiency

  use ModGITM, only: nLons, nLats, nAlts, nBlocks, Altitude_GB
  use ModEUV, only: HeatingEfficiency_CB, eHeatingEfficiency_CB
  use ModInputs, only: NeutralHeatingEfficiency

  implicit none

  integer :: iLon, iLat, iAlt
  !------------------------------------------------------------------

  HeatingEfficiency_CB(:,:,:,1:nBlocks) = NeutralHeatingEfficiency
!  max(0.1, &
!       0.40 - &
!       5.56e-5*(Altitude_GB(1:nLons,1:nLats,1:nAlts,1:nBlocks)/1000 - 165)**2)

  where(Altitude_GB(1:nLons,1:nLats,1:nAlts,1:nBlocks)/1000. > 150.)
     eHeatingEfficiency_CB(:,:,:,1:nBlocks) = 0.04
!!! min(0.4, &
!!!     0.04 + &
!!!     0.05*(Altitude_GB(1:nLons,1:nLats,1:nAlts,1:nBlocks)/1000 - 150)/100)
  elsewhere        
     eHeatingEfficiency_CB(:,:,:,1:nBlocks) = max(0.000001, &
          0.05 + &
          0.07*(Altitude_GB(1:nLons,1:nLats,1:nAlts,1:nBlocks)/1000 - 200)/100)
  end where

end subroutine init_heating_efficiency

!---------------------------------------------------------------------
! Calculate Eddy Diffusion Coefficient
!---------------------------------------------------------------------

subroutine calc_eddy_diffusion_coefficient(iBlock)

  use ModSizeGITM
  use ModGITM, only: pressure
  use ModInputs, only: EddyDiffusionPressure0,EddyDiffusionPressure1, &
       EddyDiffusionCoef
  use ModSources, only: KappaEddyDiffusion

  implicit none

  integer, intent(in) :: iBlock
  integer :: iAlt, iLat, iLon

  KappaEddyDiffusion=0.
  do iAlt = -1, nAlts+2

     do iLat = 1, nLats
        do iLon = 1, nLons

           if (pressure(iLon,iLat,iAlt,iBlock) >EddyDiffusionPressure0) then
              KappaEddyDiffusion(iLon,iLat,iAlt,iBlock) = EddyDiffusionCoef
              
           else if (pressure(iLon,iLat,iAlt,iBlock) > &
                EddyDiffusionPressure1) then

              KappaEddyDiffusion(iLon,iLat,iAlt,iBlock) = EddyDiffusionCoef * &
                   (pressure(iLon,iLat,iAlt,iBlock) - &
                   EddyDiffusionPressure1)/&
                   (EddyDiffusionPressure0 - EddyDiffusionPressure1)

           endif
        enddo
     enddo
  enddo

end subroutine calc_eddy_diffusion_coefficient

subroutine set_planet_defaults

  use ModPlanet
  use ModInputs

  implicit none
  
  iNeutralDensityOutputList(iN_4S_)=.false.
  iNeutralDensityOutputList(iHe_)=.false.
  iNeutralDensityOutputList(iN_2D_)=.false.
  iNeutralDensityOutputList(iN_2P_)=.false.
  iNeutralDensityOutputList(iH_)=.false.
  iNeutralDensityOutputList(iCO2_)=.false.
  iNeutralDensityOutputList(iO_1D_)=.false.

  iIonDensityOutputList(iNP_)=.false.
  iIonDensityOutputList(iO_2DP_)=.false.
  iIonDensityOutputList(iO_2PP_)=.false.
  iIonDensityOutputList(iHP_)=.false.
  iIonDensityOutputList(iHeP_)=.false.

  iTemperatureOutputList(2)=.false.
  iTemperatureOutputList(3)=.false.
  
  return

end subroutine set_planet_defaults


subroutine calc_co2_cooling(iBlock)

  use ModSources
  use ModEUV
  use ModGITM
  use ModTime
  
  implicit none

  integer, intent(in) :: iBlock

  integer :: iAlt, iError, iDir, iLat, iLon

  real :: tmp2(nLons, nLats, nAlts)
  real :: tmp3(nLons, nLats, nAlts)
  real :: Omega(nLons, nLats, nAlts)
  ! Updated 2022 by Jared Bell (JMB):
  ! Added more readily understood variables here
  real :: CO2_collisional_excitation(nLons, nLats, nAlts)
  real :: CO2_collisional_deexcitation(nLons, nLats, nAlts)
  real :: CO2_spontaneous_emission(nLons,nLats,nAlts)
  real :: true_temp(1:nLons,1:nLats,1:nAlts)
  real :: nCO2_vibrational(1:nLons,1:nLats,1:nAlts)
  real :: energy_CO2_emission
  real :: emission_wavelength
  ! Column Depth Efficiency
  real :: co2_column_density(1:nLons,1:nLats,1:nAlts)
  real :: co2_emission_cross_section
  real :: epsilon_co2(1:nLons,1:nLats,1:nAlts)

  
  ! [CO2] cooling 
  ! Use real temperature (in K)
  true_temp(1:nLons,1:nLats,1:nAlts) = &
        Temperature(1:nLons,1:nLats,1:nAlts,iBlock)*&
           TempUnit(1:nLons,1:nLats,1:nAlts)

!  ! First, determine what optical thickness we are in
  co2_emission_cross_section = 6.43e-19 ! m^2 15-micron absorption cross-section
!  ! Step1: Calculate the column of CO2 above a given altitude
!  ! Begin at exobase and integrate down:
!  !----
!  ! Assume that the column above the exobase is ~ n(CO2)*H(co2) in m^2
  co2_column_density(1:nLons,1:nLats,nAlts) = &
       NDensityS(1:nLons,1:nLats,nAlts,iCO2_,iBlock)*&
       (-1.0)*&
       (Boltzmanns_Constant*true_temp(1:nLons,1:nLats,nAlts)/&
       ( Mass(iCO2_)*Gravity_GB(1:nLons,1:nLats, nAlts,iBlock)))
!  !----
!  ! Next integrate downward to a specific altitude and add up the CO2 above you
!  ! Note that I just use the mean value between two cells:
  do iAlt = nAlts-1, 1, -1
     co2_column_density(1:nLons,1:nLats,iAlt) = &
     co2_column_density(1:nLons,1:nLats,iAlt+1) + &
       (0.5*(NDensityS(1:nLons,1:nLats,iAlt+1,iCO2_,iBlock) + &
             NDensityS(1:nLons,1:nLats,iAlt  ,iCO2_,iBlock)))*&
             dAlt_GB(1:nLons,1:nLats,iAlt,iBlock)
    
  enddo 

  ! Set emission efficiency factor: epsilon_co2
  ! Use Tabulated values quoted in Kumar and JAmes
  ! Adapted by Johnstone, C.P. [2018]
  where(co2_emission_cross_section*co2_column_density .ge. 2.0) 
     epsilon_co2=&
            0.7202*(co2_emission_cross_section*co2_column_density)**-0.613
  else where 
     ! Note, limit the lower value here
     epsilon_co2=&
             0.4732*(max(1.0e-10, &
             co2_emission_cross_section*co2_column_density))**-0.0069
  endwhere
  !-------
  ! Set the 15-micron cooling
  emission_wavelength = 15.0e-06  ! emission wavelength (IR) in meter
  energy_CO2_emission = Planck_Constant*Speed_Light/emission_wavelength ! Energy in J
  CO2_spontaneous_emission(1:nLons,1:nLats,1:nAlts) = 0.46  ! Einstien coef 12.54 s^-1

  ! Sum over all species [M] = O, O2, N2, CO2
  ! Based upon numerical fits to Siddles et al. [1994] and
  ! reported by Johnstone et al. [2018]
  CO2_collisional_deexcitation(1:nLons,1:nLats,1:nAlts) = &
         ( &
           ! CO2-O
           (5.10e-17)*(true_temp**-0.59)*&
           NDensityS(1:nLons,1:nLats,1:nAlts,iO_3P_,iBlock) + &
           ! CO2-O2
           (4.97e-28)*(true_temp**2.83)*&
           NDensityS(1:nLons,1:nLats,1:nAlts,iO2_,iBlock) + &
           ! CO2-N2
           (6.43e-27)*(true_temp**2.30)*&
           NDensityS(1:nLons,1:nLats,1:nAlts,iN2_,iBlock) + &
           ! CO2-N2
           (4.21e-23)*(true_temp**0.85)*&
           NDensityS(1:nLons,1:nLats,1:nAlts,iCO2_,iBlock) )

  ! From Castle et al. [2006]
  CO2_collisional_excitation(1:nLons,1:nLats,1:nAlts) = &
           2.0*CO2_collisional_deexcitation(1:nLons,1:nLats,1:nAlts)*&
            exp(-667.0/true_temp) 

  Omega = CO2_collisional_excitation/&
         (CO2_collisional_deexcitation + &
          CO2_collisional_excitation + &
          CO2_spontaneous_emission*epsilon_co2 )

  CO2Cooling = epsilon_co2*&
      ( energy_co2_emission*Omega*CO2_spontaneous_emission)*&
       NDensityS(1:nLons,1:nLats,1:nAlts,iCO2_,iBlock)

  CO2Cooling2d = 0.0
  do iAlt=1,nAlts
     RadiativeCooling2d(1:nLons, 1:nLats) = &
          RadiativeCooling2d(1:nLons, 1:nLats) + &
          CO2Cooling(1:nLons,1:nLats,iAlt) * &
          dAlt_GB(1:nLons,1:nLats,iAlt,iBlock)
     CO2Cooling2d = CO2Cooling2d + &
          CO2Cooling(1:nLons,1:nLats,iAlt) * &
          dAlt_GB(1:nLons,1:nLats,iAlt,iBlock)
  enddo
  CO2Cooling = CO2Cooling / TempUnit(1:nLons,1:nLats,1:nAlts) / &
       (Rho(1:nLons,1:nLats,1:nAlts,iBlock)*cp(:,:,1:nAlts,iBlock))

endsubroutine calc_co2_cooling

! This version can be enabled to engage the fomichev
! version of the CO2 Cooling--used in WACCM-X
! Note: This produces higher cooling rates
! Most likely due to the fact that we are not accounting
! for re-absorption of photons/NIR heating  in the formulation
! below--> it is inherently cool-to-space
subroutine calc_co2_cooling_fomichev(iBlock)

  use ModSources
  use ModEUV
  use ModGITM
  use ModTime
  
  implicit none

  integer, intent(in) :: iBlock

  integer :: iAlt, iError, iDir, iLat, iLon

  real :: tmp2(nLons, nLats, nAlts)
  real :: tmp3(nLons, nLats, nAlts)
  real :: Omega(nLons, nLats, nAlts)
  ! Updated 2022 by Jared Bell (JMB):
  ! Added more readily understood variables here
  real :: CO2_collisional_excitation(nLons, nLats, nAlts)
  real :: CO2_collisional_deexcitation(nLons, nLats, nAlts)
  real :: CO2_spontaneous_emission(nLons,nLats,nAlts)
  real :: true_temp(1:nLons,1:nLats,1:nAlts)
  real :: nCO2_vibrational(1:nLons,1:nLats,1:nAlts)
  real :: energy_CO2_emission
  real :: emission_wavelength
  ! Column Depth Efficiency
  real :: co2_column_density(1:nLons,1:nLats,1:nAlts)
  real :: co2_emission_cross_section
  real :: epsilon_co2(1:nLons,1:nLats,1:nAlts)

  
  ! [CO2] cooling 
  ! Use real temperature (in K)
  true_temp(1:nLons,1:nLats,1:nAlts) = &
        Temperature(1:nLons,1:nLats,1:nAlts,iBlock)*&
           TempUnit(1:nLons,1:nLats,1:nAlts)

!  ! First, determine what optical thickness we are in
  co2_emission_cross_section = 6.43e-19 ! m^2 15-micron absorption cross-section
!  ! Step1: Calculate the column of CO2 above a given altitude
!  ! Begin at exobase and integrate down:
!  !----
!  ! Assume that the column above the exobase is ~ n(CO2)*H(co2) in m^2
  co2_column_density(1:nLons,1:nLats,nAlts) = &
       NDensityS(1:nLons,1:nLats,nAlts,iCO2_,iBlock)*&
       (-1.0)*&
       (Boltzmanns_Constant*true_temp(1:nLons,1:nLats,nAlts)/&
       ( Mass(iCO2_)*Gravity_GB(1:nLons,1:nLats, nAlts,iBlock)))
!  !----
!  ! Next integrate downward to a specific altitude and add up the CO2 above you
!  ! Note that I just use the mean value between two cells:
  do iAlt = nAlts-1, 1, -1
     co2_column_density(1:nLons,1:nLats,iAlt) = &
     co2_column_density(1:nLons,1:nLats,iAlt+1) + &
       (0.5*(NDensityS(1:nLons,1:nLats,iAlt+1,iCO2_,iBlock) + &
             NDensityS(1:nLons,1:nLats,iAlt  ,iCO2_,iBlock)))*&
             dAlt_GB(1:nLons,1:nLats,iAlt,iBlock)
    
     ! Set emission efficiency factor: epsilon_co2
     ! if very thin, epsilon -> 
     epsilon_co2(1:nLons,1:nLats,iAlt) = &
             exp(-1.0*co2_column_density(1:nLons,1:nLats,iAlt)*&
                      co2_emission_cross_section)
  enddo 
  ! next, define epsilon = efficiency of photon escape
  ! can vary from 0 (super thick atmosphere) to 1.0 (completely thin)
  
  !-------
  ! Set the 15-micron cooling
  emission_wavelength = 15.0e-06  ! emission wavelength (IR) in meter
  energy_CO2_emission = Planck_Constant*Speed_Light/emission_wavelength ! Energy in J
  CO2_spontaneous_emission(1:nLons,1:nLats,1:nAlts) = 1.28  ! Einstien coef 12.54 s^-1

  ! Use formulation by Nischal et al. [2019]: SABER observations JGR, Space Physics, 124 2338-2356
  ! Consider only O-CO2 collisional cooling: Others are much, much less
  CO2_collisional_deexcitation(1:nLons,1:nLats,1:nAlts) = &
         ( &
           (3.5e-19)*sqrt(true_temp) + &
           (2.32e-15)*exp(-76.25/(true_temp**(1.0/3.0)))&
         )*NDensityS(1:nLons,1:nLats,1:nAlts,iO_3P_,iBlock)

  CO2_collisional_excitation(1:nLons,1:nLats,1:nAlts) = &
           2.0*CO2_collisional_deexcitation(1:nLons,1:nLats,1:nAlts)*&
            exp(-960.0/true_temp) 

  Omega = CO2_collisional_excitation/&
         (CO2_collisional_deexcitation + &
          CO2_collisional_excitation + &
          CO2_spontaneous_emission )

  CO2Cooling = epsilon_co2(1:nLons,1:nLats,1:nAlts)*energy_co2_emission*&
       Omega * CO2_spontaneous_emission(1:nLons,1:nLats,1:nAlts) *  &
       NDensityS(1:nLons,1:nLats,1:nAlts,iCO2_,iBlock)

  CO2Cooling2d = 0.0
  do iAlt=1,nAlts
     RadiativeCooling2d(1:nLons, 1:nLats) = &
          RadiativeCooling2d(1:nLons, 1:nLats) + &
          CO2Cooling(1:nLons,1:nLats,iAlt) * &
          dAlt_GB(1:nLons,1:nLats,iAlt,iBlock)
     CO2Cooling2d = CO2Cooling2d + &
          CO2Cooling(1:nLons,1:nLats,iAlt) * &
          dAlt_GB(1:nLons,1:nLats,iAlt,iBlock)
  enddo
  CO2Cooling = CO2Cooling / TempUnit(1:nLons,1:nLats,1:nAlts) / &
       (Rho(1:nLons,1:nLats,1:nAlts,iBlock)*cp(:,:,1:nAlts,iBlock))

endsubroutine calc_co2_cooling_fomichev

subroutine calc_o3p_cooling(iBlock)

  use ModSources
  use ModEUV
  use ModGITM
  use ModTime
  
  implicit none

  integer, intent(in) :: iBlock

  integer :: iAlt, iError, iDir, iLat, iLon

  real :: tmp2(nLons, nLats, nAlts)
  real :: tmp3(nLons, nLats, nAlts)
  real :: Omega(nLons, nLats, nAlts)

  ! [O] cooling 
  ! Reference: Kockarts, G., Plant. Space Sci., Vol. 18, pp. 271-285, 1970
  ! We reduce the LTE 63-um cooling rate by a factor of 2 for 
  ! the non-LTE effects.[Roble,1987]         

  tmp2 = exp(-228./(Temperature(1:nLons,1:nLats,1:nAlts,iBlock)*&
       TempUnit(1:nLons,1:nLats,1:nAlts)))
  tmp3 = exp(-326./(Temperature(1:nLons,1:nLats,1:nAlts,iBlock)*&
       TempUnit(1:nLons,1:nLats,1:nAlts)))

  ! In erg/cm3/s
  ! JMB: Updated to account for cool-to-space. Divide emission by 2.0
  !      As only 1/2 is emitted outward and cools. 
  !      Assumption: Downward propagating photons are re-absorbed.
  !      Old form has 1.67e-18 for first coefficient, assuming all photons escape:
  !      Also Note: The form below is consistent with current WACCM/WACCM-X  
  ! Formula below calculates cooling rates in ergs/cm^3/s
  OCooling = (0.8375e-18*tmp2 + 2.2545e-20*tmp3) * &
       (NDensityS(1:nLons,1:nLats,1:nAlts,iO_3P_,iBlock)/1.0e6) / &
       (1.0 + 0.6*tmp2 + 0.2*tmp3)
  ! Convert ergs/cm^3/s -> J/m^3/s
  OCooling = OCooling/10.0

  OCooling2d = 0.0
  do iAlt=1,nAlts
     RadiativeCooling2d(1:nLons, 1:nLats) = &
          RadiativeCooling2d(1:nLons, 1:nLats) + &
          OCooling(1:nLons,1:nLats,iAlt) * &
          dAlt_GB(1:nLons,1:nLats,iAlt,iBlock)
     OCooling2d = OCooling2d + &
          OCooling(1:nLons,1:nLats,iAlt) * &
          dAlt_GB(1:nLons,1:nLats,iAlt,iBlock)
  enddo

 ! Convert energy rates into [GITM Temp]/s
 OCooling = OCooling/ TempUnit(1:nLons,1:nLats,1:nAlts) / &
       (Rho(1:nLons,1:nLats,1:nAlts,iBlock)*cp(:,:,1:nAlts,iBlock))

endsubroutine calc_o3p_cooling

subroutine calc_no_cooling(iBlock)

  use ModSources
  use ModEUV
  use ModGITM
  use ModTime
  
  implicit none

  integer, intent(in) :: iBlock

  integer :: iAlt, iError, iDir, iLat, iLon

  real :: tmp2(nLons, nLats, nAlts)
  real :: tmp3(nLons, nLats, nAlts)
  real :: Omega(nLons, nLats, nAlts)
  ! Updated 2022 by Jared Bell (JMB):
  ! Added more readily understood variables here
  real :: NO_collisional_excitation(nLons, nLats, nAlts)
  real :: NO_collisional_deexcitation(nLons, nLats, nAlts)
  real :: NO_spontaneous_emission(nLons,nLats,nAlts)
  real :: Earthshine_NO_excitation(nLons,nLats,nAlts)
  real :: true_temp(1:nLons,1:nLats,1:nAlts)
  real :: nNO_excited(1:nLons,1:nLats,1:nAlts)
  real :: energy_NO_emission
  real :: emission_wavelength

  ! [NO] cooling 
  ! [Original Reference]: Kockarts,G., G.R.L.,VOL.7, PP.137-140,Feberary 1980 ]
  ! [Updated  Reference]: Oberhide et al. JGR, VOL.118, PP.7283-7293, 2013
  !---------
  ! Updated 2022 by Jared M. Bell (JMB)
  ! JMB:  Implement Update(s) from Oberhide et al. [2013] 
  !       Incorporate changes currently active in WACCM-X
  !---------
  ! Note the constants for collisional de-excit/excit are in m^3/s
  ! Thus, the net rate k_rxn * [O/O2] = s^-1 units
  ! collision excitation rate due to O-NO and O2-NO 
  !-----------
  ! Use real temperature (in K)
  true_temp(1:nLons,1:nLats,1:nAlts) = &
        Temperature(1:nLons,1:nLats,1:nAlts,iBlock)*&
           TempUnit(1:nLons,1:nLats,1:nAlts)
  !-------

  !NOCooling = Planck_Constant * Speed_Light / &
  !     5.3e-6 * &
  !     Omega * 13.3 *  &
  !     exp(- Planck_Constant * Speed_Light / &
  !     (5.3e-6 * Boltzmanns_Constant * &
  !     Temperature(1:nLons,1:nLats,1:nAlts,iBlock)* &
  !     TempUnit(1:nLons,1:nLats,1:nAlts))) * &
  !     NDensityS(1:nLons,1:nLats,1:nAlts,iNO_,iBlock)

  emission_wavelength = 5.3e-06  ! emission wavelength (IR) in meter
  energy_NO_emission = Planck_Constant*Speed_Light/emission_wavelength ! Energy in J
  NO_spontaneous_emission(1:nLons,1:nLats,1:nAlts) = 12.54  ! Einstien coef 12.54 s^-1
  Earthshine_NO_excitation(1:nLons,1:nLats,1:nAlts) = 1.06e-04 ! Earthshine excitation

  !real :: NO_collisional_deexcitation(nLons, nLats, nAlts)
  NO_collisional_deexcitation(1:nLons,1:nLats,1:nAlts) = &
          (2.8e-17)*NDensityS(1:nLons,1:nLats,1:nAlts,iO_3P_,iBlock)

  NO_collisional_excitation(1:nLons,1:nLats,1:nAlts) = &
          (2.8e-17)*NDensityS(1:nLons,1:nLats,1:nAlts,iO_3P_,iBlock)*&
           exp(-energy_NO_emission/Boltzmanns_Constant/true_temp(1:nLons,1:nLats,1:nAlts))

  Omega = (NO_collisional_excitation + Earthshine_NO_excitation)/&
       (NO_collisional_deexcitation + &
        NO_collisional_excitation + Earthshine_NO_excitation + &
        NO_spontaneous_emission )

  NOCooling = energy_no_emission*&
       Omega * &
       NO_spontaneous_emission(1:nLons,1:nLats,1:nAlts) * &
       NDensityS(1:nLons,1:nLats,1:nAlts,iNO_,iBlock)

  NOCooling2d = 0.0
  do iAlt=1,nAlts
     RadiativeCooling2d(1:nLons, 1:nLats) = &
          RadiativeCooling2d(1:nLons, 1:nLats) + &
          NOCooling(1:nLons,1:nLats,iAlt) * &
          dAlt_GB(1:nLons,1:nLats,iAlt,iBlock)
     NOCooling2d = NOCooling2d + &
          NOCooling(1:nLons,1:nLats,iAlt) * &
          dAlt_GB(1:nLons,1:nLats,iAlt,iBlock)
  enddo

  NOCooling = NOCooling / TempUnit(1:nLons,1:nLats,1:nAlts) / &
       (Rho(1:nLons,1:nLats,1:nAlts,iBlock)*cp(:,:,1:nAlts,iBlock))

endsubroutine calc_no_cooling
