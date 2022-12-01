! Copyright 2021, the GITM Development Team (see srcDoc/dev_team.md for members)
! Full license can be found in LICENSE

subroutine init_get_aurora

  use ModGITM
  use ModTime
  use ModIndicesInterfaces
  use ModInputs
  use ModUserGITM
  use ModNewell
  use ModOvationSME
  use ModAeAuroralModel
  use ModFtaModel
  use ModEIE_Interface, only: EIEr3_HaveLats, EIEr3_HaveMLTs

  implicit none

  character (len=iCharLen_), dimension(100) :: Lines
  character (len=iCharLen_) :: TimeLine
  real    :: bz

  logical :: IsFirstTime = .true.
  integer :: iError

  integer :: nAmieLats, nAmieMlts, nAmieBlocks

  iError = 0

  if (.not.IsFirstTime .or. IsFramework) return
  call report("init_get_aurora",2)

  IsFirstTime = .false.

  if (UseNewellAurora) then
     call init_newell
     UseIMF = .true.
  endif

  if (UseOvationSME) then
     call read_ovationsm_files
  endif

  if (UseAeModel) then
     call read_ae_model_files(iError)
  endif

  if (UseFtaModel) then
     call initialize_fta
  endif

end subroutine init_get_aurora

!--------------------------------------------------------------------
! get_aurora
!--------------------------------------------------------------------
subroutine get_aurora(iBlock)

  use ModGITM
  use ModTime
  use ModIndicesInterfaces
  use ModInputs
  use ModUserGITM
  use ModNewell
  use ModOvationSME, only: run_ovationsme
  use ModAeAuroralModel, only: run_ae_model
  use ModFtaModel, only: run_fta_model
  use ModEIE_Interface, only: UAl_UseGridBasedEIE
  use ModMpi

  implicit none

  integer, intent(in) :: iBlock

  integer :: iError, iLat, iLon, iAlt, iPot, nPot, iDir, nDir=1
  logical :: IsFirstTime = .true.
  logical :: IsFirstAurora(nBlocksMax) = .true.

  real, dimension(-1:nLons+2, -1:nLats+2) :: lats, mlts, EFlux
  real :: by, bz, CuspLat, CuspMlt

  call start_timing("get_aurora")
  call report("get_aurora",2)

  iError = 0

  if (index(cPlanet,"Earth") == 0) then 
     ElectronAverageEnergy = 0.1
     ElectronEnergyFlux = 0.0001
     return
  endif

  ! -----------------------------------------------------
  ! Now get the aurora.
  ! This assumes that the field lines are basically
  ! vertical starting at the top of the model.
  ! -----------------------------------------------------
  if (floor((tSimulation-dt)/DtAurora) /= &
       floor((tsimulation)/DtAurora) .or. IsFirstAurora(iBlock)) then

     call report("Getting Aurora",1)
     call init_get_aurora 

     iAlt = nAlts+1

     if (UseNewellAurora) then
        call run_newell(iBlock)
     elseif (UseOvationSME) then 
        call run_ovationsme(StartTime, CurrentTime, iBlock)
     elseif (UseAeModel) then
        call run_ae_model(CurrentTime, iBlock)
     elseif (UseFtaModel) then
        call run_fta_model(CurrentTime, iBlock)
     else

        call UA_SetGrid(                    &
             MLT(-1:nLons+2,-1:nLats+2,iAlt), &
             MLatitude(-1:nLons+2,-1:nLats+2,iAlt,iBlock), iError)

        if (iError /= 0) then
           write(*,*) "Error in routine get_potential (UA_SetGrid):"
           write(*,*) iError
           call stop_gitm("Stopping in get_potential")
        endif

        call UA_GetAveE(ElectronAverageEnergy, iError)
        if (iError /= 0) then
           write(*,*) "Error in get_potential (UA_GetAveE):"
           write(*,*) iError
           ElectronAverageEnergy = 1.0
        endif

        ! Sometimes, in AMIE, things get messed up in the
        ! Average energy, so go through and fix some of these.
        
        do iLat=-1,nLats+2
           do iLon=-1,nLons+2
              if (ElectronAverageEnergy(iLon,iLat) < 0.0) then
                 ElectronAverageEnergy(iLon,iLat) = 0.1
                 write(*,*) "ave e i,j Negative : ",iLon,iLat,&
                      ElectronAverageEnergy(iLon,iLat)
              endif
              if (ElectronAverageEnergy(iLon,iLat) > 100.0) then
                 write(*,*) "ave e i,j Positive : ",iLon,iLat,&
                      ElectronAverageEnergy(iLon,iLat)
                 ElectronAverageEnergy(iLon,iLat) = 0.1
              endif
           enddo
        enddo

        call UA_GetEFlux(ElectronEnergyFlux, iError)
        if (iError /= 0) then
           write(*,*) "Error in get_potential (UA_GetEFlux):"
           write(*,*) iError
           ElectronEnergyFlux = 0.1
        endif

        ! -----------------------------------------------------
        ! Get Ion Precipitation if desired
        ! -----------------------------------------------------

        if (UseIonPrecipitation) then

           call UA_GetIonAveE(IonAverageEnergy, iError)
           if (iError /= 0) then
              write(*,*) "Error in get_potential (UA_GetAveE):"
              write(*,*) iError
              IonAverageEnergy = 1.0
           endif

           call UA_GetIonEFlux(IonEnergyFlux, iError)
           if (iError /= 0) then
              write(*,*) "Error in get_potential (UA_GetEFlux):"
              write(*,*) iError
              IonEnergyFlux = 0.1
           endif

        endif
        
     endif

    if (UseCusp) then

       lats = abs(MLatitude(-1:nLons+2,-1:nLats+2,iAlt,iBlock))

       if (maxval(lats) > 50) then

          mlts = mod(MLT(-1:nLons+2,-1:nLats+2,iAlt)+24.0,24.0)

          call get_IMF_Bz(CurrentTime+TimeDelayHighLat, bz, iError)
          call get_IMF_By(CurrentTime+TimeDelayHighLat, by, iError)

          ! If we are in the southern hemisphere, reverse by:
          if (lats(nLons/2, nLats/2) < 0.0) by = -by

          if (bz > 0) then
             ! Newell et al., 1988:
             CuspLat = 77.2 + 0.11 * bz
             ! Asai et al., Earth Planets Space, 2005:
             CuspMlt = 11.755 + 0.169 * by
          else
             ! Asai et al., Earth Planets Space, 2005:
             CuspMlt = 11.949 + 0.0826 * by
             ! Zhang et al., JGR, 2005:
             if (Bz > -10) then
                CuspLat = 77.2 + 1.1 * bz
             else
                CuspLat = 21.7 * exp(0.1 * bz) + 58.2
             endif
          endif

          EFlux = CuspEFlux * &
               exp(-abs(lats - CuspLat)/CuspLatHalfWidth) *  &
               exp(-abs(mlts - CuspMlt)/CuspMltHalfWidth)

          do iLat=-1,nLats+2
             do iLon=-1,nLons+2
                if (EFlux(iLon,iLat) > 0.1) then
                   ElectronEnergyFlux(iLon,iLat) = EFlux(iLon,iLat)
                   ElectronAverageEnergy(iLon,iLat) = CuspAveE
                endif
             enddo
          enddo

       endif

    endif

     if (iDebugLevel > 2) &
          write(*,*) "==> Max, electron_ave_ene : ", &
          maxval(ElectronAverageEnergy), &
          maxval(ElectronEnergyFlux)

     IsFirstAurora(iBlock) = .false.

  endif

  call end_timing("get_aurora")

end subroutine get_aurora


