! Copyright 2021, the GITM Development Team (see srcDoc/dev_team.md for members)
! Full license can be found in LICENSE

!\
! MODULE ------------------------------------------------
!/

module ModVertical

  use ModSizeGitm, only: nAlts
  use ModPlanet, only: nSpecies, nIonsAdvect, nSpeciesTotal, nIons

  implicit none

  real, dimension(-1:nAlts+2)   :: dAlt_F, InvDAlt_F, Altitude_G, Gravity_G
  real, dimension(-1:nAlts+2)   :: InvRadialDistance_C, dAlt_C,Cv_1D
  real, dimension(-1:nAlts+2)   :: LogRho, Temp,MeanMajorMass_1d,Gamma_1d
  real, dimension(-1:nAlts+2,3) :: Vel_GD
  real, dimension(-1:nAlts+2,3) :: IVel
  
  real, dimension(-1:nAlts+2)             :: NewLogRho, NewTemp
  real, dimension(-1:nAlts+2,3)           :: NewVel_GD
  real, dimension(-1:nAlts+2,nSpecies) :: NewLogNS, NewVertVel
!  real, dimension(-1:nAlts+2,nIonsAdvect) :: NewLogINS
!  real, dimension(-1:nAlts+2,nIonsAdvect) :: LogINS

  real, dimension(-1:nAlts+2,nIons) :: NewLogINS
  real, dimension(-1:nAlts+2,nIons) :: LogINS
!  real, dimension(1:2,nIons) :: OldLogINSc


  real, dimension(-1:nAlts+2, nSpecies) :: LogNS, VertVel
  real, dimension(nAlts, nSpecies) :: NDensityS_1D

  real, dimension(0:nAlts+1) :: cMax
  
  real :: SZAVertical
  real :: MLatVertical
 
  integer :: iLon1D, iLat1D, iBlock1D

  real :: Heating(nAlts)
  real, dimension(-1:nAlts+2) :: EddyCoef_1d

  real :: Centrifugal, Coriolis, Lat, Lon, dAltdlon_1D, dAltdLat_1D

! JMB 01/2022
  real, dimension(1:nAlts  ) :: CellVol1D
  real, dimension(1:nAlts  ) :: Area_P12, Area_M12
  real, dimension(0:nAlts+1, 1:3) :: &
                   Mesh_ULP120, Mesh_ULP121, Mesh_ULP122, &
                   Mesh_URM120, Mesh_URM121, Mesh_URM122
  real, dimension(0:nAlts+1) :: &
                   Mesh_CLP120, Mesh_CLP121, Mesh_CLP122, &
                   Mesh_CRM120, Mesh_CRM121, Mesh_CRM122
  real, dimension(0:nAlts+1,1:3) :: &
                   Mesh_IS0, Mesh_IS1, Mesh_IS2
  real, dimension(1:nAlts,1:5) :: &
                   CD_MeshCoefs
  real, dimension(1:3,1:5) :: &
                   UB_MeshCoefs, LB_MeshCoefs
!-- 


end module ModVertical

