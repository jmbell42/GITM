!  Copyright (C) 2002 Regents of the University of Michigan, portions used with permission 
!  For more information, see http://csem.engin.umich.edu/tools/swmf

! JMB: 03/2021:
! These initialize the Mesh Coefficients for the 5th order WENO Scheme
module ModMesh

  use ModGITM, only: nBlocksmax,nLons,nLats,nAlts

  implicit none

  !LatitudeMeshCoefs
  ! We need 0:nLats+1 in the Latitude Mesh Coefs
  ! Left Faces
  real :: LatMesh_ULP120(1:nLons,0:nLats+1,1:nAlts,1:3,nBlocksMax)=0.0
  real :: LatMesh_ULP121(1:nLons,0:nLats+1,1:nAlts,1:3,nBlocksMax)=0.0
  real :: LatMesh_ULP122(1:nLons,0:nLats+1,1:nAlts,1:3,nBlocksMax)=0.0
  !
  real :: LatMesh_CLP120(1:nLons,0:nLats+1,1:nAlts,nBlocksMax)=0.0
  real :: LatMesh_CLP121(1:nLons,0:nLats+1,1:nAlts,nBlocksMax)=0.0
  real :: LatMesh_CLP122(1:nLons,0:nLats+1,1:nAlts,nBlocksMax)=0.0
  !
  ! Right Faces
  real :: LatMesh_URM120(1:nLons,0:nLats+1,1:nAlts,1:3,nBlocksMax)=0.0
  real :: LatMesh_URM121(1:nLons,0:nLats+1,1:nAlts,1:3,nBlocksMax)=0.0
  real :: LatMesh_URM122(1:nLons,0:nLats+1,1:nAlts,1:3,nBlocksMax)=0.0

  real :: LatMesh_CRM120(1:nLons,0:nLats+1,1:nAlts,nBlocksMax)=0.0
  real :: LatMesh_CRM121(1:nLons,0:nLats+1,1:nAlts,nBlocksMax)=0.0
  real :: LatMesh_CRM122(1:nLons,0:nLats+1,1:nAlts,nBlocksMax)=0.0
  !------
  real :: LatMesh_IS0(1:nLons,0:nLats+1,1:nAlts,1:3,nBlocksMax)=0.0
  real :: LatMesh_IS1(1:nLons,0:nLats+1,1:nAlts,1:3,nBlocksMax)=0.0
  real :: LatMesh_IS2(1:nLons,0:nLats+1,1:nAlts,1:3,nBlocksMax)=0.0
  !
  !LatitudeMeshCoefs
  ! We need 0:nLats+1 in the Latitude Mesh Coefs
  ! Left Faces
  real :: LonMesh_ULP120(0:nLons+1,1:nLats,1:nAlts,1:3,nBlocksMax)=0.0
  real :: LonMesh_ULP121(0:nLons+1,1:nLats,1:nAlts,1:3,nBlocksMax)=0.0
  real :: LonMesh_ULP122(0:nLons+1,1:nLats,1:nAlts,1:3,nBlocksMax)=0.0
  !
  real :: LonMesh_CLP120(0:nLons+1,1:nLats,1:nAlts,nBlocksMax)=0.0
  real :: LonMesh_CLP121(0:nLons+1,1:nLats,1:nAlts,nBlocksMax)=0.0
  real :: LonMesh_CLP122(0:nLons+1,1:nLats,1:nAlts,nBlocksMax)=0.0
  !
  ! Right Faces
  real :: LonMesh_URM120(0:nLons+1,1:nLats,1:nAlts,1:3,nBlocksMax)=0.0
  real :: LonMesh_URM121(0:nLons+1,1:nLats,1:nAlts,1:3,nBlocksMax)=0.0
  real :: LonMesh_URM122(0:nLons+1,1:nLats,1:nAlts,1:3,nBlocksMax)=0.0

  real :: LonMesh_CRM120(0:nLons+1,1:nLats,1:nAlts,nBlocksMax)=0.0
  real :: LonMesh_CRM121(0:nLons+1,1:nLats,1:nAlts,nBlocksMax)=0.0
  real :: LonMesh_CRM122(0:nLons+1,1:nLats,1:nAlts,nBlocksMax)=0.0
  !------
  real :: LonMesh_IS0(0:nLons+1,1:nLats,1:nAlts,1:3,nBlocksMax)=0.0
  real :: LonMesh_IS1(0:nLons+1,1:nLats,1:nAlts,1:3,nBlocksMax)=0.0
  real :: LonMesh_IS2(0:nLons+1,1:nLats,1:nAlts,1:3,nBlocksMax)=0.0
  !
  ! We need 0:nAlts+1 in the Altitude Mesh Coefs
  ! Left Faces
  real :: AltMesh_ULP120(1:nLons,1:nLats,0:nAlts+1,1:3,nBlocksMax)=0.0
  real :: AltMesh_ULP121(1:nLons,1:nLats,0:nAlts+1,1:3,nBlocksMax)=0.0
  real :: AltMesh_ULP122(1:nLons,1:nLats,0:nAlts+1,1:3,nBlocksMax)=0.0
  !
  real :: AltMesh_CLP120(1:nLons,1:nLats,0:nAlts+1,nBlocksMax)=0.0
  real :: AltMesh_CLP121(1:nLons,1:nLats,0:nAlts+1,nBlocksMax)=0.0
  real :: AltMesh_CLP122(1:nLons,1:nLats,0:nAlts+1,nBlocksMax)=0.0
  ! Right Faces
  real :: AltMesh_URM120(1:nLons,1:nLats,0:nAlts+1,1:3,nBlocksMax)=0.0
  real :: AltMesh_URM121(1:nLons,1:nLats,0:nAlts+1,1:3,nBlocksMax)=0.0
  real :: AltMesh_URM122(1:nLons,1:nLats,0:nAlts+1,1:3,nBlocksMax)=0.0

  real :: AltMesh_CRM120(1:nLons,1:nLats,0:nAlts+1,nBlocksMax)=0.0
  real :: AltMesh_CRM121(1:nLons,1:nLats,0:nAlts+1,nBlocksMax)=0.0
  real :: AltMesh_CRM122(1:nLons,1:nLats,0:nAlts+1,nBlocksMax)=0.0
  !------
  real :: AltMesh_IS0(1:nLons,1:nLats,0:nAlts+1,1:3,nBlocksMax)=0.0
  real :: AltMesh_IS1(1:nLons,1:nLats,0:nAlts+1,1:3,nBlocksMax)=0.0
  real :: AltMesh_IS2(1:nLons,1:nLats,0:nAlts+1,1:3,nBlocksMax)=0.0

  ! Use these for Central difference mesh coefficients
  real ::     AltCDMeshCoef(1:nLons,1:nLats,1:nAlts,1:5,nBlocksMax)=0.0
  ! These are the One-Sided Mesh Coefficients
  ! LB is the Lower boundary (cells [-1, 0, 1])
  ! UB is the Upper boundary (cells [nAlts,nAlts+1,nAlts+2])
  real :: OS_LBMeshCoef(1:nLons,1:nLats,1:3,1:5,nBlocksMax)=0.0
  real :: OS_UBMeshCoef(1:nLons,1:nLats,1:3,1:5,nBlocksMax)=0.0

end module ModMesh
