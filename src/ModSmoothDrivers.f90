!  Copyright (C) 2002 Regents of the University of Michigan, portions used with permission 
!  For more information, see http://csem.engin.umich.edu/tools/swmf

! JMB: 11/2022:
! Use these variables to smooth out drivers
module ModSmoothDrivers

  use ModGITM, only: nBlocksmax,nLons,nLats,nAlts
  use ModKind, only: Real8_

  implicit none

  !LatitudeMeshCoefs
  ! We need 0:nLats+1 in the Latitude Mesh Coefs
  ! Left Faces
  real :: PreviousPotential(-1:nLons+2,-1:nLats+2,-1:nAlts+2,1:nBlocksMax)=0.0
  real ::     NextPotential(-1:nLons+2,-1:nLats+2,-1:nAlts+2,1:nBlocksMax)=0.0
  real :: PreviousPotentialY(-1:nLons+2,-1:nLats+2,-1:nAlts+2,1:nBlocksMax)=0.0
  real ::     NextPotentialY(-1:nLons+2,-1:nLats+2,-1:nAlts+2,1:nBlocksMax)=0.0
  real(Real8_)          :: LastPotentialTime
end module ModSmoothDrivers
