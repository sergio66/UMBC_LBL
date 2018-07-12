c same as cntnm_progr.f except modified to only dump out coeffs for given input T
c ie the pressure, path length and gas amount can be defaulted to dummy values
c
c ifort -u -extend-source 132 -r8 -i8 -w -Vaxlib driver_cntnm_progr_sergio.f
c cd ../build; gmake -f make_cntnm_sergio_run8 linuxINTELdbl

c %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

C     path:      $Source$
C     author:    $Author: malvarad $
C     revision:  $Revision: 16459 $
C     created:   $Date: 2012-10-22 13:21:52 -0400 (Mon, 22 Oct 2012) $
C 
C
C  --------------------------------------------------------------------------
C |  Copyright ©, Atmospheric and Environmental Research, Inc., 2011         |
C |                                                                          |
C |  All rights reserved. This source code is part of the MT_CKD continuum   |
C |  software and is designed for scientific and research purposes.          |
C |  Atmospheric and Environmental Research, Inc. (AER) grants USER          |
C |  the right to download, install, use and copy this software              |
C |  for scientific and research purposes only. This software may be         |
C |  redistributed as long as this copyright notice is reproduced on any     |
C |  copy made and appropriate acknowledgment is given to AER. This          |
C |  software or any modified version of this software may not be            |
C |  incorporated into proprietary software or commercial software           |
C |  offered for sale.                                                       |
C |                                                                          |
C |  This software is provided as is without any express or implied          |
C |  warranties.                                                             |
C |                       (http://www.rtweb.aer.com/)                        |
C  --------------------------------------------------------------------------
C
C                  PROGRAM DRCNTNM

C>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>> NEW >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
C>>>>>>>>>>>>>>>>>>>>>>>>>>>> begin junk1.f >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>

      integer iGasID,kMaxPtsDVS,kMaxPts,iNumPtsDVS,iNumPts
      PARAMETER (kMaxPtsDVS = 1000)     !!! at coarser spacing
      PARAMETER (kMaxPts    = 100000)   !!! at user spacing
      real num_kmolesIN,ppaveIN,taveIN,paveIN
      real v1absIN,v2absIN,dvabsIN
      real raFreqDVS(kMaxPtsDVS),raSelfDVS(kMaxPtsDVS),
     $                           raFornDVS(kMaxPtsDVS)
      real raFreq(kMaxPts),raAbs(kMaxPts)

      print *,'Units as per run8 : gasID atm atm K kmoles/cm2'
      print *,'Enter [iGasID paveIN ppaveIN taveIN  num_kmolesIN] : '
      read *,iGasID,paveIN,ppaveIN,taveIN,num_kmolesIN
      paveIN  = paveIN * 1013.25
      ppaveIN = ppaveIN * 1013.25

      print *,'Enter [v1 v2 dv] : '
      read *,v1absIN,v2absIN,dvabsIN

      !! [Self Forn   o3 o2 n2]      
      IF (iGasID .EQ. 1) THEN
        call CALCON_MTCKD_32_loc(paveIN,ppaveIN,taveIN,num_kmolesIN,
     $                         v1absIN,v2absIN,dvabsIN,
     $                         raFreqDVS,raSelfDVS,raFornDVS,iNumPtsDVS,
     $                         raFreq,raAbs,iNumPts,
     $                         iGasID, 1.0, 1.0, 0.0, 0.0, 0.0)
      ELSEIF (iGasID .EQ. 3) THEN
        call CALCON_MTCKD_32_loc(paveIN,ppaveIN,taveIN,num_kmolesIN,
     $                         v1absIN,v2absIN,dvabsIN,
     $                         raFreqDVS,raSelfDVS,raFornDVS,iNumPtsDVS,
     $                         raFreq,raAbs,iNumPts,
     $                         iGasID, 0.0, 0.0, 1.0, 0.0, 0.0)
      ELSEIF (iGasID .EQ. 7) THEN
        call CALCON_MTCKD_32_loc(paveIN,ppaveIN,taveIN,num_kmolesIN,
     $                         v1absIN,v2absIN,dvabsIN,
     $                         raFreqDVS,raSelfDVS,raFornDVS,iNumPtsDVS,
     $                         raFreq,raAbs,iNumPts,
     $                         iGasID, 0.0, 0.0, 0.0, 1.0, 0.0)
      ELSEIF (iGasID .EQ. 22) THEN
        call CALCON_MTCKD_32_loc(paveIN,ppaveIN,taveIN,num_kmolesIN,
     $                         v1absIN,v2absIN,dvabsIN,
     $                         raFreqDVS,raSelfDVS,raFornDVS,iNumPtsDVS,
     $                         raFreq,raAbs,iNumPts,
     $                         iGasID, 0.0, 0.0, 0.0, 0.0, 1.0)
      ELSE
        print *,'Error : need iGasID = 1,3,7,22 for WV,O3,O3,N2, not',iGasID
        STOP
      END IF

      END

c%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
      include 'cntnm_progr_sergio.f'
