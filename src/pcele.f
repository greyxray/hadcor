*---------------------------------------------------------
      REAL FUNCTION PCele(Imod,Emes,Emip,The,Phi)

*  May 2003
*********************************************************
*
*   CELE modified for prompt photon studies 
*                                                       *
*    gives the value of corrected electron in BCAL using BPRE 
*    Jan 08, 2001, by S.Chekanov  
*    Nov 20, 2001, corrected  
*
*    written by S.Chekanov (chekanov@mail.desy.de)
*  
*    Input: 
*  
*    Imod <0  : do not correct electron  
*    Imod =0  : 5%  correction according to CALCORR  
*    Imod =1  : correct data after EM non-uniformity correction 
*    Imod =10 : correct MC after non-uniformity correction
*    Imod =2  : correct data if non-uniformity correction is not used
*    Imod =20 : correct MC if non-uniformity correction is not used 
*
*    Emes - Energy measured with Electron method
*    Emip - number of BPRE mips for the electron
*    The and Phi  - theta and phi of electron (in rad)  
*

      Integer Imod
      Real fac, Emes,Emip,The,Phi  
      Real PHImin(9), PHImax(9),THEmin(9), THEmax(9) 
      Real a0(9),a1(9),a1f(9),b0(9),b1(9) 
      real pi2 
      parameter (pi2=6.28315)
      real dphi, phi0, dthe, the0
      parameter (dphi=0.39269, phi0=0.0, dthe=0.1238, the0=0.64)

      LOGICAL FIRST
      DATA FIRST/.TRUE./
      SAVE FIRST,PHImin,PHImax,THEmin,THEmax 
      SAVE a0,a1,b0,b1

 33   FORMAT(A37,I2)   
      IF(FIRST) THEN
      FIRST = .FALSE.
      write(6,*) ' '
      write(6,33)'CALL TO CELE SUBROUTINE with option ',imod  
      write(6,*)'======================================'
      IF (Imod.lt.0) write(6,*) 'BPRE Correction is not applied' 
      IF (Imod.eq.0) write(6,*) 'CALCORR 5% correction' 
      IF (Imod.eq.1)write(6,*) 'BPRE corrects data after non-uniform'
      IF (Imod.eq.10)write(6,*) 'BPRE corrects MC  after non-uniform'
      IF (Imod.eq.2)write(6,*) 'BPRE corrects data without non-uniform'
      IF (Imod.eq.20)write(6,*) 'BPRE corrects MC  without non-uniform'
      write(6,*)'======================================'
* bin 1 ----------------------
       PHImin(1)=phi0
       PHImax(1)=pi2 
       THEmin(1)=the0
       THEmax(1)=the0+4*dthe
* bin 2----------------------
       PHImin(2)=phi0
       PHImax(2)=phi0+4*dphi
       THEmin(2)=the0+4*dthe
       THEmax(2)=the0+11*dthe
* bin 3----------------------
       PHImin(3)=phi0+4*dphi
       PHImax(3)=phi0+8*dphi
       THEmin(3)=the0+4*dthe
       THEmax(3)=the0+11*dthe
* bin 4----------------------
       PHImin(4)=phi0+8*dphi
       PHImax(4)=phi0+12*dphi
       THEmin(4)=the0+4*dthe
       THEmax(4)=the0+11*dthe
* bin 5----------------------
       PHImin(5)=phi0+12*dphi
       PHImax(5)=phi0+16*dphi
       THEmin(5)=the0+4*dthe
       THEmax(5)=the0+11*dthe
* bin 6----------------------
       PHImin(6)=phi0
       PHImax(6)=phi0+4*dphi
       THEmin(6)=the0+11*dthe
       THEmax(6)=the0+13*dthe
* bin 7----------------------
       PHImin(7)=phi0+4*dphi
       PHImax(7)=phi0+8*dphi
       THEmin(7)=the0+11*dthe
       THEmax(7)=the0+13*dthe
* bin 8----------------------
       PHImin(8)=phi0+8*dphi
       PHImax(8)=phi0+12*dphi
       THEmin(8)=the0+11*dthe
       THEmax(8)=the0+13*dthe
* bin 9----------------------
       PHImin(9)=phi0+12*dphi
       PHImax(9)=phi0+16*dphi
       THEmin(9)=the0+11*dthe
       THEmax(9)=the0+13*dthe

       
* bin1 ---------------------------------

       IF (imod.eq.1) then ! data   
       a0(1)=-0.01    
       a1(1)=-0.08                         
       b0(1)=0.999                
       b1(1)=0.0  
       elseif (imod.eq.10) then ! MC  
       a0(1)=-0.01 
       a1(1)=-0.06             
       b0(1)=0.999     
       b1(1)=0.0  
       elseif (imod.eq.2) then  ! data before non-unif  
       a0(1)=-0.02 
       a1(1)=-0.09  
       b0(1)=0.998  
       b1(1)=0.001   
       elseif (imod.eq.20) then ! MC before non-unif  
       a0(1)=-0.02 
       a1(1)=-0.06    
       b0(1)=0.998 
       b1(1)=0.001   
       endif 


* bin 2----------------------
      if (imod.eq.1) then ! data
       a0(2)=-0.02  
       a1(2)=-0.06          
       b0(2)=0.990        
       b1(2)=0.0 
      elseif (imod.eq.10) then ! MC
       a0(2)=-0.02   
       a1(2)=-0.04     
       b0(2)=0.996     
       b1(2)=0.0 
      elseif (imod.eq.2) then  ! data before non-unif
       a0(2)=-0.02 
       a1(2)=-0.06 
       b0(2)=0.980  
       b1(2)=0.002 
      elseif (imod.eq.20) then ! MC before non-unif 
       a0(2)=-0.01  
       a1(2)=-0.04 
       b0(2)=0.993   
       b1(2)=0.002 
      endif 

* bin 3----------------------

       if (imod.eq.1) then ! data
       a0(3)=-0.01 
       a1(3)=-0.05        
       b0(3)= 0.996        
       b1(3)= 0.0 
       elseif (imod.eq.10) then ! MC 
       a0(3)=-0.01  
       a1(3)=-0.04    
       b0(3)=0.999
       b1(3)=0.0 
       elseif (imod.eq.2) then  ! data before non-unif
       a0(3)=-0.03 
       a1(3)=-0.05 
       b0(3)=0.989 
       b1(3)=0.001 
       elseif (imod.eq.20) then ! MC before non-unif 
       a0(3)=-0.02   
       a1(3)=-0.04  
       b0(3)=0.990  
       b1(3)=0.001 
       endif

* bin 4----------------------

       if (imod.eq.1) then ! data
       a0(4)=-0.01   
       a1(4)=-0.07    
       b0(4)= 0.987     
       b1(4)= 0.0 
       elseif (imod.eq.10) then ! MC
       a0(4)=-0.01   
       a1(4)=-0.04   
       b0(4)= 0.99     
       b1(4)= 0.0 
       elseif (imod.eq.2) then  ! data before non-unif 
       a0(4)=-0.01
       a1(4)=-0.07 
       b0(4)=0.977   
       b1(4)=0.001  
       elseif (imod.eq.20) then ! MC before non-unif
       a0(4)=-0.01  
       a1(4)=-0.04  
       b0(4)=0.989   
       b1(4)=0.001 
       endif

* bin 5----------------------

       if (imod.eq.1) then ! data 
       a0(5)=-0.01       
       a1(5)=-0.05      
       b0(5)= 0.996            
       b1(5)= 0.0 
       elseif (imod.eq.10) then ! MC
       a0(5)=-0.01   
       a1(5)=-0.035       
       b0(5)= 0.997     
       b1(5)= 0.0 
       elseif (imod.eq.2) then  ! data before non-unif
       a0(5)=-0.01 
       a1(5)=-0.05 
       b0(5)=0.988   
       b1(5)=0.001 
       elseif (imod.eq.20) then ! MC before non-unif
       a0(5)=-0.01  
       a1(5)=-0.035  
       b0(5)=0.992    
       b1(5)=0.001 
       endif 

* bin 6----------------------

       if (imod.eq.1) then ! data
       a0(6)=-0.02    
       a1(6)=-0.07    
       b0(6)= 0.987          
       b1(6)= 0.0 
       elseif (imod.eq.10) then ! MC
       a0(6)=-0.01   
       a1(6)=-0.04   
       b0(6)= 0.99    
       b1(6)= 0.0 
       elseif (imod.eq.2) then  ! data before non-unif
       a0(6)=-0.02
       a1(6)=-0.07  
       b0(6)=0.980  
       b1(6)=0.001 
       elseif (imod.eq.20) then ! MC before non-unif
       a0(6)=-0.02  
       a1(6)=-0.04 
       b0(6)=0.989    
       b1(6)=0.001 
       endif

* bin 7----------------------

       if (imod.eq.1) then ! data
       a0(7)=-0.01   
       a1(7)=-0.06        
       b0(7)= 0.99     
       b1(7)= 0.0 
       elseif (imod.eq.10) then ! MC
       a0(7)=-0.01  
       a1(7)=-0.055    
       b0(7)= 0.99    
       b1(7)= 0.0 
       elseif (imod.eq.2) then  ! data before non-unif
       a0(7)=-0.01 
       a1(7)=-0.07
       b0(7)=0.986 
       b1(7)=0.001   
       elseif (imod.eq.20) then ! MC before non-unif
       a0(7)=-0.01 
       a1(7)=-0.06 
       b0(7)=0.996  
       b1(7)=0.001  
       endif

* bin 8----------------------

       if (imod.eq.1) then ! data
       a0(8)=-0.01   
       a1(8)=-0.05   
       b0(8)= 0.99        
       b1(8)= 0.0 
       elseif (imod.eq.10) then ! MC
       a0(8)=-0.01  
       a1(8)=-0.05     
       b0(8)= 0.99      
       b1(8)= 0.0 
       elseif (imod.eq.2) then  ! data before non-unif
       a0(8)=-0.02 
       a1(8)=-0.05
       b0(8)=0.979  
       b1(8)=0.002 
       elseif (imod.eq.20) then ! MC before non-unif
       a0(8)=-0.01 
       a1(8)=-0.03  
       b0(8)=0.985   
       b1(8)=0.001  
       endif

* bin 9----------------------

       if (imod.eq.1) then ! data 
       a0(9)=-0.01    
       a1(9)=-0.06     
       b0(9)= 0.998           
       b1(9)= 0.0
       elseif (imod.eq.10) then ! MC
       a0(9)=-0.02   
       a1(9)=-0.04   
       b0(9)= 0.99    
       b1(9)= 0.0 
       elseif (imod.eq.2) then  ! data before non-unif
       a0(9)=-0.02 
       a1(9)=-0.06 
       b0(9)=0.980 
       b1(9)=0.001 
       elseif (imod.eq.20) then ! MC before non-unif
       a0(9)=-0.02  
       a1(9)=-0.04  
       b0(9)=0.985  
       b1(9)=0.001  
       endif 

      ENDIF 

      PCele=Emes
      IF (imod.lt.0) RETURN 


      IF (imod.eq.0) THEN
         PCele=1.05*Emes
         RETURN
      ENDIF  

*- here BPRE correction goes 

      if (Phi.lt.0.0) Phi=pi2+Phi 

*- do corrections

      IF (imod.eq.2)  PCele=1.03*Emes
      IF (imod.eq.20) PCele=1.03*Emes
      IF (imod.eq.1)  PCele=1.03*Emes
      IF (imod.eq.10) PCele=1.04*Emes

      IF (PCele.gt.3.and.Emip.gt.0.1.and.Emip.lt.35) THEN
      DO J=1,9

      fac=1.0 
      IF (PCele.le.15) then
      fac=1.07
      ELSEIF (PCele.gt.15.and.PCele.le.20) then
      fac=1.04       
      ELSEIF (PCele.gt.20.and.PCele.le.25) then
      fac=1.02     
      ELSEIF (PCele.gt.25.and.PCele.le.30) then
      fac=1.0  
      ENDIF

       a1f(J)=fac*a1(J)

       if (THE.gt.THEmin(J).and.THE.lt.THEmax(J)) then
       if (PHI.gt.PHImin(J).and.PHI.lt.PHImax(J)) then
       PCele=(Emes-a0(J)-a1f(J)*Emip)/(b0(J)+b1(J)*Emip)
       endif
       endif
      ENDDO
      ENDIF
c      write (*,*)'Pcele=',Pcele

*
      RETURN
      END
