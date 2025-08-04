      Program alfgas
      include 'alfgasfin_ortho.inc'
* Var directory
*   pp   array   Wx (plasma x vel)
*   qq   array   Wz (plasma z vel)
*   mm   array   bx (mag field x)
*   nn   array   bz (mag field z)
*   gg   array   vx (gas x vel)
*   hh   array   vz (gas z vel)
*   dd   array   rho (gas density)
*   ee   array   plasma density
*   mu  real    mesh ratio = dt/dr
*   r   real    (plasma sound speed)/(alfven speed)
*   s   real    (gas speed)/(alfven speed)
*   al  real    coupling constant
*   be  real    coupling constant
*   gm  real    coupling constant
*   nt  int     total number of timesteps
*   lt  int     total no of x-increments
*   mt  int     total no of z-increments
*   l   int     index for x-dirn
*   m   int     index for z-dirn (mag field dirn)
*   n   int     index for t stepping
*   dr  real    incremental change in r - now exploited as a frequency change
cto use this code for gas-plasma momentum coupling no setu-p is required, the code can be *used dimensionless with parameters read in from alfgasfin_ortho.in.
ccccccccccccccccccccccccccccccc
cto use this code for ionisation then some setup is required, first you must choose *real * densities for both the gas and the plasma and a timestep (microsecond is a good start *as this a fluid approximation code). The ratio of these densities is then used to *calculate a temperature (from the saha equation) this temperature becomes the starting *temperature of the code.
c AT THE MOMENT THIS TEMPERATURE CALCULATION IS DONE EXTERNALLY IN MATLAB "saha.m" IF A CODE FOR THE LAMBERTW FUNCTION CAN BE FOUND IN FORTRAN THEN IT CAN BE INCORPORATED.

      do nr=1,1
       call startup
       call solve
      end do
      end
*-----------------------------------------------------------------
      Subroutine Startup
* opens all i/o files, reads in initial data and parameters
      include 'alfgasfin_ortho.inc'


      if (nr.eq.1) then

       filnam(1)='agwx .out'
       filnam(2)='agwz .out'
       filnam(3)='agvx .out'
       filnam(4)='agvz .out'
       filnam(5)='agbx .out'
       filnam(6)='agbz .out'
       filnam(7)='debg .out'
       filnam(8)='agro .out'
       filnam(9)='agpd .out'
       open(10,file='alfgasfin_ortho.in', status='old')
       rewind(10)
       read(10,datin)
       write(6,datin)
       write(6,200) lt,mt,nt,mu,nprint,fr,m1,m2,noff



       end if
* now continue to open new files below:
      do i=1,9
       write(filnam(i)(5:5),'(I1.1)') nr
       open(20+i-1,file=filnam(i), status='unknown')
       rewind(20+i-1)
      end do
200   format('lt= ',i3,' mt= ',i3,' nt= ',i3,' mu= ',f5.3,
     &       ' nprint= ',i3,' freq= ',f5.3/'m1= ',i3,' m2= ',i3,
     &        'noff =', i4)
* consistency check
      if(lt.ge.lmax) then
       lt=lmax-1
       print*,'Warning: lt reset to ',lt
      end if
       if(mt.ge.mmax) then
       mt=mmax-1
      print*,'Warning: mt reset to ',mt
      end if
* set convenient values once:
      pi2 = 2.d0*dacos(-1.d0)
      rs=r**2
      muhaf=mu/2.d0
      musqh=mu*muhaf
      musqhprs=musqh*(1.d0+rs)
      ommusq=1.d0-2.d0*musqh
      ommusqrs=1.d0-2.d0*musqh*rs
      musq8=musqh/4.d0
      musq8rs=musq8*rs
      om2musq=1.d0-4.d0*musqh
      om2musqe=om2musq-2.d0*musqh*rs
      om2musrs=1.0d0-4.d0*musqh*rs
      ss=s**2
      muhafss=muhaf*ss
      muhafrs=muhaf*rs
      musqhss=musqh*ss
      musqhrs=musqh*rs
      musq8ss=musq8*ss
      ommss=1.d0-mu**2*ss
      om2mss=1.d0-2.d0*mu**2*ss
* initialise the arrays:
      do l=0,lt+1
       do m=0,mt+1
        do n=0,1
         pp(l,m,n)=0.d0
         qq(l,m,n)=0.d0
         mm(l,m,n)=0.d0
         nn(l,m,n)=0.d0
         gg(l,m,n)=0.d0
         hh(l,m,n)=0.d0
         dd(l,m,n)=1d0
         ee(l,m,n)=1d0
	 eei(l,m)= 0.d0
	 ddi(l,m)= 0.d0
        end do
       end do
      end do


*deprecated method below 
*insert initial values: nb assign to time point 1 - SOLVE will swap
*      do l=0,lt+1
*       pp(l,0,1)=p0
*       qq(l,0,1)=q0
*       mm(l,0,1)=0.d0
*       nn(l,0,1)=0.d0
*       gg(l,0,1)=g0
*       hh(l,0,1)=h0
*       dd(l,0,1)=d0
*      end do

C Initialise temperature
 	tt=temp
        ir = ionrt(tt)
	rr = ir
C work out the constant for the adiabatic law
	chi = temp/(1**(2.0/3.0))
	print*,'ir=',ir,'rr=',rr
	print*,'chi =',chi
	
       print*, '<<<<nr = ', nr, ' fr = ',fr,' >>>>>'
       print*, 'calling initial value ritout t=0'
       call ritout(0)
      end

      Subroutine solve
* implement the 2nd order lax wendroff scheme with additional coupling
      include 'alfgasfin_ortho.inc'
*     call ritout(0)
* loop over the time integration

      do n=1,nt

       if (n.gt.1) then
* copy previous new timestep to position 0:
       do l=0,lt+1
        do m=0,mt+1
         pp(l,m,0)=pp(l,m,1)
         qq(l,m,0)=qq(l,m,1)
         mm(l,m,0)=mm(l,m,1)
         nn(l,m,0)=nn(l,m,1)
         gg(l,m,0)=gg(l,m,1)
         hh(l,m,0)=hh(l,m,1)
         dd(l,m,0)=dd(l,m,1)
         ee(l,m,0)=ee(l,m,1)
        end do
       end do
       end if
* update the dummy points:

       call extrap

       call driver(n)	


** end of boundary updates
* now the main algorithm

       do m=1,mt
        do l=1,lt
        lp1=l+1
         lm1=l-1
         mp1=m+1
         mm1=m-1

         mm(l,m,1)=ommusq*mm(l,m,0)
     &   +muhaf*(pp(l,mp1,0)-pp(l,mm1,0))
     &   +musqh*(mm(l,mp1,0)+mm(l,mm1,0))
     &   -musq8*(nn(lp1,mp1,0)+nn(lm1,mm1,0)-nn(lm1,mp1,0)
     &    -nn(lp1,mm1,0))
     &   -musq8rs*(ee(lp1,mp1,0)+ee(lm1,mm1,0)-ee(lm1,mp1,0)
     &    -ee(lp1,mm1,0) )

         nn(l,m,1)=ommusq*nn(l,m,0)
     &   -muhaf*(pp(lp1,m,0)-pp(lm1,m,0))
     &   +musqh*(nn(lp1,m,0)+nn(lm1,m,0))
     &   +musqhrs*(ee(lp1,m,0)+ee(lm1,m,0)-2.d0*ee(l,m,0))
     &   -musq8*(mm(lp1,mp1,0)+mm(lm1,mm1,0)-mm(lm1,mp1,0)
     &     -mm(lp1,mm1,0))

         pp(l,m,1)=om2musqe*pp(l,m,0)
     &   +musqh*( pp(l,mp1,0)+pp(l,mm1,0))
     &   +musqhprs*( pp(lp1,m,0)+pp(lm1,m,0) )
     &   +muhaf*( nn(lm1,m,0)+mm(l,mp1,0)-nn(lp1,m,0)-mm(l,mm1,0) )
     &   -muhafrs*( ee(lp1,m,0)-ee(lm1,m,0) )
     &   +musq8rs*( qq(lp1,mp1,0)+qq(lm1,mm1,0)
     &                -qq(lp1,mm1,0)-qq(lm1,mp1,0) )
*        + coupling terms:
     &   + gm*(al*gg(l,m,0)-be*pp(l,m,0))

          qq(l,m,1)=ommusqrs*qq(l,m,0)
     &    +musqhrs*( qq(l,mp1,0)+qq(l,mm1,0) )
     &    -muhafrs*( ee(l,mp1,0)-ee(l,mm1,0))
     &    +musq8rs*(pp(lp1,mp1,0)+pp(lm1,mm1,0)
     &                -pp(lm1,mp1,0)-pp(lp1,mm1,0) )
     &       +gm*(al*hh(l,m,0)-be*qq(l,m,0))

         gg(l,m,1)=ommss*gg(l,m,0)
     &   +musqhss*(gg(lp1,m,0)+gg(lm1,m,0))
     &   +musq8ss*(hh(lp1,mp1,0)+hh(lm1,mm1,0)-hh(lp1,mm1,0)
     &     -hh(lm1,mp1,0))
     &   -muhafss*(dd(lp1,m,0)-dd(lm1,m,0))
*        + coupling terms:
     &   + gm*(be*pp(l,m,0)-al*gg(l,m,0))

         hh(l,m,1)=ommss*hh(l,m,0)
     &   +musqhss*(hh(l,mp1,0)+hh(l,mm1,0))
     &   +musq8ss*(gg(lp1,mp1,0)+gg(lm1,mm1,0)-gg(lp1,mm1,0)
     &     -gg(lm1,mp1,0))
     &   -muhafss*(dd(l,mp1,0)-dd(l,mm1,0))
*        + coupling terms:
     &   + gm*(be*qq(l,m,0)-al*hh(l,m,0))

         dd(l,m,1)=om2mss*dd(l,m,0)
     &   +musqhss*(dd(lp1,m,0)+dd(lm1,m,0)+dd(l,mp1,0)+dd(l,mm1,0))
     &   +muhaf*(gg(lm1,m,0)+hh(l,mm1,0)-gg(lp1,m,0)-hh(l,mp1,0))

         ee(l,m,1)=om2musrs*ee(l,m,0)
     &   + musqhrs*(ee(lp1,m,0)+ee(lm1,m,0)+ee(l,mp1,0)+ee(l,mm1,0))
     &   + musqh*(nn(lp1,m,0)+nn(lm1,m,0)-2.d0*nn(l,m,0))
     &   -muhaf*( pp(lp1,m,0)-pp(lm1,m,0)+qq(l,mp1,0)-qq(l,mm1,0) )
     &   -musq8*(mm(lp1,mp1,0)+mm(lm1,mm1,0)-mm(lp1,mm1,0)
     &      -mm(lm1,mp1,0))
        end do
       end do

C calculates a new temperature based on density fluctuations.
      if(ionise.eq.1) then

       do l = 1,lt
	do m =  1,mt


	tt = chi*(((ee(l,m,1)+dd(l,m,1))*0.5)**(2.0/3.0))

	ir = ionrt(tt)
        
	ioni = ir*(pdense)*(ee(l,m,1))*gdense*dd(l,m,1)*tstep
	reco = rr*(pdense**2)*(ee(l,m,1)**2)*tstep
	if((l.eq.100).and.(m.eq.110))print*,ee(l,m,1),dd(l,m,1),
     & (ioni-reco)/pdense,tt


C add alfven ionisation stuff here, for now NOT ENERGY CONSERVING
C	ener = (ioni)/pdense*ei
c	kine = 0.5*(pp(l,m,1)**2+qq(l,m,1)**2)
C	if(ener.lt.kine) then
	 ee(l,m,1)= ee(l,m,1)+(ioni-reco)/pdense
	 dd(l,m,1)= dd(l,m,1)+(-ioni+reco)/gdense
         
C	 kine = kine-ener
C	 dir = datan(pp(l,m,1)/qq(l,m,1)+0.01d-10)
C	if((l.eq.102).and.(m.eq.102))print*,pp(l,m,1),
C     &(ioni-reco)/pdense,dir
C	 pp(l,m,1) = sqrt(kine)*dsin(dir)
C	 qq(l,m,1) = sqrt(kine)*dcos(dir)
C
C	 end if
	 


 	 end do
        end do
       end if





       if(mod(n,nprint).eq.0) call ritout(n)

* end of main loop
      end do
*      call finish
* update the parameter(s) that change...
      fr=fr+dr
      end
*-------------------------------------------
      Subroutine driver(n)
      include 'alfgasfin_ortho.inc'
      
*
*this routine introduces any changes to the array that are not as a result of the
*finite difference algorithm, ie for driving waves etc.
*note that disturbance magnitudes should be read in as d0,e0,g0 etc for consistency
*for setting a line of gas in motion in the direction of the magnetic field
       if (n.lt.noff) then
        do m=95,105
         do l = 95,105
	  
	  dd(l,m,0) = 1.0d0+0.1d0*exp(-1.0d0*abs(dble(m-100)))*
     & dexp(-1.0d0*abs(dble(l-100)))*dsin(0.02d0*3.14159*n)
	
 	
	 end do
        end do
       end if



      end
*---------------------------------------

*** this code segment is for a driver that acts on the top driving gas density waves *against the *magnetic field direction
* overwrite any time-dependent boundary conditions:
*       do m=m1,m2
*
*        mdiff=m2-m1
*        msum=m2+m1
* first the harmonic driver term..now we're driving perturbations against the magnetic
* field direction, that is, perturbations down the edge
***>><<
*        zterm=dcos(dble(n-1)*pi2*fr/dble(1+nt))
* next the spatial structure in the z (or m) direction, with gaussian
* envelope to smooth gradients
*** we have a set of possibilities: either discrete sources, or a standing wave
*        xterm=(dsin(dble(m2-m)*pi2*kx/dble(mdiff))
*     &   +dabs(dsin(dble(m2-m)*pi2*kx/dble(mdiff))))
*     &   *dexp(-decx*(m-msum/2)**2)
*       xterm=1.d0

*        if (n.ge.noff) then
*         pp(0,m,1)=p0*zterm*xterm*dexp(-dect*(n-noff))
*         qq(0,m,1)=q0*zterm*xterm*dexp(-dect*(n-noff))
*         mm(0,m,1)=m0*zterm*xterm*dexp(-dect*(n-noff))
*         nn(0,m,1)=n0*zterm*xterm*dexp(-dect*(n-noff))
*         gg(0,m,1)=g0*zterm*xterm*dexp(-dect*(n-noff))
*         hh(0,m,1)=h0*zterm*xterm*dexp(-dect*(n-noff))
*         dd(0,m,1)=1.d0+d0*zterm*xterm*dexp(-dect*(n-noff))
*         ee(0,m,1)=1.d0+e0*zterm*xterm*dexp(-dect*(n-noff))
*        else
*         pp(0,m,1)=p0*zterm*xterm
*         qq(0,m,1)=q0*zterm*xterm
*         mm(0,m,1)=m0*zterm*xterm
*         nn(0,m,1)=n0*zterm*xterm
*         gg(0,m,1)=g0*zterm*xterm
*         hh(0,m,1)=h0*zterm*xterm
*         dd(0,m,1)=1.d0+d0*xterm*zterm
*         ee(0,m,1)=1.d0+e0*xterm*zterm
*        end if
*       end do	 
*-------------------------------------------
      Subroutine extrap
      include 'alfgasfin_ortho.inc'
Cupdate the dummy points at the l=lt+1,m=0, m=mt+1 edges:
C use quadratic extrapolant
C remember that the time-dept bc are set by driving term
       do m=0,mt
Cl=lt+1 case:
        pp(lt+1,m,0)=pp(lt-2,m,0)+3.d0*pp(lt,m,0)-3.d0*pp(lt-1,m,0)
        qq(lt+1,m,0)=qq(lt-2,m,0)+3.d0*qq(lt,m,0)-3.d0*qq(lt-1,m,0)
        mm(lt+1,m,0)=mm(lt-2,m,0)+3.d0*mm(lt,m,0)-3.d0*mm(lt-1,m,0)
        nn(lt+1,m,0)=nn(lt-2,m,0)+3.d0*nn(lt,m,0)-3.d0*nn(lt-1,m,0)
        gg(lt+1,m,0)=gg(lt-2,m,0)+3.d0*gg(lt,m,0)-3.d0*gg(lt-1,m,0)
        hh(lt+1,m,0)=hh(lt-2,m,0)+3.d0*hh(lt,m,0)-3.d0*hh(lt-1,m,0)
        dd(lt+1,m,0)=dd(lt-2,m,0)+3.d0*dd(lt,m,0)-3.d0*dd(lt-1,m,0)
        ee(lt+1,m,0)=ee(lt-2,m,0)+3.d0*ee(lt,m,0)-3.d0*ee(lt-1,m,0)
       end do
C m=0 case:
       do l=1,lt
        pp(l,0,0)=pp(l,3,0)+3.d0*pp(l,1,0)-3.d0*pp(l,2,0)
        qq(l,0,0)=qq(l,3,0)+3.d0*qq(l,1,0)-3.d0*qq(l,2,0)
        mm(l,0,0)=mm(l,3,0)+3.d0*mm(l,1,0)-3.d0*mm(l,2,0)
        nn(l,0,0)=nn(l,3,0)+3.d0*nn(l,1,0)-3.d0*nn(l,2,0)
        gg(l,0,0)=gg(l,3,0)+3.d0*gg(l,1,0)-3.d0*gg(l,2,0)
        hh(l,0,0)=hh(l,3,0)+3.d0*hh(l,1,0)-3.d0*hh(l,2,0)
        dd(l,0,0)=dd(l,3,0)+3.d0*dd(l,1,0)-3.d0*dd(l,2,0)
        ee(l,0,0)=ee(l,3,0)+3.d0*ee(l,1,0)-3.d0*ee(l,2,0)
       end do

C finally, the dummy points at the far end: m=mt+1
       do l=1,lt
        pp(l,mt+1,0)=pp(l,mt-2,0)+3.d0*pp(l,mt,0)-3.d0*pp(l,mt-1,0)
        qq(l,mt+1,0)=qq(l,mt-2,0)+3.d0*qq(l,mt,0)-3.d0*qq(l,mt-1,0)
        mm(l,mt+1,0)=mm(l,mt-2,0)+3.d0*mm(l,mt,0)-3.d0*mm(l,mt-1,0)
        nn(l,mt+1,0)=nn(l,mt-2,0)+3.d0*nn(l,mt,0)-3.d0*nn(l,mt-1,0)
        gg(l,mt+1,0)=gg(l,mt-2,0)+3.d0*gg(l,mt,0)-3.d0*gg(l,mt-1,0)
        hh(l,mt+1,0)=hh(l,mt-2,0)+3.d0*hh(l,mt,0)-3.d0*hh(l,mt-1,0)
        dd(l,mt+1,0)=dd(l,mt-2,0)+3.d0*dd(l,mt,0)-3.d0*dd(l,mt-1,0)
        ee(l,mt+1,0)=ee(l,mt-2,0)+3.d0*ee(l,mt,0)-3.d0*ee(l,mt-1,0)
       end do

C for corner points (l=lt+1,m=mt+1) and (l=lt+1,m=0) need to extrapolate using the
C last row of dummy points at l=lt+1:

        pp(lt+1,0,0)=pp(lt+1,3,0)+3.d0*pp(lt+1,1,0)-3.d0*pp(lt+1,2,0)
        qq(lt+1,0,0)=qq(lt+1,3,0)+3.d0*qq(lt+1,1,0)-3.d0*qq(lt+1,2,0)
        mm(lt+1,0,0)=mm(lt+1,3,0)+3.d0*mm(lt+1,1,0)-3.d0*mm(lt+1,2,0)
        nn(lt+1,0,0)=nn(lt+1,3,0)+3.d0*nn(lt+1,1,0)-3.d0*nn(lt+1,2,0)
        gg(lt+1,0,0)=gg(lt+1,3,0)+3.d0*gg(lt+1,1,0)-3.d0*gg(lt+1,2,0)
        hh(lt+1,0,0)=hh(lt+1,3,0)+3.d0*hh(lt+1,1,0)-3.d0*hh(lt+1,2,0)
        dd(lt+1,0,0)=dd(lt+1,3,0)+3.d0*dd(lt+1,1,0)-3.d0*dd(lt+1,2,0)
        ee(lt+1,0,0)=ee(lt+1,3,0)+3.d0*ee(lt+1,1,0)-3.d0*ee(lt+1,2,0)

        pp(lt+1,mt+1,0)=pp(lt+1,mt-2,0)+3.d0*pp(lt+1,mt,0)
     &   -3.d0*pp(lt+1,mt-1,0)
        qq(lt+1,mt+1,0)=qq(lt+1, mt-2,0)+3.d0*qq(lt+1,mt,0)
     &   -3.d0*qq(lt+1,mt-1,0)
        mm(lt+1,mt+1,0)=mm(lt+1, mt-2,0)+3.d0*mm(lt+1,mt,0)
     &   -3.d0*mm(lt+1,mt-1,0)
        nn(lt+1,mt+1,0)=nn(lt+1, mt-2,0)+3.d0*nn(lt+1,mt,0)
     &   -3.d0*nn(lt+1,mt-1,0)
        gg(lt+1,mt+1,0)=gg(lt+1, mt-2,0)+3.d0*gg(lt+1,mt,0)
     &   -3.d0*gg(lt+1,mt-1,0)
        hh(lt+1,mt+1,0)=hh(lt+1, mt-2,0)+3.d0*hh(lt+1,mt,0)
     &   -3.d0*hh(lt+1,mt-1,0)
        dd(lt+1,mt+1,0)=dd(lt+1, mt-2,0)+3.d0*dd(lt+1,mt,0)
     &   -3.d0*dd(lt+1,mt-1,0)
        ee(lt+1,mt+1,0)=ee(lt+1, mt-2,0)+3.d0*ee(lt+1,mt,0)
     &   -3.d0*ee(lt+1,mt-1,0)





        end
*---------------------------------------------------------
      subroutine ritout(n)
      include 'alfgasfin_ortho.inc'

      write(20,*)'# timestep = ',n
      do l=1,lt
      write(20,100)(pp(l,m,0),m=1,mt)
      end do
      write(20,*)
      write(20,*)

      write(21,*)'# timestep = ',n
      do l=1,lt
       write(21,100)(qq(l,m,0),m=1,mt)
      end do
      write(21,*)
      write(21,*)

      write(22,*)'# timestep = ',n
      do l=1,lt
      write(22,100)(gg(l,m,0),m=1,mt)
      end do
      write(22,*)
      write(22,*)

      write(23,*)'# timestep = ',n
      do l=1,lt
      write(23,100)(hh(l,m,0),m=1,mt)
      end do
      write(23,*)
      write(23,*)

      write(24,*)'# timestep = ',n
      do l=1,lt
      write(24,100)(mm(l,m,0),m=1,mt)
      end do
      write(24,*)
      write(24,*)

      write(25,*)'# timestep = ',n
      do l=1,lt
      write(25,100)(nn(l,m,0),m=1,mt)
      end do
      write(25,*)
      write(25,*)

*      write(26,*)'# timestep = ',n
*      do m=1,lt
** this is the debug file: output the driver
       write(26,101)(dd(0,m,0),m=1,mt)
*      write(26,101)(dd(mt/2+1,m,0)-1.d0,l=0,mt)
*      write(26,100)(ee(lt/2,m,0),m=0,mt)
*      end do
*      write(26,*)
*      write(26,*)


      write(27,*)'# timestep = ',n
      do l=1,lt
      write(27,101)(dd(l,m,0)+ddi(l,m),m=1,mt) 
      end do
      write(27,*)'# timestep = ',n
      write(27,*)
      write(27,*)
  
      write(28,*)'# timestep = ',n
      do l=1,lt
      write(28,101)(ee(l,m,0)+eei(l,m),m=1,mt)
      end do
      write(28,*)
      write(28,*)
101   format(220g13.5)
100   format(300f8.4)
      print*,'timestep = ',n
      end
c---------------------------------------------------------
      subroutine finish
      include 'alfgasfin_ortho.inc'
      do i=0,8
      write(20+i,*)'========================'
      write(20+i,datin)
      close(20+i)
      end do
      end
c-------------------------------------------------------
        function ionrt(T)
	INCLUDE 'alfgasfin_ortho.inc'
*This function computes the ionisation rate of hydrogen gas as a function of temperature. *it takes T in eV and outputs the rate coefficient in cm^3/s. This can be sued to *calculate ionrates if provided with time and densities. ionrate = ratecoeff*dd*gg
        DOUBLE PRECISION A,X,GIN,GIM,GIP,ttev
C change to eV from kelvin. 
	ttev = T/(11604.5)
	A = -1.

	X = 20.9D0/ttev
*GIM is the incomplete gamma function for use in the rate coefficient equation.
        CALL INCOG(A,X,GIN,GIM,GIP)

CThis is the rate coefficient equation, it is arrhenius-like r =  A*T^b*gammainc(-1,E/T)
C A = 8.067e-6, b = -1.41, E = 20.9eV
	IONRT = (8.067e-6)*(ttev**(-1.41))*GIM
10      FORMAT(1X,F5.1,1X,F5.1,3D15.8)
	RETURN
        END


        SUBROUTINE INCOG(A,X,GIN,GIM,GIP)
C
C       ===================================================
C       Purpose: Compute the incomplete gamma function
C                r(a,x), â(a,x) and P(a,x)
C       Input :  a   --- Parameter ( a ó 170 )
C                x   --- Argument 
C       Output:  GIN --- r(a,x)
C                GIM --- â(a,x)
C                GIP --- P(a,x)
C       Routine called: GAMMA for computing â(x)
C       ===================================================
C
        IMPLICIT DOUBLE PRECISION (A-H,O-Z)
        XAM=-X+A*DLOG(X)
        IF (XAM.GT.700.0.OR.A.GT.170.0) THEN
           WRITE(*,*)'a and/or x too large'
           STOP
        ENDIF
        IF (X.EQ.0.0) THEN
           GIN=0.0
           CALL GAMMA(A,GA)
           GIM=GA
           GIP=0.0
        ELSE IF (X.LE.1.0+A) THEN
           S=1.0D0/A
           R=S
           DO 10 K=1,60
              R=R*X/(A+K)
              S=S+R
              IF (DABS(R/S).LT.1.0D-15) GO TO 15
10         CONTINUE
15         GIN=DEXP(XAM)*S
           CALL GAMMA(A,GA)
           GIP=GIN/GA
           GIM=GA-GIN
        ELSE IF (X.GT.1.0+A) THEN
           T0=0.0D0
           DO 20 K=60,1,-1
              T0=(K-A)/(1.0D0+K/(X+T0))
20         CONTINUE
           GIM=DEXP(XAM)/(X+T0)
           CALL GAMMA(A,GA)
           GIN=GA-GIM
           GIP=1.0D0-GIM/GA
        ENDIF
        END


        SUBROUTINE GAMMA(X,GA)
C
C       ==================================================
C       Purpose: Compute gamma function â(x)
C       Input :  x  --- Argument of â(x)
C                       ( x is not equal to 0,-1,-2,úúú)
C       Output:  GA --- â(x)
C       ==================================================
	
	IMPLICIT DOUBLE PRECISION (A-H,O-Z)
        
        DIMENSION G(26)
        PI=3.141592653589793D0
        IF (X.EQ.INT(X)) THEN
           IF (X.GT.0.0D0) THEN
              GA=1.0D0
              M1=X-1
              DO 10 K=2,M1
10               GA=GA*K
           ELSE
              GA=1.0D+300
           ENDIF
        ELSE
           IF (DABS(X).GT.1.0D0) THEN
              Z=DABS(X)
              M=INT(Z)
              R=1.0D0
              DO 15 K=1,M
15               R=R*(Z-K)
              Z=Z-M
           ELSE
              Z=X
           ENDIF
           DATA G/1.0D0,0.5772156649015329D0,
     &          -0.6558780715202538D0, -0.420026350340952D-1,
     &          0.1665386113822915D0,-.421977345555443D-1,
     &          -.96219715278770D-2, .72189432466630D-2,
     &          -.11651675918591D-2, -.2152416741149D-3,
     &          .1280502823882D-3, -.201348547807D-4,
     &          -.12504934821D-5, .11330272320D-5,
     &          -.2056338417D-6, .61160950D-8,
     &          .50020075D-8, -.11812746D-8,
     &          .1043427D-9, .77823D-11,
     &          -.36968D-11, .51D-12,
     &          -.206D-13, -.54D-14, .14D-14, .1D-15/
           GR=G(26)
           DO 20 K=25,1,-1
20            GR=GR*Z+G(K)
           GA=1.0D0/(GR*Z)
           IF (DABS(X).GT.1.0D0) THEN
              GA=GA*R
              IF (X.LT.0.0D0) GA=-PI/(X*GA*DSIN(PI*X))
           ENDIF
        ENDIF
        RETURN
        END
