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
*   ga  real    coupling constant
*   nt  int     total number of timesteps
*   lt  int     total no of x-increments
*   mt  int     total no of z-increments
*   l   int     index for x-dirn
*   m   int     index for z-dirn (mag field dirn)
*   n   int     index for t stepping
*   dr  real    incremental change in r - now exploited as a frequency change

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
     &   + ga*(al*gg(l,m,0)-be*pp(l,m,0))

          qq(l,m,1)=ommusqrs*qq(l,m,0)
     &    +musqhrs*( qq(l,mp1,0)+qq(l,mm1,0) )
     &    -muhafrs*( ee(l,mp1,0)-ee(l,mm1,0))
     &    +musq8rs*(pp(lp1,mp1,0)+pp(lm1,mm1,0)
     &                -pp(lm1,mp1,0)-pp(lp1,mm1,0) )
     &       +ga*(al*hh(l,m,0)-be*qq(l,m,0))

         gg(l,m,1)=ommss*gg(l,m,0)
     &   +musqhss*(gg(lp1,m,0)+gg(lm1,m,0))
     &   +musq8ss*(hh(lp1,mp1,0)+hh(lm1,mm1,0)-hh(lp1,mm1,0)
     &     -hh(lm1,mp1,0))
     &   -muhafss*(dd(lp1,m,0)-dd(lm1,m,0))
*        + coupling terms:
     &   + ga*(be*pp(l,m,0)-al*gg(l,m,0))

         hh(l,m,1)=ommss*hh(l,m,0)
     &   +musqhss*(hh(l,mp1,0)+hh(l,mm1,0))
     &   +musq8ss*(gg(lp1,mp1,0)+gg(lm1,mm1,0)-gg(lp1,mm1,0)
     &     -gg(lm1,mp1,0))
     &   -muhafss*(dd(l,mp1,0)-dd(l,mm1,0))
*        + coupling terms:
     &   + ga*(be*qq(l,m,0)-al*hh(l,m,0))

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
*       if (n.lt.noff) then
*        do m=95,105
*         do l = 95,105
*	  gg(l,m,0) = d0
*	 end do
*        end do
*       end if

*make a circular velocity disturbance moving outwards from the centre of the grid.
*want a circle centered on 100,100 will use a radius of 3 initially in order to be symmetric
*velocity will be pointed outwards from the centre with all the same mag ie vx^2+vz^2 = 1

*have the disturbance last time noff can later add some kind of decay term both in x and in t
	if (n.lt.noff) then
*decay can be used for an envelope which evolves wrt t, x,z, etc.
*offset is the time between the central time and the current time
	 offset = abs(dble(n)-ncentre)
	 decay = exp(-decx/dble(offset))
*use prescribed points ie manually read in values for each point that we have decided lies on the *circle
*first points at 90 degrees have no mix of velocities, either all vx or vz
	 gg(97,100,0)=-g0*decay
	 hh(100,103,0)=h0*decay
	 gg(103,100,0)=g0*decay
	 hh(100,97,0)=-h0*decay
*points at other angles are more complicated involve max of changes on vx and vz
	 gg(97,101) = -0.9487*g0*decay
	 hh(97,101) = 0.3163*h0*decay
 
	 gg(97,99) = -0.9487*g0*decay
	 hh(97,99) = -0.3163*h0*decay

	 gg(99,103) = -0.3163*g0*decay
	 hh(99,103) = 0.9487*h0*decay

	 gg(101,103) = 0.3163*g0*decay
	 hh(101,103) = 0.9487*h0*decay

	 gg(103,101) = 0.9487*g0*decay
	 hh(103,101) = 0.3163*h0*decay

	 gg(103,99)  = 0.9487*g0*decay
	 hh(103,99)  =-0.3163*h0*decay

	 gg(101,97)  = 0.3163*g0*decay
	 hh(101,97)  =-0.9487*h0*decay

	 gg(99,97)   =-0.3163*g0*decay
	 hh(99,97)   =-0.9487*h0*decay
*finally the 45 degree points these mix equally between x and z
	 gg(98,102)  = -0.7071*g0*decay
	 hh(98,102)  =  0.7071*h0*decay
	 
	 gg(102,102)  = 0.7071*g0*decay
	 hh(102,102)  = 0.7071*h0*decay
	 
	 gg(102,98)  =  0.7071*g0*decay
	 hh(102,98)  = -0.7071*h0*decay
	 
	 gg(98,98)  = -0.7071*g0*decay
	 hh(98,98)  = -0.7071*h0*decay
	end if
      end
*---------------------------------------

*** this code segment is for a driver that acts on the top driving gas density waves against the *magnetic field direction
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
* update the dummy points at the l=lt+1,m=0, m=mt+1 edges:
* use quadratic extrapolant
* remember that the time-dept bc are set by driving term
       do m=0,mt
*l=lt+1 case:
        pp(lt+1,m,0)=pp(lt-2,m,0)+3.d0*pp(lt,m,0)-3.d0*pp(lt-1,m,0)
        qq(lt+1,m,0)=qq(lt-2,m,0)+3.d0*qq(lt,m,0)-3.d0*qq(lt-1,m,0)
        mm(lt+1,m,0)=mm(lt-2,m,0)+3.d0*mm(lt,m,0)-3.d0*mm(lt-1,m,0)
        nn(lt+1,m,0)=nn(lt-2,m,0)+3.d0*nn(lt,m,0)-3.d0*nn(lt-1,m,0)
        gg(lt+1,m,0)=gg(lt-2,m,0)+3.d0*gg(lt,m,0)-3.d0*gg(lt-1,m,0)
        hh(lt+1,m,0)=hh(lt-2,m,0)+3.d0*hh(lt,m,0)-3.d0*hh(lt-1,m,0)
        dd(lt+1,m,0)=dd(lt-2,m,0)+3.d0*dd(lt,m,0)-3.d0*dd(lt-1,m,0)
        ee(lt+1,m,0)=ee(lt-2,m,0)+3.d0*ee(lt,m,0)-3.d0*ee(lt-1,m,0)
       end do
* m=0 case:
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

* finally, the dummy points at the far end: m=mt+1
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

*** for corner points (l=lt+1,m=mt+1) and (l=lt+1,m=0) need to extrapolate using the
*** last row of dummy points at l=lt+1:

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
      write(27,101)(dd(l,m,0),m=1,mt) 
      end do
      write(27,*)'# timestep = ',n
      write(27,*)
      write(27,*)
  
      write(28,*)'# timestep = ',n
      do l=1,lt
      write(28,101)(ee(l,m,0),m=1,mt)
      end do
      write(28,*)
      write(28,*)
101   format(220g13.5)
100   format(300f8.4)
      print*,'timestep = ',n
      end
*---------------------------------------------------------
      subroutine finish
      include 'alfgasfin_ortho.inc'
      do i=0,8
      write(20+i,*)'========================'
      write(20+i,datin)
      close(20+i)
      end do
      end
*-------------------------------------------------------
      function gamma_inc ( p, x )

!*****************************************************************************80
!
!! GAMMA_INC computes the incomplete Gamma function.
!
!  Discussion:
!
!    GAMMA_INC(P,       0) = 0,
!    GAMMA_INC(P,Infinity) = 1.
!
!    GAMMA_INC(P,X) = Integral ( 0 <= T <= X ) T**(P-1) EXP(-T) DT / GAMMA(P).
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    01 May 2001
!
!  Author:
!
!    Original FORTRAN77 version by B L Shea.
!    FORTRAN90 version by John Burkardt
!
!  Reference:
!
!    BL Shea,
!    Chi-squared and Incomplete Gamma Integral,
!    Algorithm AS239,
!    Applied Statistics,
!    Volume 37, Number 3, 1988, pages 466-473.
!
!  Parameters:
!
!    Input, real ( kind = 8 ) P, the exponent parameter.
!    0.0D+00 < P.
!
!    Input, real ( kind = 8 ) X, the integral limit parameter.
!    If X is less than or equal to 0, GAMMA_INC is returned as 0.
!
!    Output, real ( kind = 8 ) GAMMA_INC, the value of the function.
!
      implicit none

      real ( kind = 8 ) a
      real ( kind = 8 ) arg
      real ( kind = 8 ) b
      real ( kind = 8 ) c
      real ( kind = 8 ) cdf
      real ( kind = 8 ), parameter :: exp_arg_min = -88.0D+00
      real ( kind = 8 ) gamma_inc
      real ( kind = 8 ) gamma_log
      real ( kind = 8 ), parameter :: overflow = 1.0D+37
      real ( kind = 8 ) p
      real ( kind = 8 ), parameter :: plimit = 1000.0D+00
      real ( kind = 8 ) pn1
      real ( kind = 8 ) pn2
      real ( kind = 8 ) pn3
      real ( kind = 8 ) pn4
      real ( kind = 8 ) pn5
      real ( kind = 8 ) pn6
      real ( kind = 8 ) rn
      real ( kind = 8 ), parameter :: tol = 1.0D-07
      real ( kind = 8 ) x
      real ( kind = 8 ), parameter :: xbig = 1.0D+08

      gamma_inc = 0.0D+00

      if ( p <= 0.0D+00 ) then
         write ( *, '(a)' ) ' '
         write ( *, '(a)' ) 'GAMMA_INC - Fatal error!'
         write ( *, '(a)' ) '  Parameter P <= 0.'
         stop
      end if

      if ( x <= 0.0D+00 ) then
         gamma_inc = 0.0D+00
         return
      end if
!
!  Use a normal approximation if PLIMIT < P.
!
      if ( plimit < p ) then
      pn1 = 3.0D+00 * sqrt ( p ) * ( ( x / p ) ** (1.0D+00 / 3.0D+00)   &
      + 1.0D+00 / ( 9.0D+00 * p ) - 1.0D+00 )
 
         call normal_01_cdf ( pn1, cdf )
         gamma_inc = cdf
         return
      end if
!
!  Is X extremely large compared to P?
!
      if ( xbig < x ) then
         gamma_inc = 1.0D+00
         return
      end if
!
!  Use Pearson's series expansion.
!  (P is not large enough to force overflow in the log of Gamma.
!
      if ( x <= 1.0D+00 .or. x < p ) then

         arg = p * log ( x ) - x - gamma_log ( p + 1.0D+00 )
         c = 1.0D+00
         gamma_inc = 1.0D+00
         a = p

         do

            a = a + 1.0D+00
            c = c * x / a
            gamma_inc = gamma_inc + c
            
            if ( c <= tol ) then
               exit
            end if

         end do

         arg = arg + log ( gamma_inc )

         if ( exp_arg_min <= arg ) then
            gamma_inc = exp ( arg )
         else
            gamma_inc = 0.0D+00
      end if

      else
!
!  Use a continued fraction expansion.
!
         arg = p * log ( x ) - x - gamma_log ( p )
         a = 1.0D+00 - p
         b = a + x + 1.0D+00
         c = 0.0D+00
         pn1 = 1.0D+00
         pn2 = x
         pn3 = x + 1.0D+00
         pn4 = x * b
         gamma_inc = pn3 / pn4

         do

            a = a + 1.0D+00
            b = b + 2.0D+00
            c = c + 1.0D+00
            pn5 = b * pn3 - a * c * pn1
            pn6 = b * pn4 - a * c * pn2

            if ( 0.0D+00 < abs ( pn6 ) ) then

               rn = pn5 / pn6

               if ( abs ( gamma_inc - rn ) <= min ( tol, tol * rn ) ) then

                  arg = arg + log ( gamma_inc )

                  if ( exp_arg_min <= arg ) then
             gamma_inc = 1.0D+00 - exp ( arg )
          else
             gamma_inc = 1.0D+00
          end if

          return

       end if

       gamma_inc = rn
       
       end if

       pn1 = pn3
       pn2 = pn4
       pn3 = pn5
       pn4 = pn6
!
!  Rescale terms in continued fraction if terms are large.
!
       if ( overflow <= abs ( pn5 ) ) then
         pn1 = pn1 / overflow
         pn2 = pn2 / overflow
         pn3 = pn3 / overflow
         pn4 = pn4 / overflow
       end if

      end do

      end if

      return
      end



