	PROGRAM gamma
	gammainc = gammp(-0.9,10)
	print*, gammainc


	end

*----------------------	

	FUNCTION gammp(a,x)
	REAL a,gammp,x
*	USE gcf,gser
**Returns the incomplete gamma function P (a, x).
	REAL gammcf,gamser,gln
	if(x.lt.0..or.a.le.0.) then
	print*  ,'error in arguments'
	endif
	if(x.lt.a+1.)then
*Use the series representation.
	call gser(gamser,a,x,gln)
	gammp=gamser
	else
*Use the continued fraction representation
	call gcf(gammcf,a,x,gln)
	gammp=1.-gammcf
*and take its complement.
	endif
	return
	END

*----------------------------------------------------------
	FUNCTION gammq(a,x)
	REAL a,gammq,x
*	USE gcf,gser
*Returns the incomplete gamma function Q(a, x) ≡ 1 − P (a, x).
	REAL gammcf,gamser,gln
	if(x.lt.0..or.a.le.0.)then
	print*, 'bad arguments in gammq'
	end if
	if(x.lt.a+1.)then
*Use the series representation
	call gser(gamser,a,x,gln)
	gammq=1.-gamser
*and take its complement.
	else
*use the continued fraction representation.
	call gcf(gammcf,a,x,gln)
	gammq=gammcf
	endif
	return
	END
	
*---------------------------------------------------------
	SUBROUTINE gser(gamser,a,x,gln)
*	USE gammln
	INTEGER ITMAX
	REAL a,gamser,gln,x,EPS
	PARAMETER (ITMAX=1000000,EPS=3.e-7)

*Returns the incomplete gamma function P (a, x) evaluated by its 
*series representation as
*gamser. Also returns ln Γ(a) as gln.
	INTEGER n
	REAL ap,del,sum,gammln
	gln=gammln(a)
	if(x.le.0.)then
	if(x.lt.0.)then 
	print*, 'x < 0 in gser'
	endif
	gamser=0.
	return
	endif
	ap=a
	sum=1./a
	del=sum
	do n=1,ITMAX
	ap=ap+1.
	del=del*x/ap
	sum=sum+del
	if(abs(del).lt.abs(sum)*EPS)then 
	goto 1
	end if 
	end do
	print*, 'a too large, ITMAX too small in gser'
 1	gamser=sum*exp(-x+a*log(x)-gln)
	return
	END

*-------------------------------------------------------------
	SUBROUTINE gcf(gammcf,a,x,gln)
*	USE gammln
	INTEGER ITMAX
	REAL a,gammcf,gln,x,EPS,FPMIN
	PARAMETER (ITMAX=1000000,EPS=3.e-7,FPMIN=1.e-30)

*	Returns the incomplete gamma function Q(a, x) evaluated by 
*its continued fraction *repre-
*	sentation as gammcf. Also returns ln Γ(a) as gln.
*	Parameters: ITMAX is the maximum allowed number of iterations; EPS is the *relative accu-
*		racy; FPMIN is a number near the smallest representable floating-point *number.
	INTEGER i
	REAL an,b,c,d,del,h,gammln
	gln=gammln(a)
	b=x+1.-a
*	Set up for evaluating continued fraction by modified
*	Lentz’s method (§5.2) with b0 = 0.
	c=1./FPMIN
	d=1./b
	h=d
	do i=1,ITMAX
*Iterate to convergence.
		an=-i*(i-a)
		b=b+2.
		d=an*d+b
		if(abs(d).lt.FPMIN)d=FPMIN
		c=b+an/c
		if(abs(c).lt.FPMIN)c=FPMIN
		d=1./d
		del=d*c
		h=h*del
		if(abs(del-1.).lt.EPS)goto 1
	end do
	print*, 'a too large, ITMAX too small in gcf'
 1	gammcf=exp(-x+a*log(x)-gln)*h
*Put factors in front.
	return
	END
	FUNCTION gammln(xx)
	REAL gammln,xx
*Returns the value ln[Γ(xx)] for xx > 0.
	INTEGER j
	DOUBLE PRECISION ser,stp,tmp,x,y,cof(6)
*Internal arithmetic will be done in double precision, a nicety that you can omit if five-*figure
*accuracy is good enough.

     	cof(1) = 76.1800917294715 
	cof(2) =-86.5053203294168
        cof(3)=24.0140982408309 
	cof(4)=-1.23173957245015    
        cof(5)=1.208650973866179d-003 
	cof(6)= -5.395239384953000d-006
	stp= 2.5066282746310005d0
	x=xx
	y=x
	tmp=x+5.5d0
	tmp=(x+0.5d0)*log(tmp)-tmp
	ser=1.000000000190015d0
	do j=1,6
	y=y+1.d0
	ser=ser+cof(j)/y
	end do 
	gammln=tmp+log(stp*ser/x)
	return
	END


