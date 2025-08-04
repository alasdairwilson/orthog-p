      Implicit none
      Integer
     & nt,lt,mt,l,m,n,nprint,lmax,mmax,lp1,mp1,lm1,mm1,m1,m2,noff,
     & mdiff,msum,nr,i,ionise
      Parameter(lmax=250, mmax=250)
      Character filnam(10)*9
      Double Precision
     & pp(0:lmax,0:mmax,0:1),qq(0:lmax,0:mmax,0:1),
     & mm(0:lmax,0:mmax,0:1),nn(0:lmax,0:mmax,0:1),
     & gg(0:lmax,0:mmax,0:1),hh(0:lmax,0:mmax,0:1),
     & dd(0:lmax,0:mmax,0:1),ee(0:lmax,0:mmax,0:1),
     & ddi(0:lmax,0:mmax,0:1),eei(0:lmax,0:mmax,0:1),
     & mu,s,ss,al,be,gm,p0,q0,m0,n0,kx,kt,g0,h0,d0,e0,muhaf,musqh,
     & musq8,dr,decay,offset,ncentre,gdense,pdense,temp,ir,rr,tt,
     & om2musq,muhafss,musqhss,musq8ss,ommusq,om2mss,pi2,fr,ionrt,
     & xterm,zterm,xtermd,psav,ommss,decx,dect,musqhprs,ommusqrs,rs,
     & om2musqe,om2musrs,muhafrs,musq8rs,musqhrs,r,T,chi,tstep,dedt,
     & dddt, ioni,reco,emin,dir,ener,frac,vrel(2),theta,vreln(2),
     & vchange(2),vmag,totali
     
      namelist/datin/lt,mt,nt,mu,nprint,s,p0,q0,m0,n0,kx,g0,h0,d0,e0,
     & al,be,gm,fr,m1,m2,noff, decx,dect,r,dr,ncentre,gdense,pdense
     & ,temp,tstep,ionise,frac,emin
      common/vars/pp,qq,mm,nn,gg,hh,dd,ee,tt,ir,rr,chi,eei,ddi,dir,
     & emin,ener,totali
      common/handy/mu,s,al,be,gm,p0,q0,m0,n0,kx,kt,g0,h0,d0,e0,r,dr,
     & ncentre,tstep,gdense,pdense,frac,ionise,
     & muhaf,musqh,musq8,om2musq,muhafss,musqhss,musq8ss,musqhprs,
     & ommusq,om2mss,pi2,fr,nprint,lt,mt,nt,m1,m2,noff,psav,ommss,nr,
     & ommusqrs,decx,dect,rs,om2musqe, om2musrs,muhafrs,musq8rs,
     & musqhrs, filnam
