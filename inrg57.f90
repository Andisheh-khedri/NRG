!========================================================================
!*In comparison to last version (21), the fixed-number-of-states 
! Scheme has been added! 
! *********************************************************************
! Implementation of z-averaging 
! *********************************************************************
! parallel version :D (for arbitrary T grid)
! izmax=number of processors we are going to use
! *********************************************************************
! compard to inrg39.f90 mpi has been used instead of mpi-send-recieve
! *********************************************************************
! compard to inrg42.f90 saving the informations for each chain-lenght adelta
!        and then calculating thermodynamical quantities!
!========================================================================
!========================================================================
program IRLMNRG

use mpmodule

implicit none

!=============================================================================
!----------------------------------------------------------------------
!--VARIABLE DECLARATION:-----------------------------------------------
!----------------------------------------------------------------------
integer, parameter:: mmax=40,izmax=4,itmaxdim=1000,nchannel=1,echarge=7,iw1maxdim=2*100,Nbtotal=40,ie0max=100&
                     ,tnumstates=1300,numgap=0
double precision,PARAMETER ::pi = 3.1415926535897931d+0
integer*8::m,n,i,j,s,k,ss,kk,openstatus,nmax,d0,dbc,dmax,dof,icut,icdint&
          ,nmin,imin,drt,eint,mincheck,grss,maxsort&
          ,n1,i1,nexcitedprint,maxsorting,deltam,maxdeltam&
          ,cdint,nexprint,mmin,rnmax,izcount,iw1,iw1prime,nb&
          ,it,itcounter,start,itmaximum,drnmax,discounter,nii,nup,ndown,degeneracy,iw1max,ilambda
double precision::l,e,Vl,Vr,u,width,gama,ans1,ans2,ans3,ansn,emin,asumFcorr,sumFcorr,justsumFcorr&
                  ,ex1,tmr,tml,alpha,area1,area2,lambda,w0,frganagama,mapgef,gefapprox&
                  ,barbeta,Etot,EEtot,Ecut,error,Gre,alphatem,bfactor,dfactor,alphaw1
double precision,allocatable,dimension(:,:,:)::h,Trb,nd,ndb,ddagger,dbdagger,nboson,nbeboson,bpbdagger,bebpbdagger                                    
integer*8,allocatable,dimension(:)::d,d1,d2,d3,d4,dr,db1,db2,db3,db4
double precision,allocatable,dimension(:)::Di,shb,egap
integer*8,allocatable,dimension(:)::snarray,siarray
double precision,allocatable,dimension(:,:)::hm,hb
double precision,dimension(2,itmaxdim)::Hav,HHav,BEntropy,zpar,zparcheck
double precision,dimension(2,itmaxdim)::Netot,NEtotcheck
double precision::zavscimp,I0collector,ndcollector,ndzav,zavspectral,nddcollector,nddzav,nbocollector&
                  ,gcollector2,collector3,collector4,collector5,collector6,collector7,collector8,collector9,collector10,eta
double precision,dimension(2,itmaxdim)::scimp
double precision,dimension(itmaxdim)::ndtot,nTboson,ndtotcheck

double precision,dimension(1:ie0max,itmaxdim)::ThEke,ThES,ThZT0,ThEL,I0mz,I1mz,I2mz
double precision,dimension(0:izmax,1:ie0max,itmaxdim)::I0m,I1m,I2m
integer::dtotal
!=============================================================================
! Tem grid:
double precision,dimension(itmaxdim)::Tem,beta
integer,dimension(itmaxdim)::Tshell
integer,dimension(0:mmax)::itmax,itmin
!=============================================================================
! Storing Info
integer,dimension(1:2,0:mmax,0:((2*mmax)+3))::dstore,drstore&
                                            ,d1store,d2store,d3store,d4store
double precision,dimension(0:mmax,0:((nchannel*(mmax+1))+1),1:4000)::Hstore,nbomat,matndbpbmre,bmat,bdagmat,ndmat,bpbnurmat
double precision,allocatable,dimension(:,:,:,:)::Trstore,ddmat,bpbmat
integer,dimension(1:2)::m0
double precision,dimension(1:2,0:mmax)::eground
integer,dimension(1:2,0:mmax)::nground
integer,dimension(1:4,1:nchannel*(mmax+1)+1)::nupmat,ndownmat,ssmat
integer,dimension(1:4)::niimat
double precision,dimension(1:itmaxdim,1:iw1maxdim)::spectral,sumrule,aspectral,newspectral,Fcorr,aFcorr,anewspectral&
                                                    ,zSpectral,zaspectral,zFcorr,zaFcorr
double precision,dimension(1:iw1maxdim)::w1,deltaw1,adeltaw1
integer,dimension(1:iw1maxdim)::w1shell
double precision,dimension(0:mmax)::zpartition
double precision,dimension(0:1,0:Nbtotal)::H0local
double precision,dimension(0:1,0:Nbtotal,0:Nbtotal)::Tr0local
double precision,dimension(1:ie0max)::e0
double precision,dimension(1:ie0max,1:itmaxdim)::nde0
double precision,dimension(0:izmax,1:ie0max,1:itmaxdim)::ndde0
double precision,dimension(1:itmaxdim)::efgamma,efgammadd
double precision,allocatable,dimension(:)::Htotal

complex*16,dimension(1:itmaxdim,1:iw1maxdim)::phselfenergy,DotgR,MatDotg,MatFcorr,Matphselfenergy,aDotgR,aphselfenergy&
                                            ,zMatDotgRe,zMatFcorrRe,zMatDotgIm,zMatFcorrIm
complex*16::ic
!=============================================================================
! Storing the discarded states:
!double precision,dimension(0:mmax,0:((nchannel*(mmax+1))+1),1:4000)::
!=============================================================================
! parameters for the lapack lib (DSYEV)
!===============================================================
character::JoBz,UPLO
integer*8::LWMAX
integer*8::LDA,LWORK,INFO
double precision,allocatable,dimension(:)::WORK
!===============================================================
! parameters for the blas lib (dgemm)
!===============================================================
double precision::alphaconst,betaconst
double precision,allocatable,dimension(:,:)::Amatrix,Bmatrix,Dmatrix
!===============================================================
! parameters for Tridiagonalization SUBROUTINE (tri)
!===============================================================
type(mp_real)::lmulti,gamamulti,zmulti
type(mp_real)::tm(0:mmax),epsilonm(0:mmax)
!===============================================================
! parameters for MPI:
!===============================================================
integer::ierr,rank,size,partner,iz
include 'mpif.h'
integer status(MPI_STATUS_SIZE)
!----------------------------------------------------------------------
!---------------------------------------------------------------------
real :: start_time, stop_time
!---------------------------------------------------------------------
!allocate(Trstore(0:mmax,1:2*echarge+1,1:4000,1:1000))
allocate(ddmat(0:mmax,1:2*(mmax+1)+1,1:500,1:500),bpbmat(0:mmax,1:2*(mmax+1)+1,1:500,1:500))
!-----------------------------------------------------------------------
call cpu_time(start_time)
!-----------------------------------------------------------------------
ic=cmplx(0.0d+0,1.0d+0, kind = 16)
!-----------------------------------------------------------------------
!parameters of the model:
!-----------------------------------------------------------------------
lmulti=mpreal(4.0d0)        !scale parameter
l=dble(lmulti)
alphatem=1.0d+00
!-----------------------------------------------------------------------
!computational:
!-----------------------------------------------------------------------                                          
maxsorting=700
maxdeltam=20
nexcitedprint=5  
!barbeta=0.7    !0.46d+0  !0.727d+0    !0.46d+0  !0.88!
!-----------------------------------------------------------------------
! Initialization of Parallelization:
!-----------------------------------------------------------------------
call mpi_init(ierr)
call mpi_comm_rank(MPI_COMM_WORLD,rank,ierr)
call mpi_comm_size(MPI_COMM_WORLD,size,ierr) 
!------------------------------------------------------------------------------------------------ 
Do ilambda=1,1
!-----------------------------------------------------------------------
Do eint=51,53
!----------------------------------------------------------------------- 
!------------------------------------------------------------------------------------------------
if(rank==0)then
!------------------------------------------------------------------------------------------------
open(unit=2,file="k.txt",status="new",action="write",position="rewind",IOSTAT=openstatus)
!------------------------------------------------------------------------------------------------
endif
!------------------------------------------------------------------------------------------------
!------------------------------------------------------------------------------------------------
!iz=rank
!-----------------------------------------------------------------------
zSpectral=0.0d+00
zFcorr=0.0d+00
zaSpectral=0.0d+00
zaFcorr=0.0d+00
zMatDotgRe=0.0d+00
zMatDotgIm=0.0d+00
zMatFcorrRe=0.0d+00
zMatFcorrIm=0.0d+00
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
Do iz=0,izmax-1
!-----------------------------------------------------------------------
!print*,"Rank=",iz
!-----------------------------------------------------------------------
!zmulti=mpreal((2*(iz+1)-1)*1.0d0)/mpreal(izmax*2.0d0)
zmulti=mpreal((iz+1)*1.0d0)/mpreal(izmax*1.0d0)
print*,iz,dble(zmulti)
!-----------------------------------------------------------------------
!physical parameters:
!-----------------------------------------------------------------------
width=1.0d+00                               !Bandwidth
gamamulti=mpreal(1.0d0)/mpreal(10000.0d00)  !Hybridization function
gama=dble(gamamulti)
Vl=sqrt((2.0d+0)*(gama*width)/(pi))         !T-amplitude to left  lead
Vr=sqrt((2.0d+0)*(gama*width)/(pi))         !T-amplitude to right lead                         
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
eta=gama/(10000000000.0d+00)
u=0.0d+00!gama/10000.0d+0!gama/10000.0
!u=2*((10.0d+0)**(-3.0d+0+(eint-1)*6.0d+0/20.0d+0))
!-----------------------------------------------------------------------
w0=20.0d+00*gama!0.00025d+00!1.0d+00*gama
lambda=0.5d+00*w0!w0*((ilambda-1)*3.0d+00/40.0d+00)!0.2d+00*w0
!-------------------------------------------------------------------------------------
!-------------------------------------------------------------------------------------
! approximated effective Tunneling rate from FRG (flow eq for effective mass): 
!-------------------------------------------------------------------------------------
frganagama=gama*(exp(((-((lambda/w0)**2))*(((4.0d+00*gama/(pi*w0))*log(gama/w0))&
                                           +(1.0d+00-(gama/w0)**2)&
                                           +((2.0d+00*gama/(w0*pi))*(1.0d+00+(gama/w0)**2))))&
                        /((1.0d+00+(gama/w0)**2)**2)))
!-------------------------------------------------------------------------------------
! approximated effective Tunneling rate from IRLM-mapping: 
!-------------------------------------------------------------------------------------------------
mapgef=w0*(((gama/w0)*(exp(-((lambda/w0)**2))))**(1.0d+00/(1.0d+00+(4.0d+00*gama/(pi*w0))*((w0/gama)**2))))
!-------------------------------------------------------------------------------------------------
gefapprox=min(frganagama,gama)
!-------------------------------------------------------------------------------------------------
print*,"w0/Gamma=",w0/gama
print*,"approximated effective Tunneling rate from FRG=",frganagama/gama,"gamma"
print*,"approximated effective Tunneling rate from IRLM-maping=",mapgef/gama
!-------------------------------------------------------------------------------------------------
!e=0.0d+00!*gama  !40*gama*((( eint-51.0)*exp((abs(51.0- eint)-50.0)/23)/50.0)) !Impurity energy   
e=0.0d+00*(lambda**2)/w0+30.0d+00*0.1d+00*(gefapprox)*(((eint-51.0)*exp((abs(51.0- eint)-50.0)/23.0d+00)/50.0d+00))   !Impurity energy

if(eint==53)then
e=-1.0d+00*gama+0.0d+00*(lambda**2)/w0!
endif
!e=gama*(10.0d+00**(-2.0d+00+(4.0d+00*(ilambda-1)/40.0d+00)))
!e=0.0d+00
e0(eint)=e                           
!-----------------------------------------------------------------------
!HHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHH
!-----------------------------------------------------------------------
! Hoppings:
!-----------------------------------------------------------------------
call tri('alt',lmulti,gamamulti,zmulti,mmax,tm,epsilonm)
!-----------------------------------------------------------------------
!HHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHH
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
!IIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIII
!-----------------------------------------------------------------------
! Iterative Diagonalization:
!-----------------------------------------------------------------------
!IIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIII
!-----------------------------------------------------------------------
Do icdint=0,0
!-----------------------------------------------------------------------
cdint=2-icdint
itcounter=0
!-----------------------------------------------------------------------
if(cdint==1)then
mmin=1
else
mmin=0
!-----------------------------------------------------------------------
!  m=-1 : the local Hamiltonian:
!-----------------------------------------------------------------------
Do n=0,1
!-----------------------------------------------------------------------
allocate(hm(1:Nbtotal+1,1:Nbtotal+1),Di(1:Nbtotal+1))
!-----------------------------------------------------------------------
Do s=1,Nbtotal+1
Do k=1,Nbtotal+1
if(k==s)then
hm(s,k)=(n-1.0d+00*0.50d+00)*e+(s-1)*w0
elseif(k==s-1)then
hm(s,k)=(n-1.0d+00*0.50d+00)*lambda*sqrt(1.0d+00*(k))
elseif(k==s+1)then
hm(s,k)=(n-1.0d+00*0.50d+00)*lambda*sqrt(1.0d+00*(s))
else
hm(s,k)=0.0d+00
endif
end do !k
!print*,(hm(s,k),k=1,Nbtotal+1)
end do !s

!-----------------------------------------------------------------------
LWORK=max(1,3*(Nbtotal+1)-1)
allocate(WORK(max(1,lWORK)))
!-----------------------------------------------------------------------
LDA=Nbtotal+1

if(LDA>1)then

CALL DSYEV( 'Vectors', 'Upper', Nbtotal+1, hm, LDA, Di, WORK, LWORK, INFO )

else if(LDA==1)then
hm(1:Nbtotal+1,1:Nbtotal+1)=1.0d+00
Di(1:Nbtotal+1)=(n-0.50d+00)*e
endif

H0local(n,0:Nbtotal)=Di(1:Nbtotal+1)    
Tr0local(n,0:Nbtotal,0:Nbtotal)=hm(1:Nbtotal+1,1:Nbtotal+1)

deallocate(hm,Di,WORK)

end do !n

!-----------------------------------------------------------------------

endif
!-----------------------------------------------------------------------
discounter=0
!-----------------------------------------------------------------------
Do m=mmin,mmax-1
!-----------------------------------------------------------------------
if(nchannel==2)then
  if(cdint==1)then
  nmax=2*m+2
  !nmax=m+1  !1cband
  elseif(cdint==2)then
  nmax=2*m+3
  !nmax=m+2  !cband
  endif
elseif(nchannel==1)then
  if(cdint==1)then
  nmax=m+1  
  elseif(cdint==2)then
  nmax=m+2  
  endif
endif  
!-----------------------------------------------------------------------
if(m>0)then
! tmr=l**(-(m-1)/2.0d+0)
! tml=l**(-(m-1)/2.0d+0)
tmr=dble(tm(m-1))
tml=tmr
else
tmr=1.0d+00!l**(-(m-1)/2.0d+0)
tml=1.0d+00!l**(-(m-1)/2.0d+0)
endif
!-----------------------------------------------------------------------
drt=0  !Total reduced dimension in each steps
!-----------------------------------------------------------------------
allocate(d(0:nmax),d1(0:nmax),d2(0:nmax),d3(0:nmax),d4(0:nmax))
!-----------------------------------------------------------------------
if(m==mmin)then
!-----------------------------------------------------------------------
if(cdint==2)then
!-----------------------------------------------------------------------
allocate(dr(0:1),ndb(0:1,1:(Nbtotal+1),1:(Nbtotal+1)),dbdagger(0:0,1:Nbtotal+1,1:Nbtotal+1),nbeboson(0:1,1:Nbtotal+1,1:Nbtotal+1))
allocate(bebpbdagger(0:0,1:Nbtotal+1,1:Nbtotal+1))
dr(0)=1*(Nbtotal+1)
dr(1)=1*(Nbtotal+1)
ndb=0.0d+00
dbdagger=0.0d+00
bebpbdagger=0.0d+00
nbeboson=0.0d+00
ndb(0,1:(Nbtotal+1),1:(Nbtotal+1))=0.0d+0

Do i=1,Nbtotal+1
ndb(1,i,i)=1.0d+00
Do j=1,Nbtotal+1
Do nb=0,Nbtotal
dbdagger(0,i,j)=dbdagger(0,i,j)+Tr0local(1,nb,i-1)*Tr0local(0,nb,j-1)
nbeboson(0,i,j)=nbeboson(0,i,j)+nb*Tr0local(0,nb,i-1)*Tr0local(0,nb,j-1)
nbeboson(1,i,j)=nbeboson(1,i,j)+nb*Tr0local(1,nb,i-1)*Tr0local(1,nb,j-1)
if(nb>0)then
bebpbdagger(0,i,j)=bebpbdagger(0,i,j)+Tr0local(1,nb-1,i-1)*Tr0local(0,nb,j-1)*sqrt(nb*1.0d+00)&
                                     +Tr0local(1,nb,i-1)*Tr0local(0,nb-1,j-1)*sqrt(nb*1.0d+00)

endif
end do !nb
end do !j
end do !i
!-----------------------------------------------------------------------
else
!-----------------------------------------------------------------------
if(nchannel==2)then
  allocate(dr(0:2))
  dr(0)=1*(Nbtotal+1)
  dr(1)=2*(Nbtotal+1)
  dr(2)=1*(Nbtotal+1)
  elseif(nchannel==1)then
  allocate(dr(0:2))
  dr(0)=1*(Nbtotal+1)
  dr(1)=1*(Nbtotal+1)
endif
!-----------------------------------------------------------------------
endif  !cdints
end if !m
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
! Dimension of h:
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
Do n=0,nmax
 d(n)=0
d1(n)=0
d2(n)=0
d3(n)=0
d4(n)=0
Do i=1,4
Do j=1,4
!-----------------------------------------------------------------------
if(i==j)then
!-----------------------------------------------------------------------
if(i==1)then
dbc=dr(n)
Do s=1,d0(n,nmax-nchannel,dbc)

d1(n)=d1(n)+1
end do
dof=d1(n)
!-----------------------------------------------------------------------
elseif(i==2)then
dbc=dr(n-1)
Do s=1,d0(n-1,nmax-nchannel,dbc)

d2(n)=d2(n)+1
end do
dof=d1(n)+d2(n)
!-----------------------------------------------------------------------
elseif(i==3)then
dbc=dr(n-1)
Do s=1,d0(n-1,nmax-nchannel,dbc)

d3(n)=d3(n)+1
end do
dof=d1(n)+d2(n)+d3(n)
!-----------------------------------------------------------------------
elseif(i==4)then
dbc=dr(n-2)
Do s=1,d0(n-2,nmax-nchannel,dbc)

d4(n)=d4(n)+1
end do
!-----------------------------------------------------------------------
endif
endif
!-----------------------------------------------------------------------
end do
end do
!-----------------------------------------------------------------------
! 1band case:
!-----------------------------------------------------------------------
if(nchannel==1)then
d3(n)=0
d4(n)=0
endif
!-----------------------------------------------------------------------
d(n)=d1(n)+d2(n)+d3(n)+d4(n)
end do
!-----------------------------------------------------------------------
deallocate(dr)
!-----------------------------------------------------------------------
dmax=0
dtotal=0
Do n=0,nmax
if(d(n)>dmax)then
dmax=d(n)
endif
dtotal=dtotal+d(n)
end do
!-----------------------------------------------------------------------
!print*,m,dtotal
!-----------------------------------------------------------------------
! Definition of nonvanishing Blocks for m=0
!-----------------------------------------------------------------------
allocate(h(0:nmax,1:dmax,1:dmax))
!-----------------------------------------------------------------------
if(cdint>1)then
allocate(nd(0:nmax,1:dmax,1:dmax),ddagger(0:nmax-1,1:dmax,1:dmax),nboson(0:nmax,1:dmax,1:dmax))
allocate(bpbdagger(0:nmax-1,1:dmax,1:dmax))
nd=0.0d+0
ddagger=0.0d+00
nboson=0.0d+00
bpbdagger=0.0d+00
endif
!-----------------------------------------------------------------------
h=0.0d+0
if(m==mmin)then
!-----------------------------
if(cdint>1)then
!-----------------------------------------------------------------------
  if (nchannel==2)then

    Do i=1,3
    Do j=1,3
    h(1,1,1)=(e/2.0d+0)-(u/2.0d+0)
    h(1,2,2)=(-e/2.0d+0)
    h(1,3,3)=(-e/2.0d+0)
    h(1,1,2)=Vl
    h(1,1,3)=Vr
    h(1,2,3)=0.0d+0
    h(1,j,i)=h(1,i,j)

    h(2,1,1)=e/2.0d+0
    h(2,1,2)=0.0d+0
    h(2,1,3)=-Vr
    h(2,2,2)=(e/2.0d+0)
    h(2,2,3)=Vl
    h(2,3,3)=-e/2.0d+0-(u/2.0d+0)
    h(2,j,i)=h(2,i,j)
    end do
    end do
    h(0,1,1)=(-e/2.0d+0)+(u/2.0d+0)
    h(3,1,1)=(e/2.0d+0)+(u/2.0d+0)
    
  elseif(nchannel==1)then
    
    Do nb=0,Nbtotal
    h(0,nb+1,nb+1)=H0local(0,nb)+(u/4.0)
    h(2,nb+1,nb+1)=H0local(1,nb)+(u/4.0)
    end do
    
    Do i=1,2*(Nbtotal+1)
    Do j=1,2*(Nbtotal+1)
    if(i<Nbtotal+2 .and. j<Nbtotal+2)then
    if(i==j)then
    h(1,i,i)=H0local(1,i-1)
    endif    
    elseif(i>Nbtotal+1 .and. j>Nbtotal+1)then
    if(i==j)then
    h(1,i,i)=H0local(0,i-Nbtotal-2)
    endif 
    elseif(i<Nbtotal+2 .and. j>Nbtotal+1)then
    h(1,i,j)=Vl*dbdagger(0,i,j-Nbtotal-1)
    else
    h(1,i,j)=h(1,j,i)
    endif
    
    end do
!print*,m,i,(h(1,i,j),j=1,2*(Nbtotal+1))
    end do

  endif
!-----------------------------
!-----------------------------------------------------------------------
elseif(cdint==1)then
  if(nchannel==2)then

    h(0,1,1)=0.0d+0
    h(4,1,1)=0.0d+0

    h(1,1,3)=1.0d+0
    h(1,2,4)=1.0d+0

    h(2,1,3)=-1.0d+0
    h(2,1,4)=1.0d+0
    h(2,3,6)=-1.0d+0
    h(2,4,6)=1.0d+0

    h(3,1,3)=-1.0d+0
    h(3,2,4)=-1.0d+0

  elseif(nchannel==1)then

    Do i=1,Nbtotal+1
    h(0,i,i)=(i-1)*w0
    h(2,i,i)=(i-1)*w0
    end do
    
    Do i=1,2*(Nbtotal+1)
    Do j=1,2*(Nbtotal+1)
    if(i<Nbtotal+2 .and. j<Nbtotal+2)then
    if(i==j)then
    h(1,i,i)=(i-1)*w0
    endif    
    elseif(i>Nbtotal+1 .and. j>Nbtotal+1)then
    if(i==j)then
    h(1,i,i)=(i-Nbtotal-2)*w0
    endif 
    elseif(i==j-Nbtotal-1 .or. i==j+Nbtotal+1)then
    h(1,i,j)=tml
    endif
    end do
!print*,m,i,(h(1,i,j),j=1,2*(Nbtotal+1))
    end do

  endif
!-----------------------------------------------------------------------
endif !cdint
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
else   !(for m>mmin)
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
! Definition of h(n)
!-----------------------------------------------------------------------
Do n=0,nmax
Do i=1,4
Do j=1,4
!-----------------------------------------------------------------------
! Diagonal elements:
!-----------------------------------------------------------------------
if(i==j)then
!-----------------------------------------------------------------------
if(i==1)then
Do s=1,d1(n)
h(n,s,s)=hb(n,s)!*(l**(0.5))
end do
!-----------------------------------------------------------------------
elseif(i==2)then
Do s=1,d2(n)
h(n,d1(n)+s,d1(n)+s)=hb(n-1,s)!*(l**(0.5))
end do
!-----------------------------------------------------------------------
elseif(i==3)then
Do s=1,d3(n)
h(n,d1(n)+d2(n)+s,d1(n)+d2(n)+s)=hb(n-1,s)!*(l**(0.5))
end do
!-----------------------------------------------------------------------
elseif(i==4)then
Do s=1,d4(n)
h(n,d1(n)+d2(n)+d3(n)+s,d1(n)+d2(n)+d3(n)+s)=hb(n-2,s)!*(l**(0.5))
end do
endif
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
! off-Diagonal elements:
!-----------------------------------------------------------------------
else if(i==1 .and.j==2)then
Do s=1,d1(n)
Do k=1,d2(n)
ans1=0.0d+0
Do ss=1,db1(n-1)
ans1=ans1+Trb(n-1,ss,k)*Trb(n,db1(n)+ss,s)
end do
ans3=0.0d+0
Do kk=1,db3(n-1)
ans3=ans3+Trb(n-1,kk+db1(n-1)+db2(n-1),k)*Trb(n,kk+db1(n)+db2(n)+db3(n),s)
end do
h(n,s,d1(n)+k)=(ans1-ans3)*(tml)
end do
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
end do
!-----------------------------------------------------------------------
else if(i==1 .and.j==3)then
Do s=1,d1(n)
Do k=1,d3(n)
ans1=0.0d+0
Do ss=1,db1(n-1)
ans1=ans1+Trb(n-1,ss,k)*Trb(n,ss+db1(n)+db2(n),s)
end do
ans2=0.0d+0
Do kk=1,db2(n-1)
ans2=ans2+Trb(n-1,db1(n-1)+kk,k)*Trb(n,db1(n)+db2(n)+db3(n)+kk,s)
end do
h(n,s,d1(n)+d2(n)+k)=(ans1+ans2)*(tmr)
end do
end do
!-----------------------------------------------------------------------
else if(i==2 .and.j==4)then
Do s=1,d2(n)
Do k=1,d4(n)
ans1=0.0d+0
Do ss=1,db1(n-2)
ans1=ans1+Trb(n-2,ss,k)*Trb(n-1,db1(n-1)+db2(n-1)+ss,s)
end do
ans2=0.0d+0
Do kk=1,db2(n-2)
ans2=ans2+Trb(n-2,db1(n-2)+kk,k)*Trb(n-1,db1(n-1)+db2(n-1)+db3(n-1)+kk,s)
end do
h(n,d1(n)+s,d1(n)+d2(n)+d3(n)+k)=-(ans1+ans2)*(tmr)
end do
end do
!-----------------------------------------------------------------------
elseif(i==3 .and.j==4)then
Do s=1,d3(n)
Do k=1,d4(n)
ans1=0.0d+0
Do ss=1,db1(n-2)
ans1=ans1+Trb(n-2,ss,k)*Trb(n-1,db1(n-1)+ss,s)
end do
ans3=0.0d+0
Do kk=1,db3(n-2)
ans3=ans3+Trb(n-2,db1(n-2)+db2(n-2)+kk,k)*Trb(n-1,db1(n-1)+db2(n-1)+db3(n-1)+kk,s)
end do
h(n,d1(n)+d2(n)+s,d1(n)+d2(n)+d3(n)+k)=(ans1-ans3)*(tml)
end do
end do
!-----------------------------------------------------------------------
endif
end do
end do
!-----------------------------------------------------------------------
! cause h is symmetric ...
!-----------------------------------------------------------------------
!Do s=1,d(n)
!Do k=1,d(n)
!h(n,k,s)=h(n,s,k)
!end do
!end do
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
end do !n
!-----------------------------------------------------------------------
deallocate(db1,db2,db3,db4)
!-----------------------------------------------------------------------
deallocate(hb,Trb)
!-----------------------------------------------------------------------
endif     !m>0
!-----------------------------------------------------------------------
allocate(hb(0:nmax,1:dmax),Trb(0:nmax,1:dmax,1:dmax))
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
! Diagonalization 
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
do n=0,nmax
!-----------------------------------------------------------------------
allocate(hm(1:d(n),1:d(n)),Di(1:d(n)))
!-----------------------------------------------------------------------
hm(1:d(n),1:d(n))=h(n,1:d(n),1:d(n))
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
LWORK=max(1,3*d(n)-1)
allocate(WORK(max(1,lWORK)))
!-----------------------------------------------------------------------
if(d(n)>0)then

  LDA=d(n)

  CALL DSYEV( 'Vectors', 'Upper', d(n), hm, LDA, Di, WORK, LWORK, INFO )
    
      
endif
!-----------------------------------------------------------------------
!----------------------------------------------------------------------- 
Trb(n,1:d(n),1:d(n))=(hm(1:d(n),1:d(n)))
hb(n,1:d(n))=(Di(1:d(n)))
!-----------------------------------------------------------------------
deallocate(hm,Di,WORK)
!-----------------------------------------------------------------------
end do     !n
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
allocate(dr(0:nmax))
!-----------------------------------------------------------------------
allocate(Htotal(1:dtotal))
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
!Ground-state Energy:
!-----------------------------------------------------------------------
emin=100.0!hb(0,1)   !10
nmin=0
imin=1
j=0
Do n=0,nmax
Do i=1,d(n)
if(hb(n,i)<emin)then
emin=hb(n,i)
nmin=n
imin=i
endif
j=j+1
Htotal(j)=hb(n,i)
end do
end do

if(m==14)then
print*,nmin
endif

eground(cdint,m)=emin
nground(cdint,m)=nmin
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
!1st excited Energy:
!-----------------------------------------------------------------------
ex1=hb(0,1)   
Do n=0,nmax
Do i=1,d(n)
if(n==nmin .and. i==imin)then
ex1=ex1
else
if(hb(n,i)<ex1)then
ex1=hb(n,i)
endif
endif

end do
end do
!-----------------------------------------------------------------------
! Sorting the states
!-----------------------------------------------------------------------
CALL  Sort(Htotal, dtotal)
!-----------------------------------------------------------------------
if(dtotal<tnumstates+numgap+2)then  !1702
drnmax=dmax
drt=dtotal
Do n=0,nmax
dr(n)=d(n)
end do
else
ex1=Htotal(tnumstates+1)-Htotal(tnumstates)  !1501,1500
icut=tnumstates+1  !1501
Do i=1,numgap !200
if(Htotal(tnumstates+2+i)-Htotal(tnumstates+1+i)>ex1)then  !1502,1501
ex1=Htotal(tnumstates+2+i)-Htotal(tnumstates+1+i)  !1502,1501
icut=tnumstates+2+i    !1502+i
endif
end do
Ecut=Htotal(icut)-emin!Htotal(1001)-emin
!print*,m,Htotal(1)-emin,ex1/tml
!endif
!-----------------------------------------------------------------------
!-------------------------------------------------------------------------------------------
! Keeping low-Energy states only!( Cut-off Scheme! )
!-------------------------------------------------------------------------------------------
! note: Ecut depends on nchannel and l (l=2.5 ecut=11*tml for 2c and ecut=60*tml for 1c)
!                                      (l=4.0                    and ecut=10*tml for 1c) 
!-------------------------------------------------------------------------------------------
        !Ecut=400.0d+00*tml!200.0d+00*tml   !

 	drt=0
 	Do n=0,nmax
 	dr(n)=0
 	Do j=1,d(n)
 	if((hb(n,j)-emin)<Ecut)then   !(2.0*(l**(-(m-1)/2.0))))then
 	dr(n)=dr(n)+1
 	!hb(n,dr(n))=hb(n,j)
 	
 	hb(n,dr(n))=hb(n,j)
 	
 	!CCCCCCCCCCCCCCCCCCCCCCCCCCC
 	if(hb(n,j)-emin<0)then
 	!hb(n,j)=emin
 	print*,'error'
 	endif
 	!CCCCCCCCCCCCCCCCCCCCCCCCCCC
 	
 	Do i=1,d(n)
 	Trb(n,i,dr(n))=Trb(n,i,j)
 	end do
 	endif 	
 	end do
 	drt=drt+dr(n)
        if(n==0)then
        drnmax=dr(n)
        elseif(dr(n)>drnmax)then
        drnmax=dr(n)
        endif
 	end do
        print*,ilambda,eint,iz,m,dtotal,drt,dmax,drnmax,Ecut,emin-Htotal(1)
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
endif
!-----------------------------------------------------------------------
deallocate(Htotal)
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
	!-----------------------------------------------------------------------
	! printing the first few excited energies:
	!-----------------------------------------------------------------------
	if(2*int(m/2.0)==m)then
	Do ss=1,min(4,dr(nmin))
  	!write(2,*),m,nmin,(hb(nmin,ss)-emin)/tml,(hb(nmin+1,ss)-emin)/tml&
  	!          ,(hb(nmin+2,ss)-emin)/tml&
  	!          ,(hb(nmin+3,ss)-emin)/tml&
  	!          ,(hb(nmin-1,ss)-emin)/tml,(ex1-emin)/tml,tml
	end do
	else
	Do ss=1,min(10,dr(nmin))
!  	write(2,*),m,nmin,(hb(m,ss)-emin)/(l**(-(m-1)/2.0)),(hb(m+1,ss)-emin)/(l**(-(m-1)/2.0))&
!  	          ,(hb(m+2,ss)-emin)/(l**(-(m-1)/2.0))&
!  	          ,(hb(m-1,ss)-emin)/(l**(-(m-1)/2.0)),(ex1-emin)/tml,tml
	end do
	endif
	!-----------------------------------------------------------------------
!---------------------------------------------------------------------------
   Do n=0,nmax
   !print*,'finding effective charge',m,n,dr(n),nmin
   if(nmax>15)then
   if(dr(nmin+echarge+1)>0.or.dr(nmin-echarge-1)>1)then
   !print*,"Echarge is larger"
   endif
   endif
   if(d(n)>dr(n).and.itcounter==0)then
   itcounter=itcounter+1
   m0(cdint)=m
   print*,'Starting discarding the states!',m
   endif

   end do
   !end do
!-----------------------------------------------------------------------
!dgemm parameters:
!-----------------------------------------------------------------------
alphaconst=1.0d+0
betaconst=0.0d+0
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
allocate(db1(0:nmax),db2(0:nmax),db3(0:nmax),db4(0:nmax))
Do n=0,nmax
db1(n)=d1(n)
db2(n)=d2(n)
db3(n)=d3(n)
db4(n)=d4(n)
!-----------------------------------------------------------------------
! Storing Infos:
!-----------------------------------------------------------------------
if(drnmax>500)then
print*,"NOT ENOUGHT MEMORY!!!!",d(n),dr(n),drt
endif
d1store(cdint,m,n)=d1(n)
d2store(cdint,m,n)=d2(n)
d3store(cdint,m,n)=d3(n)
d4store(cdint,m,n)=d4(n)
drstore(cdint,m,n)=dr(n)
dstore(cdint,m,n)=d(n)
Hstore(m,n,1:dr(n))=hb(n,1:dr(n))
!-----------------------------------------------------------------------
if(abs(n-nmin)<8)then
!Trstore(m,n,1:d(n),1:dr(n))=Trb(n,1:d(n),1:dr(n))
endif
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
if(cdint>1)then
!-----------------------------------------------------------------------
Do i=1,4
if(i==1)then
nii=0
ss=d1(n)
ndown=0
nup=d1(n)
else if(i==2)then
nii=1
ss=d2(n)
ndown=d1(n)
nup=d1(n)+d2(n)
elseif(i==3)then
nii=1
ss=d3(n)
ndown=d1(n)+d2(n)
nup=d1(n)+d2(n)+d3(n)
elseif(i==4)then
nii=2
ss=d4(n)
ndown=d1(n)+d2(n)+d3(n)
nup=d1(n)+d2(n)+d3(n)+d4(n)
end if


ndownmat(i,n)=ndown
nupmat(i,n)=nup
ssmat(i,n)=ss
niimat(i)=nii

if(ss>0 .and. dr(n)>0 .and. n-nii>=0 )then

!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
! Occupation Number:
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
allocate(Amatrix(1:ss,1:dr(n)),Bmatrix(1:ss,1:ss),Dmatrix(1:dr(n),1:ss))

Amatrix(1:ss,1:dr(n))=Trb(n,(ndown+1):(nup),1:dr(n))
Bmatrix(1:ss,1:ss)=ndb(n-nii,:,:)

CALL dgemm ('T','N',dr(n),ss,ss,alphaconst,Amatrix,ss,Bmatrix,ss,betaconst,Dmatrix,dr(n)) 

deallocate(Amatrix,Bmatrix)

allocate(Amatrix(1:dr(n),1:ss),Bmatrix(1:ss,1:dr(n)))

Amatrix=Dmatrix
Bmatrix=Trb(n,(ndown+1):nup,1:dr(n))

deallocate(Dmatrix)

allocate(Dmatrix(1:dr(n),1:dr(n)))

CALL dgemm ('N','N',dr(n),dr(n),ss,alphaconst,Amatrix,dr(n),Bmatrix,ss,betaconst,Dmatrix,dr(n)) 


nd(n,1:dr(n),1:dr(n))=nd(n,1:dr(n),1:dr(n))+Dmatrix


deallocate(Amatrix,Bmatrix,Dmatrix)
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
! Average number of Phonons
!-----------------------------------------------------------------------
allocate(Amatrix(1:ss,1:dr(n)),Bmatrix(1:ss,1:ss),Dmatrix(1:dr(n),1:ss))

Amatrix(1:ss,1:dr(n))=Trb(n,(ndown+1):(nup),1:dr(n))
Bmatrix(1:ss,1:ss)=nbeboson(n-nii,:,:)

CALL dgemm ('T','N',dr(n),ss,ss,alphaconst,Amatrix,ss,Bmatrix,ss,betaconst,Dmatrix,dr(n)) 

deallocate(Amatrix,Bmatrix)

allocate(Amatrix(1:dr(n),1:ss),Bmatrix(1:ss,1:dr(n)))

Amatrix=Dmatrix
Bmatrix=Trb(n,(ndown+1):nup,1:dr(n))

deallocate(Dmatrix)

allocate(Dmatrix(1:dr(n),1:dr(n)))

CALL dgemm ('N','N',dr(n),dr(n),ss,alphaconst,Amatrix,dr(n),Bmatrix,ss,betaconst,Dmatrix,dr(n)) 


nboson(n,1:dr(n),1:dr(n))=nboson(n,1:dr(n),1:dr(n))+Dmatrix


deallocate(Amatrix,Bmatrix,Dmatrix)
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
end if    !(ss>0 condition)

end do  !i
!-----------------------------------------------------------------------
endif !cdint>1
!-----------------------------------------------------------------------
end do  !n
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
! Checking:
!-----------------------------------------------------------------------
if(m>0 .and. m<mmax-1)then
zparcheck(cdint,m+1)=0.0d+00
Netotcheck(cdint,m+1)=0.0d+00
if(cdint>1)then
ndtotcheck(m+1)=0.0d+00
endif
Do n=0,nmax
Do s=1,dr(n)
zparcheck(cdint,m+1)=zparcheck(cdint,m+1)+exp(-(1.0d+00/tml)*(hb(n,s)-eground(cdint,m)))
Netotcheck(cdint,m+1)=Netotcheck(cdint,m+1)+n*exp(-(1.0d+00/tml)*(hb(n,s)-eground(cdint,m)))
!-----------------------------------------------------------------------
if(cdint>1)then
ndtotcheck(m+1)=ndtotcheck(m+1)+nd(n,s,s)*exp(-(1.0d+00/tml)*(hb(n,s)-eground(cdint,m)))
endif
!-----------------------------------------------------------------------
end do   !s
end do   !n
if(cdint>1)then
ndtotcheck(m+1)=ndtotcheck(m+1)/zparcheck(2,m+1)
Netotcheck(2,m+1)=(Netotcheck(2,m+1)/zparcheck(2,m+1))
else
Netotcheck(1,m+1)=(Netotcheck(1,m+1)/zparcheck(1,m+1))
endif
!-----------------------------------------------------------------------
!write(2,*)m,ndtotcheck(m+1)/zparcheckcdint,m+1)
!print*,m,dtotal,drt,ndtotcheck(m+1),(Netotcheck(2,m+1)-Netotcheck(1,m+1))
endif
!-----------------------------------------------------------------------
! end checking
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
! ddagger Matrix elements:
!-----------------------------------------------------------------------
if(cdint>1)then
do n=0,nmax-1
if(dr(n+1)>0 .and. d(n)>0)then
Do i=1,4
!-----------------------------------------------------------------------
if(ssmat(i,n)>0 .and. ssmat(i,n+1)>0 .and. n-niimat(i)>-1)then

allocate(Amatrix(1:ssmat(i,n+1),1:dr(n+1)),Bmatrix(1:ssmat(i,n+1),1:ssmat(i,n)),Dmatrix(1:dr(n+1),1:ssmat(i,n)))

Amatrix(1:ssmat(i,n+1),1:dr(n+1))=Trb(n+1,(ndownmat(i,n+1)+1):(nupmat(i,n+1)),1:dr(n+1))
Bmatrix(1:ssmat(i,n+1),1:ssmat(i,n))=dbdagger(n-niimat(i),1:ssmat(i,n+1),1:ssmat(i,n))

CALL dgemm ('T','N',dr(n+1),ssmat(i,n),ssmat(i,n+1),alphaconst,Amatrix,ssmat(i,n+1),Bmatrix,ssmat(i,n+1),betaconst,Dmatrix,dr(n+1))


deallocate(Amatrix,Bmatrix)
allocate(Amatrix(1:dr(n+1),1:ssmat(i,n)),Bmatrix(1:ssmat(i,n),1:dr(n)))

Amatrix=Dmatrix
Bmatrix=Trb(n,(ndownmat(i,n)+1):nupmat(i,n),1:dr(n))
deallocate(Dmatrix)
allocate(Dmatrix(1:dr(n+1),1:dr(n)))

CALL dgemm ('N','N',dr(n+1),dr(n),ssmat(i,n),alphaconst,Amatrix,dr(n+1),Bmatrix,ssmat(i,n),betaconst,Dmatrix,dr(n+1)) 

ddagger(n,1:dr(n+1),1:dr(n))=ddagger(n,1:dr(n+1),1:dr(n))+Dmatrix(1:dr(n+1),1:dr(n))
!-----------------------------------------------------------------------
deallocate(Amatrix,Bmatrix,Dmatrix)
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
! ddagger(b+bdagger) Matrix elements:
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
allocate(Amatrix(1:ssmat(i,n+1),1:dr(n+1)),Bmatrix(1:ssmat(i,n+1),1:ssmat(i,n)),Dmatrix(1:dr(n+1),1:ssmat(i,n)))

Amatrix(1:ssmat(i,n+1),1:dr(n+1))=Trb(n+1,(ndownmat(i,n+1)+1):(nupmat(i,n+1)),1:dr(n+1))
Bmatrix(1:ssmat(i,n+1),1:ssmat(i,n))=bebpbdagger(n-niimat(i),1:ssmat(i,n+1),1:ssmat(i,n))

CALL dgemm ('T','N',dr(n+1),ssmat(i,n),ssmat(i,n+1),alphaconst,Amatrix,ssmat(i,n+1),Bmatrix,ssmat(i,n+1),betaconst,Dmatrix,dr(n+1))


deallocate(Amatrix,Bmatrix)
allocate(Amatrix(1:dr(n+1),1:ssmat(i,n)),Bmatrix(1:ssmat(i,n),1:dr(n)))

Amatrix=Dmatrix
Bmatrix=Trb(n,(ndownmat(i,n)+1):nupmat(i,n),1:dr(n))
deallocate(Dmatrix)
allocate(Dmatrix(1:dr(n+1),1:dr(n)))

CALL dgemm ('N','N',dr(n+1),dr(n),ssmat(i,n),alphaconst,Amatrix,dr(n+1),Bmatrix,ssmat(i,n),betaconst,Dmatrix,dr(n+1)) 

bpbdagger(n,1:dr(n+1),1:dr(n))=bpbdagger(n,1:dr(n+1),1:dr(n))+Dmatrix(1:dr(n+1),1:dr(n))
!-----------------------------------------------------------------------
deallocate(Amatrix,Bmatrix,Dmatrix)
!-----------------------------------------------------------------------
endif
!-----------------------------------------------------------------------

end do !i

end if
end do !n
endif
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
if(cdint>1)then
!-----------------------------------------------------------------------
deallocate(ndb,dbdagger,nbeboson,bebpbdagger)
allocate(ndb(0:nmax,1:dmax,1:dmax),dbdagger(0:nmax-1,1:dmax,1:dmax),nbeboson(0:nmax,1:dmax,1:dmax))
allocate(bebpbdagger(0:nmax-1,1:dmax,1:dmax))
Do n=0,nmax
Do s=1,dr(n)
Do k=1,dr(n)
ndb(n,s,k)=nd(n,s,k)
nbeboson(n,s,k)=nboson(n,s,k)
end do
ndmat(m,n,s)=nd(n,s,s)
nbomat(m,n,s)=nboson(n,s,s)
end do

if(n<nmax)then
dbdagger(n,1:dr(n+1),1:dr(n))=ddagger(n,1:dr(n+1),1:dr(n))
ddmat(m,n,1:dr(n+1),1:dr(n))=ddagger(n,1:dr(n+1),1:dr(n))

bebpbdagger(n,1:dr(n+1),1:dr(n))=bpbdagger(n,1:dr(n+1),1:dr(n))
bpbmat(m,n,1:dr(n+1),1:dr(n))=bpbdagger(n,1:dr(n+1),1:dr(n))
endif
end do
!-----------------------------------------------------------------------
deallocate(nd,ddagger,nboson,bpbdagger)
!-----------------------------------------------------------------------
endif !cdint>1
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
!print*,m,dtotal,drt,discounter,Ecut/tml,dmax,drnmax
!-----------------------------------------------------------------------
deallocate(d1,d2,d3,d4)
!-----------------------------------------------------------------------
deallocate(h,d)
!-----------------------------------------------------------------------
end do !m
!-----------------------------------------------------------------------
if(cdint>1)then
deallocate(ndb,dbdagger,nbeboson,bebpbdagger)
endif
!-----------------------------------------------------------------------
!IIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIII
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
deallocate(dr,hb,trb,db1,db2,db3,db4)
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
! Thermodynamical Quantities:
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
if(cdint==1)then
m0(1)=min(m0(1),m0(2))
endif
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
!TTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTT
! Temerature Grid:
!-----------------------------------------------------------------------
itcounter=0
Do it=1,itmaxdim!1!itmaxdim!itmaxdim-15

!Tem(it)=1.0d+00*gama*(10.0d+00**(-4.0d+00))!1.0d+00*(l**(-19.d+00*(it-1)/((itmaxdim-1)*1.0d+00)))!l**(-(mmax-3)/2.0d+00)   !10.0d+00*(l**(-19.d+00*(it-1)/((itmaxdim-1)*1.0d+00)))
!Tem(it)=1.0d+00*(l**(-19.d+00*(it-1)/((itmaxdim-1)*1.0d+00)))!l**(-(mmax-3)/2.0d+00)
Tem(it)=(gefapprox)*(10.0d+00**(6.0d+00+(-10.0d+00*(it-1)/itmaxdim)))
!Tem(it)=gama*0.01d+00*it
if(Tem(it)<dble(tm(mmax-2)))then
print*,it,"Tem too small!!!"
else
itmaximum=it
end if
end do
!-----------------------------------------------------------------------
! Best Shell:
!-----------------------------------------------------------------------
Do it=1,itmaximum

if(alphatem*Tem(it).ge. dble(tm(m0(cdint)-2)))then
Tshell(it)=m0(cdint)-1
else
Do m=m0(cdint),mmax-1
if(alphatem*Tem(it).lt.dble(tm(m-2)).and.alphatem*Tem(it).ge.dble(tm(m-1)))then
Tshell(it)=m
endif
end do  !m
endif

if(rank==0)then
!write(1,*)m0(cdint),it,Tem(it),Tshell(it),alphatem*Tem(it)/dble(tm(Tshell(it)-1))
endif
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
beta(it)=1.0d+00/Tem(it)
!-----------------------------------------------------------------------
! Initialization
!-----------------------------------------------------------------------
zpar(cdint,it)=0.0d+00

Hav(cdint,it)=0.0d+00
HHav(cdint,it)=0.0d+00
Netot(cdint,it)=0.0d+00
if(cdint>1)then
ndtot(it)=0.0d+00
nTboson(it)=0.0d+00
endif
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
Do n=0,(nchannel*(Tshell(it)+1)+cdint-1)
Do s=1,drstore(cdint,Tshell(it),n)
zpar(cdint,it)=zpar(cdint,it)+exp(-beta(it)*(Hstore(Tshell(it),n,s)-eground(cdint,Tshell(it))))

Hav(cdint,it)=Hav(cdint,it)+((Hstore(Tshell(it),n,s)-eground(cdint,Tshell(it))))&
                            *exp(-beta(it)*(Hstore(Tshell(it),n,s)-eground(cdint,Tshell(it))))

HHav(cdint,it)=HHav(cdint,it)+(((Hstore(Tshell(it),n,s)-eground(cdint,Tshell(it))))**2.0)&
                              *exp(-beta(it)*(Hstore(Tshell(it),n,s)-eground(cdint,Tshell(it))))

Netot(cdint,it)=Netot(cdint,it)+n*exp(-beta(it)*(Hstore(Tshell(it),n,s)-eground(cdint,Tshell(it))))

!-----------------------------------------------------------------------
if(cdint>1)then
ndtot(it)=ndtot(it)+ndmat(Tshell(it),n,s)*exp(-beta(it)*(Hstore(Tshell(it),n,s)-eground(cdint,Tshell(it))))
nTboson(it)=nTboson(it)+nbomat(Tshell(it),n,s)*exp(-beta(it)*(Hstore(Tshell(it),n,s)-eground(cdint,Tshell(it))))
endif
!-----------------------------------------------------------------------
end do   !s
end do   !n
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
if(cdint>1)then
ndtot(it)=ndtot(it)/zpar(cdint,it)
nTboson(it)=nTboson(it)/zpar(cdint,it)
endif

Netot(cdint,it)=Netot(cdint,it)/zpar(cdint,it)
Hav(cdint,it)=Hav(cdint,it)/zpar(cdint,it)
HHav(cdint,it)=HHav(cdint,it)/zpar(cdint,it)

BEntropy(cdint,it)=beta(it)*Hav(cdint,it)+log(zpar(cdint,it))

!-----------------------------------------------------------------------
! specific heat:
!-----------------------------------------------------------------------
scimp(cdint,it)=(beta(it)**2.0d+00)*((HHav(cdint,it))-((Hav(cdint,it))**2.0d+00))
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
! Effective Tunneling:
!-----------------------------------------------------------------------
if(cdint>1)then
!ndde0(eint,it)=ndtot(it)
!if(eint==52)then
!efgammadd(it)=-1.0d+00*(e0(52)-e0(51))/(pi*(ndde0(52,it)-ndde0(51,it)))  
!call MPI_Reduce(efgammadd(it),collector2,1,MPI_double_precision,MPI_SUM,0, MPI_COMM_WORLD, ierr)
!efgammadd(it)=collector2/(1.0d+00*izmax)
!if(rank==0)then
!write(2,*)it,Tem(it)/gama,efgammadd(it)/gama,e0(52)/gama
!endif
!endif
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------        
! Moments of spectral functions without broadening!
!-----------------------------------------------------------------------
I0m(iz,eint,it)=0.0d+00
I1m(iz,eint,it)=0.0d+00
I2m(iz,eint,it)=0.0d+00
!-----------------------------------------------------------------------
!if(cdint>1)then
Do n=0,nchannel*(Tshell(it)+1)
Do s=1,drstore(2,Tshell(it),n+1)
Do k=1,drstore(2,Tshell(it),n)

I0m(iz,eint,it)=I0m(iz,eint,it)+((ddmat(Tshell(it),n,s,k)*ddmat(Tshell(it),n,s,k)))&
                                 /((exp(beta(it)*(Hstore(Tshell(it),n+1,s)-eground(2,Tshell(it)))))&
                                  +(exp(beta(it)*(Hstore(Tshell(it),n,k)-eground(2,Tshell(it))))))


I1m(iz,eint,it)=I1m(iz,eint,it)+((ddmat(Tshell(it),n,s,k)*ddmat(Tshell(it),n,s,k)))&
                                 *(Hstore(Tshell(it),n+1,s)-Hstore(Tshell(it),n,k))&
                                 /((exp(beta(it)*(Hstore(Tshell(it),n+1,s)-eground(2,Tshell(it)))))&
                                  +(exp(beta(it)*(Hstore(Tshell(it),n,k)-eground(2,Tshell(it))))))

I2m(iz,eint,it)=I2m(iz,eint,it)+((ddmat(Tshell(it),n,s,k)*ddmat(Tshell(it),n,s,k)))&
                                 *((Hstore(Tshell(it),n+1,s)-Hstore(Tshell(it),n,k))**2)&
                                 /((exp(beta(it)*(Hstore(Tshell(it),n+1,s)-eground(2,Tshell(it)))))&
                                  +(exp(beta(it)*(Hstore(Tshell(it),n,k)-eground(2,Tshell(it))))))

end do
end do
end do

I0m(iz,eint,it)=I0m(iz,eint,it)/(Tem(it)*zpar(cdint,it))
I1m(iz,eint,it)=I1m(iz,eint,it)/(Tem(it)*zpar(cdint,it))
I2m(iz,eint,it)=I2m(iz,eint,it)/(Tem(it)*zpar(cdint,it))

endif
!-----------------------------------------------------------------------        
!-----------------------------------------------------------------------
if(cdint>1)then
ndde0(iz,eint,it)=ndtot(it)
!-----------------------------------------------------------------------
!call MPI_Reduce(nTboson(it),nbocollector,1,MPI_double_precision,MPI_SUM,0, MPI_COMM_WORLD, ierr)
!-----------------------------------------------------------------------
if(rank==0) then
!write(2,*)it,Tem(it)/gama,nbocollector/(1.0d+00*izmax)
end if
!-----------------------------------------------------------------------
elseif(cdint==1)then
!-----------------------------------------------------------------------
nde0(eint,it)=Netot(2,it)-Netot(1,it)
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
! z-averaging
!-----------------------------------------------------------------------
zavscimp=scimp(2,it)-scimp(1,it)
!call MPI_Reduce(zavscimp,collector,1,MPI_double_precision,MPI_SUM,0, MPI_COMM_WORLD, ierr)
!-----------------------------------------------------------------------
!call MPI_Reduce(nde0(eint,it),ndcollector,1,MPI_double_precision,MPI_SUM,0, MPI_COMM_WORLD, ierr)
!-----------------------------------------------------------------------
!call MPI_Reduce(ndde0(eint,it),nddcollector,1,MPI_double_precision,MPI_SUM,0, MPI_COMM_WORLD, ierr)
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
if(rank==0) then
!write(2,*)it,Tem(it)/gama,collector/(1.0d+00*izmax),ndcollector/(1.0d+00*izmax)&
!          ,Tshell(it),Tem(it)/dble(tm(Tshell(it)-1)),nddcollector/(1.0d+00*izmax)&
!          ,nbocollector/(1.0d+00*izmax),1.0d+00/((exp(beta(it)*w0))-1.0d+00)
end if
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
endif
!-----------------------------------------------------------------------
if(cdint>1)then
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
! checking:
!-----------------------------------------------------------------------
Do m=0,mmax-1
spectral(1,1)=0.0d+00
zpartition(m)=0.0d+00
degeneracy=0 


Do n=0,nchannel*(m+1)+1

Do s=1,drstore(2,m,n)
zpartition(m)=zpartition(m)+exp(-beta(it)*(Hstore(m,n,s)-eground(2,m)))
end do
end do

!print*,m,spectral(1,1)/zpartition(m)

!spectral(1,1)=spectral(1,1)
!write(1,*)it,Tem(it)/gama,m,zpartition(m),spectral(1,1)/zpartition(m),Tshell(it)
end do
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
! Spectral Function:
!-----------------------------------------------------------------------
bfactor=0.3d+00
alphaw1=1.0d+00
area1=0.0d+00
area2=0.0d+00
sumFcorr=0.0d+00
asumFcorr=0.0d+00
!-----------------------------------------------------------------------
! Incoming-Energy Grid:
!-----------------------------------------------------------------------
Do iw1=1,1!iw1maxdim/2
!-----------------------------------------------------------------------
w1(iw1)=1.0d+00*(l**(-24.d+00*(iw1-1)/(((iw1maxdim/2)-1)*1.0d+00)))
if(alphaw1*w1(iw1)<l**(-(mmax-3)/2.0d+00))then  !dble(tm(mmax-2)))then
!print*,iw1,"Incoming-energy grid too small!!!"
else
iw1max=iw1
end if
end do !iw1
!-----------------------------------------------------------------------
Do iw1=1,1!2*iw1max
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
Spectral(it,iw1)=0.0d+00
sumrule(it,iw1)=0.0d+00
Fcorr(it,iw1)=0.0d+00
phselfenergy=0.0d+00
DotgR(it,iw1)=0.0d+00
newSpectral(it,iw1)=0.0d+00

aSpectral(it,iw1)=0.0d+00
aFcorr(it,iw1)=0.0d+00
aphselfenergy=0.0d+00
aDotgR(it,iw1)=0.0d+00
anewSpectral(it,iw1)=0.0d+00

MatFcorr(it,iw1)=0.0d+00
MatDotg(it,iw1)=0.0d+00
Matphselfenergy(it,iw1)=0.0d+00
!-----------------------------------------------------------------------
if(iw1<iw1max+1)then
!-----------------------------------------------------------------------
! Best Shell:
!-----------------------------------------------------------------------
if(alphaw1*w1(iw1).ge. dble(tm(m0(cdint)-2)))then
W1shell(iw1)=m0(cdint)-1
else
Do m=m0(cdint),mmax-1
if(alphaw1*w1(iw1).lt.dble(tm(m-2)).and.alphaw1*w1(iw1).ge.dble(tm(m-1)))then
w1shell(iw1)=m
endif
end do  !m
endif
!-----------------------------------------------------------------------
else
w1(iw1)=-w1(2*iw1max-iw1+1)
w1shell(iw1)=w1shell(2*iw1max-iw1+1)
endif
!-----------------------------------------------------------------------
dfactor=1.0d+00*(dble(tm(w1shell(iw1)-1)))  !1.9
!-----------------------------------------------------------------------
Do n=0,nchannel*(w1shell(iw1)+1)
Do s=1,drstore(2,w1shell(iw1),n+1)
do k=1,drstore(2,w1shell(iw1),n)
!------------------------------------------------------------------------------------------------------------------------------
if(SIGN(1.0d+00,w1(iw1))*(Hstore(w1shell(iw1),n+1,s)-Hstore(w1shell(iw1),n,k)).LE.0.0d+00)then
adeltaw1(iw1)=0.0d+00
else
adeltaw1(iw1)=(exp(-bfactor**2/4.0d+00)/(bfactor*sqrt(pi)*abs(Hstore(w1shell(iw1),n+1,s)-Hstore(w1shell(iw1),n,k))))&
              *exp(-((log(abs((w1(iw1))/(Hstore(w1shell(iw1),n+1,s)-Hstore(w1shell(iw1),n,k)))))/bfactor)**2)
endif
!------------------------------------------------------------------------------------------------------------------------------
deltaw1(iw1)=(exp(-(((w1(iw1)-(Hstore(w1shell(iw1),n+1,s)-Hstore(w1shell(iw1),n,k)))/(dfactor))**2)))&
              /((sqrt(pi))*dfactor)

sumrule(it,iw1)=sumrule(it,iw1)+((ddmat(w1shell(iw1),n,s,k))**2)&
                                 *((exp(-beta(it)*(Hstore(w1shell(iw1),n+1,s)-eground(2,w1shell(iw1)))))&
                                  +(exp(-beta(it)*(Hstore(w1shell(iw1),n,k)-eground(2,w1shell(iw1))))))


Spectral(it,iw1)=Spectral(it,iw1)+((ddmat(w1shell(iw1),n,s,k)*ddmat(w1shell(iw1),n,s,k)))&
                                 *((exp(-beta(it)*(Hstore(w1shell(iw1),n+1,s)-eground(2,w1shell(iw1)))))&
                                  +(exp(-beta(it)*(Hstore(w1shell(iw1),n,k)-eground(2,w1shell(iw1))))))*deltaw1(iw1)


aSpectral(it,iw1)=aSpectral(it,iw1)+((ddmat(w1shell(iw1),n,s,k)*ddmat(w1shell(iw1),n,s,k)))&
                                 *((exp(-beta(it)*(Hstore(w1shell(iw1),n+1,s)-eground(2,w1shell(iw1)))))&
                                  +(exp(-beta(it)*(Hstore(w1shell(iw1),n,k)-eground(2,w1shell(iw1))))))*adeltaw1(iw1)

Fcorr(it,iw1)=Fcorr(it,iw1)+ddmat(w1shell(iw1),n,s,k)*bpbmat(w1shell(iw1),n,s,k)&
                            *(exp(-beta(it)*(Hstore(w1shell(iw1),n,k)-eground(2,w1shell(iw1))))&
                             +exp(-beta(it)*(Hstore(w1shell(iw1),n+1,s)-eground(2,w1shell(iw1)))))*deltaw1(iw1)

aFcorr(it,iw1)=aFcorr(it,iw1)+ddmat(w1shell(iw1),n,s,k)*bpbmat(w1shell(iw1),n,s,k)&
                            *(exp(-beta(it)*(Hstore(w1shell(iw1),n,k)-eground(2,w1shell(iw1))))&
                             +exp(-beta(it)*(Hstore(w1shell(iw1),n+1,s)-eground(2,w1shell(iw1)))))*adeltaw1(iw1)
!-------------------------------------------------------------------------------------------------------------------------
! Dot Green's function and Self-energy in imaginary space:
!-------------------------------------------------------------------------------------------------------------------------
MatFcorr(it,iw1)=MatFcorr(it,iw1)+ddmat(w1shell(iw1),n,s,k)*bpbmat(w1shell(iw1),n,s,k)&
                            *(exp(-beta(it)*(Hstore(w1shell(iw1),n,k)-eground(2,w1shell(iw1))))&
                             +exp(-beta(it)*(Hstore(w1shell(iw1),n+1,s)-eground(2,w1shell(iw1)))))&!*deltaw1(iw1)!&
                             /(ic*w1(iw1)+Hstore(w1shell(iw1),n,k)-Hstore(w1shell(iw1),n+1,s))


MatDotg(it,iw1)=MatDotg(it,iw1)+((ddmat(w1shell(iw1),n,s,k)*ddmat(w1shell(iw1),n,s,k)))&
                                 *((exp(-beta(it)*(Hstore(w1shell(iw1),n+1,s)-eground(2,w1shell(iw1)))))&
                                  +(exp(-beta(it)*(Hstore(w1shell(iw1),n,k)-eground(2,w1shell(iw1))))))&
                             /(ic*w1(iw1)+Hstore(w1shell(iw1),n,k)-Hstore(w1shell(iw1),n+1,s))

!-------------------------------------------------------------------------------------------------------------------------
!-------------------------------------------------------------------------------------------------------------------------
end do !k
end do !s
end do !n

Spectral(it,iw1)=Spectral(it,iw1)/zpartition(w1shell(iw1))   !zpar(2,it)
sumrule(it,iw1)=sumrule(it,iw1)/zpartition(w1shell(iw1))
Fcorr(it,iw1)=Fcorr(it,iw1)/zpartition(w1shell(iw1))
DotgR(it,iw1)=DotgR(it,iw1)/zpartition(w1shell(iw1))

aSpectral(it,iw1)=aSpectral(it,iw1)/zpartition(w1shell(iw1))
aFcorr(it,iw1)=aFcorr(it,iw1)/zpartition(w1shell(iw1))
aDotgR(it,iw1)=aDotgR(it,iw1)/zpartition(w1shell(iw1))

MatFcorr(it,iw1)=MatFcorr(it,iw1)/zpartition(w1shell(iw1))
MatDotg(it,iw1)=MatDotg(it,iw1)/zpartition(w1shell(iw1))

!-------------------------------------------------------------------------------------------------------------------------
!-----------------------------------------------------------------------
! checking the normalization of spectral function:
!-----------------------------------------------------------------------
if(iw1>1)then
area1=area1-(w1(iw1)-w1(iw1-1))*(Spectral(it,iw1)+Spectral(it,iw1-1))*0.50d+00
area2=area2-(w1(iw1)-w1(iw1-1))*(aSpectral(it,iw1)+aSpectral(it,iw1-1))*0.50d+00
sumFcorr=sumFcorr-(w1(iw1)-w1(iw1-1))*(Fcorr(it,iw1)+aSpectral(it,iw1-1))*0.50d+00
asumFcorr=asumFcorr-(w1(iw1)-w1(iw1-1))*(aFcorr(it,iw1)+aSpectral(it,iw1-1))*0.50d+00
endif
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
! z-averaging
!-----------------------------------------------------------------------
!zSpectral(it,iw1)=zSpectral(it,iw1)+(Spectral(it,iw1)/(1.0d+00*izmax))
!call MPI_Reduce(Spectral(it,iw1),zSpectral(it,iw1),1,MPI_double_precision,MPI_SUM,0, MPI_COMM_WORLD, ierr)
zSpectral(it,iw1)=zSpectral(it,iw1)/(1.0d+00*izmax)

!zFcorr(it,iw1)=zFcorr(it,iw1)+(Fcorr(it,iw1)/(1.0d+00*izmax))
!call MPI_Reduce(Fcorr(it,iw1),zFcorr(it,iw1),1,MPI_double_precision,MPI_SUM,0, MPI_COMM_WORLD, ierr)
zFcorr(it,iw1)=zFcorr(it,iw1)/(1.0d+00*izmax)

!zaSpectral(it,iw1)=zaSpectral(it,iw1)+(aSpectral(it,iw1)/(1.0d+00*izmax))
!call MPI_Reduce(aSpectral(it,iw1),zaSpectral(it,iw1),1,MPI_double_precision,MPI_SUM,0, MPI_COMM_WORLD, ierr)
zaSpectral(it,iw1)=zaSpectral(it,iw1)/(1.0d+00*izmax)

!zaFcorr(it,iw1)=zaFcorr(it,iw1)+(aFcorr(it,iw1)/(1.0d+00*izmax))
!call MPI_Reduce(aFcorr(it,iw1),zaFcorr(it,iw1),1,MPI_double_precision,MPI_SUM,0, MPI_COMM_WORLD, ierr)
zaFcorr(it,iw1)=zaFcorr(it,iw1)/(1.0d+00*izmax)

!zMatDotgRe(it,iw1)=zMatDotgRe(it,iw1)+(real(MatDotg(it,iw1))/(1.0d+00*izmax))
!zMatDotgIm(it,iw1)=zMatDotgIm(it,iw1)+(aimag(MatDotg(it,iw1))/(1.0d+00*izmax))
!call MPI_Reduce(real(MatDotg(it,iw1)),zMatDotgRe(it,iw1),1,MPI_double_precision,MPI_SUM,0, MPI_COMM_WORLD, ierr)
zMatDotgRe(it,iw1)=zMatDotgRe(it,iw1)/(1.0d+00*izmax)
!call MPI_Reduce(aimag(MatDotg(it,iw1)),zMatDotgIm(it,iw1),1,MPI_double_precision,MPI_SUM,0, MPI_COMM_WORLD, ierr)
zMatDotgIm(it,iw1)=zMatDotgIm(it,iw1)/(1.0d+00*izmax)

!zMatFcorrRe(it,iw1)=zMatFcorrRe(it,iw1)+(real(MatFcorr(it,iw1))/(1.0d+00*izmax))
!zMatFcorrIm(it,iw1)=zMatFcorrIm(it,iw1)+(aimag(MatFcorr(it,iw1))/(1.0d+00*izmax))
!call MPI_Reduce(real(MatFcorr(it,iw1)),zMatFcorrRe(it,iw1),1,MPI_double_precision,MPI_SUM,0, MPI_COMM_WORLD, ierr)
zMatFcorrRe(it,iw1)=zMatFcorrRe(it,iw1)/(1.0d+00*izmax)
!call MPI_Reduce(aimag(MatFcorr(it,iw1)),zMatFcorrIm(it,iw1),1,MPI_double_precision,MPI_SUM,0, MPI_COMM_WORLD, ierr)
zMatFcorrIm(it,iw1)=zMatFcorrIm(it,iw1)/(1.0d+00*izmax)
!-----------------------------------------------------------------------
! checking:
!write(2,*)iz,iw1,w1(iw1)/gama,w1shell(iw1),Spectral(it,iw1),aSpectral(it,iw1),zSpectral(it,iw1)
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
end do !iw1
!-----------------------------------------------------------------------
endif  !cdint>1
!-----------------------------------------------------------------------
end do   !it
!-----------------------------------------------------------------------
!TTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTT
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
end do !cdint
!-----------------------------------------------------------------------
end do !iz
!-----------------------------------------------------------------------
Do it=itmaximum-1,itmaximum-1!1,itmaximum!
!-----------------------------------------------------------------------
!I0m(it)=0.0d+00
!I1m(it)=0.0d+00
!I2m(it)=0.0d+00
!-----------------------------------------------------------------------
Do iw1=1,1!2*iw1max
!-----------------------------------------------------------------------
Do iw1prime=2,2*iw1max-1
if(iw1prime .ne. iw1 )then
if(w1(iw1prime)>0)then
DotgR(it,iw1)=DotgR(it,iw1)+abs(w1(iw1prime)-w1(iw1prime-1))*(zSpectral(it,iw1prime)-zSpectral(it,iw1))/(w1(iw1prime)-w1(iw1))
phselfenergy(it,iw1)=phselfenergy(it,iw1)+abs(w1(iw1prime)-w1(iw1prime-1))&
                                         *(zFcorr(it,iw1prime)-zFcorr(it,iw1))/(w1(iw1prime)-w1(iw1))

aDotgR(it,iw1)=aDotgR(it,iw1)+abs(w1(iw1prime)-w1(iw1prime-1))*(zaSpectral(it,iw1prime)-zaSpectral(it,iw1))/(w1(iw1prime)-w1(iw1))
aphselfenergy(it,iw1)=aphselfenergy(it,iw1)+abs(w1(iw1prime)-w1(iw1prime-1))&
                      *(zaFcorr(it,iw1prime)-zaFcorr(it,iw1))/(w1(iw1prime)-w1(iw1))
else
DotgR(it,iw1)=DotgR(it,iw1)+abs(w1(iw1prime+1)-w1(iw1prime))*(zSpectral(it,iw1prime)-zSpectral(it,iw1))/(w1(iw1prime)-w1(iw1))
phselfenergy(it,iw1)=phselfenergy(it,iw1)+abs(w1(iw1prime+1)-w1(iw1prime))&
                                         *(zFcorr(it,iw1prime)-zFcorr(it,iw1))/(w1(iw1prime)-w1(iw1))

aDotgR(it,iw1)=aDotgR(it,iw1)+abs(w1(iw1prime+1)-w1(iw1prime))*(zaSpectral(it,iw1prime)-zaSpectral(it,iw1))/(w1(iw1prime)-w1(iw1))
aphselfenergy(it,iw1)=aphselfenergy(it,iw1)+abs(w1(iw1prime+1)-w1(iw1prime))&
                      *(zaFcorr(it,iw1prime)-zaFcorr(it,iw1))/(w1(iw1prime)-w1(iw1))

endif
endif
end do
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------

DotgR(it,iw1)=-DotgR(it,iw1)/(pi*2.0d+00)-(zSpectral(it,iw1)/(2.0d+00))*((1.0d+00/pi)*log(abs((width-w1(iw1))/(width+w1(iw1))))+ic)
phselfenergy(it,iw1)=-(phselfenergy(it,iw1)/(pi*2.0d+00))&
                     -(zFcorr(it,iw1)/2.0d+00)*((1.0d+00/pi)*log(abs((width-w1(iw1))/(width+w1(iw1))))+ic)
phselfenergy(it,iw1)=lambda*phselfenergy(it,iw1)/DotgR(it,iw1)
newspectral(it,iw1)=-aimag(1.0d+00/(w1(iw1)-e0(eint)+ic*gama-phselfenergy(it,iw1)))


aDotgR(it,iw1)=-aDotgR(it,iw1)/(pi*2.0d+00)-(zaSpectral(it,iw1)/(2.0d+00))*((1.0d+00/pi)&
               *log(abs((width-w1(iw1))/(width+w1(iw1))))+ic)
aphselfenergy(it,iw1)=-(aphselfenergy(it,iw1)/(pi*2.0d+00))&
                     -(zaFcorr(it,iw1)/2.0d+00)*((1.0d+00/pi)*log(abs((width-w1(iw1))/(width+w1(iw1))))+ic)
aphselfenergy(it,iw1)=lambda*aphselfenergy(it,iw1)/aDotgR(it,iw1)
anewspectral(it,iw1)=-aimag(1.0d+00/(w1(iw1)-e0(eint)+ic*gama-aphselfenergy(it,iw1)))
!-----------------------------------------------------------------------
Matphselfenergy(it,iw1)=lambda*(zMatFcorrRe(it,iw1)+ic*zMatFcorrIm(it,iw1))/(zMatDotgRe(it,iw1)+ic*zMatDotgIm(it,iw1))
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
! Moments of Spectral Function:
!-----------------------------------------------------------------------
if(iw1>2 .and. iw1<2*iw1max-1 )then
if(w1(iw1prime)>0)then
!I0m(it)=I0m(it)+abs(w1(iw1)-w1(iw1-1))*(anewspectral(it,iw1))&
!               *((-exp(w1(it)/Tem(it)))/(Tem(it)*(1.0d+00+exp(w1(iw1)/Tem(it)))**2))

!I0m(it)=I0m(it)+((1.0d+00-TANH(w1(iw1)/(2.0d+00*Tem(it))))/2.0d+00)*(anewspectral(it,iw1-1)-anewspectral(it,iw1))

!I1m(it)=I1m(it)+abs(w1(iw1)-w1(iw1-1))*(anewspectral(it,iw1))*(w1(iw1))&
!               *((-exp(w1(it)/Tem(it)))/(Tem(it)*(1.0d+00+exp(w1(iw1)/Tem(it)))**2))

!I1m(it)=I1m(it)+((1.0d+00-TANH(w1(iw1)/(2.0d+00*Tem(it))))/2.0d+00)*(w1(iw1-1)*anewspectral(it,iw1-1)-w1(iw1)*anewspectral(it,iw1))

!I2m(it)=I2m(it)+abs(w1(iw1)-w1(iw1-1))*(anewspectral(it,iw1))*(w1(iw1)**2)&
!               *((-exp(w1(it)/Tem(it)))/(Tem(it)*(1.0d+00+exp(w1(iw1)/Tem(it)))**2))

!I2m(it)=I2m(it)+((1.0d+00-TANH(w1(iw1)/(2.0d+00*Tem(it))))/2.0d+00)&
!               *((w1(iw1-1)*2)*anewspectral(it,iw1-1)-(w1(iw1)**2)*anewspectral(it,iw1))

else
!I0m(it)=I0m(it)+abs(w1(iw1+1)-w1(iw1))*(anewspectral(it,iw1))&
!               *((-exp(w1(it)/Tem(it)))/(Tem(it)*(1.0d+00+exp(w1(iw1)/Tem(it)))**2))

!I0m(it)=I0m(it)+((1.0d+00-TANH(w1(iw1)/(2.0d+00*Tem(it))))/2.0d+00)*(anewspectral(it,iw1)-anewspectral(it,iw1+1))

!I1m(it)=I1m(it)+abs(w1(iw1+1)-w1(iw1))*(anewspectral(it,iw1))*(w1(iw1))&
!               *((-exp(w1(it)/Tem(it)))/(Tem(it)*(1.0d+00+exp(w1(iw1)/Tem(it)))**2))

!I1m(it)=I1m(it)+((1.0d+00-TANH(w1(iw1)/(2.0d+00*Tem(it))))/2.0d+00)*(w1(iw1)*anewspectral(it,iw1)-w1(iw1+1)*anewspectral(it,iw1+1))

!I2m(it)=I2m(it)+abs(w1(iw1+1)-w1(iw1))*(anewspectral(it,iw1))*(w1(iw1)**2)&
!               *((-exp(w1(it)/Tem(it)))/(Tem(it)*(1.0d+00+exp(w1(iw1)/Tem(it)))**2))

!I2m(it)=I2m(it)+((1.0d+00-TANH(w1(iw1)/(2.0d+00*Tem(it))))/2.0d+00)&
!               *((w1(iw1)*2)*anewspectral(it,iw1)-(w1(iw1+1)**2)*anewspectral(it,iw1+1))

endif
endif
!-----------------------------------------------------------------------
!write(2,*)it,Tem(it)/gama,iw1,w1(iw1)/gama,w1shell(iw1),gama*anewspectral(it,iw1)!,I0m(it)/anewspectral(it,iw1)&
!          ,(anewspectral(it,iw1-1)-anewspectral(it,iw1)),((1.0d+00-TANH(w1(iw1)/(2*Tem(it))))/2.0d+00)
!-----------------------------------------------------------------------
!write(2,*)it,Tem(it)/gama,iw1,w1(iw1)/gama,w1shell(iw1),pi*gama*zSpectral(it,iw1)&
!         ,pi*gama*zaSpectral(it,iw1),gama*newSpectral(it,iw1),gama*anewSpectral(it,iw1)&
!         ,real(phselfenergy(it,iw1))/((lambda**2/w0)),aimag(phselfenergy(it,iw1))/(gama)&
!         ,real(aphselfenergy(it,iw1))/((lambda**2/w0)),aimag(aphselfenergy(it,iw1))/(gama)&
!         ,(1.0d+00/(lambda**2/w0))*real(Matphselfenergy(it,iw1))&
!         ,(1.0d+00/gama)*aimag(Matphselfenergy(it,iw1))&
!         ,w1(iw1)/dble(tm(w1shell(iw1)-1)),area1,area2
!-----------------------------------------------------------------------
if(iw1==iw1max-1)then
!write(2,*)it,Tem(it)/gama,iw1,w1(iw1)/gama,w1shell(iw1),pi*gama*zSpectral(it,iw1)&
!         ,pi*gama*zaSpectral(it,iw1),gama*newSpectral(it,iw1),gama*anewSpectral(it,iw1),eint,e0(eint)/gama
endif         
!-----------------------------------------------------------------------
end do !iw1
!-----------------------------------------------------------------------
!write(2,*)it,Tem(it)/gama,I0m(it)/anewSpectral(it,iw1max-1),I1m(it),I2m(it),eint,e0(eint)/gama&
!            ,((1.0d+00-TANH(w1(iw1max-1)/(2.0d+00*Tem(it))))/2.0d+00) 
!-----------------------------------------------------------------------
end do !it
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
end do !eint
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
! Effective Tunneling:
!-----------------------------------------------------------------------
!efgammadd(itmaximum-5)=-1.0d+00*(e0(52)-e0(51))/(pi*(ndde0(52,itmaximum-5)-ndde0(51,itmaximum-5)))  
!call MPI_Reduce(efgammadd(it),collector3,1,MPI_double_precision,MPI_SUM,0, MPI_COMM_WORLD, ierr)
!efgammadd(itmaximum-5)=collector3/(1.0d+00*izmax)
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
Do eint=53,53
!-----------------------------------------------------------------------
Do it=1,itmaximum

!-----------------------------------------------------------------------
! Effective Tunneling:
!-----------------------------------------------------------------------
gcollector2=0.0d+00
I0collector=0.0d+00
ndcollector=0.0d+00
nddcollector=0.0d+00
!-----------------------------------------------------------------------
Do iz=0,izmax-1
!-----------------------------------------------------------------------
efgammadd(it)=-1.0d+00*(e0(52)-e0(51))/(pi*(ndde0(iz,52,it)-ndde0(iz,51,it)))  
!call MPI_Reduce(efgammadd(it),collector2,1,MPI_double_precision,MPI_SUM,0, MPI_COMM_WORLD, ierr)
!efgammadd(it)=collector2/(1.0d+00*izmax)
gcollector2=gcollector2+efgammadd(it)
!-----------------------------------------------------------------------
! z-averaging
!-----------------------------------------------------------------------
!call MPI_Reduce(I0m(eint,it),collector,1,MPI_double_precision,MPI_SUM,0, MPI_COMM_WORLD, ierr)
I0collector=I0collector+I0m(iz,eint,it)
!-----------------------------------------------------------------------
!call MPI_Reduce(I1m(eint,it),ndcollector,1,MPI_double_precision,MPI_SUM,0, MPI_COMM_WORLD, ierr)
ndcollector=ndcollector+I1m(iz,eint,it)
!-----------------------------------------------------------------------
!call MPI_Reduce(I2m(eint,it),nddcollector,1,MPI_double_precision,MPI_SUM,0, MPI_COMM_WORLD, ierr)
nddcollector=nddcollector+I2m(iz,eint,it)
!-----------------------------------------------------------------------
end do  !iz
!-----------------------------------------------------------------------
efgammadd(it)=gcollector2/(1.0d+00*izmax)
I0mz(eint,it)=pi*gama*I0collector/(1.0d+00*izmax)
I1mz(eint,it)=pi*gama*ndcollector/(1.0d+00*izmax)
I2mz(eint,it)=pi*gama*nddcollector/(1.0d+00*izmax)
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
end do  !it
!-----------------------------------------------------------------------
Do it=1,itmaximum
!-----------------------------------------------------------------------
ThEke(eint,it)=(I2mz(eint,it)-(((I1mz(eint,it))**2)/(I0mz(eint,it))))/Tem(it)
ThES(eint,it)=-I1mz(eint,it)/(Tem(it)*I0mz(eint,it))
ThZT0(eint,it)=I0mz(eint,it)*(ThES(eint,it)**2)*Tem(it)/ThEke(eint,it)
ThEL(eint,it)=ThEke(eint,it)/(I0mz(eint,it)*Tem(it))
!-----------------------------------------------------------------------
if(rank==0 )then
write(2,*)ilambda,lambda/w0&
          ,Tem(it)/gama,I0mz(eint,it),I1mz(eint,it),I2mz(eint,it),ThEke(eint,it),ThES(eint,it),ThZT0(eint,it),ThEL(eint,it)&
         ,Tem(it)/efgammadd(itmaximum-1),e0(53)/efgammadd(itmaximum-5),e0(53)/gama,efgammadd(it)/gama,e0(52)/gama&
         ,Tem(it)/dble(tm(mmax-2))&
         ,(ndde0(0,51,it)+ndde0(1,51,it)+ndde0(2,51,it)+ndde0(3,51,it))/(1.0d+00*izmax)&
         ,(ndde0(0,52,it)+ndde0(1,52,it)+ndde0(2,52,it)+ndde0(3,52,it))/(1.0d+00*izmax)
!write(2,*)it,Tem(it)/gama,collector/(1.0d+00*izmax),ndcollector/(1.0d+00*izmax),nddcollector/(1.0d+00*izmax)

!write(2,*)it,Tem(it)/gama,ndde0(0,52,it),e0(52)/gama,efgammadd(it)/gama

endif
!-----------------------------------------------------------------------
end do
!-----------------------------------------------------------------------
end do  ! eint
!-----------------------------------------------------------------------
print*,ilambda,gama/w0,efgammadd(itmaximum-5)/gama,frganagama/gama
!-----------------------------------------------------------------------
end do ! ilambda
!-----------------------------------------------------------------------
call mpi_finalize(ierr)
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
call cpu_time(stop_time)
print*,"total time",stop_time - start_time
!-----------------------------------------------------------------------
close(1)
!-----------------------------------------------------------------------
close(2)
!-----------------------------------------------------------------------
! ---------------------------------------------------------------------
CONTAINS
! ---------------------------------------------------------------------
! INTEGER FUNCTION  FindMinimum():
!    This function returns the location of the minimum in the section
!    between Start and End.
! ---------------------------------------------------------------------

   INTEGER FUNCTION  FindMinimum(x, Start, Endd)
      IMPLICIT  NONE
      Double Precision, DIMENSION(1:), INTENT(IN) :: x
      INTEGER, INTENT(IN)                :: Start, Endd
      Double precision                   :: Minimum
      INTEGER                            :: Location
      INTEGER                            :: i

      Minimum  = x(Start)		! assume the first is the min
      Location = Start			! record its position
      DO i = Start+1, Endd		! start with next elements
         IF (x(i) < Minimum) THEN	!   if x(i) less than the min?
            Minimum  = x(i)		!      Yes, a new minimum found
            Location = i                !      record its position
         END IF
      END DO
      FindMinimum = Location        	! return the position
   END FUNCTION  FindMinimum

! --------------------------------------------------------------------
! SUBROUTINE  Swap():
!    This subroutine swaps the values of its two formal arguments.
! --------------------------------------------------------------------

   SUBROUTINE  Swap(a, b)
      IMPLICIT  NONE
      Double Precision, INTENT(INOUT) :: a, b
      Double Precision                :: Tempor

      Tempor = a
      a    = b
      b    = Tempor
   END SUBROUTINE  Swap

! --------------------------------------------------------------------
! SUBROUTINE  Sort():
!    This subroutine receives an array x() and sorts it into ascending
! order.
! --------------------------------------------------------------------

   SUBROUTINE  Sort(x, Sizee)
      IMPLICIT  NONE
      Double Precision, DIMENSION(1:), INTENT(INOUT) :: x
      INTEGER, INTENT(IN)                   :: Sizee
      INTEGER                               :: i
      INTEGER                               :: Location

      DO i = 1, Sizee-1			! except for the last
         Location = FindMinimum(x, i, Sizee)	! find min from this to last
         CALL  Swap(x(i), x(Location))	! swap this and the minimum
      END DO
   END SUBROUTINE  Sort
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
end program IRLMNRG
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
!-------------------------------------------------------------
!*************************************************************
!   dimension check
!*************************************************************
Integer*8 function d0(n,nmax,dbc)
implicit none
integer*8,intent(in)::n,nmax
integer*8,intent(in)::dbc
if(n<0 .or.  n> nmax)then
d0=0
else
d0=dbc
endif
end function d0
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
!----------------------------------------------------------------------
!ssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssss
! Tridiagonalization:
!----------------------------------------------------------------------
SUBROUTINE tri(dtype,l,gama,z,mmax,t,epsilonn)

use mpmodule

implicit none
!----------------------------------------------------------------------
!--VARIABLE DECLARATION:-----------------------------------------------
!----------------------------------------------------------------------
character*3, INTENT(IN) :: dtype
integer, INTENT(IN) :: mmax
type(mp_real), INTENT(IN) ::l,gama,z
type(mp_real), INTENT(OUT)::t(0:mmax),epsilonn(0:mmax)

type(mp_real)::xi0
type(mp_real)::xip(0:mmax),xin(0:mmax),gmp(0:mmax),gmn(0:mmax)&
               ,anst(0:mmax),anse(0:mmax)
type(mp_real)::umat(0:mmax,0:mmax),vmat(0:mmax,0:mmax)
integer::m,mi,iz
!----------------------------------------------------------------------
umat=mpreal(0.0d0)
vmat=mpreal(0.0d0)
!----------------------------------------------------------------------

Do mi=0,mmax
anst(mi)=mpreal(0.0d0)
anse(mi)=mpreal(0.0d0)
Do m=0,mmax
if(m==0)then
  if(dtype=='log')then
    xip(m)=(l**(1-z-m))*(mpreal(1.0d0)+ (l**(-mpreal(1.0d0))))/mpreal(2.0d0) 
    elseif(dtype=='alt')then
    xip(0)=(mpreal(1.0d0)-(l**(-z)))/(z*log(l))                            
  endif
xin(0)=-xip(0)

gmp(0)=sqrt(gama*(1-(l**(-z))))
gmn(0)=gmp(0)
else

  if(dtype=='log')then
    xip(m)=(l**(mpreal(1.0d0)-z-m))*(mpreal(1.0d0)+(l**(-mpreal(1.0d0))))/mpreal(2.0d0)
    elseif(dtype=='alt')then
    xip(m)=(mpreal(1.0d0)-(l**(-mpreal(1.0d0))))*(l**(1-m-z))/log(l)                             
  endif
  
xin(m)=-xip(m)
gmp(m)=sqrt(gama*(l**(mpreal(1.0d0)-m-z))*(1-(l**(-mpreal(1.0d0)))))
gmn(m)=gmp(m)
endif

xi0=2d+0*gama
!----------------------------------------------------------------------
! Initialization:
!----------------------------------------------------------------------
if(mi==0)then

umat(mi,m)=gmp(m)/sqrt(xi0)
vmat(mi,m)=gmn(m)/sqrt(xi0)
anse(mi)=anse(mi)+(xip(m)*((umat(mi,m))**2))+(xin(m)*((vmat(mi,m))**2))
anst(mi)=anst(mi)+(((gmp(m)**2)*((xip(m))**2))+((gmn(m)**2)*((xin(m))**2)))/Xi0

elseif(mi==1)then

umat(mi,m)=umat(0,m)*xip(m)/t(0)
vmat(mi,m)=vmat(0,m)*xin(m)/t(0)
anse(mi)=anse(mi)+(xip(m)*((umat(mi,m))**2))+(xin(m)*((vmat(mi,m))**2))
anst(mi)=anst(mi)+(xip(m)**2)*((umat(mi,m))**2)+(xin(m)**2)*((vmat(mi,m))**2) 

else

umat(mi,m)=(umat(mi-1,m)*(xip(m)-epsilonn(mi-1))-t(mi-2)*umat(mi-2,m))/t(mi-1)
vmat(mi,m)=(vmat(mi-1,m)*(xin(m)-epsilonn(mi-1))-t(mi-2)*vmat(mi-2,m))/t(mi-1)
anse(mi)=anse(mi)+(xip(m)*((umat(mi,m))**2))+(xin(m)*((vmat(mi,m))**2))
anst(mi)=anst(mi)+(xip(m)**2)*((umat(mi,m))**2)+(xin(m)**2)*((vmat(mi,m))**2)  

endif
end do !m
!----------------------------------------------------------------------
epsilonn(mi)=anse(mi)
if(mi==0)then
t(mi)=sqrt(anst(mi))
else
t(mi)=sqrt(anst(mi)-(t(mi-1)**2)-(epsilonn(mi))**2)
endif
!----------------------------------------------------------------------
!----------------------------------------------------------------------
end do !mi
!----------------------------------------------------------------------
END SUBROUTINE tri
! ---------------------------------------------------------------------
