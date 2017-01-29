! ==========================================================
! Programa ARM - [A]mazonas [R]odrigo [M]elissa
! Analise estatistica da conformacao em cadeias oligomericas
! ==========================================================

PROGRAM ARMinterPTC60

USE vectors
IMPLICIT NONE

integer::i,j,k,l,n,m,o
integer::atm1,atm2,atm5,atmS,s

integer::nat,natol,ncadeiaspt,ncadeiasc60
integer::natmon,naneisolig  

integer, PARAMETER::natCARBmon = 8

REAL(8),dimension(:,:),allocatable::xm, ym, zm 
REAL(8),dimension(:,:,:,:),allocatable:: R
REAL(8),dimension(:,:),allocatable::xvm, yvm, zvm 
real(8)::xmv2, ymv2, zmv2 
real(8)::xmv4, ymv4, zmv4
real(8)::xmv6, ymv6, zmv6 

TYPE (vector), DIMENSION (:,:),allocatable:: v_normal, vnn
TYPE (vector), DIMENSION (:,:),allocatable:: v_eixo, v_eixo_pvm
TYPE (vector):: v_eixo_pv1, v_eixo_pv2
TYPE (vector), DIMENSION (:,:,:,:),allocatable:: vrr
TYPE (vector)::v1,v2

TYPE (vector), DIMENSION (:,:),allocatable:: vv1_normal,vv2_normal

character(len=25)::nimp1
character(len=45)::arquivo,arquivo1,arquivo2

real(8), PARAMETER::PI= 3.141592654

integer,dimension(:,:),allocatable::num 
integer,dimension(:,:),allocatable::atom
character(len=10):: chainfile   
real(8),dimension(:,:),allocatable::cx,cy,cz 

real(8):: alpha, beta, gamma, theta1,theta2, phi, tau, lambda

real(8) A(3,3),B(3,3)         !matriz de vetores da celula unitaria e inversa
integer(8) IPBX,IPBY,IPBZ

real(8),dimension(:,:),    allocatable:: sx,sy,sz          ! ...para posicoes
real(8):: ssx,ssy,ssz                                      ! ...para distancias
real(8):: xxm,yym,zzm   

real(8),dimension(:,:),    allocatable:: svx,svy,svz          ! ...para posicoes

real(8),dimension(:,:),    allocatable:: vnnx,vnny,vnnz

double precision:: dotp,dotp1,dotp2,fco1,fco2

TYPE (vector) anelanel
TYPE (vector) v_comp, v_plano
real(8):: Lplanar, Lp
real(8):: Estlm, Estno

real(8) limphi, limtheta,limR,limhpp,limtau,limdelta,limlambda
integer flgphi,flgtheta1,flgtheta2,flgR,flghpp,flgtau,flgdelta,flglambda

integer calpha, dC1,dC2,dC3 
integer cvinil, dCv1,dCv2,dCv3 

print *, '==========================================='
print *, 'Statistical Analisis'
print *, 'ARM pack- [A]mazonas, [R]amos, [M]elissa' 
print *, '==========================================='
print *, ''

!lendo matriz de vetores da celula unitaria (A) - sim eh transposto assim mesmo =P
open(12,file='MAT.in') 
read(12,*) A(1,1),A(2,1),A(3,1)
read(12,*) A(1,2),A(2,2),A(3,2)
read(12,*) A(1,3),A(2,3),A(3,3)
close(12)

! invertendo matriz vetores da celula unitaria (B)
do i = 1, 3
   do j =1,3
      B(i,j)=A(i,j)
   end do
end do
call INV(B,3)

open(14,file='MAT.inv')
write(14,*) B(1,1),B(1,2),B(1,3)
write(14,*) B(2,1),B(2,2),B(2,3)
write(14,*) B(3,1),B(3,2),B(3,3)
close(14)

!===================================
!Structural information 
!===================================

print *, 'Collecting Structural information'
print *, '-------------------------------------------------------'
print *, 'Notes:'
print *, '--> C60 + P3HT'
print *, '-------------------------------------------------------'

print*, '>>processando p3ht'

print *, 'Reading info.stat'
open(10,file='info.stat')
read(10,*) natol          !numero de atomos de um oligomero
read(10,*) ncadeiaspt      !numero de cadeias 
read(10,*) naneisolig      !numero de aneis por oligomero
read(10,*) limR           ! raio de corte para avaliar estatisticas sinter-cadeia
read(10,*) calpha         ! indice do calpha da estrutura
read(10,*) dC1            ! dC1: do calpha para o calpha oposto (define centro do anel)
read(10,*) s              ! sulfur atom
read(10,*) dC2            ! passo ateh o prox. primeiro calpha
close(10)

!fixando
calpha = 1
dc1    = 3
s      = 5 
dC2    = 11
print *, calpha, dc1,s,dC2
print *, ncadeiaspt,natol
Lp = 3.9175
if (naneisolig.eq.1) then
  Lplanar = 0.0
  else
    Lplanar = Lp*(naneisolig - 1)
end if
print *, 'Lp=', Lp, 'Lplanar=', Lplanar

allocate(cx(ncadeiaspt,natol))
allocate(cy(ncadeiaspt,natol))
allocate(cz(ncadeiaspt,natol))
cx = 0.0
cy = 0.0
cz = 0.0

!READING
open(15,file='chains.in')   ! colocar arquivos do dir ./chains/
                            ! na lista 'chains.in'
do i=1,ncadeiaspt
   read(15,'(a7)') chainfile   ! formato NNN.xyz
   open(16, file='./chainsout/'//chainfile)
   read(16,*) nimp1
   do j=1,natol
      read(16,*) nimp1, cx(i,j), cy(i,j), cz(i,j)
   end do
   close(16)
end do
close(15)

allocate(v_normal(ncadeiaspt,naneisolig))
allocate(vnn(ncadeiaspt,naneisolig))
allocate(xm(ncadeiaspt,naneisolig))
allocate(ym(ncadeiaspt,naneisolig))
allocate(zm(ncadeiaspt,naneisolig))
xm = 0.0
ym = 0.0
zm = 0.0

open(18,file='all.vec.rg-p3ht.dat')                  !coordenadas dos centros
do l = 1, ncadeiaspt
   atm2 = calpha
   atm5 = calpha+dC1
   atmS = s
   do m = 1, naneisolig
      xm(l,m)=( cx(l,atm2)+cx(l,atm5) )/2   !\vec centro
      ym(l,m)=( cy(l,atm2)+cy(l,atm5) )/2
      zm(l,m)=( cz(l,atm2)+cz(l,atm5) )/2

      xmv4=cx(l,atmS)-xm(l,m)        !\vec centro para 4 (at)
      ymv4=cy(l,atmS)-ym(l,m) 
      zmv4=cz(l,atmS)-zm(l,m)

      v_normal(l,m) = (/ xmv4, ymv4, zmv4 /)
      v_normal(l,m) = v_normal(l,m)/sqrt( v_normal(l,m) .DOT. v_normal(l,m) )
	
      write(18,'(2i3,9f12.5)')l,m, xm(l,m), ym(l,m), zm(l,m),v_normal(l,m)

      atm2 =  atm2   +  dC2
      atm5 =  atm5   +  dC2
      atmS =  atmS   +  dC2         
   end do
end do
close(18)

print*, '=========================================================================='
print*, '>>processando c60'

print *, 'Reading infoc60.stat'
open(10,file='infoc60.stat')
read(10,*) natol          !numero de atomos de um oligomero
read(10,*) ncadeiasc60    !numero de cadeias 
read(10,*) limR           ! raio de corte para avaliar estatisticas sinter-cadeia
read(10,*) calpha         ! indice do calpha da estrutura
read(10,*) dC1            ! dC1: do calpha para o calpha oposto (define centro do anel)
close(10)

!fixando pra garantir
calpha = 1
dc1    = 56
print *, calpha, calpha+dc1
print*, natol,ncadeiasc60

deallocate(cx)
deallocate(cy)
deallocate(cz)

allocate(cx(ncadeiasc60,natol))
allocate(cy(ncadeiasc60,natol))
allocate(cz(ncadeiasc60,natol))
cx = 0.0
cy = 0.0
cz = 0.0

allocate(xvm(ncadeiasc60,1))
allocate(yvm(ncadeiasc60,1))
allocate(zvm(ncadeiasc60,1))
xvm = 0.0
yvm = 0.0
zvm = 0.0

open(15,file='c60.in')   ! colocar arquivos do dir ./chains/
                            ! na lista 'chains.in'
do i=1,ncadeiasc60
   read(15,'(a7)') chainfile   ! formato NNN.xyz
   open(16, file='./chainsout/'//chainfile)
   read(16,*) nimp1
   do j=1,natol
      read(16,*) nimp1, cx(i,j), cy(i,j), cz(i,j)
   end do
   close(16)
end do
close(15)

open(18,file='all.vec.rg-c60.dat')                  !coordenadas dos centros
m=1
do l = 1, ncadeiasc60
   atm2 = calpha
   atm5 = calpha+dC1
      xvm(l,m)=( cx(l,atm2)+cx(l,atm5) )/2   !\vec centro
      yvm(l,m)=( cy(l,atm2)+cy(l,atm5) )/2
      zvm(l,m)=( cz(l,atm2)+cz(l,atm5) )/2
      write(18,'(2i3,9f12.5)')l,m, xvm(l,m), yvm(l,m), zvm(l,m)
end do
close(18)

!PBC
print *, 'inserting periodic boundary condictions....'

allocate(sx(ncadeiaspt,naneisolig))
allocate(sy(ncadeiaspt,naneisolig))
allocate(sz(ncadeiaspt,naneisolig))

allocate(svx(ncadeiasc60,1))
allocate(svy(ncadeiasc60,1))
allocate(svz(ncadeiasc60,1))
sx = 0.0
sy = 0.0
sz = 0.0
svx = 0.0
svy = 0.0
svz = 0.0

open(18,file='all.vec.pbc-p3ht.dat')                  !coord. centros pbc
do l = 1, ncadeiaspt
   do m = 1, naneisolig

   sx(l,m) = B(1,1)*xm(l,m)+B(1,2)*ym(l,m)+B(1,3)*zm(l,m)
   sy(l,m) = B(2,1)*xm(l,m)+B(2,2)*ym(l,m)+B(2,3)*zm(l,m)
   sz(l,m) = B(3,1)*xm(l,m)+B(3,2)*ym(l,m)+B(3,3)*zm(l,m)

   IPBX = dint(sx(l,m)+1.d0)-1
   IPBY = dint(sy(l,m)+1.d0)-1
   IPBZ = dint(sz(l,m)+1.d0)-1

   sx(l,m)= sx(l,m) - IPBX
   sy(l,m)= sy(l,m) - IPBY
   sz(l,m)= sz(l,m) - IPBZ

   xm(l,m) = A(1,1)*sx(l,m)+A(1,2)*sy(l,m)+A(1,3)*sz(l,m)
   ym(l,m) = A(2,1)*sx(l,m)+A(2,2)*sy(l,m)+A(2,3)*sz(l,m)
   zm(l,m) = A(3,1)*sx(l,m)+A(3,2)*sy(l,m)+A(3,3)*sz(l,m)

   write(18,'(2i3,6f12.5)')l,m, xm(l,m), ym(l,m), zm(l,m),v_normal(l,m)
   end do 
end do
close(18)

print*, 'p3ht --> ok'

!C60
m=1
open(18,file='all.vec.pbc-c60.dat')
do l = 1, ncadeiasc60

   svx(l,m) = B(1,1)*xvm(l,m)+B(1,2)*yvm(l,m)+B(1,3)*zvm(l,m)
   svy(l,m) = B(2,1)*xvm(l,m)+B(2,2)*yvm(l,m)+B(2,3)*zvm(l,m)
   svz(l,m) = B(3,1)*xvm(l,m)+B(3,2)*yvm(l,m)+B(3,3)*zvm(l,m)

   IPBX = dint(svx(l,m)+1.d0)-1
   IPBY = dint(svy(l,m)+1.d0)-1
   IPBZ = dint(svz(l,m)+1.d0)-1

   svx(l,m)= svx(l,m) - IPBX
   svy(l,m)= svy(l,m) - IPBY
   svz(l,m)= svz(l,m) - IPBZ

   xvm(l,m) = A(1,1)*svx(l,m)+A(1,2)*svy(l,m)+A(1,3)*svz(l,m)
   yvm(l,m) = A(2,1)*svx(l,m)+A(2,2)*svy(l,m)+A(2,3)*svz(l,m)
   zvm(l,m) = A(3,1)*svx(l,m)+A(3,2)*svy(l,m)+A(3,3)*svz(l,m)

   write(18,'(i3,9f12.5)')l,xvm(l,m), yvm(l,m), zvm(l,m)
end do
close(18)
print*,'c60 ---> ok!'
print*,'pbc ok!'

allocate(vrr(ncadeiasc60,1,ncadeiaspt,naneisolig))
allocate(  R(ncadeiasc60,1,ncadeiaspt,naneisolig))
R   = 0.0d0

print*, 'calculando correlacao c60 x p3ht......'

do l=1,ncadeiasc60
  m=1
  do n = 1, ncadeiaspt
   do o = 1, naneisolig

   ssx = sx(n,o) - svx(l,m)
   ssx = ssx  - (dint((2.d0*ssx +3.d0)/2.d0)-1.d0)
   ssy = sy(n,o) - svy(l,m)
   ssy = ssy  - (dint((2.d0*ssy +3.d0)/2.d0)-1.d0)
   ssz = sz(n,o) - svz(l,m)
   ssz = ssz  - (dint((2.d0*ssz +3.d0)/2.d0)-1.d0)

   xxm  = A(1,1)*ssx +A(1,2)*ssy +A(1,3)*ssz 
   yym  = A(2,1)*ssx +A(2,2)*ssy +A(2,3)*ssz 
   zzm  = A(3,1)*ssx +A(3,2)*ssy +A(3,3)*ssz 

   vrr(l,m,n,o) = (/ xxm, yym, zzm /)
   R(l,m,n,o)= sqrt( vrr(l,m,n,o).DOT.vrr(l,m,n,o) )
!   print*, 'R:',l,n,o,R(l,m,n,o), vrr(l,m,n,o)
   end do
  end do 
end do

print*, 'vrr e R calculados ....'

deallocate(sx)
deallocate(sy)
deallocate(sz)
deallocate(svx)
deallocate(svy)
deallocate(svz)

!B.2. Angulos dentro de um raio de corte
!Estabelecendo relacao com orientacao relativa dos aneis

print*, 'limR:', limR

open(22,file='inter.pbc-p3ht+c60.dat')
do l=1,ncadeiasc60
      m=1
      do n = 1, ncadeiaspt
         do o = 1, naneisolig
	    if (R(l,m,n,o) < limR) then

            anelanel = vrr(l,m,n,o)/R(l,m,n,o)
            dotp1 = anelanel.DOT.v_normal(n,o)
	
            write(22,*) l+20,n,o,R(l,m,n,o),dotp1!,vrr(l,m,n,o)
            end if
         end do
      end do
end do
close(22)

	
print *, 'ARM-INTER-P3HTC60: concluido!'

END PROGRAM ARMinterPTC60

SUBROUTINE INV(ARRAY,N)
	implicit double precision (a-h,o-z)
        dimension ARRAY(N,N)
		!AMAX(n,n),SAV
!  troquei os dimensions ik e jk de 4 p/ 16 .
        DIMENSION IK(N),JK(N)
11      DO 100 K=1,N
        AMAX=0.0d0
21      DO 30 I=K,N
        DO 30 J=K,N
23      IF(dABS(AMAX)-dABS(ARRAY(I,J)))24,24,30
24      AMAX=ARRAY(I,J)
        IK(K)=I
        JK(K)=J
30      CONTINUE
41      I=IK(K)
        IF(I-K)21,51,43
43      DO 50 J=1,N
        SAV=ARRAY(K,J)
        ARRAY(K,J)=ARRAY(I,J)
50      ARRAY(I,J)=-SAV
51      J=JK(K)
        IF(J-K)21,61,53
53      DO 60 I=1,N
        SAV=ARRAY(I,K)
        ARRAY(I,K)=ARRAY(I,J)
60      ARRAY(I,J)=-SAV
61      DO 70 I=1,N
        IF(I-K)63,70,63
63      ARRAY(I,K)=-ARRAY(I,K)/AMAX
70      CONTINUE
71      DO 80 I=1,N
        DO 80 J=1,N
        IF(I-K)74,80,74
74      IF(J-K)75,80,75
75      ARRAY(I,J)=ARRAY(I,J)+ARRAY(I,K)*ARRAY(K,J)
80      CONTINUE
81      DO 90 J=1,N
        IF(J-K)83,90,83
83      ARRAY(K,J)=ARRAY(K,J)/AMAX
90      CONTINUE
        ARRAY(K,K)=1./AMAX
100     CONTINUE
101     DO 130 L=1,N
        K=N-L+1
        J=IK(K)
        IF(J-K)111,111,105
105     DO 110 I=1,N
        SAV=ARRAY(I,K)
        ARRAY(I,K)=-ARRAY(I,J)
110     ARRAY(I,J)=SAV
111     I=JK(K)
        IF(I-K)130,130,113
113     DO 120 J=1,N
        SAV=ARRAY(K,J)
        ARRAY(K,J)=-ARRAY(I,J)
120     ARRAY(I,J)=SAV
130     CONTINUE
140     RETURN
        END

