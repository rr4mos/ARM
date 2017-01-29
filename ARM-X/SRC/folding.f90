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
print *, '--> folding inside the unitcell'
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

allocate(cx(ncadeiaspt,natol))
allocate(cy(ncadeiaspt,natol))
allocate(cz(ncadeiaspt,natol))
cx = 0.0
cy = 0.0
cz = 0.0

allocate(sx(ncadeiaspt,natol))
allocate(sy(ncadeiaspt,natol))
allocate(sz(ncadeiaspt,natol))
sx = 0.0
sy = 0.0
sz = 0.0

open(15,file='chains.in')   ! colocar arquivos do dir ./chains/
open(98,file='allp3htpbc.xyz')
write(98,'(i10)')natol*ncadeiaspt
write(98,*)

do l = 1, ncadeiaspt
   read(15,'(a7)') chainfile   ! formato NNN.xyz
   open(16, file='./chainsout/'//chainfile)
   read(16,*) nimp1

   open(99, file='./cpbc/'//chainfile)
   write(99, '(a3)'),'330'
   write(99,*)

   do m = 1, natol
   read(16,*) nimp1, cx(l,m), cy(l,m), cz(l,m)

   sx(l,m) = B(1,1)*cx(l,m)+B(1,2)*cy(l,m)+B(1,3)*cz(l,m)
   sy(l,m) = B(2,1)*cx(l,m)+B(2,2)*cy(l,m)+B(2,3)*cz(l,m)
   sz(l,m) = B(3,1)*cx(l,m)+B(3,2)*cy(l,m)+B(3,3)*cz(l,m)

   IPBX = dint(sx(l,m)+1.d0)-1
   IPBY = dint(sy(l,m)+1.d0)-1
   IPBZ = dint(sz(l,m)+1.d0)-1

   sx(l,m)= sx(l,m) - IPBX
   sy(l,m)= sy(l,m) - IPBY
   sz(l,m)= sz(l,m) - IPBZ

   cx(l,m) = A(1,1)*sx(l,m)+A(1,2)*sy(l,m)+A(1,3)*sz(l,m)
   cy(l,m) = A(2,1)*sx(l,m)+A(2,2)*sy(l,m)+A(2,3)*sz(l,m)
   cz(l,m) = A(3,1)*sx(l,m)+A(3,2)*sy(l,m)+A(3,3)*sz(l,m)

   write(99, '(a1,3f10.3)'),nimp1, cx(l,m), cy(l,m), cz(l,m)
   write(98, '(a1,3f10.3)'),nimp1, cx(l,m), cy(l,m), cz(l,m)
   end do 
end do
close(15)
close(16)
close(99)
close(98)

deallocate(cx)
deallocate(cy)
deallocate(cz)

print *, 'Reading infoc60.stat'
open(10,file='infoc60.stat')
read(10,*) natol          !numero de atomos de um oligomero
read(10,*) ncadeiasc60    !numero de cadeias 
read(10,*) limR           ! raio de corte para avaliar estatisticas sinter-cadeia
read(10,*) calpha         ! indice do calpha da estrutura
read(10,*) dC1            ! dC1: do calpha para o calpha oposto (define centro do anel)
close(10)

allocate(cx(ncadeiasc60,natol))
allocate(cy(ncadeiasc60,natol))
allocate(cz(ncadeiasc60,natol))
cx = 0.0
cy = 0.0
cz = 0.0

allocate(svx(ncadeiasc60,natol))
allocate(svy(ncadeiasc60,natol))
allocate(svz(ncadeiasc60,natol))
svx=0.0
svy=0.0
svz=0.0

open(15,file='c60.in')
open(98,file='allc60pbc.xyz')
write(98,'(i10)')natol*ncadeiasc60
write(98,*)

open(15,file='c60.in')   ! colocar arquivos do dir ./chains/
do l = 1, ncadeiasc60

   read(15,'(a7)') chainfile   ! formato NNN.xyz
   open(16, file='./chainsout/'//chainfile)
   read(16,*) nimp1

   open(99, file='./cpbc/'//chainfile)
   write(99, '(a2)'),'60'
   write(99,*)

   do m=1, natol
   read(16,*) nimp1, cx(l,m), cy(l,m), cz(l,m)

   svx(l,m) = B(1,1)*cx(l,m)+B(1,2)*cy(l,m)+B(1,3)*cz(l,m)
   svy(l,m) = B(2,1)*cx(l,m)+B(2,2)*cy(l,m)+B(2,3)*cz(l,m)
   svz(l,m) = B(3,1)*cx(l,m)+B(3,2)*cy(l,m)+B(3,3)*cz(l,m)

   IPBX = dint(svx(l,m)+1.d0)-1
   IPBY = dint(svy(l,m)+1.d0)-1
   IPBZ = dint(svz(l,m)+1.d0)-1

   svx(l,m)= svx(l,m) - IPBX
   svy(l,m)= svy(l,m) - IPBY
   svz(l,m)= svz(l,m) - IPBZ

   cx(l,m) = A(1,1)*svx(l,m)+A(1,2)*svy(l,m)+A(1,3)*svz(l,m)
   cy(l,m) = A(2,1)*svx(l,m)+A(2,2)*svy(l,m)+A(2,3)*svz(l,m)
   cz(l,m) = A(3,1)*svx(l,m)+A(3,2)*svy(l,m)+A(3,3)*svz(l,m)

   write(99, '(a1,3f10.3)'),nimp1, cx(l,m), cy(l,m), cz(l,m)
   write(98, '(a1,3f10.3)'),nimp1, cx(l,m), cy(l,m), cz(l,m)

   end do
end do

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

