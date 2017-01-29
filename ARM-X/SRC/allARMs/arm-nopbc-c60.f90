! ==========================================================
! Programa ARM - [A]mazonas [R]odrigo [M]elissa
! Analise estatistica da conformacao em cadeias oligomericas
! ==========================================================

PROGRAM ARMintraC60
USE vectors
IMPLICIT NONE

integer::i,j,k,l,n,m,o
integer::atm1,atm2,atm5,atmS
integer::count
integer::nat,natol,ncadeias 
integer::natmon,naneisolig  

real(8),dimension(:,:)    ,allocatable::xm, ym, zm 
real(8),dimension(:,:,:,:),allocatable:: R
real(8)::xmv2, ymv2, zmv2 
real(8)::xmv4, ymv4, zmv4
real(8)::xmv6, ymv6, zmv6 

TYPE (vector), dimension (:,:)    ,allocatable:: v_normal, v_eixo, vnn
TYPE (vector), dimension (:,:,:,:),allocatable:: vrr
TYPE (vector)::v1,v2,vstret
TYPE (vector) v_comp, v_plano
real(8) A(3,3)
integer sdir

integer,dimension(:,:),allocatable::num 
integer,dimension(:,:),allocatable::atom
real(8),dimension(:,:),allocatable::cx,cy,cz 

real(8),dimension(:,:),    allocatable:: vnnx,vnny,vnnz
real(8):: Lplanar, Lp

integer,dimension(:,:),allocatable:: sitio
integer,dimension(:,:),allocatable:: conjuga
integer conjlgth,mant

double precision:: dotp,dotp1,dotp2
real(8) gamma, tau, lambda
real(8) limR,limtau,limlambda
integer flgR,flgtau,flglambda
real(8) sumRg2, sumRa2

integer calpha, dC1,dC2,dC3,s

real(8), parameter::         PI = 3.141592654
integer, parameter:: natCARBmon = 8
character(len=25)::nimp1
character(len=10):: chainfile

print *, '==========================================='
print *, 'Statistical Analisis'
print *, 'ARM pack- [A]mazonas, [R]amos, [M]elissa' 
print *, '==========================================='
print *, ''

print *, 'Reading info.stat'
open(10,file='infoc60.stat')
read(10,*) natol          !numero de atomos de um oligomero
read(10,*) ncadeias       !numero de cadeias 
read(10,*) limR           ! raio de corte para avaliar estatisticas sinter-cadeia
read(10,*) calpha         ! indice do calpha da estrutura
read(10,*) dC1            ! dC1: do calpha para o calpha oposto (define centro do anel)
close(10)

!fixando
calpha = 1
dc1    = 56
print *, calpha, calpha+dc1
print*, natol,ncadeias,naneisolig

open(12,file='MAT.in') 

!lendo matriz de vetores da celula unitaria (A) - sim eh transposto assim mesmo =P
read(12,*) A(1,1),A(2,1),A(3,1)
read(12,*) A(1,2),A(2,2),A(3,2)
read(12,*) A(1,3),A(2,3),A(3,3)
close(12)

allocate(cx(ncadeias,natol))
allocate(cy(ncadeias,natol))
allocate(cz(ncadeias,natol))

cx = 0.0
cy = 0.0
cz = 0.0

allocate(num(ncadeias,natol))    !numero do atomo no filme
num  = 0.0

nat = 0
open(15,file='c60.in')   ! colocar arquivos do dir ./chains/
                            ! na lista 'chains.in'
do i=1,ncadeias
   read(15,'(a7)') chainfile   ! formato NNN.xyz
   open(16, file='./chainsout/'//chainfile)
   read(16,*) nimp1
   do j=1,natol
      read(16,*) nimp1, cx(i,j), cy(i,j), cz(i,j)
      nat = nat + 1
      num(i,j) = nat
   end do
   close(16)
end do
close(15)

allocate(v_normal(ncadeias,naneisolig))
allocate(v_eixo(ncadeias,naneisolig))
allocate(vnn(ncadeias,naneisolig))

allocate(xm(ncadeias,naneisolig))
allocate(ym(ncadeias,naneisolig))
allocate(zm(ncadeias,naneisolig))
xm = 0.0
ym = 0.0
zm = 0.0

allocate(R(ncadeias,naneisolig,ncadeias,naneisolig))
R = 0.0

allocate(vnnx(ncadeias,naneisolig))
allocate(vnny(ncadeias,naneisolig))
allocate(vnnz(ncadeias,naneisolig))
vnnx = 0.0
vnny = 0.0
vnnz = 0.0

allocate(sitio(ncadeias,naneisolig))
allocate(conjuga(ncadeias,naneisolig))

!===================================
!Structural information 
!===================================

print *, 'Collecting Structural information'
print *, '-------------------------------------------------------'
print *, 'Notes:'
print *, '--> C60'
print *, '-------------------------------------------------------'

!calculating centers and normal vectors

open(18,file='all.vec.rg-c60.dat')                  !coordenadas dos centros
m=1
do l = 1, ncadeias
   atm2 = calpha
   atm5 = calpha+dC1
      xm(l,m)=( cx(l,atm2)+cx(l,atm5) )/2   !\vec centro
      ym(l,m)=( cy(l,atm2)+cy(l,atm5) )/2
      zm(l,m)=( cz(l,atm2)+cz(l,atm5) )/2
      write(18,'(2i3,9f12.5)')l,m, xm(l,m), ym(l,m), zm(l,m)
end do
close(18)

!Distancia entre centros de aneis intra (todos com todos)

open(18,file='dists-c60.dat')                  

m=1
o=1
do l=1, ncadeias
do n=l+1, ncadeias-1
         v1 = (/ xm(l,m), ym(l,m), zm(l,m) /)
         v2 = (/ xm(n,o), ym(n,o), zm(n,o) /)
         R(l,m,n,o) = sqrt((v2 - v1).DOT.(v2 - v1))
         write(18,*)l,n,R(l,m,n,o)
end do
end do

print *, 'ARM-INTRA-C60: concluido!'
print *, ''
END PROGRAM ARMintraC60


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

