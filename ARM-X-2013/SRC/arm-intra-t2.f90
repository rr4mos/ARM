! ==========================================================
! Programa ARM - [A]mazonas [R]odrigo [M]elissa
! Analise estatistica da conformacao em cadeias oligomericas
! ==========================================================

PROGRAM ARMintraP3HT
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
open(10,file='info.stat')
read(10,*) natol          !numero de atomos de um oligomero
read(10,*) ncadeias       !numero de cadeias 
read(10,*) naneisolig     !numero de aneis por oligomero
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
dC2    = 5
print *, calpha, dc1,s,dC2

Lp = 3.9175
if (naneisolig.eq.1) then
  Lplanar = 0.0
  else
    Lplanar = Lp*(naneisolig - 1)
end if
print *, 'Lp=', Lp, 'Lplanar=', Lplanar

open(12,file='MAT.in') 

!lendo matriz de vetores da celula unitaria (A) - sim eh transposto assim mesmo =P

read(12,*) A(1,1),A(2,1),A(3,1)
read(12,*) A(1,2),A(2,2),A(3,2)
read(12,*) A(1,3),A(2,3),A(3,3)

close(12)
sdir = 1
vstret = (/ A(1,sdir), A(2,sdir), A(3,sdir) /)
vstret = vstret/sqrt(vstret.DOT.vstret)

allocate(cx(ncadeias,natol))
allocate(cy(ncadeias,natol))
allocate(cz(ncadeias,natol))

cx = 0.0
cy = 0.0
cz = 0.0

allocate(num(ncadeias,natol))    !numero do atomo no filme
num  = 0.0

nat = 0
open(15,file='chains.in')   ! colocar arquivos do dir ./chains/
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
print *, '--> P3HT, only backbone'
print *, '--> hydrogen atoms excluded' 
print *, '-------------------------------------------------------'

!calculating centers and normal vectors

open(18,file='all.vec.rg.dat')                  !coordenadas dos centros
do l = 1, ncadeias
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

      v1 = (/ xmv4, ymv4, zmv4 /)

      v_normal(l,m) = v1 
      v_normal(l,m) = v_normal(l,m)/sqrt( v_normal(l,m) .DOT. v_normal(l,m) )
    
!      v_eixo(l,m) = v1
!      v_eixo(l,m) = v_eixo(l,m)/sqrt( v_eixo(l,m).DOT.v_eixo(l,m) )
	
      write(18,'(2i3,9f12.5)')l,m, xm(l,m), ym(l,m), zm(l,m),v_normal(l,m)!,v_eixo(l,m)

         atm2 =  atm2   +  dC2
         atm5 =  atm5   +  dC2
         atmS =  atmS   +  dC2         
   end do
end do
close(18)

open(18,file='intra.vec.1nb.dat')
do l=1, ncadeias
   do m = 1, naneisolig-1     
      v1 = (/ xm(l,m)  , ym(l,m)  , zm(l,m)   /)
      v2 = (/ xm(l,m+1), ym(l,m+1), zm(l,m+1) /)
      vnn(l,m) = v2 - v1 
      write(18,'(2i3, 4f10.5)')  l,m, vnn(l,m)
   end do
   write(18,'(2i3, 4f10.5)')  l,m, 0.0, 0.0, 0.0
end do

close(18)

!Distancia entre centros de aneis intra (todos com todos)
do l=1, ncadeias
   do m = 1, naneisolig-1    
      do o = m+1, naneisolig

         v1 = (/ xm(l,m), ym(l,m), zm(l,m) /)
         v2 = (/ xm(l,o), ym(l,o), zm(l,o) /)
         R(l,m,l,o) = sqrt((v2 - v1).DOT.(v2 - v1))
         R(l,o,l,m) = R(l,m,l,o)

      end do
   end do
end do

!Comprimento end-to-end (Kuhn factor) (relativa ao tam. planar)
open(18, file='R.ete.dat')
do l = 1, ncadeias
   write(18,*) l, R(l,1,l,naneisolig)/Lplanar
end do
close(18)

!Radius of gyration, Rg2, e analogia para ângulo, Ra2.

open(18,file='intra-TT.dat')
open(19,file='Rg2.dat')
open(20,file='Ra2.dat')

do l=1, ncadeias
   sumRg2 = 0.0d0
   sumRa2 = 0.0d0
   do m = 1, naneisolig-1    
      do o = m+1, naneisolig

      dotp  = v_normal(l,m) .DOT. v_normal(l,o)
      dotp  = dotp*0.9999999d0 ! gambiarra para instabilidade do acos. nota (*)
      lambda= dacos( dotp )*(180.d0/PI)

      sumRg2 = sumRg2 + R(l,m,l,o)**2  ! somando para calcular o Rg
      sumRa2 = sumRa2 + dotp**2
      write(18,*) l,m,l,o,R(l,m,l,o),lambda,abs(dotp)

      end do
   end do
      write(19,*) l,sumRg2/(naneisolig**2)
      write(20,*) l,sumRa2/(naneisolig**2)
end do
close(18)
close(19)
close(20)

!Vetores unitarios entre aneis. -> anisotropia.
!Correlacao de angulo entre aneis vizinhos

open(18,file='intra.anis.dat') 
open(19,file='ang-corr.dat')
open(20,file='normal-normal.dat')
do l=1, ncadeias
   sumRa2=0.0d0
   do m = 1, naneisolig-1          
      vnn(l,m) = vnn(l,m)/sqrt( vnn(l,m) .DOT. vnn(l,m) )
 
      write(18,*) l,m, vnn(l,m)
      write(18,*) m,l, real_times_vector (-1.0d0, vnn(l,m))

      dotp  = v_normal(l,m) .DOT. v_normal(l,m+1)
      write(20,*) l,m, acos(dotp*0.99999d0)*(180.d0/PI)
      
      sumRa2 = sumRa2 + abs(dotp)
   end do
   write(19,*) l,sumRa2/naneisolig
end do
close(18)
close(19)
close(20)
	
!Angulo entre dois segmentos de cadeia adjacentes // 3 aneis // tambem em desuso atualmente
open(18,file='intra.ang-sg.dat')
do l=1, ncadeias
   do m = 1, naneisolig-2
      gamma = acos( (vnn(l,m) .DOT. vnn(l,m+1))*0.9999999d0 )*(180.0/PI) 
      write(18,'(3i5, 1f10.3)') l, m, m+1, gamma      
   end do
end do
close(18)

print *, 'ARM-INTRA-P3HT: concluido!'
print *, ''
END PROGRAM ARMintraP3HT


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

