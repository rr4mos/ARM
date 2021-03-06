! ==========================================================
! Programa ARM - [A]mazonas [R]odrigo [M]elissa
! Analise estatistica da conformacao em cadeias oligomericas
! ==========================================================

PROGRAM ARMinter

USE vectors
IMPLICIT NONE

integer::i,j,k,l,n,m,o
integer::atm1,atm2,atm5

integer::nat,natol,ncadeias 
integer::natmon,naneisolig  

integer, PARAMETER::natCARBmon = 8

REAL(8),dimension(:,:),allocatable::xm, ym, zm 
REAL(8),dimension(:,:,:,:),allocatable:: R
real(8)::xmv2, ymv2, zmv2 
real(8)::xmv4, ymv4, zmv4
real(8)::xmv6, ymv6, zmv6 

TYPE (vector), DIMENSION (:,:),allocatable:: v_normal, v_eixo, vnn
TYPE (vector), DIMENSION (:,:,:,:),allocatable:: vrr
TYPE (vector)::v1,v2

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

real(8),dimension(:,:),    allocatable:: vnnx,vnny,vnnz

double precision:: dotp,dotp1,dotp2

TYPE (vector) anelanel
TYPE (vector) v_comp, v_plano
real(8):: Lplanar, Lp
real(8):: Estlm, Estno

real(8) limphi, limtheta,limR,limhpp,limtau,limdelta,limlambda
integer flgphi,flgtheta1,flgtheta2,flgR,flghpp,flgtau,flgdelta,flglambda

integer calpha, dC1,dC2,dC3 

print *, '==========================================='
print *, 'Statistical Analisis'
print *, 'ARM pack- [A]mazonas, [R]amos, [M]elissa' 
print *, '==========================================='
print *, ''

print *, '>> info.stat'
open(10,file='info.stat')
read(10,*) natol          !numero de atomos de um oligomero
read(10,*) ncadeias       !numero de cadeias 
read(10,*) naneisolig     !numero de aneis por oligomero
read(10,*) limR           ! raio de corte para avaliar estatisticas sinter-cadeia
read(10,*) calpha         ! indice do calpha da estrutura
read(10,*) dC1            ! dC1: do calpha para o calpha oposto (define centro do anel)
read(10,*) dC2            ! dC2,dC3: do calpha para outros Cs do anel (define normal do anel)
read(10,*) dC3 
close(10)

Lp = 6.667
if (naneisolig.eq.1) then
  Lplanar = 2.80
  else
    Lplanar = Lp*(naneisolig - 1)
end if
print *, 'Lp=', Lp, 'Lplanar=', Lplanar

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

allocate(cx(ncadeias,natol))
allocate(cy(ncadeias,natol))
allocate(cz(ncadeias,natol))
cx = 0.0
cy = 0.0
cz = 0.0

allocate(num(ncadeias,natol))    !numero do atomo no filme
num  = 0.0

!READING
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

!===================================
!Structural information 
!===================================

print *, 'Collecting Structural information'
print *, '-------------------------------------------------------'
print *, 'Notes:'
print *, '--> PPV phenil capped structure'
print *, '--> Each ring followed by vinil structures in xyz file' 
print *, '--> 1st ring *without* its vinil structure'
print *, '--> hydrogen atoms excluded' 
print *, '-------------------------------------------------------'

!calculating centers and normal vectors
open(18,file='all.vec.rg.dat')                  !coordenadas dos centros
do l = 1, ncadeias
   atm2 = calpha
   atm5 = calpha+dC1
   do m = 1, naneisolig

      xm(l,m)=( cx(l,atm2)+cx(l,atm5) )/2   !\vec centro
      ym(l,m)=( cy(l,atm2)+cy(l,atm5) )/2
      zm(l,m)=( cz(l,atm2)+cz(l,atm5) )/2

      xmv4=xm(l,m)-cx(l,atm2+dC2)        !\vec centro para 4 (at)
      ymv4=ym(l,m)-cy(l,atm2+dC2) 
      zmv4=zm(l,m)-cz(l,atm2+dC2)

      xmv6=xm(l,m)-cx(l,atm2+dC3)        !\vec centro para 6 (at)
      ymv6=ym(l,m)-cy(l,atm2+dC3)
      zmv6=zm(l,m)-cz(l,atm2+dC3)

      v1 = (/ xmv4, ymv4, zmv4 /)
      v2 = (/ xmv6, ymv6, zmv6 /)

      v_normal(l,m) = ( v1*v2 ) 
      v_normal(l,m) = v_normal(l,m)/sqrt( v_normal(l,m) .DOT. v_normal(l,m) )

      xmv2=cx(l,atm2)-xm(l,m)          !do centro para o carbono alpha (atm2)
      ymv2=cy(l,atm2)-ym(l,m)
      zmv2=cz(l,atm2)-zm(l,m)
     
      v_eixo(l,m) = (/ xmv2, ymv2, zmv2 /)
      v_eixo(l,m) = v_eixo(l,m)/sqrt( v_eixo(l,m).DOT.v_eixo(l,m) )

! impondo a orientacao simetrica das normais, pode pois: [0,180] = [0,90]
      dotp = v_normal(l,m-1).DOT.v_normal(l,m)
      if (dotp <= 0.0d0) then
         v_normal(l,m) = real_times_vector(-1.0d00,v_normal(l,m))
      end if

      write(18,'(2i3,9f12.5)')l,m, xm(l,m), ym(l,m), zm(l,m),v_normal(l,m)!,v_eixo(l,m)
!
!WARNING: Isso (prox. linhas) pode ser arbitrario 
!         dependendo da numeracao da cadeia.
!         A estrutura phenil capped tem esse problema de nao ter o 
!         primeiro monomero completo.
     if (atm2 .eq. calpha) then
	 atm2 = atm2 + natCARBmon - 2
	 atm5 = atm5 + natCARBmon - 2
     else
         atm2 =  atm2   +  natCARBmon
         atm5 =  atm5   +  natCARBmon
     end if
   end do
end do
close(18)

!----------------------
!B. INTER-CADEIAS
!----------------------
print *, 'interchain DATA'

!inserting periodic boundary condictions 
!(important to relevant to inter chain distances).

allocate(sx(ncadeias,naneisolig))
allocate(sy(ncadeias,naneisolig))
allocate(sz(ncadeias,naneisolig))
sx = 0.0
sy = 0.0
sz = 0.0

open(18,file='all.vec.rg-pbc.dat')                  !coord. centros pbc
do l = 1, ncadeias
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

   write(18,'(2i3,9f12.5)')l,m, xm(l,m), ym(l,m), zm(l,m),v_normal(l,m)
   end do 
end do

close(18)

allocate(vrr(ncadeias,naneisolig,ncadeias,naneisolig))

!B.1. DISTANCIA ENTRE OS CENTROS DOS ANEIS

do l=1,ncadeias
 do m = 1, naneisolig
  do n = 1, ncadeias
   do o = 1, naneisolig

   ssx = sx(n,o) - sx(l,m)
   ssx = ssx  - (dint((2.d0*ssx +3.d0)/2.d0)-1.d0)
   ssy = sy(n,o) - sy(l,m)
   ssy = ssy  - (dint((2.d0*ssy +3.d0)/2.d0)-1.d0)
   ssz = sz(n,o) - sz(l,m)
   ssz = ssz  - (dint((2.d0*ssz +3.d0)/2.d0)-1.d0)

   xxm  = A(1,1)*ssx +A(1,2)*ssy +A(1,3)*ssz 
   yym  = A(2,1)*ssx +A(2,2)*ssy +A(2,3)*ssz 
   zzm  = A(3,1)*ssx +A(3,2)*ssy +A(3,3)*ssz 

   vrr(l,m,n,o) = (/ xxm, yym, zzm /)    
   R(l,m,n,o)= sqrt( vrr(l,m,n,o).DOT.vrr(l,m,n,o) )
   end do
  end do 
 end do
end do

deallocate(sx)
deallocate(sy)
deallocate(sz)

!Estabelecendo relacao com orientacao relativa dos aneis
open(22,file='inter.pbc-PP.dat')
do l=1,ncadeias-1
   do m = 1, naneisolig
      do n = l+1, ncadeias
         do o = 1, naneisolig
	    if (R(l,m,n,o) < limR) then

            !dotp: funcao de correlacao orientacional = (|rij.ni| + 1)*(|rij.nj| + 1)/4
            anelanel = vrr(l,m,n,o)/R(l,m,n,o)
            dotp = abs(anelanel.DOT.v_normal(l,m)) + 1
            dotp = dotp*(abs(anelanel.DOT.v_normal(n,o)) +1)/4.0d0

            write(22,*) l,m,n,o,R(l,m,n,o),dotp
            end if
         end do
      end do
   end do
end do
close(22)

!TERMOS INTRA
open(22,file='intra.pbc-PP.dat')
do l=1,ncadeias
   do m = 1, naneisolig-1
      do o = m+1, naneisolig
	    n = l
            !dotp: funcao de correlacao orientacional = (|rij.ni| + 1)*(|rij.nj| + 1)/4
            anelanel = vrr(l,m,n,o)/R(l,m,n,o)
            dotp = abs(anelanel.DOT.v_normal(l,m)) + 1
            dotp = dotp*(abs(anelanel.DOT.v_normal(n,o)) +1)/4.0d0
 	    if (R(l,m,n,o) < limR) write(22,*) l,m,n,o,R(l,m,n,o),dotp

         end do
      end do
   end do

print *, 'ARM-INTER: concluido!'
print *, ''

END PROGRAM ARMinter

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

