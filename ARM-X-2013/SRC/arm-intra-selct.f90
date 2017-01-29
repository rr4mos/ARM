! ==========================================================
! Programa ARM - [A]mazonas [R]odrigo [M]elissa
! Analise estatistica da conformacao em cadeias oligomericas
! ==========================================================

PROGRAM ARMintrasg
USE vectors
IMPLICIT NONE

integer::i,j,k,l,n,m,o
integer::atm1,atm2,atm5
integer::count
integer::nat,natol,ncadeias 
integer::natmon,naneisolig  

real(8),dimension(:,:)    ,allocatable::xm, ym, zm 
real(8),dimension(:,:,:),allocatable:: R
real(8)::xmv2, ymv2, zmv2 
real(8)::xmv4, ymv4, zmv4
real(8)::xmv6, ymv6, zmv6 

TYPE (vector), dimension (:,:)    ,allocatable:: v_normal, v_eixo, vnn
TYPE (vector), dimension (:,:,:),allocatable:: vrr
TYPE (vector)::v1,v2
TYPE (vector) v_comp, v_plano

integer,dimension(:,:),allocatable::num 
integer,dimension(:,:),allocatable::atom
real(8),dimension(:,:),allocatable::cx,cy,cz 

real(8),dimension(:,:),allocatable:: vnnx,vnny,vnnz
real(8):: Lplanar, Lp

integer,dimension(:,:),allocatable:: sitio
integer,dimension(:,:),allocatable:: conjuga
integer conjlgth,mant

double precision:: dotp,dotp1,dotp2
real(8) gamma, tau, lambda
real(8) limR,limtau,limlambda
integer flgR,flgtau,flglambda
real(8) sumRg2, sumRa2

real(8), parameter::         PI = 3.141592654
integer, parameter:: natCARBmon = 8
character(len=25)::nimp1
character(len=10):: chainfile

integer nsegments,ch,fst,lst
integer,dimension(:),allocatable::chlist
integer,dimension(:),allocatable::chlistfst
integer,dimension(:),allocatable::chlistlst


print *, '==========================================='
print *, 'Statistical Analisis'
print *, '==========================================='
print *, ''

print *, 'Reading info.stat'
open(10,file='info.stat')

read(10,*) natol          !numero de atomos de um oligomero
read(10,*) ncadeias       !numero de cadeias 
read(10,*) naneisolig     !numero de aneis por oligomero
read(10,*) limR           ! raio de corte para avaliar angulos entre-aneis e hopping e etc... inter-cadeia
read(10,*) limtau         ! corte do sitio na torcao
read(10,*) limlambda      ! corte do sitio no kink
close(10)

print *, 'new: Reading chain list // segments'
open(10,file='chainlist.in')

read(10,*) nsegments      !numero de segmentos
print *, nsegments

allocate(chlist(nsegments))       !cadeias selecionadas
allocate(chlistfst(nsegments))    !primeiro anel do segmento dentro da cadeia
allocate(chlistlst(nsegments))    !ultimo   "

do i=1,nsegments
	read(10,*) chlist(i),chlistfst(i),chlistlst(i)
	print*, chlist(i),chlistfst(i),chlistlst(i)
end do
close(10)

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
allocate(vrr(ncadeias,naneisolig,naneisolig))

allocate(xm(ncadeias,naneisolig))
allocate(ym(ncadeias,naneisolig))
allocate(zm(ncadeias,naneisolig))
xm = 0.0
ym = 0.0
zm = 0.0

allocate(R(ncadeias,naneisolig,naneisolig))
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
print *, '--> PPV phenil capped structure'
print *, '--> Each ring followed by vinil structures in xyz file' 
print *, '--> 1st ring *without* its vinil structure'
print *, '--> hydrogen atoms excluded' 
print *, '-------------------------------------------------------'

!calculating centers and normal vectors
!checado (vpython) / ok.
open(18,file='all.vec.rg.selct.dat')                  !coordenadas dos centros

do l = 1, ncadeias
   atm2 = 2
   atm5 = 5
   do m = 1, naneisolig

      xm(l,m)=( cx(l,atm2)+cx(l,atm5) )/2   !\vec centro
      ym(l,m)=( cy(l,atm2)+cy(l,atm5) )/2
      zm(l,m)=( cz(l,atm2)+cz(l,atm5) )/2

      xmv4=xm(l,m)-cx(l,atm2+5)        !\vec centro para 4 (at)
      ymv4=ym(l,m)-cy(l,atm2+5) 
      zmv4=zm(l,m)-cz(l,atm2+5)

      xmv6=xm(l,m)-cx(l,atm5-2)        !\vec centro para 6 (at)
      ymv6=ym(l,m)-cy(l,atm5-2)
      zmv6=zm(l,m)-cz(l,atm5-2)

      v1 = (/ xmv4, ymv4, zmv4 /)
      v2 = (/ xmv6, ymv6, zmv6 /)

      v_normal(l,m) = ( v1*v2 ) 
      v_normal(l,m) = v_normal(l,m)/sqrt( v_normal(l,m) .DOT. v_normal(l,m) )

! diversas possibilidades de definir o vetor do eixo do anel. 
! *nota: a opcao usando os vetores 4 e 6 tem a vantagem de preservar c/ mais rigor a ortogonalidade
!        com o vetor normal. >> melhor para determinar as rotacoes // precisa avaliar a imprecisao.
!
!      xmv2=( cx(l,atm2)-cx(l,atm5) )   !do atm5 para o atm2
!      ymv2=( cy(l,atm2)-cy(l,atm5) )
!      zmv2=( cz(l,atm2)-cz(l,atm5) )
!
!      usando os vetores que sao perpend. a normal.
!
!      v_eixo(l,m)= v1+v2                ! usando a simetria dos vetores 4 e 6
!

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
     if (atm2 .eq. 2) then
	 atm2 = atm2 + 6
	 atm5 = atm5 + 6
     else
         atm2 =  atm2   +  natCARBmon
         atm5 =  atm5   +  natCARBmon         
     end if
   end do
end do
close(18)
close(19)

!open(18,file='intra.vec.selct.dat')

do i=1,nsegments
	l=chlist(i)
	do m=chlistfst(i),chlistlst(i)
	do n=chlistfst(i),chlistlst(i)

        v1 = (/ xm(l,m), ym(l,m), zm(l,m) /)
        v2 = (/ xm(l,n), ym(l,n), zm(l,n) /)
      
        vrr(l,m,n) = v2 - v1 
        R(l,m,n) = sqrt((v2 - v1).DOT.(v2 - v1))
        R(l,n,m) = R(l,m,n)
!        write(18,'(2i3, 3f10.5)') l,m, vrr(l,m,n)
	enddo
        enddo
end do
!close(18)

open(19,file='centerpos.xyz')
do i=1,nsegments
   l=chlist(i)
   do m=chlistfst(i),chlistlst(i)
   write(19,'(a3,3f12.5)') 'C',xm(l,m), ym(l,m), zm(l,m)
   end do
end do
close(19)

!correlacao da orientacao relativa
open(18,file='intra.rg.selct.dat')
do i=1,nsegments
   l=chlist(i)
   do m=chlistfst(i),chlistlst(i)
   do n=chlistfst(i),chlistlst(i)
   if ((n .gt. m) .and. (R(l,m,n) < 8.0)) then
     !LAMBDA // angulo entre eixos de aneis, atualmente normal/normal
      dotp1  = v_normal(l,m) .DOT. ( vrr(l,m,n)/R(l,m,n) )
      dotp2  = v_normal(l,n) .DOT. ( vrr(l,m,n)/R(l,m,n) )

      write(18,'(3i5, 2f10.3)') l, m, n,R(l,m,n),(abs(dotp1)+1)*(abs(dotp2)+1)/4
   end if
   enddo
   enddo
end do
close(18)

print *, 'ARM-INTRA com restricao de segmento: concluido!'
print *, ''
END PROGRAM ARMintrasg


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

