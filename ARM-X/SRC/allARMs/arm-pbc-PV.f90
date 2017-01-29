! ==========================================================
! Programa ARM - [A]mazonas [R]odrigo [M]elissa
! Analise estatistica da conformacao em cadeias oligomericas
! ==========================================================

PROGRAM ARMinterPV

USE vectors
IMPLICIT NONE

integer::i,j,k,l,n,m,o
integer::atm1,atm2,atm5

integer::nat,natol,ncadeias 
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
TYPE (vector)::v1,v2,vstret

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

integer calpha, dC1,dC2,dC3 ,sdir
integer cvinil, dCv1,dCv2,dCv3 

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

sdir = 1
vstret = (/ A(1,sdir), A(2,sdir), A(3,sdir) /)
vstret = vstret/sqrt(vstret.DOT.vstret)

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
allocate(v_eixo_pvm(ncadeias,naneisolig-1))

allocate(vnn(ncadeias,naneisolig))
allocate(xm(ncadeias,naneisolig))
allocate(ym(ncadeias,naneisolig))
allocate(zm(ncadeias,naneisolig))
xm = 0.0
ym = 0.0
zm = 0.0

allocate(vv1_normal(ncadeias,naneisolig-1))
allocate(vv2_normal(ncadeias,naneisolig-1))
allocate(xvm(ncadeias,naneisolig-1))
allocate(yvm(ncadeias,naneisolig-1))
allocate(zvm(ncadeias,naneisolig-1))
xvm = 0.0
yvm = 0.0
zvm = 0.0

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

      xmv2=cx(l,atm2)-xm(l,m)            !do centro para o carbono alpha (atm2)
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

print*, '>>processando vinileno'


! nao generalizado pra leitura ainda
cvinil =13
dCv1 = 1
dCv2 = -2
dCv3 = -14

open(18,file='all.vec.vinil.dat')
open(19,file='dihedral.dat')
!VINIL 
do l = 1, ncadeias
   atm2 = cvinil
   atm5 = cvinil+dCv1                        
   do m = 1, naneisolig-1                    !1 vinil a menos que aneis
!      print*, atm2, atm5
      xvm(l,m)=( cx(l,atm2)+cx(l,atm5) )/2   !\vec centro
      yvm(l,m)=( cy(l,atm2)+cy(l,atm5) )/2
      zvm(l,m)=( cz(l,atm2)+cz(l,atm5) )/2

!      print*, atm2, '->', atm5      !cvinil-f -> cvinil-i
      xmv4= cx(l,atm5)-cx(l,atm2) 
      ymv4= cy(l,atm5)-cy(l,atm2)
      zmv4= cz(l,atm5)-cz(l,atm2)

!      print*,  atm2,'->', atm2+dCv2  !cvinil-f -> calpha-f
      xmv6= cx(l,atm2+dCv2)-cx(l,atm2)        
      ymv6= cy(l,atm2+dCv2)-cy(l,atm2)
      zmv6= cz(l,atm2+dCv2)-cz(l,atm2)

      v1 = (/ xmv4, ymv4, zmv4 /)
      v2 = (/ xmv6, ymv6, zmv6 /)

      vv1_normal(l,m) = ( v1*v2 ) 
      vv1_normal(l,m) = vv1_normal(l,m)/sqrt( vv1_normal(l,m) .DOT. vv1_normal(l,m) )

!      print*, atm5,'->',atm2  ! cvinil-i -> cvinil-f
      xmv4= -xmv4
      ymv4= -ymv4
      zmv4= -zmv4

!      print*, atm5,'->',atm5+dCv3 ! cvinil-i -> calpha
      if (atm2 .ne. cvinil) then
      xmv6= cx(l,atm5+dCv3)-cx(l,atm5) 
      ymv6= cy(l,atm5+dCv3)-cy(l,atm5)
      zmv6= cz(l,atm5+dCv3)-cz(l,atm5)
      else
      xmv6= cx(l,atm5+dCv3+2)-cx(l,atm5) 
      ymv6= cy(l,atm5+dCv3+2)-cy(l,atm5)
      zmv6= cz(l,atm5+dCv3+2)-cz(l,atm5)
      end if

      v1 = (/ xmv4, ymv4, zmv4 /)
      v2 = (/ xmv6, ymv6, zmv6 /)

      vv2_normal(l,m) = ( v2*v1 ) 
      vv2_normal(l,m) = vv2_normal(l,m)/sqrt( vv2_normal(l,m) .DOT. vv2_normal(l,m) )

      dotp = vv1_normal(l,m).DOT.vv2_normal(l,m)

      write(19,*) l,m,dacos(dotp*0.9999999d0)*180.d0/PI

      if (dotp < 0.0d0) then
!         vv2_normal(l,m) = real_times_vector(-1.0d00,vv2_normal(l,m))
      end if

      dotp = vv1_normal(l,m).DOT.v_normal(l,m)
      if (dotp < 0.0d0) then
 !        vv2_normal(l,m) = real_times_vector(-1.0d00,vv2_normal(l,m))
 !        vv1_normal(l,m) = real_times_vector(-1.0d00,vv1_normal(l,m))
      end if
      write(18,'(2i3,10f12.5)')l,m, xvm(l,m), yvm(l,m), zvm(l,m),vv1_normal(l,m),vv2_normal(l,m),abs(dotp)

      atm2 =  atm2   +  natCARBmon
      atm5 =  atm5   +  natCARBmon         
   end do
   write(18,'(2i3,a20)') l,m,'x x x x x x x x x'       ! completamento para posproc.
end do
close(18)
close(19)

!calculando vetores  fenil-vinil, para analisar o paralelismo local.
do l=1,ncadeias
 do m = 1, naneisolig-1
   xxm =  xvm(l,m) - xm(l,m) 
   yym =  yvm(l,m) - ym(l,m) 
   zzm =  zvm(l,m) - zm(l,m) 

   v_eixo_pv1 = (/ xxm, yym, zzm /) 
   v_eixo_pv1 = v_eixo_pv1/sqrt(v_eixo_pv1.DOT.v_eixo_pv1)

   xxm =  xm(l,m+1) - xvm(l,m)  
   yym =  ym(l,m+1) - yvm(l,m) 
   zzm =  zm(l,m+1) - zvm(l,m) 

   v_eixo_pv2 = (/ xxm, yym, zzm /) 
   v_eixo_pv2 = v_eixo_pv2/sqrt(v_eixo_pv2.DOT.v_eixo_pv2)

   v_eixo_pvm(l,m) = (v_eixo_pv1+v_eixo_pv2)/2.d0
   v_eixo_pvm(l,m) =  v_eixo_pvm(l,m)/sqrt(v_eixo_pvm(l,m).DOT.v_eixo_pvm(l,m))
end do
end do

!----------------------
!B. INTER-CADEIAS
!----------------------
print *, 'interchain DATA'

!inserting periodic boundary condictions 
!(important to relevant to inter chain distances).

allocate(sx(ncadeias,naneisolig))
allocate(sy(ncadeias,naneisolig))
allocate(sz(ncadeias,naneisolig))

allocate(svx(ncadeias,naneisolig))
allocate(svy(ncadeias,naneisolig))
allocate(svz(ncadeias,naneisolig))
sx = 0.0
sy = 0.0
sz = 0.0
svx = 0.0
svy = 0.0
svz = 0.0

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

   write(18,'(2i3,6f12.5)')l,m, xm(l,m), ym(l,m), zm(l,m),v_normal(l,m)
   end do 
end do
close(18)
 
!VINIL
open(18,file='all.vec.vinil-pbc.dat')
do l = 1, ncadeias
   do m = 1, naneisolig-1

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

   write(18,'(2i3,9f12.5)')l,m, xvm(l,m), yvm(l,m), zvm(l,m),vv1_normal(l,m),vv2_normal(l,m)
   end do 
   write(18,'(2i3,a20)') l,m,'x x x x x x x x x'
end do
close(18)

print*,'pbc ok!'

allocate(vrr(ncadeias,naneisolig,ncadeias,naneisolig-1))

!B.1. DISTANCIA ENTRE OS CENTROS DOS ANEIS x VINIS
do l=1,ncadeias
 do m = 1, naneisolig
  do n = 1, ncadeias
   do o = 1, naneisolig-1

   ssx = sx(l,m) - svx(n,o)
   ssx = ssx  - (dint((2.d0*ssx +3.d0)/2.d0)-1.d0)
   ssy = sy(l,m) - svy(n,o)
   ssy = ssy  - (dint((2.d0*ssy +3.d0)/2.d0)-1.d0)
   ssz = sz(l,m) - svz(n,o)
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
deallocate(svx)
deallocate(svy)
deallocate(svz)

!B.2. Angulos inter_cadeia dentro de um raio de corte
!Estabelecendo relacao com orientacao relativa dos aneis

open(22,file='inter.pbc-PV.dat')
open(21,file='pqfasiso.dat')
do l=1,ncadeias-1
   do m = 1, naneisolig
      do n = l+1, ncadeias
         do o = 1, naneisolig-1
	    if (R(l,m,n,o) < limR) then

            dotp2 = abs(v_normal(l,m).DOT.vstret)
            dotp1 = abs(vv1_normal(n,o).DOT.vstret)
            write(21,*) dotp2,dotp1

            if ((dotp1 <= 0.5) .and. (dotp2 <= 0.5)) then


            !funcao de correlacao orientacional = (|rij.ni| + 1)*(|rij.nj| + 1)/4
            anelanel = vrr(l,m,n,o)/R(l,m,n,o)
            dotp1 = abs(anelanel.DOT.vv1_normal(n,o))+1
            dotp2 = abs(anelanel.DOT.  v_normal(l,m))+1
            fco1 = dotp1*dotp2/4.0d0

            anelanel = vrr(l,m,n,o)/R(l,m,n,o)
            dotp1 = abs(anelanel.DOT.vv2_normal(n,o))+1
            dotp2 = abs(anelanel.DOT.  v_normal(l,m))+1
            fco2 = dotp1*dotp2/4.0d0

            ! correlacao com eixos fenil-vinil (avaliar paralelismo)
            dotp = abs(v_eixo(l,m).DOT.v_eixo_pvm(n,o))

            write(22,*) l,m,n,o,R(l,m,n,o),fco1,fco2,dotp
	    end if
            end if
         end do
      end do
   end do
end do
close(22)
close(21)

!TERMOS INTRA
open(22,file='intra.pbc-PV.dat')
do l=1,ncadeias
   n = l
   do m = 1, naneisolig
      do o = 1, naneisolig-1
	    if (R(l,m,n,o) < limR) then

            !funcao de correlacao orientacional = (|rij.ni| + 1)*(|rij.nj| + 1)/4
            anelanel = vrr(l,m,n,o)/R(l,m,n,o)
            dotp1 = abs(anelanel.DOT.vv1_normal(n,o))+1
            dotp2 = abs(anelanel.DOT.  v_normal(l,m))+1
            fco1 = dotp1*dotp2/4.0d0

            anelanel = vrr(l,m,n,o)/R(l,m,n,o)
            dotp1 = abs(anelanel.DOT.vv2_normal(n,o))+1
            dotp2 = abs(anelanel.DOT.  v_normal(l,m))+1
            fco2 = dotp1*dotp2/4.0d0

            write(22,*) l,abs(m-o),R(l,m,n,o),fco1,fco2
            end if
      end do
   end do
end do
close(22)

print *, 'ARM-INTER-PV: concluido!'

END PROGRAM ARMinterPV

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

