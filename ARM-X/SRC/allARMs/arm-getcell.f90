! ==========================================================
! Programa ARM - [A]mazonas [R]odrigo [M]elissa
! Analise estatistica da conformacao em cadeias oligomericas
! ==========================================================

PROGRAM ARMinterPTC60

USE vectors
IMPLICIT NONE

integer::i,j,k,l,n,m,o,ii
integer::atm1,atm2,atm5,atmS,s

integer::nat,natol,natolT,ncadeiaspt,ncadeiasc60
integer::natmon,naneisolig  

integer, PARAMETER::natCARBmon = 8

INTEGER,dimension(:),allocatable:: Cn,Rn

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
real(8)::Tx,Ty,Tz,limBall
integer::nC60

integer,dimension(:,:),allocatable::num 
integer,dimension(:,:),allocatable::atom
character(len=10):: chainfile   
real(8),dimension(:,:),allocatable::cx,cy,cz 
real(8),dimension(:,:),allocatable::c6x,c6y,c6z 

real(8):: alpha, beta, gamma, theta1,theta2, phi, tau, lambda

real(8) A(3,3),B(3,3)         !matriz de vetores da celula unitaria e inversa
integer(8) IPBX,IPBY,IPBZ

real(8),dimension(:,:),    allocatable:: sx,sy,sz          ! ...para posicoes
real(8):: ssx,ssy,ssz                                      ! ...para distancias
real(8):: xxm,yym,zzm   

real(8),dimension(:,:),    allocatable:: svx,svy,svz          ! ...para posicoes
real(8),dimension(:),      allocatable:: fsvx,fsvy,fsvz       ! ...para posicoes
real(8),dimension(:),      allocatable:: fxv,fyv,fzv          ! ...para posicoes

real(8),dimension(:,:),    allocatable:: sxT,syT,szT          ! ...para posicoes
real(8),dimension(:,:),    allocatable:: xxT,yyT,zzT          ! ...para posicoes

real(8),dimension(:,:),    allocatable:: vnnx,vnny,vnnz

double precision:: dotp,dotp1,dotp2,fco1,fco2

TYPE (vector) anelanel
TYPE (vector) v_comp, v_plano
real(8):: Lplanar, Lp
real(8):: Estlm, Estno

real(8) limphi, limtheta,limRC,limRT,limhpp,limtau,limdelta,limlambda
integer flgphi,flgtheta1,flgtheta2,flgR,flghpp,flgtau,flgdelta,flglambda

integer calphaT, dC1T,dC2T,dC3T 
integer calphaC, dC1C,dC2C,dC3C 
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
print *, '--> C60 + P3HT:: getting representative cell'
print *, '-------------------------------------------------------'

print*, '>>processando p3ht'

print *, 'Reading info.stat'
open(10,file='info.stat')
read(10,*) natolT          !numero de atomos de um oligomero
read(10,*) ncadeiaspt      !numero de cadeias 
read(10,*) naneisolig      !numero de aneis por oligomero
read(10,*) limRT           ! raio de corte para avaliar estatisticas sinter-cadeia
read(10,*) calphaT        ! indice do calpha da estrutura
read(10,*) dC1T            ! dC1: do calpha para o calpha oposto (define centro do anel)
read(10,*) s              ! sulfur atom
read(10,*) dC2T            ! passo ateh o prox. primeiro calpha
close(10)

!fixando
calphaT = 1
dC1T    = 3
s      = 5 
dC2T    = 11
print *, calphaT, dC1T,s,dC2T
print *, ncadeiaspt,natolT

allocate(cx(ncadeiaspt,natolT))
allocate(cy(ncadeiaspt,natolT))
allocate(cz(ncadeiaspt,natolT))
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
   do j=1,natolT
      read(16,*) nimp1, cx(i,j), cy(i,j), cz(i,j)
   end do
   close(16)
end do
close(15)


print*, '=========================================================================='
print*, '>>processando c60'

print *, 'Reading infoc60.stat'
open(10,file='infoc60.stat')
read(10,*) natol          !numero de atomos de um oligomero
read(10,*) ncadeiasc60    !numero de cadeias 
read(10,*) limRC           ! raio de corte para avaliar estatisticas sinter-cadeia
read(10,*) calphaC         ! indice do calpha da estrutura
read(10,*) dC1C            ! dC1: do calpha para o calpha oposto (define centro do anel)
close(10)


print *, 'Reading infoGet.stat'
open(10,file='infoGet.stat')
read(10,*) limBall          !raio de corte da bola
read(10,*) nC60             !C60 escolhido 


!fixando pra garantir
calphaC = 1
dC1C    = 56
print *, calphaC, calphaC+dc1C
print*, natol,ncadeiasc60

allocate(c6x(ncadeiasc60,natol))
allocate(c6y(ncadeiasc60,natol))
allocate(c6z(ncadeiasc60,natol))

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
      read(16,*) nimp1, c6x(i,j), c6y(i,j), c6z(i,j)
   end do
   close(16)
end do
close(15)

!Translation of c60-1 to the center of box
m=1
do l = 1, ncadeiasc60
   atm2 = calphaC
   atm5 = calphaC+dC1C
      xvm(l,m)=( c6x(l,atm2)+c6x(l,atm5) )/2   !\vec centro
      yvm(l,m)=( c6y(l,atm2)+c6y(l,atm5) )/2
      zvm(l,m)=( c6z(l,atm2)+c6z(l,atm5) )/2
end do

Tx = A(1,1)/2.0d0 + A(1,2)/2.0d0 + A(1,3)/2.0d0  - xvm(nC60,1)
Ty = A(2,1)/2.0d0 + A(2,2)/2.0d0 + A(2,3)/2.0d0  - yvm(nC60,1)
Tz = A(3,1)/2.0d0 + A(3,2)/2.0d0 + A(3,3)/2.0d0  - zvm(nC60,1)

print *, Tx,Ty,Tz

do i=1, ncadeiasc60
do j=1, natol
   write(18,*) 'C ', c6x(i,j),c6y(i,j),c6z(i,j)
   c6x(i,j) = c6x(i,j)+Tx
   c6y(i,j) = c6y(i,j)+Ty
   c6z(i,j) = c6z(i,j)+Tz
   write(17,*) 'C ', c6x(i,j),c6y(i,j),c6z(i,j)
end do
end do

do i=1,ncadeiaspt
do j=1, natolT
   cx(i,j) = cx(i,j)+Tx
   cy(i,j) = cy(i,j)+Ty
   cz(i,j) = cz(i,j)+Tz
end do
end do
!end of translation


allocate(v_normal(ncadeiaspt,naneisolig))
allocate(vnn(ncadeiaspt,naneisolig))
allocate(xm(ncadeiaspt,naneisolig))
allocate(ym(ncadeiaspt,naneisolig))
allocate(zm(ncadeiaspt,naneisolig))
xm = 0.0
ym = 0.0
zm = 0.0

do l = 1, ncadeiaspt
   atm2 = calphaT
   atm5 = calphaT+dC1T
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
	
      atm2 =  atm2   +  dC2T
      atm5 =  atm5   +  dC2T
      atmS =  atmS   +  dC2T         
   end do
end do

m=1
do l = 1, ncadeiasc60
   atm2 = calphaC
   atm5 = calphaC+dC1C
      xvm(l,m)=( c6x(l,atm2)+c6x(l,atm5) )/2   !\vec centro
      yvm(l,m)=( c6y(l,atm2)+c6y(l,atm5) )/2
      zvm(l,m)=( c6z(l,atm2)+c6z(l,atm5) )/2
end do

!PBC
print *, 'inserting periodic boundary condictions....'

allocate(sx(ncadeiaspt,naneisolig))
allocate(sy(ncadeiaspt,naneisolig))
allocate(sz(ncadeiaspt,naneisolig))
allocate(svx(ncadeiasc60,1))
allocate(svy(ncadeiasc60,1))
allocate(svz(ncadeiasc60,1))
sx  = 0.0
sy  = 0.0
sz  = 0.0
svx = 0.0
svy = 0.0
svz = 0.0

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

   end do 
end do

!C60
m=1
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
end do


allocate(vrr(ncadeiasc60,1,ncadeiaspt,naneisolig))
allocate(  R(ncadeiasc60,1,ncadeiaspt,naneisolig))
R   = 0.0d0


!folding C60's

allocate(fsvx(60))
allocate(fsvy(60))
allocate(fsvz(60))
allocate(fxv(60))
allocate(fyv(60))
allocate(fzv(60))

l=nC60
!write(16,*) '60 '
!write(16,*)
 do j = 1, 60

   fsvx(j) = B(1,1)*c6x(l,j)+B(1,2)*c6y(l,j)+B(1,3)*c6z(l,j)
   fsvy(j) = B(2,1)*c6x(l,j)+B(2,2)*c6y(l,j)+B(2,3)*c6z(l,j)
   fsvz(j) = B(3,1)*c6x(l,j)+B(3,2)*c6y(l,j)+B(3,3)*c6z(l,j)

   IPBX = dint(fsvx(j)+1.d0)-1
   IPBY = dint(fsvy(j)+1.d0)-1
   IPBZ = dint(fsvz(j)+1.d0)-1

   fsvx(j)= fsvx(j) - IPBX
   fsvy(j)= fsvy(j) - IPBY
   fsvz(j)= fsvz(j) - IPBZ

   fxv(j) = A(1,1)*fsvx(j)+A(1,2)*fsvy(j)+A(1,3)*fsvz(j)
   fyv(j) = A(2,1)*fsvx(j)+A(2,2)*fsvy(j)+A(2,3)*fsvz(j)
   fzv(j) = A(3,1)*fsvx(j)+A(3,2)*fsvy(j)+A(3,3)*fsvz(j)

   write(16,*) 'C ', fxv(j), fyv(j), fzv(j)
!   write(16,*) 'C ', c6x(l,j), c6y(l,j), c6z(l,j)
 end do


!folding p3ht's

allocate(sxT(ncadeiaspt,natolT))
allocate(syT(ncadeiaspt,natolT))
allocate(szT(ncadeiaspt,natolT))
allocate(xxT(ncadeiaspt,natolT))
allocate(yyT(ncadeiaspt,natolT))
allocate(zzT(ncadeiaspt,natolT))

do l = 1, ncadeiaspt
   do m = 1, natolT

   sxT(l,m) = B(1,1)*cx(l,m)+B(1,2)*cy(l,m)+B(1,3)*cz(l,m)
   syT(l,m) = B(2,1)*cx(l,m)+B(2,2)*cy(l,m)+B(2,3)*cz(l,m)
   szT(l,m) = B(3,1)*cx(l,m)+B(3,2)*cy(l,m)+B(3,3)*cz(l,m)

   IPBX = dint(sxT(l,m)+1.d0)-1
   IPBY = dint(syT(l,m)+1.d0)-1
   IPBZ = dint(szT(l,m)+1.d0)-1

   sxT(l,m)= sxT(l,m) - IPBX
   syT(l,m)= syT(l,m) - IPBY
   szT(l,m)= szT(l,m) - IPBZ

   xxT(l,m) = A(1,1)*sxT(l,m)+A(1,2)*syT(l,m)+A(1,3)*szT(l,m)
   yyT(l,m) = A(2,1)*sxT(l,m)+A(2,2)*syT(l,m)+A(2,3)*szT(l,m)
   zzT(l,m) = A(3,1)*sxT(l,m)+A(3,2)*syT(l,m)+A(3,3)*szT(l,m)

   write (99,*) 'C', xxT(l,m),yyT(l,m),zzT(l,m)
   end do 
end do
allocate(Cn(natolT))
allocate(Rn(natolT))

  l=nC60
  m=1
  k = 1
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

   if (R(l,m,n,o)<= limBall) then
    write(11,'(3f12.5)') xxm,yym,zzm
    Cn(k)=n
    Rn(k)=o
    k=k+1
   endif
   end do
  end do 

do i=1,k-1
	print *, Cn(i),Rn(i)
	do ii=1,5
	n = (Rn(i)-1)*11 + ii
	print*, n
	write(16,*) 'C', xxT(Cn(i),n),yyT(Cn(i),n),zzT(Cn(i),n)
	end do
end do

print*, 'calculando correlacao c60 x p3ht...... gerando lista'

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

   if (R(l,m,n,o)<= limBall) then
    print*, 'Rcr:',ncadeiaspt+l,n,o,R(l,m,n,o)
   endif

  end do
  end do 
  print*, 'Rcr:'
end do

deallocate(sx)
deallocate(sy)
deallocate(sz)
deallocate(svx)
deallocate(svy)
deallocate(svz)
	
print *, 'ARM-GETCELL: concluido!'

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

