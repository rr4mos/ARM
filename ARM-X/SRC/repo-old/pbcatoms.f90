! ==========================================================
! Programa ARM - [A]mazonas [R]odrigo [M]elissa
! Analise estatistica da conformacao em cadeias oligomericas
!
! Versao 1  : Jarlesson Amazonas
! Versao 2-7: Rodrigo Ramos e Melissa F. S. Pinto 
! Versao 8  : Rodrigo Ramos e Jarlesson Amazonas 
! Versao 9  : Rodrigo Ramos 
! Versao 10 : Rodrigo Ramos 
! Versao 10.1 : J. Amazonas & R. Ramos.
! Versao 10.1 + vetores de conexao de sitios: R. Ramos
! ==========================================================
!
! notas:
! * http://www.megasolutions.net/fortran/Epsilon,-Precision-or-Tiny_-64109.aspx

PROGRAM ang_anel

IMPLICIT NONE

integer::i,j

integer::natol,ncadeias 
integer::naneisolig  

character(len=25)::nimp1
character(len=10):: chainfile   
real(8),dimension(:,:),allocatable::cx,cy,cz 

							  !tratamento das condicoes periodicas
real(8) A(3,3),B(3,3)         !matriz de vetores da celula unitaria e inversa
integer(8) IPBX,IPBY,IPBZ

real(8),dimension(:,:),allocatable:: sx,sy,sz              ! ...para posicoes

print *, '==========================================='
print *, 'PBC nos atomos'
print *, '==========================================='

print *, 'Reading info.stat'
open(10,file='info.stat')

 read(10,*) natol          !numero de atomos de um oligomero
 read(10,*) ncadeias       !numero de cadeias 
 read(10,*) naneisolig     !numero de aneis por oligomero
 close(10)

open(12,file='MAT.in') 

!lendo matriz de vetores da celula unitaria (A) - sim eh transposto assim mesmo =P

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

allocate(sx(ncadeias,natol))
allocate(sy(ncadeias,natol))
allocate(sz(ncadeias,natol))

cx = 0.0
cy = 0.0
cz = 0.0

sx = 0.0
sy = 0.0
sz = 0.0

open(15,file='chains.in')   ! colocar arquivos do dir ./chainsout/
                            ! na lista 'chains.in'
							! gerar o diretorio ./chainsoutpbc/
do i=1,ncadeias

read(15,'(a7)') chainfile   ! formato NNN.xyz
open(16, file='./chainsout/'//chainfile)
open(17, file='./chainsoutpbc/'//chainfile)

read(16,*)  nimp1
write(17,*) nimp1
write(17,*) 
print *, 'cadeia: ', i

do j=1,natol

    read(16,*) nimp1, cx(i,j), cy(i,j), cz(i,j)
	print *, 'atomo: ', j

	sx(i,j) = -1.0
	sy(i,j) = -1.0
	sz(i,j) = -1.0

	do while ((sx(i,j) < 0.0).or.(sy(i,j)<0.0).or.(sz(i,j)<0.0))

		sx(i,j) = B(1,1)*cx(i,j)+B(1,2)*cy(i,j)+B(1,3)*cz(i,j)
		sy(i,j) = B(2,1)*cx(i,j)+B(2,2)*cy(i,j)+B(2,3)*cz(i,j)
		sz(i,j) = B(3,1)*cx(i,j)+B(3,2)*cy(i,j)+B(3,3)*cz(i,j)

		IPBX = dint(sx(i,j)+1.d0)-1
		IPBY = dint(sy(i,j)+1.d0)-1
		IPBZ = dint(sz(i,j)+1.d0)-1

		sx(i,j)= sx(i,j) - IPBX
		sy(i,j)= sy(i,j) - IPBY
		sz(i,j)= sz(i,j) - IPBZ

		cx(i,j) = A(1,1)*sx(i,j)+A(1,2)*sy(i,j)+A(1,3)*sz(i,j)
		cy(i,j) = A(2,1)*sx(i,j)+A(2,2)*sy(i,j)+A(2,3)*sz(i,j)
		cz(i,j) = A(3,1)*sx(i,j)+A(3,2)*sy(i,j)+A(3,3)*sz(i,j)
		
		print *, 's',sx(i,j), sy(i,j), sz(i,j)
		print *, 'p',IPBX,IPBY,IPBZ
		print *, 'c',cx(i,j), cy(i,j), cz(i,j)
!		print *, IPBX,IPBY,IPBZ

!	stop
	end do
	write(17,'(a6, 3f15.5)') nimp1, cx(i,j), cy(i,j),cz(i,j)
	print *, 'atomo: ', j, 'ok!'

end do ! atomos
close(16)
close(17)
print *, 'cadeia: ', i, 'ok!'

end do ! cadeias

end program

! inversao da matriz - subrotina auxiliar

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


