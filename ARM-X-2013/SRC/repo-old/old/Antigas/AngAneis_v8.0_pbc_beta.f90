! ==========================================================
! Programa ARM - nome provisorio
! Analise estatistica da conformacao em cadeias oligomericas
!
! Versao 8 / set.2010 (under dev.)
! Versao 1  : Jarlesson Amazonas
! Versao 2-7: Rodrigo Ramos e Melissa F. S. Pinto 
! Versao 8  : Rodrigo Ramos e Jarlesson Amazonas 
! ==========================================================

PROGRAM ang_anel

USE vectors
IMPLICIT NONE

integer::i,j,k,l,n,m,o
integer::atm1,atm2,atm5
integer::count

integer::nat,natol,ncadeias 
integer::natmon,naneisolig  

integer::uni,unit,unid,unida

integer, PARAMETER::natCARBmon = 8

REAL(8),dimension(:,:),allocatable::xm, ym, zm 
REAL(8),dimension(:,:,:,:),allocatable:: R
real(8)::xmv4, ymv4, zmv4
real(8)::xmv6, ymv6, zmv6 

TYPE (vector), DIMENSION (:,:),allocatable:: v_normal, vnn

TYPE (vector), DIMENSION (:,:,:,:),allocatable:: vrr

TYPE (vector), DIMENSION (:,:),allocatable::vdesl

TYPE (vector)::v1,v2

character(len=25)::nimp1, nimp2, nimp3, nimp4, nimp5, nimp6, nimp10
character(len=45)::arquivo,arquivo1,arquivo2

real(8), PARAMETER::PI= 3.141592654
REAL(8), DIMENSION (3) :: array_out


integer,dimension(:,:),allocatable::num 
integer,dimension(:,:),allocatable::atom
character(len=10):: chainfile   
real(8),dimension(:,:),allocatable::cx,cy,cz 

REAL(8):: DELTAxm, DELTAym, DELTAzm

real(8):: alpha, beta
real(8):: gamma

! tratamento das condicoes periodicas

real(8) A(3,3),B(3,3)         !matriz de vetores da celula unitaria e inversa
integer(8) IPBX,IPBY,IPBZ
real(8) D11,D22,D33,D21,D32
real(8) D13,D31,D12,D23
real(8) C

real(8),dimension(:,:),    allocatable:: sx,sy,sz          ! ...para posicoes
real(8):: ssx,ssy,ssz                                      ! ...para distancias
real(8):: xxm,yym,zzm   

real(8),dimension(:,:),    allocatable:: vnnx,vnny,vnnz

double precision:: dotp

TYPE (vector)anelanel
real(8) Lplanar


print *, '============================='
print *, 'Statistical Analisis'
print *, 'ARM pack- [A]mazonas, [R]amos, [M]elissa' 
print *, 'Version 8 / set.2010'
print *, '============================='
print *, ''

print *, 'Reading info.stat'
open(10,file='info.stat')

read(10,*) natol          !numero de atomos de um oligomero  
read(10,*) ncadeias       !numero de cadeias 
read(10,*) naneisolig     !numero de aneis por oligomero
read(10,*) Lplanar        !numero de aneis por oligomero

close(10)

print *,'\/\ ok'

open(12,file='MAT.in') 

!lendo matriz de vetores da celula unitaria (A)

read(12,*) A(1,1),A(2,1),A(3,1)
read(12,*) A(1,2),A(2,2),A(3,2)
read(12,*) A(1,3),A(2,3),A(3,3)

close(12)

! invertendo matriz vetores da celula unitaria (B)

D11=A(2,2)*A(3,3)-A(2,3)*A(3,2)
D22=A(3,3)*A(1,1)-A(3,1)*A(1,3)
D33=A(1,1)*A(2,2)-A(1,2)*A(2,1)
D21=A(2,3)*A(3,1)-A(2,1)*A(3,3)
D32=A(3,1)*A(1,2)-A(3,2)*A(1,1)
D13=A(1,2)*A(2,3)-A(1,3)*A(2,2)
D31=A(2,1)*A(3,2)-A(3,1)*A(2,2)
D12=A(3,2)*A(1,3)-A(1,2)*A(3,3)
D23=A(1,3)*A(2,1)-A(2,3)*A(1,1)
C=A(1,1)*D11+A(1,2)*D21+A(1,3)*D31

B(1,1)=D11/C
B(2,2)=D22/C
B(3,3)=D33/C
B(1,2)=D12/C
B(2,3)=D23/C
B(3,1)=D31/C
B(2,1)=D21/C
B(3,2)=D32/C
B(1,3)=D13/C

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

allocate(atom(ncadeias,natol))   !tipo de atomo
allocate(num(ncadeias,natol))    !numero do atomo no filme

num  = 0.0

!READING

nat = 0
open(15,file='chains.in')   ! colocar arquivos do dir ./chains/
                            ! na lista 'chains.in'


!WARNING: isso pode ser arbitrario jah que o pdb pode ter campos extras
!         rode cautelosamente os scripts de separacao

do i=1,ncadeias

   read(15,'(a7)') chainfile   ! formato NNN.pdb
   open(16, file='./chainsout/'//chainfile)

   do j=1,natol
      read(16,*) nimp1, atom(i,j), nimp2, nimp2, cx(i,j), cy(i,j), cz(i,j), nimp4, nimp5, nimp6
      nat = nat + 1
      num(i,j) = nat
   end do

   close(16)

end do

close(15)

allocate(v_normal(ncadeias,naneisolig))
allocate(vnn(ncadeias,naneisolig))
allocate(vdesl(ncadeias,naneisolig))

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

print *, 'Colecting Structural information'
print *, ''


print *, '-------------------------------------------------------'
print *, 'Notes:'
print *, '--> PPV phenil capped structure'
print *, '--> Each ring followed by vinil structures in pdb file' 
print *, '--> 1st ring *without* its vinil structure'
print *, '--> hydrogen atoms excluded' 
print *, '-------------------------------------------------------'

!calculating centers and normal vectors
!checado (vpython) / ok.
open(18,file='all.vec.rg.dat')                  !coordenadas dos centros
do l = 1, ncadeias
   atm2 = 2
   atm5 = 5
   do m = 1, naneisolig

      xm(l,m)=( cx(l,atm2)+cx(l,atm5) )/2   !\vec centro
      ym(l,m)=( cy(l,atm2)+cy(l,atm5) )/2
      zm(l,m)=( cz(l,atm2)+cz(l,atm5) )/2
      
      xmv4=xm(l,m)-cx(l,atm2+2)        !\vec centro para 4 (at)
      ymv4=ym(l,m)-cy(l,atm2+2) 
      zmv4=zm(l,m)-cz(l,atm2+2)

      xmv6=xm(l,m)-cx(l,atm5+1)        !\vec centro para 6 (at)
      ymv6=ym(l,m)-cy(l,atm5+1)              
      zmv6=zm(l,m)-cz(l,atm5+1)
   
      v1 = (/ xmv4, ymv4, zmv4 /)
      v2 = (/ xmv6, ymv6, zmv6 /)

      v_normal(l,m) = ( v1*v2 ) 
      v_normal(l,m) = v_normal(l,m)/sqrt( v_normal(l,m) .DOT. v_normal(l,m) )
       
      write(18,'(2i3,9f10.5)')l,m, xm(l,m), ym(l,m), zm(l,m),v_normal(l,m)
!
!WARNING: Isso (prox. linhas) pode ser arbitrario 
!         dependendo da numeracao da cadeia.
!         A estrutura phenil capped tem esse problema de nao ter o 
!         primeiro monomero completo.// VER PRINT* ANTERIOR.
      
      if (m .eq. 1) then 
         atm2 = atm2 + 6
         atm5 = atm5 + 6         
      else         
         atm2 =  atm2   +  natCARBmon
         atm5 =  atm5   +  natCARBmon         
      end if
   end do
end do
close(18)

open(18,file='intra.vec.1nb.dat')

do l=1, ncadeias
   do m = 1, naneisolig-1
      
      v1 = (/ xm(l,m)  , ym(l,m)  , zm(l,m)   /)
      v2 = (/ xm(l,m+1), ym(l,m+1), zm(l,m+1) /)
      
      vnn(l,m) = v2 - v1 
      
      write(18,'(2i3, 3f10.5)') l,m, vnn(l,m)
      
   end do
   write(18,'(2i3, 3f10.5)') l,m, 0.0,0.0,0.0
end do
close(18)

!=======================================
! Calculating statistical distributions
!=======================================

!-------------------
!A. INTRACADEIAS
!-------------------
print *, 'intrachain DATA'

!A.1. Distancia entre centros de aneis (todos com todos)

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

!S. Todos com todos
open(18,file='intra.dist-rg.dat')  

do l=1, ncadeias
   do m = 1, naneisolig-1    
      do o = 1, naneisolig

         write(18,*) l,m,o, R(l,m,l,o)
         
      end do
   end do
end do

close(18)

!S. Primeiros vizinhos intracadeia 
!open(18,file='intra.dist1st.dat') 
!
!do l=1, ncadeias
!   do m = 1, naneisolig-1    
!         write(18,*) m,R(l,m,l,m+1)
!   end do
!end do
!
!close(18)

!A.2. Calculando linearidade end-to-end (relativa ao tam. planar)
!(*Kuhn factor)

open(18, file='intra.lin-ete.dat')

do l = 1, ncadeias
   write(18,*) l, R(l,1,l,naneisolig)/Lplanar
end do

close(18)

!A.4. Vetores unitarios entre aneis. -> anisotropia.

open(18,file='intra.anis.dat') 

do l=1, ncadeias
   do m = 1, naneisolig-1          

      vnn(l,m) = vnn(l,m)/sqrt( vnn(l,m) .DOT. vnn(l,m) )
      write(18,*) l,m, vnn(l,m)      

   end do
end do

close(18)

!A.5. Angulo entre dois segmentos de cadeia adjacentes
open(18,file='intra.ang-sg.dat')

! jah conferido (vpython) / ok.
do l=1, ncadeias
   do m = 1, naneisolig-2

      gamma = acos( (vnn(l,m) .DOT. vnn(l,m+1)) )*(180.0/PI) 
      write(18,'(4i5, 1f10.3)') l, m, l, m+1, gamma      

   end do
end do

close(18)

!A.6. Angulos entre aneis 
open(18,file='intra.ang-rg.dat')      

do l=1,ncadeias
   count = 0                              !mapping conjugation breaking
   do m = 1, naneisolig-1
 
!alpha: ang entre normais dos aneis
!beta : ang. entre normal do primeiro anel com vetor que liga ao segundo

      dotp = v_normal(l,m) .DOT. v_normal(l,m+1) 

      dotp = dotp*0.9999999d0 ! gambiarra para superar instabilidade do acos... 
                              ! ler http://www.megasolutions.net/fortran/Epsilon,-Precision-or-Tiny_-64109.aspx

      alpha = dacos( dotp )*(180.d0/PI)

      dotp = v_normal(l,m) .DOT. vnn(l,m)        
      beta =  dacos( dotp ) *(180.d0/PI)

      if ((abs(alpha) >=  60.0) .or. (abs(beta) <= 45.0)) count = count + 1 ! quebra de conjugacao  

      write(18,'(3i5, 2f10.3, 1i5)') l, m, m+1, alpha, beta, count

      
   end do
end do

close(18)

print*, '\/\/\/ok'

! talvez seja bom desalocar a memoria que nao for mais ser usada.
! como tratar as imagens periodicas considerando cadeias distintas?

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

sx(l,m) = -1.0
sy(l,m) = -1.0
sz(l,m) = -1.0

! estah com essa gambiarra para colocar todo mundo na celula 
! [0,a1]x[0,a2]x[0,a3] // fazer a coisa mais bonitinha com algebra melhor
! essa pbc funcionava antes (tbmd) porque nunca os atomos se afastavam
! para alem das celulas periodicas vizinhas. Essa autoconsistencia
! resolve a coisa, reaplicando a transformacao de periodicidade.

!print *, l,m

do while ((sx(l,m) < 0.0).or.(sy(l,m)<0.0).or.(sz(l,m)<0.0))

sx(l,m) = B(1,1)*xm(l,m)+B(1,2)*ym(l,m)+B(1,3)*zm(l,m)
sy(l,m) = B(2,1)*xm(l,m)+B(2,2)*ym(l,m)+B(2,3)*zm(l,m)
sz(l,m) = B(3,1)*xm(l,m)+B(3,2)*ym(l,m)+B(3,3)*zm(l,m)

IPBX = dint(sx(l,m)+1.d0)-1
IPBY = dint(sy(l,m)+1.d0)-1
IPBZ = dint(sz(l,m)+1.d0)-1

!print *, xm(l,m),ym(l,m),zm(l,m)
!print *, sx(l,m),sy(l,m),sz(l,m)
!print *, IPBX,IPBY,IPBZ

sx(l,m)= sx(l,m) - IPBX
sy(l,m)= sy(l,m) - IPBY
sz(l,m)= sz(l,m) - IPBZ

xm(l,m) = A(1,1)*sx(l,m)+A(1,2)*sy(l,m)+A(1,3)*sz(l,m)
ym(l,m) = A(2,1)*sx(l,m)+A(2,2)*sy(l,m)+A(2,3)*sz(l,m)
zm(l,m) = A(3,1)*sx(l,m)+A(3,2)*sy(l,m)+A(3,3)*sz(l,m)

end do

!print *, xm(l,m),ym(l,m),zm(l,m)
!print *, '---'

write(18,'(2i3,9f10.5)')l,m, xm(l,m), ym(l,m), zm(l,m),v_normal(l,m)

   end do 
end do

close(18)

allocate(vrr(ncadeias,naneisolig,ncadeias,naneisolig))

!B.1. DISTANCIA ENTRE OS CENTROS DOS ANEIS

open(18,file='inter.dist-rg.dat')      

! isso aqui eh a parte mais lenta... 
! talvez convenha partir para lista de vizinhos... (wish list versao 9)

print *, 'entrando na etapa mais lenta... / inter.dist-rg'

do l=1,ncadeias-1
   do m = 1, naneisolig              
      do n = l+1, ncadeias
         do o = 1, naneisolig
                     
ssx = sx(l,m) - sx(n,o)
ssx = ssx  - (dint((2.d0*ssx +3.d0)/2.d0)-1.d0)
ssy = sy(l,m) - sy(n,o)
ssy = ssy  - (dint((2.d0*ssy +3.d0)/2.d0)-1.d0)
ssz = sz(l,m) - sz(n,o)
ssz = ssz  - (dint((2.d0*ssz +3.d0)/2.d0)-1.d0)

xxm  = A(1,1)*ssx +A(1,2)*ssy +A(1,3)*ssz 
yym  = A(2,1)*ssx +A(2,2)*ssy +A(2,3)*ssz 
zzm  = A(3,1)*ssx +A(3,2)*ssy +A(3,3)*ssz 

vrr(l,m,n,o) = (/ xxm, yym, zzm /)    
vrr(n,o,l,m) = vrr(l,m,n,o) - 2.0d0*vrr(l,m,n,o)

R(l,m,n,o)= sqrt( vrr(l,m,n,o).DOT.vrr(l,m,n,o) )
R(n,o,l,m) = R(l,m,n,o)

write(18,*) l,m,n,o, R(l,m,n,o)

         end do
      end do      
   end do
end do

close(18)

print *, 'concluido.'

deallocate(sx)
deallocate(sy)
deallocate(sz)

!B.2. Angulos inter_cadeia dentro de um raio de corte
!Estabelecendo relacao com orientacao relativa dos aneis.

! passar o raio de corte para leitura externa no info.stat

open(18,file='inter.ang-rgsg.dat')      !coordenadas

do l=1,ncadeias-1
   do m = 1, naneisolig              
      do n = l+1, ncadeias               ! corta as intra-cadeias       
         do o = 1, naneisolig  
 
            if ( R(l,m,n,o).le.(4.0) ) then
               
               anelanel = vrr(l,m,n,o)/sqrt( R(l,m,n,o))
               
               alpha  = (v_normal(l,m) .DOT. v_normal(n,o))
               beta = (v_normal(l,m) .DOT. anelanel)*(180.0/PI) 
               
               write(18,'(4i5, 2f10.3)') l,m,n,o, alpha, beta
            end if

         end do
      end do
   end do
end do

close(18)
print *, '\/\/\/\ ok'

print *, ''
print *, 'Successful termination. ;^)'
print *, ''


END PROGRAM ANG_ANEL


MODULE vectors

IMPLICIT NONE

!Declare vector data type:
TYPE :: vector
REAL(8) :: x
REAL(8) :: y
REAL(8) :: z
END TYPE

!Declare interface operators
INTERFACE ASSIGNMENT (=)
        MODULE PROCEDURE array_to_vector
        MODULE PROCEDURE vector_to_array
END INTERFACE

INTERFACE OPERATOR (+)
        MODULE PROCEDURE vector_add
END INTERFACE

INTERFACE OPERATOR (-)
        MODULE PROCEDURE vector_subtract
END INTERFACE

INTERFACE OPERATOR (*)
        MODULE PROCEDURE vector_times_real
        MODULE PROCEDURE real_times_vector
        MODULE PROCEDURE vector_times_int
        MODULE PROCEDURE int_times_vector
        MODULE PROCEDURE cross_product
END INTERFACE

INTERFACE OPERATOR (/)
        MODULE PROCEDURE vector_div_real
        MODULE PROCEDURE vector_div_int
END INTERFACE

INTERFACE OPERATOR (.DOT.)
        MODULE PROCEDURE dot_product
END INTERFACE

!       Now define the implementing functions.

CONTAINS
        SUBROUTINE array_to_vector (vec_result, array)
        TYPE (vector), INTENT (OUT) :: vec_result
        REAL(8), DIMENSION (3), INTENT (IN) :: array
        vec_result%x = array (1)
        vec_result%y = array (2)
        vec_result%z = array (3)
        END SUBROUTINE array_to_vector

        SUBROUTINE vector_to_array(array_result, vec_1)
        REAL(8), DIMENSION (3), INTENT (OUT) :: array_result
        TYPE (vector), INTENT (IN) :: vec_1
        array_result (1) = vec_1%x
        array_result (2) = vec_1%y
        array_result (3) = vec_1%z
        END SUBROUTINE vector_to_array

        FUNCTION vector_add (vec_1, vec_2)
        TYPE (vector) :: vector_add
        TYPE (vector), INTENT (IN) :: vec_1, vec_2
        vector_add%x = vec_1%x + vec_2%x
        vector_add%y = vec_1%y + vec_2%y
        vector_add%z = vec_1%z + vec_2%z
        END  FUNCTION vector_add

        FUNCTION vector_subtract (vec_1, vec_2)
        TYPE (vector) :: vector_subtract
        TYPE (vector), INTENT (IN) :: vec_1, vec_2
        vector_subtract%x = vec_1%x - vec_2%x
        vector_subtract%y = vec_1%y - vec_2%y
        vector_subtract%z = vec_1%z - vec_2%z
        END  FUNCTION vector_subtract

        FUNCTION vector_times_real (vec_1, real_2)
        TYPE (vector) ::  vector_times_real
        TYPE (vector), INTENT (IN) :: vec_1
        REAL(8), INTENT (IN) :: real_2
        vector_times_real%x = vec_1%x * real_2
        vector_times_real%y = vec_1%y * real_2
        vector_times_real%z = vec_1%z * real_2
        END  FUNCTION  vector_times_real

        FUNCTION real_times_vector (real_1, vec_2)
        TYPE (vector) ::  real_times_vector
        REAL(8), INTENT (IN) :: real_1
        TYPE (vector), INTENT (IN) :: vec_2
        real_times_vector%x = real_1 * vec_2%x
        real_times_vector%y = real_1 * vec_2%y
        real_times_vector%z = real_1 * vec_2%z
        END  FUNCTION real_times_vector

        FUNCTION vector_times_int (vec_1, int_2)
        TYPE (vector) :: vector_times_int
        TYPE (vector), INTENT (IN) :: vec_1
        INTEGER, INTENT (IN) :: int_2
        vector_times_int%x = vec_1%x * REAL(int_2)
        vector_times_int%y = vec_1%y * REAL(int_2)
        vector_times_int%z = vec_1%z * REAL(int_2)
        END  FUNCTION vector_times_int

        FUNCTION int_times_vector (int_1, vec_2)
        TYPE (vector) :: int_times_vector
        INTEGER, INTENT (IN) :: int_1
        TYPE (vector), INTENT (IN) :: vec_2
        int_times_vector%x = REAL(int_1) * vec_2%x
        int_times_vector%y = REAL(int_1) * vec_2%y
        int_times_vector%z = REAL(int_1) * vec_2%z
        END  FUNCTION  int_times_vector

        FUNCTION vector_div_real(vec_1, real_2)
        TYPE (vector) ::  vector_div_real
        TYPE (vector), INTENT(IN) :: vec_1
        REAL(8), INTENT(IN) :: real_2
        vector_div_real%x = vec_1%x / real_2
        vector_div_real%y = vec_1%y / real_2
        vector_div_real%z = vec_1%z / real_2
        END  FUNCTION  vector_div_real

        FUNCTION vector_div_int(vec_1, int_2)
        TYPE (vector) :: vector_div_int
        TYPE (vector), INTENT(IN) :: vec_1
        INTEGER, INTENT(IN) :: int_2
        vector_div_int%x = vec_1%x / REAL(int_2)
        vector_div_int%y = vec_1%y / REAL(int_2)
        vector_div_int%z = vec_1%z / REAL(int_2)
        END FUNCTION vector_div_int

        FUNCTION dot_product (vec_1, vec_2)
        REAL(8) :: dot_product
        TYPE (vector), INTENT (IN) :: vec_1, vec_2
        dot_product = vec_1%x*vec_2%x + vec_1%y*vec_2%y + vec_1%z*vec_2%z
        END FUNCTION dot_product

        FUNCTION cross_product (vec_1, vec_2)
                TYPE (vector) :: cross_product
                TYPE (vector), INTENT (IN) :: vec_1, vec_2
                cross_product%x = vec_1%y*vec_2%z - vec_1%z*vec_2%y
                cross_product%y = vec_1%z*vec_2%x - vec_1%x*vec_2%z
                cross_product%z = vec_1%x*vec_2%y - vec_1%y*vec_2%x
        END FUNCTION cross_product

        END MODULE vectors
