program MINIMIZACION
implicit none
double precision, external:: chichinadze, beale, fibonacci, invent1,invent2                         
double precision, parameter:: eps=10.d-20, l=10.d-11                                  
integer, parameter:: n2=2
double precision x2(n2),x2p(n2),norma, grad2(n2), hess2(n2,n2),invhess2(n2,n2)
double precision aa, bb, landa, mu,a,b
external:: gradchichinadze, hesschichinadze,inversamatriz2x2, modulo         
integer i,ene,k

open(11,file='ResultadosPracticaAlgoritmos.txt')
write(11,'(A)') '-------------------------------------------------------------'
write(11,'(A)') 'Ejecucion del Algoritmo de Fibonacci con funciones de una variable:'
write(11,'(A)') '-------------------------------------------------------------'

!ALGORITMO DE FIBONACCI CON LA PRIMERA FUNCION
write(11,'(A)') 'La primera funcion para hacer el test es:'
write(11,'(A)') 'F(X)=sin(3*x/2)*exp(-2*sqrt(x))+log(x**2+3*x+0.1)-8*x*cos(x+3)'
i=0
a=30.d0    !Intervalo inicial
b=40.d0  
write(11,'(A,f10.7,a,f10.7,a)') 'Intervalo inicial [a,b]=[',a,',',b,']'
do
  if (fibonacci(i)>(b-a)/l) then
    ene=i
    goto 199   !Buscamos el N tal que Fn>(b-a)/L
  else
    i=i+1
  endif
enddo
199 continue 
aa=a; bb=b; landa=aa+fibonacci(ene-2)*(bb-a)/fibonacci(ene)
mu=aa+fibonacci(ene-1)*(bb-aa)/fibonacci(ene); k=2
22 if (invent1(landa)<=invent1(mu)) then
  bb=mu; mu=landa; landa=aa+fibonacci(ene-k-1)*(bb-aa)/fibonacci(ene-k+1)
else
  aa=landa; landa=mu; mu=aa+fibonacci(ene-k)*(bb-aa)/fibonacci(ene-k+1)     
endif
if (k<ene-1) then
  k=k+1; goto 22
else
  mu=landa+l
endif
if (invent1(landa)<=invent1(mu)) then
  bb=mu
else
  aa=landa
endif
write(11,'(a,f16.13,a,f16.13,a)') 'El valor minimo de la funcion esta restringido al intervalo [a´,b´]=[',aa,',',bb,']'
write(11,'(a,f16.13,a,f18.13)')'Tomando el punto medio del intervalo se tiene xmin=',(aa+bb)/2.d0,', f(xmin)=',invent1((aa+bb)/2.d0)
write(11,'(A,/)') '-------------------------------------------------------------'
!!!!

!!!FIBONACCI FUNCION 2
!ALGORITMO DE FIBONACCI
write(11,'(A)') 'La segunda funcion para hacer el test es:'
write(11,'(A)') 'F(X)=sin(x)+sin(10*x/3)+log(x)-0.84*x+3'
i=0
a=2.7d0   !Intervalo inicial
b=7.5d0
write(11,'(A,f9.7,a,f9.7,a)') 'Intervalo inicial [a,b]=[',a,',',b,']'  
do
  if (fibonacci(i)>(b-a)/l) then
    ene=i
    goto 1990
  else
    i=i+1
  endif
enddo
1990 continue 
aa=a; bb=b; landa=aa+fibonacci(ene-2)*(bb-a)/fibonacci(ene)
mu=aa+fibonacci(ene-1)*(bb-aa)/fibonacci(ene); k=2
220 if (invent2(landa)<=invent2(mu)) then
  bb=mu; mu=landa; landa=aa+fibonacci(ene-k-1)*(bb-aa)/fibonacci(ene-k+1)
else
  aa=landa; landa=mu; mu=aa+fibonacci(ene-k)*(bb-aa)/fibonacci(ene-k+1)     
endif
if (k<ene-1) then
  k=k+1; goto 220
else
  mu=landa+l
endif
if (invent2(landa)<=invent2(mu)) then
  bb=mu
else
  aa=landa
endif
write(11,'(a,f16.13,a,f16.13,a)') 'El valor minimo de la funcion esta restringido al intervalo [a´,b´]=[',aa,',',bb,']'
write(11,'(a,f16.13,a,f18.13)')'Tomando el punto medio del intervalo se tiene xmin=',(aa+bb)/2.d0,', f(xmin)=',invent2((aa+bb)/2.d0)
write(11,'(A,/)') '-------------------------------------------------------------'
!!!!

write(11,'(A)') '-------------------------------------------------------------'
write(11,'(A)') 'Ejecucion del Algoritmo de Newton con funciones de dos variables:'
write(11,'(A)') '-------------------------------------------------------------'
!METODO DE NEWTON CON LA FUNCION DE CHICHINADZE
write(11,'(A)') 'La primera funcion para hacer el test es la FUNCION DE CHICHINADZE:'
write(11,'(A)') 'F(x,y)=x**2-12*x+11+10*cos(pi*x/2)+8*sin(5*pi*x)-exp(-(y-0.5)**2/2)/sqrt(5)'
x2(1)=5.83   !Estimacion inicial de la solucion
x2(2)=0
write(11,'(A,f8.5,a,f8.5,a)') 'Estimacion inicial X0=(x,y)=(',x2(1),',',x2(2),')'
222 call gradchichinadze(x2,grad2)   !Calculamos el gradiente en ese punto
call modulo(grad2,n2,norma)  !Calculamos la norma del gradiente
if (norma<eps) then
  goto 111
else
 call hesschichinadze(x2,hess2)    !Calculamos la matriz HESSIANA
 call inversamatriz2x2(hess2,invhess2) !Invertimos la Hessiana
 x2p=x2-matmul(invhess2,grad2) !Calculamos la multiplicacion del paso del algoritmo
endif
call modulo(x2p-x2,n2,norma)
if (norma<eps) then
  goto 111
else
  x2=x2p
  goto 222
endif
111 continue
write(11,'(a,f16.13,a,f16.13,a,f17.13)') 'El valor minimo encontrado es Xmin=(x,y)=(',x2(1),',',x2(2),'), F(Xmin)=',chichinadze(x2)
write(11,'(A,/)') '-------------------------------------------------------------'
!FIN DEL METODO DE NEWTON CON LA FUNCION CHICHINADZE


!METODO DE NEWTON CON LA FUNCION DE BEALE
write(11,'(A)') 'La segunda funcion para hacer el test es la FUNCION DE BEALE:'
write(11,'(A)') 'F(x,y)=(1.5-x+x*y)**2 + (2.25-x+x*y**2)**2 + (2.625-x+x*y**3)**2'
x2(1)=2
x2(2)=0
write(11,'(A,f8.5,a,f8.5,a)') 'Estimacion inicial X0=(x,y)=(',x2(1),',',x2(2),')'
223 call gradbeale(x2,grad2)
call modulo(grad2,n2,norma)
if (norma<eps) then
  goto 112
else
 call hessbeale(x2,hess2)
 call inversamatriz2x2(hess2,invhess2)
 x2p=x2-matmul(invhess2,grad2)
endif
call modulo(x2p-x2,n2,norma)
if (norma<eps) then
  goto 111
else
  x2=x2p
  goto 223
endif
112 continue
write(11,'(a,f16.13,a,f16.13,a,f17.13)') 'El valor minimo encontrado es Xmin=(x,y)=(',x2(1),',',x2(2),'), F(Xmin)=',beale(x2)
write(11,'(A,/)') '-------------------------------------------------------------'

!FIN DEL METODO DE NEWTON CON LA FUNCION DE BEALE
endprogram


double precision function chichinadze(x)   !Primera funcion Test de dos variables
implicit none
integer, parameter:: n=2
double precision x(n), pi
pi=atan(1.d0)*4.d0
chichinadze=x(1)**2.d0-12.d0*x(1)+11.d0+10.d0*cos(pi*x(1)/2.d0)+8.d0*sin(5.d0*pi*x(1))-exp(-(x(2)-0.5d0)**2.d0/2.d0)/dsqrt(5.d0)
endfunction   

double precision function beale(x)   !Segunda funcion Test de dos variables
implicit none
integer, parameter:: n=2
double precision x(n)
beale=(1.5d0-x(1)+x(1)*x(2))**2.d0 + (2.25d0-x(1)+x(1)*x(2)**2.d0)**2.d0 + (2.625d0-x(1)+x(1)*x(2)**3.d0)**2.d0
endfunction  


subroutine gradchichinadze(x,g)   !Calculo del gradiente de la Funcion 1 de dos variables
implicit none
integer, parameter:: n=2
double precision x(n),g(n),pi
pi=atan(1.d0)*4.d0
g(1)=40.d0*pi*cos(5.d0*pi*x(1))-5.d0*pi*sin(pi*x(1)/2.d0)+2.d0*x(1)-12.d0
g(2)=(x(2)-0.5d0)*exp(-((x(2)-0.50d0)**2.d0)/2.d0)/dsqrt(5.d0)
endsubroutine                                                   


subroutine hesschichinadze(x,h)   !Calculo de la matriz Hessiana de la Funcion 1 de dos variables
implicit none
integer, parameter:: n=2
double precision x(n), h(n,n),pi
pi=atan(1.d0)*4.d0
h(1,1)=-5.d0*pi**2.d0*cos(pi*x(1)/2.d0)/2.d0-200.d0*pi**2.d0*sin(5.d0*pi*x(1))+2.d0
h(2,2)=-((x(2)-0.5d0)**2.d0)*exp(-((x(2)-0.5d0)**2.d0)/2.d0)/dsqrt(5.d0)+exp(-((x(2)-0.5d0)**2.d0)/2.d0)/dsqrt(5.d0)       
h(1,2)=0.d0
h(2,1)=h(1,2)
endsubroutine


subroutine inversamatriz2x2(mat,inv)   !Procedimiento que calcula la inversa de una matriz de 2x2
implicit none
integer, parameter:: n=2
double precision mat(n,n),inv(n,n)
inv(1,1)=1.d0/mat(1,1)-mat(1,2)*mat(2,1)/(mat(1,1)**2.d0*(mat(1,2)*mat(2,1)/mat(1,1)-mat(2,2)))
inv(1,2)=mat(1,2)/(mat(1,1)*(mat(1,2)*mat(2,1)/mat(1,1)-mat(2,2)))
inv(2,1)=mat(2,1)/(mat(1,1)*(mat(1,2)*mat(2,1)/mat(1,1)-mat(2,2)))
inv(2,2)=-1.d0/(mat(1,2)*mat(2,1)/mat(1,1)-mat(2,2))
endsubroutine

subroutine modulo(vector,dim,norma)  !Procedimiento que calcula la norma de un vector de tamaño n
implicit none
integer dim,i
double precision vector(dim),norma
norma=0.d0
do i=1,dim
  norma=norma+vector(i)**2.d0
enddo
norma=dsqrt(norma)
endsubroutine 


subroutine gradbeale(x,g)  !Gradiente de la segunda funcion de 2 variables
implicit none
integer, parameter:: n=2
double precision x(n),g(n)
g(1)=2.d0*x(1)*x(2)**6.d0+2*x(1)*x(2)**4.d0+(-4.d0*x(1)+5.25d0)*x(2)**3.d0+(-2.d0*x(1)+4.5d0)*x(2)**2.d0+x(2)*(-4.d0*x(1)+3.d0)&
& +6.d0*x(1)-12.75d0
g(2)=6.d0*x(1)**2.d0*x(2)**5.d0 + 4.d0*x(1)**2.d0*x(2)**3.d0 + (-6.d0*x(1)**2.d0 + 15.75d0*x(1))*x(2)**2.d0 - 2.d0*x(1)**2.d0 +&
&(-2.d0*x(1)**2.d0 + 9.d0*x(1))*x(2) + 3.d0*x(1)
endsubroutine 


subroutine hessbeale(x,h)   !Hessiana de la segunda funcion de dos variables
implicit none
integer, parameter:: n=2
double precision x(n), h(n,n)
h(1,1)=2.d0*x(2)**6.d0 + 2.d0*x(2)**4.d0 - 4.d0*x(2)**3.d0 - 2.d0*x(2)**2.d0 - 4.d0*x(2) + 6.d0
h(2,2)=30.d0*x(1)**2.d0*x(2)**4.d0+12.d0*x(1)**2.d0*x(2)**2.d0-2.d0*x(1)**2.d0+(-12.d0*x(1)**2.d0+31.5d0*x(1))*x(2)+9.d0*x(1)       
h(1,2)=12.d0*x(1)*x(2)**5.d0+8.d0*x(1)*x(2)**3.d0+(-12.d0*x(1)+15.75d0)*x(2)**2.d0+(-4.d0*x(1)+9.d0)*x(2)-4.d0*x(1)+3.d0
h(2,1)=h(1,2)
endsubroutine                   


double precision function invent1(x)   !Primera funcion de una variable INVENTADA
implicit none
double precision x
invent1=sin(3.d0*x/2.d0)*exp(-2.d0*dsqrt(x))+log(x**2.d0+3.d0*x+0.1d0)-8.d0*x*cos(x+3.d0)
invent1=-invent1
endfunction                                                                               

double precision function fibonacci(i)    !Funcion que calcula los numeros de Fibonacci
implicit none
integer i
fibonacci=(((1.d0+dsqrt(5.d0))/2.d0)**(i+1)-(((1.d0-dsqrt(5.d0))/2.d0)**(i+1)))/dsqrt(5.d0)
endfunction

double precision function invent2(x)  !Segunda funcion de una variable INVENTADA
implicit none
double precision x
invent2=sin(x)+sin(10.d0*x/3.d0)+log(x)-0.84d0*x+3.d0
endfunction