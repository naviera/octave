function [optval,sol,solarray] = FOM(A,b,tau,NITER)

###########
#
#  Esta funci�n recibe como input una matriz A y un vecor b y busca resolver el problema
#  de regularizaci�n L-1 (LASSO) con un par�metro tau y haciendo a lo m�s NITER iteraciones
#  de un m�todo b�sico de subgradiente.
#
#  La funci�n  devuelve el valor �ptimo, optval, el vector soluci�n sol y una matriz 
#  solarray donde cada fila representa una iteraci�n. 
#  Las columnas de solarray corresponden a:
#  - el n�mero de iteraci�n
#  - el valor objetivo actual en esa iteraci�n
#  - la norma 1 de la soluci�n actual.
#  - el error de estimaci�n: ||Ax-b||/||b||  (con normas infinito, valor m�ximo de desviaci�n)
#
#

cpuini = cputime;

[m,n] = size(A)
#
# xk es el punto actual xk1 es x(k-1) y xk2 es x(k-2)
# Se usan para evaluar el avance de la trayectoria del algoritmo
#
xk = zeros(n,1);
xk1 = xk;
xk2 = xk1;
zig = 0;
thetak = 1;
L2 = 0;
angle = 0;

#
# Ac� se estiman R y L seg�n como se desarrolla aprox. en clases
#
R = sqrt(norm((A*A')^-1))*norm(b)
L = tau*sqrt(n) + norm(A*A')*R + norm(A'*b)

optval = tau*norm(xk,1) + norm(A*xk-b)**2;

#
# Este es el n�mero de iteraciones a realizar, seg�n 
# el teorema de convergencia
#
thetak = R/(L*sqrt(NITER+1));

for t = 1:NITER,
  #
  # Aqu� se calcula un subgradiente de la funci�n objetivo
  #
  g = tau*sign(xk) + 2*A'*(A*xk - b);
  xk2 = xk1;
  xk1 = xk;
  #
  # Aqu� se hace el avance en la direcci�n
  #
  xk = xk - thetak*g;
  #
  # Ac� se eval�a el error de ajuste de Ax = b
  #
  aux = max(abs(A*xk-b))/max(abs(b));
  optval = tau*norm(xk,1) + norm(A*xk-b)**2;
  #
  # Aqu� se calcula el producto interno entre los vectores de avance
  # en iteraciones sucesivas
  #
  if t >= 3
    angle = ((xk2-xk1)'*(xk-xk1)/(norm(xk2-xk1)*norm(xk-xk1)));
  endif
  #
  # La matriz solarray almacena datos de las iteraciones
  #
  solarray(t,1) = t;
  solarray(t,2) = optval;
  solarray(t,3) = thetak;
  solarray(t,4) = norm(xk,1);
  solarray(t,5) = aux;
  solarray(t,5) = angle;
  printf("%6i %10.6f %10.6f %10.6f %10.6f  %10.6f \n",t,optval,thetak,norm(xk,1),aux,angle);
  fflush(stdout);
end
 
cpufinal = cputime-cpuini
sol = xk;


endfunction