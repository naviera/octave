m=30;
n=300;
A=rand(m,n)*50;
b=rand(m,1)*50;

tau=500;
NITER=8000;
FOM(A,b,tau,NITER)

