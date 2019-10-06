% CIR approximation using implicit euler scheme
%
% dX(t) = (a + b * X(t)) dt + sigma * sqrt(X(t)) dW(t)
% with the substitution Y(t) = sqrt(X(t)) the eq. becomes
% dY(t) = ((a - sigma^2/4) / (2Y(t)) - b/2 Y(t)) dt + sigma/2 dW(t)

% equation parameters
a = 1;
b = 1;
sigma = 1;

t0 = 0;
T = 1;
x0 = 0.04;

n = 200;

% generating gaussian noise
noise = randn(1,n+1);

% method setup
tol = 1e-2;
maxiter = 50;

drift = @(t,x) 0.5 * (a - sigma^2 / 4) / x +  0.5 * b * x;
ddrift = @(t,x) -0.5 * (a - sigma^2 / 4) / (x*x) + 0.5 * b;
diffusion = @(t,x) sigma / 2 + 0*x;

% solving
[t,y] =  eul_mar_imp (drift,diffusion,t0,T,x0,n,noise,ddrift,tol,maxiter);
plot(t,y.*y);