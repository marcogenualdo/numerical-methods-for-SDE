% implicit-euler for SDE: dX(t) = b(t,X(t))*dt + sigma(t,X(t))*dW(t)
%
% INPUTS: b,sigma = handles of SDE coefficients; 
%         t0,T = starting time and time span;
%         x0 = starting value;
%         n = number of steps;
%         noise = 1xn array of numbers iid N(0,1).
%
%         db = handle of the partial x derivative of b
%         tol, maxiter = increment-based arrest and max iterations for
%         newton algorithm

function [t,u] = eul_mar_imp (b,sigma,t0,T,x0,n,noise,db,tol,maxiter)

t = linspace (t0,t0+T,n+1);
dt = T/n;
dW = sqrt(dt) * noise;
u = zeros (1,n+1);
u(1) = x0;

for k=1:n
    f = @(y) y - u(k) - b(t(k+1), y)*dt - sigma(t(k), u(k)) * dW(k);
    df = @(y) 1 - db(t(k+1), y) * dt;
    
    u(k+1) = newton(f,df,tol,maxiter,u(k));
end
return