% euler-maruyama for SDE: dX(t) = b(t,X(t))*dt + sigma(t,X(t))*dW(t)
%
% INPUTS: b,sigma = handles of SDE coefficients; 
%         t0,T = starting time and time span;
%         x0 = starting value;
%         n = number of steps;
%         noise = 1xn array of iid N(0,1) distributed numbers.

function [t,u] = eul_mar (b,sigma,t0,T,x0,n,noise)

t = linspace (t0,t0+T,n+1);
dt = T/n;
dW = sqrt(dt) * noise;
u = zeros (1,n+1);
u(1) = x0;

for i=1:n
  u(i+1) = u(i) + b(t(i), u(i))*dt + sigma(t(i), u(i)) * dW(i);
end
return