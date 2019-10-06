% CIR approximation using implicit euler scheme
%
% dX(t) = (a + b * X(t)) dt + sigma * sqrt(X(t)) dW(t)
% with the substitution Y(t) = sqrt(X(t)) the eq. becomes
% dY(t) = ((a - sigma^2/4) / (2Y(t)) + b/2 Y(t)) dt + sigma/2 dW(t)

% equation parameters
a = 1;
b = 1;
sigma = 2;

t0 = 0;
T = 10;
x0 = 1;

% discretization parameters
n = 200;
dt = T / n;

% generating gaussian noise
dW = sqrt(dt) * randn(1,n+1);
t = linspace(t0, t0 + T, n+1);
y = zeros(1,n+1);
y(1) = x0;

for k=1:n
    y(k+1) = (y(k) + sigma / 2 * dW(k) + sqrt( (y(k) + sigma / 2 * dW(k))^2 ...
        + 2 * (1 - b * dt / 2) * (a - sigma^2 / 4) * dt)) ...
        / (2 * (1 - b * dt / 2));
end
x = y .* y;

plot(t, x);