% CIR split-step approximation using euler scheme
%
% dX(t) = (a + b * X(t)) dt + sigma * sqrt(X(t)) dW(t)
% split into:
% dY(t) = a * dt + sigma * sqrt(Y(t)) dW(t)  sampled from chi_d(lambda)
% dZ(t) = b * Z(t) * dt  solved exactly/with deterministic euler

% equation parameters
a = 1;
b = 1;
sigma = 2;

t0 = 0;
T = 1;
x0 = 1;

% setting up
n = 200;
dt = T / n;
t = linspace(t0, t0+T, n+1);
z = zeros(1,n+1);
z(1) = x0;
y = x0;

% solving
d = 4 * a / sigma^2;

for k=1:n
    lambda = 4 * z(k) / (sigma^2 * dt);
    y = 0.25 * sigma^2 * dt * ncx2rnd(d, lambda);
    z(k+1) = y + dt * b * y;
end

plot(t,z);