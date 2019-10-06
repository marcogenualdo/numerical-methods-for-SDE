% CIR approximation using implicit euler scheme
%
% dX(t) = (a + b * X(t)) dt + sigma * sqrt(X(t)) dW(t)
% with the substitution Y(t) = sqrt(X(t)) the eq. becomes
% dY(t) = ((a - sigma^2/4) / (2Y(t)) + b/2 Y(t)) dt + sigma/2 dW(t)

function [dts, errors, logfit, logfit_mse] = weak_ord_cir_eulimp (a,b,sigma)

if ~((length(a) == length(b)) && (length(b) == length(sigma)))
    error('Input vectors must be of the same length.');
end
len_param = length(a);

% equation parameters
t0 = 0;
T = 1;
x0 = 1;

% order test parameters
tries = 10;
dt_max = 0.1;
dt_min = 0.005;
samples = 5e+06;

dts = logspace(log10(dt_max), log10(dt_min), tries);
errors = zeros(len_param, tries);

for n = 1:tries
    % setting up
    dt = dts(n);
    steps = floor(T / dt);
    t = t0 : dt : (t0+T);

    y = zeros(samples, len_param);
    yp = sqrt(x0) * ones(samples, len_param);

    % running the scheme
    for k=1:steps
        dW = sqrt(dt) * randn(samples,1);
        
        for j=1:len_param
            y(:,j) = (yp(:,j) + sigma(j) / 2 * dW + sqrt( (yp(:,j) + sigma(j) / 2 * dW).^2 ...
                + 2 * (1 - b(j) * dt / 2) * (a(j) - sigma(j)^2 / 4) * dt)) ...
                / (2 * (1 - b(j) * dt / 2));
        end
        yp = y;
    end
    x = y .* y;
    
    % error estimation
    true_expectation = (x0 + a ./ b) .* exp(b * t(end)) - a ./ b;
    errors(:,n) = abs(mean(x,1) - true_expectation);
end

% fitting lines
logfit = zeros(len_param, 2);
logfit_mse = zeros(len_param, 1);
for j=1:len_param
    logfit(j,:) = polyfit(log(dts), log(errors(j,:)), 1);
    logfit_mse(j) = mean((errors(j,:) - ...
        exp(logfit(j,1) * log(dts) + logfit(j,2))).^2);
end

% test plot
loglog(dts,errors(1,:), 'o');
line_x = linspace(dts(1), dts(end), 2);
line_y = exp(logfit(1,1) * log(line_x) + logfit(1,2));
hold on;
loglog(line_x, line_y);
xlabel('\Deltat');
ylabel('errore');
hold off;

return