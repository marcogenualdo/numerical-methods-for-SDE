% CIR approximation using implicit euler scheme
%
% dX(t) = (a + b * X(t)) dt + sigma * sqrt(X(t)) dW(t)
% with the substitution Y(t) = sqrt(X(t)) the eq. becomes
% dY(t) = ((a - sigma^2/4) / (2Y(t)) + b/2 Y(t)) dt + sigma/2 dW(t)

function [dts, error_mat, logfit, logfit_mse] = strong_ord_cir_eulimp(a,b,sigma)
tests_num = length(a);

% equation parameters
t0 = 0;
T = 1;
x0 = 1;

% test parameters
tries = 1 + 8;
start_steps = 2*2048;
samples = 500;
 
error_mat = zeros(tests_num, tries-1);
logfit = zeros(tests_num, 2);
logfit_mse = zeros(tests_num, 1);

% generating gaussian noise
start_dW = sqrt(T / start_steps) * randn(samples, start_steps);

% loop over coefficients
for test = 1:tests_num
    % initializing
    steps = start_steps;
    dW = start_dW;

    errors = zeros(1,tries);
    x_prev = zeros(samples, steps);

    % loop over time steps
    for run = 1:tries
        dt = T / steps;

        t = linspace(t0, t0 + T, steps);
        y = zeros(samples, steps);
        y(:,1) = sqrt(x0) * ones(samples, 1);

        for k=1:steps-1
            y(:,k+1) = (y(:,k) + sigma(test) / 2  * dW(:,k) ...
                + sqrt( (y(:,k) + sigma(test) / 2 * dW(:,k)).^2 ...
                + 2 * (1 - b(test) * dt / 2) * (a(test) - sigma(test)^2 / 4) * dt)) ...
                / (2 * (1 - b(test) * dt / 2));
        end
        x = y .* y;

        if run > 1
            %plot(t,x(1,:), t, x_prev(1, 1:2:end));
            %pause(1);
            errors(run) = mean(max((x - x_prev(:,1:2:end)) .^ 2, [], 2));
        end

        x_prev = x;
        steps = steps / 2;
        dW = dW(:, 1:2:end-1) + dW(:, 2:2:end);
    end

    % creating error table
    steps = start_steps;
    dts = zeros(1,tries);
    for run=1:tries
        dts(run) = 1 / steps;
        steps = steps / 2;
    end

    dts = dts(2:end);
    errors = errors(2:end);
    error_mat(test,:) = errors;

    % fitting line
    logfit(test,:) = polyfit(log(dts), log(errors), 1);
    logfit_mse(test) = mean((errors - ...
        exp(logfit(test,1) * log(dts) + logfit(test,2))).^2);
end

% test plot
line_x = linspace(dts(1), dts(end), 2);
line_y = exp(logfit(1,1) * log(line_x) + logfit(1,2));

loglog(dts, error_mat(1,:), 'o');
hold on
loglog(line_x, line_y);
hold off

return