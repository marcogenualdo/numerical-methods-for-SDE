% CIR split-step approximation using euler scheme
%
% dX(t) = (a + b * X(t)) dt + sigma * sqrt(X(t)) dW(t)
% split into:
% dY(t) = a * dt + sigma * sqrt(Y(t)) dW(t)  sampled from chi_d(lambda)
% dZ(t) = b * Z(t) * dt  solved exactly/with deterministic euler

function [dts, error_mat, logfit, logfit_mse] = strong_ord_cir_split(a,b,sigma)
tests_num = length(b);

% equation parameters
t0 = 0;
T = 1;
x0 = 1;

% test parameters
tries = 1 + 8;
start_steps = 2*2048;
samples = 200;

% initializing
error_mat = zeros(tests_num, tries-1);
logfit = zeros(tests_num, 2);
logfit_mse = zeros(tests_num, 1);

% loop over coefficients
for test=1:tests_num
    % initializing
    steps = start_steps;
    errors = zeros(1,tries);
    
    % generating gaussian noise
    %d = floor(4 * a(test) / sigma(test)^2);
    %dN = randn(samples, d, start_steps);
    dN = randn(samples, start_steps);
    B = zeros(samples, start_steps);
    for k=1:start_steps-1
        B(:,k+1) = B(:,k) + dN(:,k);
    end

    % loop over time steps
    for run = 1:tries
        dt = T / steps;
        Bsq = dt * B .^ 2;
        skip = 2 ^ (run - 1);

        t = linspace(t0, t0 + T, steps);
        z = zeros(samples, steps);
        z(:,1) = x0 * ones(1, samples,1);

        % running the scheme    
        for k=1:(steps-1)
            %lambda = 4 / (sigma(test)^2 * dt) * z(:,k) * ones(1,d);
            %mi = sqrt(lambda / d);
            %y = 0.25 * sigma(test)^2 * dt * sum((dN(:,:,k) + mi).^2, 2);
            
            y = z(:,k) + sigma(test)^2 / 4 * (Bsq(:,skip * (k+1)) - Bsq(:,skip * k)); 
            z(:,k+1) = y * exp(b(test) * dt);
        end

        if run > 1
            %plot(t,z(1,:), t, z_prev(1, 1:2:end));
            %pause(1);
            errors(run) = mean(max((z - z_prev(:,1:2:end)) .^ 2, [], 2));
        end

        z_prev = z;
        steps = steps / 2;
        dN = dN(:, 1:2:end-1) + dN(:, 2:2:end);
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

plot(dts, error_mat(1,:), 'o');
hold on
plot(line_x, line_y);
hold off

return