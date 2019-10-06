% CIR split-step approximation using euler scheme
%
% dX(t) = (a + b * X(t)) dt + sigma * sqrt(X(t)) dW(t)
% split into:
% dY(t) = a * dt + sigma * sqrt(Y(t)) dW(t)  sampled from chi_d(lambda)
% dZ(t) = b * Z(t) * dt  solved exactly/with deterministic euler

function [dts, error_mat, logfit, logfit_mse] = weak_ord_cir_split2(n,b,sigma)
tests_num = length(b);

% equation parameters
t0 = 0;
T = 1;
x0 = 1;

% order test parameters
n_dts = 10;
dt_max = 0.1;
dt_min = 0.005;
samples = 1e+06;

dts = logspace(log10(dt_max), log10(dt_min), n_dts);

% initializing
error_mat = zeros(tests_num, n_dts);
logfit = zeros(tests_num, 2);
logfit_mse = zeros(tests_num, 1);

% loop over coefficients
for test=1:tests_num
    errors = zeros(1, n_dts);
    
    %{
    % computing exact solution
    dt = dts(end);
    steps = floor(T / dt);
    sol = zeros(samples, steps);
    
    for ni=1:n(test)
        % generating gaussian noise
        dB = sqrt(dt) * randn(samples, steps);

        % ou = Ornstein-Ulhenbeck
        ou = zeros(samples, steps);
        ou(:,1) = sqrt(x0 / n(test)) * ones(samples,1);
        int = zeros(samples,1);

        for k=1:steps-1
            % computing ou(t(k))
            int = int + exp(- b(test) / 2 * dt * (k-1)) * dB(:,k);
            ou(:,k+1) = exp(b(test) / 2 * dt * k) * (ou(:,1) + sigma(test) / 2 * int);
        end
        sol = sol + ou .^ 2;
    end

    t = linspace(t0, t0+T, steps);
    a = n * sigma(test)^2/4;
    plot(t,mean(sol,1), t, (x0 + a(test) / b(test)) .* exp(b(test) * t) - a(test) / b(test));
    pause(4);
    %}
    
    % loop over time steps
    for run = 1:n_dts
        dt = dts(run);
        steps = floor(T / dt);

        z = zeros(samples, steps);
        z(:,1) = x0 * ones(1, samples,1);

        % running the scheme    
        for k=1:(steps-1)
            y = z(:,k) + sigma(test)^2 / 4 * dt * sum(randn(samples,n(test)) .^2, 2); 
            z(:,k+1) = y * exp(b(test) * dt);
        end

        % optional: plotting one path of z vs exact sol
        %t = linspace(t0, t0 + T, steps);
        %plot(t,z(1,:), t, sol(1, 1:skip:end));
        %pause(1);
        
        % estimating E(sup_t |z(t) - sol(t)|^2)
        t = linspace(t0, t0 + T, steps);
        a = n(test) * sigma(test) ^ 2 / 4;
        true_expectation = (x0 + a / b(test)) .* exp(b(test) * t) - a / b(test);
        errors(run) = sqrt(max((mean(z,1) - true_expectation) .^ 2));
    end
    
    % filling error table
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