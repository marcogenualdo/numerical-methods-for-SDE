% CIR split-step approximation using euler scheme
%
% dX(t) = (a + b * X(t)) dt + sigma * sqrt(X(t)) dW(t)
% split into:
% dY(t) = a * dt + sigma * sqrt(Y(t)) dW(t)  sampled from chi_d(lambda)
% dZ(t) = b * Z(t) * dt  solved exactly/with deterministic euler

function [dts, error_mat, logfit, logfit_mse] = strong_ord_cir_split2(a,b,sigma)
tests_num = length(b);

% equation parameters
t0 = 0;
T = 1;
x0 = 1;

% test parameters
n_dts = 8;
start_steps = 2*2048;
samples = 2000;

% initializing
error_mat = zeros(tests_num, n_dts-1);
logfit = zeros(tests_num, 2);
logfit_mse = zeros(tests_num, 1);

% loop over coefficients
for test=1:tests_num
    % initializing
    steps = start_steps;
    dt = T / steps;
    errors = zeros(1,n_dts);
    
    % generating gaussian noise
    dB = sqrt(dt) * randn(samples, start_steps);
    B = zeros(samples, start_steps);
    
    % ou = Ornstein-Ulhenbeck
    ou = zeros(samples, start_steps);
    ou(:,1) = sqrt(x0) * ones(samples,1);
    int = zeros(samples,1);
    
    for k=1:start_steps-1
        B(:,k+1) = B(:,k) + dB(:,k);    
        
        %computing exact solution
        int = int + exp(- b(test) / 2 * dt * (k-1)) * dB(:,k);
        ou(:,k+1) = exp(b(test) / 2 * dt * k) * (ou(:,1) + sigma(test) / 2 * int);
    end
    Bsq = B .^ 2;
    sol = ou .* ou;

    % loop over time steps
    for run = 1:n_dts
        dt = T / steps;
        skip = 2 ^ (run - 1);

        z = zeros(samples, steps);
        z(:,1) = x0 * ones(1, samples,1);

        % running the scheme    
        for k=1:(steps-1)
            y = z(:,k) + sigma(test)^2 / 4 * (Bsq(:,skip * (k+1)) - Bsq(:,skip * k)); 
            %y = ((z(:,k) + sigma(test) / 2  * dB(:,skip * k) ...
            %    + sqrt( (z(:,k) + sigma(test) / 2 * dB(:,skip * k)).^2 ...
            %    + 2 * (a(test) - sigma(test)^2 / 4) * dt))  / 2) .^ 2;
            z(:,k+1) = y * exp(b(test) * dt);
        end

        % optional: plottin g one path of z vs exact sol
        t = linspace(t0, t0 + T, steps);
        plot(t,z(1,:), t, sol(1, 1:skip:end));
        pause(1);
        
        % estimating E(sup_t |z(t) - sol(t)|^2)
        errors(run) = mean(max((z - sol(:,1:skip:end)) .^ 2, [], 2));
        
        steps = steps / 2;
    end
    
    % creating error table
    steps = start_steps;
    dts = zeros(1,n_dts);
    for run=1:n_dts
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