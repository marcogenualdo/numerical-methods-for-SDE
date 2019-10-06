% CIR approximation using implicit euler scheme
%
% dX(t) = (a + b * X(t)) dt + sigma * sqrt(X(t)) dW(t)
% with the substitution Y(t) = sqrt(X(t)) the eq. becomes
% dY(t) = ((a - sigma^2/4) / (2Y(t)) + b/2 Y(t)) dt + sigma/2 dW(t)

% equation parameters
a = 1;
b = 1;
sigma = 1;

t0 = 0;
T = 1;
x0 = 1;

% order test parameters
tries = 10;
dt_max = 0.5;
dt_min = 0.01;
samples = 500000;

dts = logspace(log10(dt_max), log10(dt_min), tries);
errors = zeros(1, tries);

for n = 1:tries
    % setting up
    dt = dts(n);
    steps = floor(T / dt);
    t = t0 : dt : t0+T;
    
    y = zeros(samples, steps+1);
    y(:,1) = sqrt(x0) * ones(samples, 1);

    % generating gaussian noise
    dW = sqrt(dt) * randn(samples,steps+1);

    % running the scheme
    for k=1:steps
        y(:,k+1) = (y(:,k) + sigma / 2 * dW(:,k) + sqrt( (y(:,k) + sigma / 2 * dW(:,k)).^2 ...
            + 2 * (1 - b * dt / 2) * (a - sigma^2 / 4) * dt)) ...
            / (2 * (1 - b * dt / 2));
    end
    x = y .* y;
    
    % error estimation
    true_expectations = (x0 + a/b) * exp(b * t) - a/b;
    errors(n) = max(abs(mean(x,1) - true_expectations));
    
    %plot(t, mean(x,1),t,true_expectations);
    %pause(1);
end

% fitting line
fit_parameters = polyfit(log(dts), log(errors), 1);
line_x = linspace(10^(-2), 1, 2);
line_y = exp(fit_parameters(1) * log(line_x) + fit_parameters(2));

% plotting
loglog(dts,errors, 'o');
hold on;
loglog(line_x, line_y);
xlabel('\Deltat');
ylabel('errore');
hold off;