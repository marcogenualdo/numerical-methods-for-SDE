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

% order test parameters
dt_max = 0.5;
dt_min = 0.01;
tries = 10;
samples = 50000;

dts = logspace(log10(dt_max), log10(dt_min), tries);
errors = zeros(1, tries);

for n = 1:tries
    % setting up
    dt = dts(n);
    steps = floor(T / dt);
    t = t0 : dt : t0+T;

    z = zeros(samples,steps+1);
    z(:,1) = x0 * ones(samples,1);
    y = x0 * ones(samples,1);

    % running the scheme
    d = 4 * a / sigma^2 * ones(samples,1);

    for k=1:steps
        lambda = 4 * z(:,k) / (sigma^2 * dt);
        y = 0.25 * sigma^2 * dt * ncx2rnd(d, lambda);
        z(:,k+1) = y * exp(b * dt);
    end
    
    % error estimation
    true_expectations = ((x0 + a/b) * exp(b * t) - a/b);
    errors(n) = max(abs(mean(z,1) - true_expectations)); 
    
    %plot(t, mean(z,1),t,true_expectations);
    %pause(5);
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