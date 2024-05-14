% The Model with disease
function dYdt = predator_prey_model(t, Y, params)
    % Unpack the parameters
    r = params.r;
    K = params.K;
    a1 = params.a1;
    a2 = params.a2;
    gamma1 = params.gamma1;
    beta = params.beta;
    delta = params.delta;
    e1 = params.e1;
    e2 = params.e2;
    
    % Unpack the state variables
    S = Y(1);
    I = Y(2);
    P1 = Y(3);
    P2 = Y(4);
    
    % Define the differential equations
    dSdt = r*S*(1 - (S + I)/K) - a1*S*P1 - a2*S*P2;
    dIdt = beta*S - gamma1*I - a1*I*P1 - a2*I*P2;
    dP1dt = e1*a1*(S + I)*P1 - delta*P1;
    dP2dt = e2*a2*(S + I)*P2 - delta*P2;
    
    dYdt = [dSdt; dIdt; dP1dt; dP2dt];
end
% Parameters
params.r = 0.5; % Growth rate of prey
params.K = 100; % Carrying capacity of prey
params.a1 = 0.01; % Predation rate of P1
params.a2 = 0.01; % Predation rate of P2
params.gamma1 = 0.1; % Recovery rate of infected prey
params.beta = 0.1; % Infection rate
params.delta = 0.1; % Death rate of predators
params.e1 = 0.1; % Conversion efficiency of P1
params.e2 = 0.1; % Conversion efficiency of P2

% Initial conditions
S0 = 100;
I0 = 20;
P1_0 = 10;
P2_0 = 2;
Y0 = [S0; I0; P1_0; P2_0];

% Time span
tspan = [0 200];
% Solve the system
[t, Y] = ode45(@(t, Y) predator_prey_model(t, Y, params), tspan, Y0);

% Plot the results
figure;
subplot(2,1,1);
plot(t, Y(:,1), 'b', 'DisplayName', 'S (Susceptible Prey)');
hold on;
plot(t, Y(:,2), 'r', 'DisplayName', 'I (Infected Prey)');
xlabel('Time');
ylabel('Population');
legend;
title('Prey Population Dynamics');

subplot(2,1,2);
plot(t, Y(:,3), 'g', 'DisplayName', 'P1 (Susceptible Predators)');
hold on;
plot(t, Y(:,4), 'm', 'DisplayName', 'P2 (Infected Predators)');
xlabel('Time');
ylabel('Population');
legend;
title('Predator Population Dynamics');