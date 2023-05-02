%% Sliding Mode Control of Inverted Pendulum on Cart

clc; clear all; close all;

%% System Model

%  States: x = [x1, x2]^T = [theta, thetaDot]^T
%
%  | x1 | = |           x2           |
%  | x2 |   | f(x1, x2) + b(x1, x2)u |
%
%  f(x1, x2) = ((M + m) * g * sin(x1) - m * l * x2^2sin(x1)cos(x1)) / ...
%              ((4 / 3) * (M + m) * l - m * g * l * cos^2(x1))
%  b(x1, x2) = cos(x1) / ((4 / 3) * (M + m) * l - m * g * l * cos^2(x1))

%% Controller Design



%% System Specifications

M     = 1; % kg
m     = 0.1; % kg
l     = 0.5; % m
g     = 9.81; % m / s

%% Control System Simulation

dt = 0.001;
t = dt:dt:10;
A = [1, dt; 0, 1];
H = eye(2);
P = 10 .* eye(2); P_obs = P; Pplus = P; P1 = P(1, 1); P2 = P(2, 2);
x = [0; 0]; xObs = x;
xHat = x; xHatPlus = x;
sigmaProcess = [deg2rad(0.1);
                deg2rad(0.002)]; % Disturbance

sigmaMeasure = [deg2rad(0.05);
                deg2rad(0.001)]; % Measurement Noise
Q = [sigmaProcess(1)^2, 0; 0, sigmaProcess(2)^2];
R = [sigmaMeasure(1)^2, 0; 0, sigmaMeasure(2)^2];

lambda = 5; k = 18;
controlSaturation = 30;

dt = 0.001;
t = dt:dt:10;

%{
x1Des = zeros(size(t));
x2Des = zeros(size(t));
x2DesDot = zeros(size(t));
%}

x1Des = (pi / 30) .* sin(t);
x2Des = (pi / 30) .* cos(t);
x2DesDot = -(pi / 30) .* sin(t);
%}
%% System with No Observer

for i = 2:length(t)

    e = x(1, i - 1) - x1Des(i);
    eDot = x(2, i - 1) - x2Des(i);

    s = lambda * e + eDot;
    f = ((M + m) * g * sin(x(1, i - 1)) - m * l * x(2, i - 1)^2 * sin(x(1, i - 1)) * cos(x(1, i - 1))) / ((4 / 3) * (M + m) * l - m * l * cos(x(1, i - 1))^2);
    b = cos(x(1, i - 1)) / ((4 / 3) * (M + m) * l - m * l * cos(x(1, i - 1))^2);
    h = lambda * (x(2, i - 1) - x2Des(i)) + (f - x2DesDot(i));
    u(:, i) = -(1 / b) * (h + k * sign(s));
    if abs((u(i))) > controlSaturation
        u(i) = sign(u(i)) * controlSaturation;
    end

    xDot1 = x(2, i - 1);
    xDot2 = f + b * u(:, i);
    x1 = x(1, i - 1) + xDot1 * dt;
    x2 = x(2, i - 1) + xDot2 * dt;
    xProcess = [x1; x2] + [sigmaProcess(1) * randn; sigmaProcess(2) * randn];
    x(:, i) = H * xProcess + [sigmaMeasure(1) * randn; sigmaMeasure(2) * randn];

end

%% System with Observer

for i = 2:length(t)

    e = xObs(1, i - 1) - x1Des(i);
    eDot = xObs(2, i - 1) - x2Des(i);

    s = lambda * e + eDot;
    f = ((M + m) * g * sin(xObs(1, i - 1)) - m * l * xObs(2, i - 1)^2 * sin(xObs(1, i - 1)) * cos(xObs(1, i - 1))) / ((4 / 3) * (M + m) * l - m * l * cos(xObs(1, i - 1))^2);
    b = cos(xObs(1, i - 1)) / ((4 / 3) * (M + m) * l - m * l * cos(xObs(1, i - 1))^2);
    h = lambda * (xObs(2, i - 1) - x2Des(i)) + (f - x2DesDot(i));
    uObs(:, i) = -(1 / b) * (h + k * sign(s));
    if abs((uObs(i))) > controlSaturation
        uObs(i) = sign(uObs(i)) * controlSaturation;
    end

    xDot1 = xObs(2, i - 1);
    xDot2 = f + b * uObs(:, i);
    x1 = xObs(1, i - 1) + xDot1 * dt;
    x2 = xObs(2, i - 1) + xDot2 * dt;
    y(:, i) = [x1; x2] + [sigmaMeasure(1) * randn; sigmaMeasure(2) * randn];

    xHatMinus = [x1; x2] + [sigmaProcess(1) * randn; sigmaProcess(2) * randn];
    Pminus = A * Pplus * A' + Q;

    L = Pminus * H' * (inv(H * Pminus * H' + R)); L1(i) = L(1, 1); L2(i) = L(2, 2);
    xHatPlus = xHatMinus + L * (y(:, i) - H * xHatMinus);
    Pplus = (eye(2) - L * H) * Pminus; P1(:, i) = Pplus(1, 1); P2(:, i) = Pplus(2, 2);
    xObs(:, i) = xHatPlus;

end

meanTrackingErrorNoObs = mean(rad2deg(x1Des - x(1, :)));
meanTrackingErrorObs = mean(rad2deg(x1Des - xObs(1, :)));
stdTrackingErrorNoObs = std(rad2deg(x1Des - x(1, :)));
stdTrackingErrorObs = std(rad2deg(x1Des - xObs(1, :)));
disp(['Mean Tracking Error without Observer: ', num2str(meanTrackingErrorNoObs), ' Degrees'])
disp(['Mean Tracking Error with Observer: ', num2str(meanTrackingErrorObs), ' Degrees'])
disp(['Standard Deviation Tracking Error without Observer: ', num2str(stdTrackingErrorNoObs), ' Degrees'])
disp(['Standard Deviation Tracking Error with Observer: ', num2str(stdTrackingErrorObs), ' Degrees'])

%% Plotting

settings.markerSize = 15;
settings.linewidth = 5;
settings.axisFont = 25;
settings.titleFont = 30;

figure()
plot(t, rad2deg(x(1, :)), '.', 'MarkerSize', settings.markerSize)
hold on
plot(t, rad2deg(xObs(1, :)), '.', 'MarkerSize', settings.markerSize)
hold on
plot(t, rad2deg(x1Des), '.', 'MarkerSize', settings.markerSize)
ax = gca; ax.FontSize = settings.axisFont; axis tight
xlabel('Time [Seconds]', 'Interpreter', 'Latex', 'FontSize', settings.axisFont)
ylabel('$\theta$ [Degrees]', 'Interpreter', 'Latex', 'FontSize', settings.axisFont)
legend('Raw Measurements', 'Kalman Filter', 'Reference', 'Interpreter', 'Latex', 'FontSize', settings.axisFont)
title('Control System Response', 'Interpreter', 'Latex', 'FontSize', settings.titleFont)

figure()
subplot(2, 1, 1)
plot(t, u, 'LineWidth', settings.linewidth)
ax = gca; ax.FontSize = settings.axisFont; axis tight
xlabel('Time [Seconds]', 'Interpreter', 'Latex', 'FontSize', settings.axisFont)
ylabel('$\tau$ [Nm]', 'Interpreter', 'Latex', 'FontSize', settings.axisFont)
subplot(2, 1, 2)
plot(t, uObs, 'LineWidth', settings.linewidth)
ax = gca; ax.FontSize = settings.axisFont; axis tight
xlabel('Time [Seconds]', 'Interpreter', 'Latex', 'FontSize', settings.axisFont)
ylabel('$\tau$ [Nm]', 'Interpreter', 'Latex', 'FontSize', settings.axisFont)
sgtitle('Controller Input', 'Interpreter', 'Latex', 'FontSize', settings.titleFont) 

figure()
plot(t, rad2deg(x1Des - x(1, :)), '.', 'MarkerSize', settings.markerSize)
hold on
plot(t, rad2deg(x1Des - xObs(1, :)), '.', 'MarkerSize', settings.markerSize)
ax = gca; ax.FontSize = settings.axisFont; axis tight
xlabel('Time [Seconds]', 'Interpreter', 'Latex', 'FontSize', settings.axisFont)
ylabel('$\theta$ [Degrees]', 'Interpreter', 'Latex', 'FontSize', settings.axisFont)
legend('Raw Measurements', 'Kalman Filter', 'Interpreter', 'Latex', 'Location', 'southeast', 'FontSize', settings.axisFont)
title('Tracking Error', 'Interpreter', 'Latex', 'FontSize', settings.titleFont)
