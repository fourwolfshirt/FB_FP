%% Sliding Mode Control of Inverted Pendulum on Cart

clc; clear all; close all;

%% System Model

%  States: x = [x1, x2]^T = [theta, thetaDot]^T
%
%  | x1 | = |           x2           |
%  | x2 |   | f(x1, x2) + b(x1, x2)u |
%
%  f(x1, x2) = (m * g * l * sin(x1)) / (2 * J) + (m * g * (l + r) *
%  sin(x1)) / J
%  b(x1, x2) = 1 / J

%% System Parameters

mRod      = 0.01814; % [kg]
mBob      = 0.038; % [kg]
rBob      = 0.013; % [m]
l         = 0.1016; % [m]
g         = 9.81; % [m / s^2]
J         = 0.1; % [Nm]
J_rod     = (1 / 3) * mRod * l^2;
J_bob     = mBob * (l + rBob)^2;
J_eff     = J_rod + J_bob;

%% Control System Simulation

dt = 0.001;
t = dt:dt:10;
A = [1, dt; 0, 1];
B = 1 / J_eff;
H = eye(2);
P = 10 .* eye(2); P_obs = P; Pplus = P; P1 = P(1, 1); P2 = P(2, 2);
x = [pi / 2; 0]; xObs = x;
xHat = x; xHatPlus = x;
sigmaProcess = [deg2rad(0.1);
                deg2rad(0.002)]; % Disturbance

sigmaMeasure = [deg2rad(0.05);
                deg2rad(0.001)]; % Measurement Noise
Q = [sigmaProcess(1)^2, 0; 0, sigmaProcess(2)^2];
R = [sigmaMeasure(1)^2, 0; 0, sigmaMeasure(2)^2] / 100;

lambda = 5; k = 18;
controlSaturation = 30; % [Nm]

%{
x1Des = zeros(size(t));
x2Des = zeros(size(t));
x2DesDot = zeros(size(t));
%}

x1Des = pi / 2 + (pi / 30) .* sin(t);
x2Des = (pi / 30) .* cos(t);
x2DesDot = -(pi / 30) .* sin(t);
%}

%% System with No Observer
for i = 2:length(t)
    e = x(1, i - 1) - x1Des(i);
    eDot = x(2, i - 1) - x2Des(i);

    s = lambda * e + eDot;
    f = (mRod * g * (l / 2) * sin(x(1, i - 1))) / J_eff + (mBob * g * (l * rBob) * sin(x(1, i - 1))) / J_eff;
    b = 1 / J_eff;
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
    f = (mRod * g * (l / 2) * sin(xObs(1, i - 1))) / J_eff + (mBob * g * (l * rBob) * sin(xObs(1, i - 1))) / J_eff;
    b = 1 / J_eff;
    h = lambda * (xObs(2, i - 1) - x2Des(i)) + (f - x2DesDot(i));
    uObs(:, i) = -(1 / b) * (h + k * sign(s));
    if abs((uObs(i))) > controlSaturation
        u(i) = sign(uObs(i)) * controlSaturation;
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

%% Tracking Error Statistics

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

