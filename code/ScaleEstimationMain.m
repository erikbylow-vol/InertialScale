% Finds the absolute scale of the visual reconstruction given the camera 
% poses and inertial measurements. The code is based on the paper: 
% 
% Mustaniemi J., Kannala J., S�rkk� S., Matas J., Heikkil� J.
% "Inertial-Based Scale Estimation for Structure from Motion on Mobile 
% Devices", International Conference on Intelligent Robots and Systems 
% (IROS), 2017

% See the README.txt for more information on how to use the code.
tic
clear all
close all
addpath quaternions

base_path = "/Users/erikbylow/Code/Reconstruction/reports/moby/imu_test/"
suffix = "2*"
% Run the scale estimation using the following dataset
datasets = dir(strcat(base_path, suffix));

for k = 1:length(datasets)
    % Read camera poses and timestamps
    dataset = sprintf('%s%s/', base_path, datasets(k).name);
    [posVis,qtVis,tVis,scaleGT] = readVisual(dataset);

    % Read inertial measurements and timestamps
    [accImu,angImu,tImu] = readInertial(dataset);

    % Kalman filtering and RTS smoothing
    accVis = kalmanRTS(posVis,tVis);

    % Temporal and spatial alignment of the camera and IMU
    [Rs,td,bg] = estimateAlignment(qtVis,tVis,angImu,tImu);
    [accVis,qtVis,accImu,t] = alignCameraIMU(accVis,qtVis,tVis,accImu,tImu,Rs,td);

    % Transform visual accelerations from world frame to local frame
    accVis = qt_rot(qtVis',accVis')';

    % Find initial estimates for the scale, gravity and bias by solving
    % a linear system of equations Ax = b
    [A,b,s0,g0,b0] = initializeEstimates(accVis,qtVis,accImu);

    % Perform final estimation in the frequency domain while enforcing
    % gravity constraint: norm(g) = 9.81
    [scale,gravity,bias] = estimateScale(A,b,s0,g0,b0,t);

    fprintf('%s', repmat('-', 1, 60));
    fprintf('\nFinal estimates\n');
    if (scaleGT > 0)
        scaleErr = 100*abs(scale-scaleGT)/scaleGT;
        fprintf('scale = %.4f (error = %.1f%%)\n',scale,scaleErr);
    else
        fprintf('scale = %.4f ', scale);
    end
    fprintf('gravity = [%.4f, %.4f, %.4f]\n', gravity);
    fprintf('bias = [%.4f, %.4f, %.4f]\n',bias);
    fprintf('td = %.4f seconds\n',td);
    fprintf('Rs = [%.2f %.2f %.2f; %.2f %.2f %.2f; %.2f %.2f %.2f]\n', Rs');
    fprintf('\n');
    toc
    close all
end
