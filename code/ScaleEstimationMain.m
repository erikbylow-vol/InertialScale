% Finds the absolute scale of the visual reconstruction given the camera 
% poses and inertial measurements. The code is based on the paper: 
% 
% Mustaniemi J., Kannala J., Särkkä S., Matas J., Heikkilä J.
% "Inertial-Based Scale Estimation for Structure from Motion on Mobile 
% Devices", International Conference on Intelligent Robots and Systems 
% (IROS), 2017

% See the README.txt for more information on how to use the code.
tic
clear all
close all
addpath quaternions
METHOD_NAMES = ["arkit"];
folder = "charuco_imu";

% base_path = "/Users/erikbylow/Code/Reconstruction/reports/moby/imu_test/"
base_path = "/home/erikbylow/Code/Reconstruction/reports/moby/cameras/" + folder + "/";
% base_path = "/home/erikbylow/Code/myFork/InertialScale/data/";

suffix = "2*";
% Run the scale estimation using the following dataset
datasets = dir(strcat(base_path, suffix));
for method = METHOD_NAMES
    fileID = fopen(method + "_" + folder + ".txt", 'w');
    fname = base_path + method + ".json";
    val = jsondecode(fileread(fname));
    for k = 1:length(datasets)
        capture = datasets(k).name;
        % Read camera poses and timestamps
        dataset = sprintf('%s%s/%s/', base_path, capture, method);
        [posVis,qtVis,tVis,scaleGT] = readVisual(dataset);
        posVis = poseToWorldFrame(posVis, qtVis);

        % Read inertial measurements and timestamps
        [accImu,angImu,tImu, userAccImu] = readInertial(dataset);

        % Kalman filtering and RTS smoothing
        accVis = kalmanRTS(posVis,tVis);

        % Temporal and spatial alignment of the camera and IMU
        [Rs,td,bg] = estimateAlignment(qtVis,tVis,angImu,tImu);
        [accVis,qtVis,userAccImu,t] = alignCameraIMU(accVis,qtVis,tVis,userAccImu,tImu,Rs,td);
    %     
    %     close all;
    %     for l = 1:3
    %         figure(l);
    %         plot(accVis(:, l), 'r-');
    %         hold on
    %         plot(userAccImu(:, l), 'b-');
    %         hold off
    %     end

        % Transform visual accelerations from world frame to local frame
        accVis = qt_rot(qtVis',accVis')';

    %     for l = 1:3
    %         figure(l+3);
    %         plot(accVis(:, l), 'r-');
    %         hold on
    %         plot(userAccImu(:, l), 'b-');
    %         hold off
    %     end

        % Find initial estimates for the scale, gravity and bias by solving
        % a linear system of equations Ax = b
        [A,b,s0,b0] = initializeEstimates(accVis,qtVis,userAccImu);

        % Perform final estimation in the frequency domain while enforcing
        % gravity constraint: norm(g) = 9.81
        [scale, bias] = estimateScale(A,b,s0,b0,t);
        fprintf('%s', repmat('-', 1, 60));
        fprintf('\nFinal estimates\n');
        if (scaleGT > 0)
            scaleErr = 100*abs(scale-scaleGT)/scaleGT;
            fprintf('scale = %.4f (error = %.1f%%)\n',scale,scaleErr);
        else
            fprintf('scale = %.4f ', scale);
        end
        fprintf('bias = [%.4f, %.4f, %.4f]\n',bias);
        fprintf('td = %.4f seconds\n',td);
        fprintf('Rs = [%.2f %.2f %.2f; %.2f %.2f %.2f; %.2f %.2f %.2f]\n', Rs');
        fprintf('\n');
        fprintf('%s', dataset);
        toc
        close all
        fprintf(fileID,'%s %f\n', capture, scale);
    end
    fprintf('bias = [%.4f, %.4f, %.4f]\n',bias);
    fprintf('td = %.4f seconds\n',td);
    fprintf('Rs = [%.2f %.2f %.2f; %.2f %.2f %.2f; %.2f %.2f %.2f]\n', Rs');
    fprintf('\n');
    toc
    close all
end
