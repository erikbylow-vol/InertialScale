function [accImu, angImu, timeImu, userAccImu] = readInertial(dataset)
%
% Read inertial measurements and timestamps from the input files. 
% See the README.txt for more information about the input format.
%
% INPUT:    dataset : Name of the dataset (string)
%
% OUTPUT:   accImu  : Inertial accelerations [m/s^2] (Mx3 matrix)
%           angImu  : Inertial angular velocities [rad/s] (Mx3 matrix)
%           timeImu : Timestamps in seconds (Mx1 vector)
%

fprintf('%s', repmat('-', 1, 60));
fprintf('\nReading inertial data (%s)\n', dataset);

% Read accelerometer readings
if exist([dataset 'accelerometer.txt'],'file')
    fileID = fopen([dataset 'accelerometer.txt']);
    acc = textscan(fileID,'%f %f %f %f');
    fclose(fileID);
    timeAcc = 1e-9*acc{1};
    accImu = [acc{2} acc{3} acc{4}];
else
   error('accelerometer.txt not found! Check inputs.'); 
end

% Read accelerometer readings
if exist([dataset 'user_accelerometer.txt'],'file')
    fileID = fopen([dataset 'user_accelerometer.txt']);
    acc = textscan(fileID,'%f %f %f %f');
    fclose(fileID);
    timeAcc = 1e-9*acc{1};
    userAccImu = [acc{2} acc{3} acc{4}];
else
    userAccImu = [];
end

% Read gyroscope readings
if exist([dataset 'gyroscope.txt'],'file')
    fileID = fopen([dataset 'gyroscope.txt']);
    gyr = textscan(fileID,'%f %f %f %f');
    fclose(fileID);
    timeGyr = 1e-9*gyr{1};
    gyrImu = [gyr{2} gyr{3} gyr{4}];
else
   error('gyroscope.txt not found! Check inputs.'); 
end

dtAcc = mean(diff(timeAcc));
dtGyr = mean(diff(timeGyr));
fprintf('Sampling rate of the accelerometer = %.2f Hz\n', 1/(dtAcc));
fprintf('Sampling rate of the gyroscope = %.2f Hz\n', 1/(dtGyr));

% Resample gyroscope readings to match the sampling of the accelerometer
t1 = min([timeAcc(1),timeGyr(1)]);
timeImu = timeAcc - t1;
angImu = interp1(timeGyr-t1,gyrImu,timeImu,'linear','extrap');

end

