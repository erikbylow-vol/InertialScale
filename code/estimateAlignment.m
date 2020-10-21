function [Rs,td,bg] = estimateAlignment(qtVis,tVis,angImu,tImu)
%
% Estimate temporal and spatial alignment between the camera and IMU.
% Gyroscope bias is also estimated in the process.
%
% INPUT:    qtVis   : Visual orientations (Nx4 matrix)
%           tVis    : Visual timestamps in seconds (Nx1 vector)
%           angImu  : Inertial angular velocities [rad/s] (Mx3 matrix)
%           tImu    : Inertial timestamps in seconds (Mx1 vector)
%
% OUTPUT:   Rs      : Rotation between the camera and IMU (3x3 matrix)
%           td      : Time offset between the camera and IMU (scalar)
%           bg      : Gyroscope bias [rad/s] (1x3 vector)
%

fprintf('%s', repmat('-', 1, 60));
fprintf('\nTemporal and spatial alignment\n');
tic;

% Use only time period which all sensors have values
timeStop = min([tVis(end) tImu(end)]);
tVis = tVis(tVis <= timeStop);
tImu = tImu(tImu <= timeStop);

qtVis = qtVis(1:length(tVis),:);
angImu = angImu(1:length(tImu),:);

% Upsample visual-data to match the sampling of the IMU
t = tImu;
dt = mean(diff(t));
qtVis = interp1(tVis,qtVis,t,'linear','extrap'); % Consider using SLERP
qtVis = qtVis./sqrt(sum(qtVis.^2,2));
% Compute visual angular velocities
qtDiffs = diff(qtVis);
qtDiffs = [qtDiffs; qtDiffs(end,:)]';
angVis = -(2/dt)*qt_mul(qtDiffs, qt_inv(qtVis'));
angVis = angVis(2:4,:)';
% Smooth angular velocities
angVis(:,1) = smooth(angVis(:,1),15);
angVis(:,2) = smooth(angVis(:,2),15);
angVis(:,3) = smooth(angVis(:,3),15);
angImu(:,1) = smooth(angImu(:,1),15);
angImu(:,2) = smooth(angImu(:,2),15);
angImu(:,3) = smooth(angImu(:,3),15);

for k = 1:3
    figure(k);
    plot(angVis(:, k), 'r');
    hold on
    plot(angImu(:, k), 'b');
    hold off
end

gRatio = (1 + sqrt(5)) / 2;
tolerance = 0.0001;

maxOffset = 0.5;
a = -maxOffset;
b = maxOffset;

c = b - (b - a) / gRatio;
d = a + (b - a) / gRatio;

iter = 0;

while abs(c - d) > tolerance
     % Evaluate function at f(c) and f(d)
    [Rsc,biasc,fc] = solveClosedForm(angVis,angImu,t,c);
    [Rsd,biasd,fd] = solveClosedForm(angVis,angImu,t,d);
    fprintf("cost %f and limit %d\n", fc, c);
    fprintf("cost %f and limit %d\n", fd, d);    
    if fc < fd
        b = d;
        Rs = Rsc;
        bg = biasc;
    else
        a = c;
        Rs = Rsd;
        bg = biasd;
    end
    
    c = b - (b - a) / gRatio;
    d = a + (b - a) / gRatio;
    
    iter = iter + 1;
end

td = (b + a) / 2;

% timeInterval = linspace(a,b, 2000);
% minCost = inf;
% td = a;
% 
% for k = 1:length(timeInterval)
%     [Rsc,biasc,fc] = solveClosedForm(angVis,angImu,t,timeInterval(k));
%     if (fc < minCost)
%         minCost = fc;
%         td = timeInterval(k);
%         Rs = Rsc;
%         bg = biasc;
%         fprintf("Cost: %d\n", fc);
%     end
% end

fprintf('Golden-section search (%.0f iterations)\n', iter);
fprintf('Finished in %.3f seconds\n', toc);
angVis = interp1(t-td,angVis,t,'linear','extrap');
angImu = angImu*Rs + repmat(bg, size(angImu,1),1);
for k = 1:3
    figure(k+3);
    plot(angVis(:, k), 'r');
    hold on
    plot(angImu(:, k), 'b');
    hold off
end
end



function [Rs,bias,f] = solveClosedForm(angVis,angImu,t,td)
%
% Finds the relative rotation between the camera and IMU when using the 
% provided time offset td. Gyroscope bias is estimated in the process.
%
% INPUT:    angVis  : Visual angular velocities [rad/s] (Mx3 matrix)
%           angImu  : Inertial angular velocities [rad/s] (Mx3 matrix)
%           t       : Timestamps in seconds
%           td      : Time offset in seconds
%
% OUTPUT:   Rs      : Rotation between the camera and IMU (3x3 matrix)
%           bias    : Gyroscope bias [rad/s] (1x3 vector)
%           f       : Function value (sum of squared differences)
%

% Adjust visual angular velocities based on current offset
angVis = interp1(t-td,angVis,t,'linear','extrap');
squaredVelocity = sum(angVis.^2, 2);
K = 150;
[~, indMax] = maxk(squaredVelocity, K);
indKeep = true(size(angVis,1),1);
indKeep(indMax) = false;

N = size(angVis, 1);

% Compute mean vectors
meanImu = repmat(mean(angImu(indKeep,:)), sum(indKeep), 1);
meanVis = repmat(mean(angVis(indKeep,:)), sum(indKeep), 1);

% Compute centralized point sets
P = angImu(indKeep,:) - meanImu;
Q = angVis(indKeep,:) - meanVis;

% Singular value decomposition
[U,S,V] = svd(P'*Q);

% Ensure a right handed coordinate system and correct if necessary
C = eye(3);
if (det(V*U') < 0)
    C(3,3) = -1;
end

Rs = V*C*U';

% Find the translation, which is the gyroscope bias
bias = mean(angVis(indKeep,:)) - mean(angImu(indKeep,:))*Rs;

% Residual
D = angVis(indKeep, :) - (angImu(indKeep,:)*Rs + repmat(bias,sum(indKeep),1));
errors = sort(sum(D.^2,2));
% f = sum(D(:).^2);
% f = sum(errors(1:end-15));
f = sum(errors);

end

