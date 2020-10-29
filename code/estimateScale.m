function [scale, bias] = ...
    estimateScale(A,b,scale0,bias0,t)
% Estimation is performed in the frequency domain

fprintf('%s', repmat('-', 1, 60));
fprintf('\nFinal estimation in the frequency domain\n');
tic;

% Select valid range of frequencies [0Hz - 1.2Hz]
N = length(t);
fmax = 1.2;
fs = 1/mean(diff(t));
f = fs*(0:(N/2))/N;
freqRange = (f <= fmax);
fprintf('Upper limit for the frequencies = %.2f Hz\n',fmax);

% Optimize while enforcing gravity constraint: norm(g) = 9.82
% options = optimset(fmincon, 'Display', 'off');
x0 = [scale0; bias0];

Av = A*x0; % Visual accelerations
Ai = b;    % Inertial accelerations
Av = [Av(1:3:end) Av(2:3:end) Av(3:3:end)];
Ai = [Ai(1:3:end) Ai(2:3:end) Ai(3:3:end)];
for k = 1:3
    figure;
    plot(Av(:, k), 'r');
    hold on
    plot(Ai(:,k), 'b');
    hold off
end
x = fmincon(@(x)minFunc(x, A, b, freqRange), ...
    x0, [],[],[],[],[],[],[]);

scale = x(1);
bias = x(2:4);
fprintf('Error: %.5f\n', minFunc(x, A, b, freqRange));

fprintf('Finished in %.3f seconds\n', toc);
Av = A*x; % Visual accelerations
Ai = b;    % Inertial accelerations
Av = [Av(1:3:end) Av(2:3:end) Av(3:3:end)];
Ai = [Ai(1:3:end) Ai(2:3:end) Ai(3:3:end)];
AiSmooth(:,1) = smooth(Ai(:,1),35);
AiSmooth(:,2) = smooth(Ai(:,2),35);
AiSmooth(:,3) = smooth(Ai(:,3),35);
close all;
for k = 1:3
    figure(k+3)
    plot(Av(:, k), 'r');
    hold on
    plot(Ai(:, k), 'b');
    hold off
    figure(k+6);
    plot(Av(:,k), 'r');
    hold on
    plot(AiSmooth(:, k), 'g');
    hold off
end
display(scale);
display(bias)
end


function [c,ceq] = gravityConstraint(x)

c = [];
ceq = norm([x(2) x(3) x(4)])-9.80;

end


function f = minFunc(x, A, b, freqRange)

Av = A*x; % Visual accelerations
Ai = b;    % Inertial accelerations

Av = [Av(1:3:end) Av(2:3:end) Av(3:3:end)];
Ai = [Ai(1:3:end) Ai(2:3:end) Ai(3:3:end)];

Fv = abs(fft(Av));
Fi = abs(fft(Ai));

Fv = Fv(freqRange,:);
Fi = Fi(freqRange,:);

f = (Fv - Fi).^2;
f = sum(f(:));

end