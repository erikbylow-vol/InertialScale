function [A,b,scale,bias] = initializeEstimates(accVis,qtVis,accImu)

fprintf('%s', repmat('-', 1, 60));
fprintf('\nFinding initial estimates by solving a linear system\n');
tic;

N = size(accVis,1);
R = qt_dircos(qtVis');

% Construct linear system
A = zeros(3*N, 4);
b = zeros(3*N, 1);

for n = 0:N-1
    A(3*n+1:3*n+3,1) = accVis(n+1,:)';
    A(3*n+1:3*n+3,2:4) = eye(3);
    b(3*n+1:3*n+3,1) = accImu(n+1,:)';
end

x = A\b;
scale = x(1);
bias = x(2:4);

fprintf('scale0 = %.4f\n', scale);
fprintf('bias0 = [%.4f, %.4f, %.4f]\n',bias);
fprintf('Finished in %.3f seconds\n', toc);

end

