function transformImuAcc(gravity, bias, qtVis, accImu, accVis, scale)
    N = size(accVis,1);
    R = qt_dircos(qtVis');
    accImuTransformed = accImu;
    for k = 1:length(accImu)
        accImuTransformed(k, :) = accImu(k,:) - (bias + R(:,:,k)*gravity)';
    end
    
    for k = 1:3
        figure;
        plot(accImuTransformed(:, k), 'b-');
        hold on
        plot(accVis(:, k), 'r-');
        hold off
        title('No Scaling');
    end
    
    for k = 1:3
        figure;
        plot(accImuTransformed(:, k), 'b-');
        hold on
        plot(abs(scale)*accVis(:, k), 'r-');
        hold off
        title('Scaling');
    end
end