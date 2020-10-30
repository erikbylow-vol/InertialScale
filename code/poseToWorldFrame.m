function posVis = poseToWorldFrame(posVis, qtVis)
for k = 1:length(posVis)
    R = qt_dircos(qtVis(k,:)');
    tmp = posVis(k,:)';
    tmp = R'*tmp;
    posVis(k,:) = tmp';
end
end