function json = getNewJsonWithScale(old_json, scale)
    json = old_json;
    for k = 1:length(json)
        json(k).scale = scale;
        similarity_transform = eye(4)*scale;
        similarity_transform(4, 4) = 1;
        world_from_camera = reshape(json(k).world_from_camera,[4, 4])';
        world_from_camera(:,4) = similarity_transform * world_from_camera(:,4);
        world_from_camera = world_from_camera';
        json(k).world_from_camera = world_from_camera(:);
        
    end
end