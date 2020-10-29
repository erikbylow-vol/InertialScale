function json = getNewJsonWithScale(old_json, scale)
    json = old_json;
    for k = 1:length(json)
        json(k).scale = scale;
    end
end