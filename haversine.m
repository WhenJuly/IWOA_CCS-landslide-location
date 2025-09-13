% 计算两点间大圆距离的 haversine 函数
function distance = haversine(lon1, lat1, lon2, lat2)
    % 将经纬度从度转换为弧度
    R = 6371; % 地球半径，单位为公里
    lon1 = deg2rad(lon1);
    lat1 = deg2rad(lat1);
    lon2 = deg2rad(lon2);
    lat2 = deg2rad(lat2);
    
    % Haversine公式
    dphi = lat2 - lat1;
    dlambda = lon2 - lon1;
    a = sin(dphi / 2).^2 + cos(lat1) .* cos(lat2) .* sin(dlambda / 2).^2;
    c = 2 * atan2(sqrt(a), sqrt(1 - a));
    
    % 计算距离
    distance = R * c;
end
