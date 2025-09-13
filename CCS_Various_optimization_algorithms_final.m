clear; clc; close all; % 清除工作区的所有变量，清空命令行窗口
%% 导入数据与计算走时表
filepath = '...\CCS_landslide location\iceland'; 
Coordfile = fullfile(filepath,'\Coordination.txt');
Coordinates = load(Coordfile);  % 台站经纬度 [Lon, Lat]
Coordinates = Coordinates([7 4 3 6 1 2],:);
Sourcefile = fullfile(filepath,'source.txt');
source = load(Sourcefile);      % 震源初始位置 [Source_Lon, Source_Lat]
Prossesed_datafile = fullfile(filepath,'\Aligned_Waveforms.mat'); 
load(Prossesed_datafile);  % 处理后的波形
SNR_data = fullfile(filepath,'\SNR_list.txt');
SNR_list = load(SNR_data);
SNR_list_w = SNR_list([7 4 3 6 1 2],:);
recvz = cell2mat(Aligned_Waveforms([7 4 3 6 1 2])');

% %****用SNR调节波形****
% 归一化 SNR 到 [0, 1]
SNR_norm = (SNR_list_w - min(SNR_list_w)) / (max(SNR_list_w) - min(SNR_list_w));
% 自定义 Sigmoid 函数
w = 1 ./ (1 + exp(-10 * (SNR_norm - 0.5)));  % 调整到 [0, 1]
% w = [1;1;1;1;1;1];
recvz = recvz .* w';
% ****Envelope****
for i = 1:size(recvz,2)
     [Svalue(:,i),~] = envelope(recvz(:,i),500,'rms');
end
recvz = Svalue;
% ****STA/LTA****
% for i = 1:size(recvz,2)
%     Svalue(:,i) = stalta(recvz(:,i),500,100,0.1); 
% end
% recvz = Svalue;
% ****extractSSTRidge****
% Fs = 100;
% recvz_sst = extractSSTRidge(recvz, Fs);
% recvz = recvz_sst;
% ****近似熵 提取近似熵特征****
% m = 2; % 嵌入维数
% r = 0.2 * std(recvz(:)); % 相似容限
% recvz_apen = extractApEn(recvz, m, r);
% ****STFT脊线****
% Fs = 100;
% recvz = extractSTFTRidge(recvz, Fs);
% 计算走时表
Resolution = 50;  % 网格分辨率
%step 1
Lon_lb = mean(Coordinates(:, 2)) - 2;
Lon_ub = mean(Coordinates(:, 2)) + 2;
Lat_lb = mean(Coordinates(:, 1)) - 2;
Lat_ub = mean(Coordinates(:, 1)) + 2;

lon_range = linspace(Lon_lb, Lon_ub, Resolution);
lat_range = linspace(Lat_lb, Lat_ub, Resolution);
[lon_grid, lat_grid] = meshgrid(lon_range, lat_range);
Vs = 2.95;  % 假设传播速度    3 km/s 2.95 2.327(阿尔卑斯山脉)
numStations = size(Coordinates, 1);

% 计算走时表
numGridPoints = numel(lon_grid);
TravelTimeTable = zeros(numStations, numGridPoints);
for i = 1:numStations
    distances = haversine(lon_grid(:), lat_grid(:), Coordinates(i, 2), Coordinates(i, 1));
    TravelTimeTable(i, :) = distances / Vs;
end

%% colormap
baseCM = {[240,240,240; 114,148,197; 131,158,192; 155,170,197; 168,184,197; ...
            181,189,189; 200,210,204; 214,225,201; 228,235,203; 241,248,207; ...
            250,245,198; 253,229,178; 247,211,157; 241,193,141; 244,175,133; ...
            248,152,111; 240,136,101; 240,114,89; 228,94,62; 227,68,54]./255};
Cmap=baseCM{1};
Ci=1:size(Cmap,1);Cq=linspace(1,size(Cmap,1),200);% 插值到200行
Cmap=[interp1(Ci,Cmap(:,1),Cq,'linear')',...
     interp1(Ci,Cmap(:,2),Cq,'linear')',...
     interp1(Ci,Cmap(:,3),Cq,'linear')'];  

Fs = 100;
dt = 1/Fs;  % 设置时间步长为 0.001秒
s0 = size(Coordinates, 1);  % 设置震源台站的数量
k0 = length(recvz);  % 数据长度
nn = Resolution;  % 设置地震模型的水平网格大小为 901
mm = Resolution;  % 设置地震模型的垂直网格大小为 301
TTPS = TravelTimeTable';  % 将传输时间表转换为 (nn*mm) x s0 矩阵
%% cross-correlation 计算互相关
recW5 = zeros(2*(k0-1)+1, s0, s0);  % 初始化互相关结果矩阵
for ii = 1:s0  % 遍历所有台站 ii
    for jj = 1:s0  % 遍历所有台站 jj
        recW5(:, jj, ii) = xcorr(recvz(:, ii), recvz(:, jj));  % 计算并存储台站 ii 和 jj 之间的互相关
    end
end
clear recvz;  % 清除接收到的地震波形数据变量 recvz

nnt = size(recW5, 1);  % 计算互相关结果的长度
recm1z = zeros(1, nn*mm);  % 初始化成像值矩阵（最终的震源位置值）

%% 参数设置
SearchAgents_no = 50;  % 种群数量
Max_iter = 200;         % 最大迭代次数
lb = [1, 1];            % 下界
ub = [nn, mm];          % 上界
dim = 2;                % 搜索空间维度
xmin = diag([1, 1]) * ones(dim, SearchAgents_no);  % 粒子位置的最小值
xmax = diag([nn, mm]) * ones(dim, SearchAgents_no);  % 粒子位置的最大值
init_method = "pwclm";  % 选择初始化方法 ('pwclm', 'sine', 'tent')
pop_init = initialization(xmin, xmax, dim, SearchAgents_no, init_method);
tic;  % 记录程序开始的时间

%% ESOA 优化过程 各种优化算法途径
[y_global_best, x_global_best, Convergence_curve, recm1z] = WOA_origin(SearchAgents_no, Max_iter, lb, ub, dim, pop_init, nn, mm, recm1z, TTPS, Fs, k0, nnt, recW5);
% [y_global_best, x_global_best, Convergence_curve, recm1z] = WOA(SearchAgents_no, Max_iter, lb, ub, dim, pop_init, nn, mm, recm1z, TTPS, Fs, k0, nnt, recW5);
% [y_global_best, x_global_best, Convergence_curve, recm1z] = POA(SearchAgents_no, Max_iter, lb, ub, dim, pop_init, nn, mm, recm1z, TTPS, Fs, k0, nnt, recW5);
% [y_global_best, x_global_best, Convergence_curve, recm1z] = ESOA(SearchAgents_no, Max_iter, lb, ub, dim, pop_init, nn, mm, recm1z, TTPS, Fs, k0, nnt, recW5);
% [y_global_best, x_global_best, Convergence_curve, recm1z] = SSA(SearchAgents_no, Max_iter, lb, ub, dim, pop_init, nn, mm, recm1z,TTPS, Fs, k0, nnt, recW5);
%% 最优位置映射到网格
% x_global_best 是二维坐标
bestI = round(x_global_best(1));  % 行索引
bestJ = round(x_global_best(2));  % 列索引
bestX = [bestI, bestJ];

toc;  % 计算并显示程序运行时间
%% 绘制成像结果
recm1z = reshape(recm1z, nn, mm);  % 将成像结果重塑为 nn x mm 矩阵

ind = find(recm1z == max(recm1z(:)));  % 找到最大成像值的位置
[xgmax, ygmax] = ind2sub([nn, mm], ind);  % 将线性索引转换为二维坐标
[source_Lon, source_Lat] = find(recm1z==max(max(recm1z)));
figure  % 显示成像结果的 2D 切片
min_val = min(recm1z(recm1z~=0));
recm1z(recm1z ~= 0) = recm1z(recm1z ~= 0) - min_val;
imagesc((recm1z./ max(recm1z(:)))');  % 绘制成像图，归一化显示
colormap(Cmap)  % 设置色图为  Turbo jet
hc = colorbar('horiz');  % 显示颜色条
% hc = colorbar('verti');
minn = 0; maxx = 1;  % 设置颜色条的范围
caxis([minn maxx])  % 设置颜色条的显示范围
set(hc, 'ytick', [minn maxx], 'yticklabel', {'0' '1'}, 'fontsize', 14)  % 设置颜色条的标签
% set(hc, 'pos', [0.918 0.145 0.02 0.78]);  % 设置颜色条的位置
set(hc, 'pos', [0.125 0.938 0.78 0.02 ]);  % 设置颜色条的位置

%%######输出缩略图时用##########
num_ticks = 2; % 保持5个刻度
xtick_values = linspace(Lon_lb, Lon_ub, num_ticks);
ytick_values = linspace(Lat_lb, Lat_ub, num_ticks);

% 四舍五入为整数
xtick_labels = round(xtick_values);
ytick_labels = round(ytick_values);

% 应用设置
set(gca, 'fontsize', 23);
set(gca, 'xtick', linspace(1, Resolution, num_ticks), ...
         'xticklabel', xtick_labels, ...
         'ytick', linspace(1, Resolution, num_ticks), ...
         'yticklabel', ytick_labels);

set(gcf, 'pos', [500 500 700 700])  % 设置图形窗口的大小
% set(gcf, 'pos', [100 100 550 700])  % 设置图形窗口的大小
xlabel('Longitude', 'fontsize', 24);
ylabel('Latitude', 'fontsize', 24);
set(gca, 'YDir', 'normal'); 
% axis([1 901 1 301])  % 设置坐标轴的显示范围
hold on
% *******网格生成导入Arcgis********* 可以选择不要这一步
[xGrid, yGrid] = meshgrid(1:nn, 1:mm);  % 生成网格坐标
% 绘制垂直网格线（偏移 0.5 个像素点）
for i = 1:nn+1  % 需要绘制 nn+1 条垂直线
    x = (i - 0.5);  % 偏移 0.5 个像素点
    plot([x, x], [0.5, mm+0.5], 'Color', [0.5, 0.5, 0.5], 'LineWidth', 0.5);  % 绘制垂直网格线
end
% 绘制水平网格线（偏移 0.5 个像素点）
for j = 1:mm+1  % 需要绘制 mm+1 条水平线
    y = (j - 0.5);  % 偏移 0.5 个像素点
    plot([0.5, nn+0.5], [y, y], 'Color', [0.5, 0.5, 0.5], 'LineWidth', 0.5);  % 绘制水平网格线
end

% 绘制真实震源位置
realSourceLonIndex = round((source(2) - Lon_lb) / (Lon_ub - Lon_lb) * (Resolution - 1)) + 1;
realSourceLatIndex = round((source(1) - Lat_lb) / (Lat_ub - Lat_lb) * (Resolution - 1)) + 1;
scatter(realSourceLonIndex, realSourceLatIndex, 200, 'r', 'p', 'filled', 'DisplayName', 'Real Source');
legend;

% 绘制估计震源位置
estimatedSourceLon = lon_range(bestI);  
estimatedSourceLat = lat_range(bestJ);  
scatter(bestI, bestJ, 200, 'r', 'p', 'filled','MarkerEdgeColor', [0 0 0], 'DisplayName', 'Estimated Source');

% 绘制检波器的位置
realStationsLonIndex = round((Coordinates(:,2) - Lon_lb) / (Lon_ub - Lon_lb) * (Resolution - 1)) + 1; % 经度映射到 X 轴
realStationsLatIndex = round((Coordinates(:,1) - Lat_lb) / (Lat_ub - Lat_lb) * (Resolution - 1)) + 1; % 纬度映射到 Y 轴
scatter(realStationsLonIndex, realStationsLatIndex, 100, '^', 'MarkerEdgeColor','k','MarkerFaceColor','k','DisplayName','Stations');

hold off;

figure;  % 绘制适应度值的收敛曲线
h = plot(max(Convergence_curve) - Convergence_curve);
% h = plot(Convergence_curve);
set(h(1),'LineWidth',2.2, 'LineStyle','-.', 'Color', [0.6350 0.0780 0.1840]);
xlabel('Iteration');
ylabel('Fitness Value');
title('Convergence Curve');
Convergence_curve = (max(Convergence_curve) - Convergence_curve)./max(Convergence_curve);
Convergence_curve_new = Convergence_curve';
distance_to_real = haversine(estimatedSourceLon, estimatedSourceLat, source(2), source(1));
fprintf('Distance to real source: %.4f km\n', distance_to_real);
fprintf('%.4f %.4f\n', estimatedSourceLon,estimatedSourceLat);

