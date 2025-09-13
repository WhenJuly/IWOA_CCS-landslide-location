function [y_global_best, x_global_best, Convergence_curve, recm1z] = POA(SearchAgents_no, Max_iter, lb, ub, dim, pop_init, nn, mm, recm1z, TTPS, Fs, k0, nnt, recW5)
    % POA 抛物线优化算法
    % 参数初始化
    lb = ones(dim, 1) .* lb';  % 变量的下限
    ub = ones(dim, 1) .* ub';  % 变量的上限

    %% 种群初始化
    X = pop_init;  % 初始种群
    fit = zeros(1, SearchAgents_no);  % 适应度值
    Convergence_curve = zeros(1, Max_iter);  % 记录每轮的最优适应度值

    % 计算初始适应度值（第一轮单独运行）
    dt = 1 / Fs;
    xx = sub2ind([nn, mm], round(X(1, :)), round(X(2, :)));  % 将位置转换为索引
%     % *****原始的成像*****
%     recm1z(xx(:)) = 0;  % 初始化成像值为零
%     for ii = 1:size(TTPS, 2)  % 遍历每个台站 ii
%         for jj = 1:size(TTPS, 2)  % 遍历每个台站 jj
%             ntp = round((TTPS(xx, ii) - TTPS(xx, jj)) / dt) + k0;  % 计算传播时间差
%             ntp = min(max(ntp, 1), nnt);  % 限制 ntp 范围
%             fit = fit + recW5(ntp, jj, ii)'.^2;  % 更新适应度值
%             recm1z(xx(:)) = recm1z(xx(:)) + recW5(ntp, jj, ii)'.^2;  % 更新成像值
%         end
%     end
    %******更改成像逻辑*******
    recm1z(xx(:)) = 0;  % recm1z 初始值为 0
    % 外层循环：组间（ii）
    for ii = 1:size(TTPS, 2)
        % 初始化组内乘积值（每组的结果初始值为 1）
        group_product = zeros(1, nn*mm);
        group_product(xx(:)) = 1; 
        % 内层循环：组内（jj）
        for jj = 1:size(TTPS, 2)
            % 计算传播时间差
            ntp = round((TTPS(xx, ii) - TTPS(xx, jj)) / dt) + k0;
            ntp = min(max(ntp, 1), nnt);  % 限制 ntp 范围
            % 组内相乘
            group_product(xx(:)) = group_product(xx(:)) .* recW5(ntp, jj, ii)';
            fit = fit + recW5(ntp, jj, ii)';
        end
        % 组间相加
        recm1z(xx(:)) = recm1z(xx(:)) + group_product(xx(:));
    end

    % 初始化全局最优
    [y_global_best, bestI] = max(fit);  % 找到全局最优适应度值
    x_global_best = X(:, bestI);  % 全局最优位置
    Convergence_curve(1) = y_global_best;  % 记录第一轮的最优适应度值

    %% 主循环（从 2 到 Max_iter）
    for t = 2:Max_iter
        % 随机选择一个粒子作为食物来源
        k = randi(SearchAgents_no);  % 随机选择一个粒子
        X_FOOD = X(:, k);  % 食物的位置
        F_FOOD = fit(k);  % 食物的适应度值

        %% PHASE 1: 捕食阶段（探索阶段）
        I = round(1 + rand());  % 随机因子
        X_new = X + rand(1, SearchAgents_no) .* (X_FOOD - I .* X);  % 向食物移动
        X_new = max(X_new, lb);  % 边界处理
        X_new = min(X_new, ub);
        fit_new = zeros(1, SearchAgents_no);  % 临时存储新适应度值

        % 计算新位置的适应度值
        xx_new = sub2ind([nn, mm], round(X_new(1, :)), round(X_new(2, :)));
        recm1z(xx_new(:)) = 0;
        for ii = 1:size(TTPS, 2)
            for jj = 1:size(TTPS, 2)
                ntp = round((TTPS(xx_new, ii) - TTPS(xx_new, jj)) / dt) + k0;
                ntp = min(max(ntp, 1), nnt);
                fit_new = fit_new + recW5(ntp, jj, ii)'.^2;
                recm1z(xx_new(:)) = recm1z(xx_new(:)) + recW5(ntp, jj, ii)'.^2;
            end
        end

        % 更新粒子和适应度值
        update_idx = fit_new >= fit;
        X(:, update_idx) = X_new(:, update_idx);
        fit(update_idx) = fit_new(update_idx);

        %% PHASE 2: 水面滑行阶段（开发阶段）
        X_new = X + 0.2 * (1 - t / Max_iter) .* (2 * rand(dim, SearchAgents_no) - 1) .* X;  % 随机扰动
        X_new = max(X_new, lb);  % 边界处理
        X_new = min(X_new, ub);

        % 计算新位置的适应度值
        xx_new = sub2ind([nn, mm], round(X_new(1, :)), round(X_new(2, :)));
        %******原始成像逻辑*******
%         recm1z(xx_new(:)) = 0;
%         fit_new = zeros(1, SearchAgents_no);  % 临时存储新适应度值
%         for ii = 1:size(TTPS, 2)
%             for jj = 1:size(TTPS, 2)
%                 ntp = round((TTPS(xx_new, ii) - TTPS(xx_new, jj)) / dt) + k0;
%                 ntp = min(max(ntp, 1), nnt);
%                 fit_new = fit_new + recW5(ntp, jj, ii)'.^2;
%                 recm1z(xx_new(:)) = recm1z(xx_new(:)) + recW5(ntp, jj, ii)'.^2;
%             end
%         end
        
        %******更改成像逻辑*******
        recm1z(xx_new(:)) = 0;  % recm1z 初始值为 0
        % 外层循环：组间（ii）
        for ii = 1:size(TTPS, 2)
            group_product = zeros(1, nn*mm);
            group_product(xx_new(:)) = 1; 
            for jj = 1:size(TTPS, 2)
                ntp = round((TTPS(xx_new, ii) - TTPS(xx_new, jj)) / dt) + k0;  % 计算传播时间差
                ntp = min(max(ntp, 1), nnt);  % 限制 ntp 范围
                fit_new = fit_new + recW5(ntp, jj, ii)';  % 更新适应度值  
                group_product(xx_new(:)) = group_product(xx_new(:)) .* recW5(ntp, jj, ii)';% 更新成像值
            end
            recm1z(xx_new(:)) = recm1z(xx_new(:)) + group_product(xx_new(:));  % 更新成像值
        end
        
        % 更新粒子和适应度值
        update_idx = fit_new >= fit;
        X(:, update_idx) = X_new(:, update_idx);
        fit(update_idx) = fit_new(update_idx);

        % 更新全局最优
        [current_best, bestI] = max(fit);
        if current_best > y_global_best
            y_global_best = current_best;
            x_global_best = X(:, bestI);
        end

        % 记录当前迭代的全局最优适应度值
        Convergence_curve(t) = y_global_best;
    
    end
end
