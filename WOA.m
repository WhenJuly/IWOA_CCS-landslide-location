function [y_global_best, x_global_best, Convergence_curve, recm1z] = WOA(SearchAgents_no, Max_iter, lb, ub, dim, pop_init, nn, mm, recm1z, TTPS, Fs, k0, nnt, recW5)
    % WOA 鲸鱼优化算法
    % 参数初始化
    lb = ones(dim, 1) .* lb';  % 变量的下限
    ub = ones(dim, 1) .* ub';  % 变量的上限

    %% 种群初始化
    X = pop_init;  % 初始种群
    X_new = pop_init;
    fit = zeros(1, SearchAgents_no);  % 适应度值
    Convergence_curve = zeros(1, Max_iter);  % 记录每轮的最优适应度值

    % 计算初始适应度值（第一轮单独运行）
    dt = 1 / Fs;
    xx = sub2ind([nn, mm], round(X(1, :)), round(X(2, :)));  % 将位置转换为索引
    xx_new = sub2ind([nn, mm], round(X_new(1, :)), round(X_new(2, :)));
%     %******原始成像逻辑*******
%     recm1z(xx(:)) = 0;  % 初始化成像值为零
%     for ii = 1:size(TTPS, 2)
%         for jj = 1:size(TTPS, 2)
%             ntp = round((TTPS(xx, ii) - TTPS(xx, jj)) / dt) + k0;  % 计算传播时间差
%             ntp = min(max(ntp, 1), nnt);  % 限制 ntp 范围
%             fit = fit + recW5(ntp, jj, ii)';  % 更新适应度值
%             recm1z(xx(:)) = recm1z(xx(:)) + recW5(ntp, jj, ii)';  % 更新成像值
%         end
%     end
        % 初始化 recm1z（假设 xx 的大小为 1x100）
%******更改成像逻辑*******
    recm1z(xx(:)) = 0;  % recm1z 初始值为 0
    recm1z_new = zeros(1, nn*mm);
    for ii = 1:size(TTPS, 2)
        group_product = zeros(1, nn*mm);
        group_product(xx(:)) = 1; 
        for jj = 1:size(TTPS, 2)
            ntp = round((TTPS(xx, ii) - TTPS(xx, jj)) / dt) + k0;
            ntp = min(max(ntp, 1), nnt);  % 限制 ntp 范围
            group_product(xx(:)) = group_product(xx(:)) .* recW5(ntp, jj, ii)';% 组内相乘
            fit = fit + recW5(ntp, jj, ii)';
        end
        % 组间相加
        recm1z(xx(:)) = recm1z(xx(:)) + group_product(xx(:));
    end
    recm1z_new(xx(:)) = recm1z(xx(:));
    % 初始化全局最优
    [y_global_best, bestI] = max(fit);  % 找到全局最优适应度值
    x_global_best = X(:, bestI);  % 全局最优位置
    Convergence_curve(1) = y_global_best;  % 记录第一轮的最优适应度值

    %% 主循环（从 2 到 Max_iter）
    for t = 2:Max_iter
        % 计算参数 a 和 a2
        a = 2 - t * (2 / Max_iter);  % 线性递减参数
        a2 = -1 + t * (-1 / Max_iter);  % 线性递减参数
%         w = cos(pi/2*(1-t/Max_iter)); % 收敛速度要变慢 可以拿来作对比
        w = 1;
        b = exp(cos(pi*(1 - t./Max_iter))); % 
%         b = 1;
        % 更新所有搜索代理的位置
        for i = 1:SearchAgents_no
            % 边界处理
            X(:, i) = min(max(X(:, i), lb), ub);

            % 随机参数
            r1 = rand();  % 随机数 [0, 1]
            r2 = rand();  % 随机数 [0, 1]
            A = 2 * a * r1 - a;  % 参数 A
            C = 2 * r2;  % 参数 C

            % 随机选择 p
            p = rand();  % 参数 p

            % 更新位置
            for j = 1:dim
                if p < 0.5
                    if abs(A) >= 1
                        % 随机选择一个搜索代理
                        rand_leader_index = randi(SearchAgents_no);
                        X_rand = X(:, rand_leader_index);
                        D_X_rand = abs(C * X_rand(j) - X(j, i));  % 计算距离
                        X(j, i) = w * X_rand(j) - A * D_X_rand;  % 更新位置
                    else
                        % 向领导者移动
                        D_Leader = abs(C * x_global_best(j) - X(j, i));  % 计算距离
                        X(j, i) = w * x_global_best(j) - A * D_Leader;  % 更新位置
                    end
                else
                    % 螺旋捕食行为
                    distance2Leader = abs(x_global_best(j) - X(j, i));
                    b = 1;  % 参数 b
                    l = (a2 - 1) * rand + 1;  % 参数 l
                    X(j, i) = b * distance2Leader * exp(b * l) * cos(l * 2 * pi) + w * x_global_best(j);  % 更新位置
                    X(j, i) = abs(X(j, i));
                end
                % 边界检查
                if X(j, i) < lb(j) || X(j, i) > ub(j)
                    X(j, i) = lb(j) + (ub(j) - lb(j)) * rand();  % 在边界内随机生成一个新位置
                end
            end
        %在这里加入小波依概率的小波变异策略
        X_new(:, i) = wavelet_mutation(X(:, i), lb(j), ub(j), t, Max_iter);
%         X_new = X(:, i);
        end

        % 计算新位置的适应度值
        fit_new = ones(1, SearchAgents_no);
        fit_group_product = ones(1, SearchAgents_no);
        fit_new_variation = ones(1, SearchAgents_no);  % 变异
        fit_group_product_variation = ones(1, SearchAgents_no);  % 变异
        xx_new = sub2ind([nn, mm], round(X(1, :)), round(X(2, :)));  % 将位置转换为索引
        xx_new_variation = sub2ind([nn, mm], round(X_new(1, :)), round(X_new(2, :)));  % 将变异位置转换为索引
        recm1z(xx_new(:)) = 0;  % 初始化成像值为零
        recm1z_new(xx_new_variation(:)) = 0;  % 变异初始化成像值为零
        % *****1.原始成像函数逻辑***********
%         for ii = 1:size(TTPS, 2)
%             for jj = 1:size(TTPS, 2)
%                 ntp = round((TTPS(xx_new, ii) - TTPS(xx_new, jj)) / dt) + k0;  % 计算传播时间差
%                 ntp = min(max(ntp, 1), nnt);  % 限制 ntp 范围
%                 fit_new = fit_new + recW5(ntp, jj, ii)';  % 更新适应度值
%                 recm1z(xx_new(:)) = recm1z(xx_new(:)) + recW5(ntp, jj, ii)';  % 更新成像值
%             end
%         end
        % *****2.更改成像函数逻辑***********
        % 原始的
        for ii = 1:size(TTPS, 2)
            group_product = zeros(1, nn*mm);
            group_product(xx_new(:)) = 1; 
            for jj = 1:size(TTPS, 2)
                ntp = round((TTPS(xx_new, ii) - TTPS(xx_new, jj)) / dt) + k0;  % 计算传播时间差
                ntp = min(max(ntp, 1), nnt);  % 限制 ntp 范围
                fit_group_product = fit_group_product .* recW5(ntp, jj, ii)';  % 更新适应度值  
                group_product(xx_new(:)) = group_product(xx_new(:)) .* recW5(ntp, jj, ii)';% 更新成像值
            end
            fit_new = fit_new +fit_group_product;
            recm1z(xx_new(:)) = recm1z(xx_new(:)) + group_product(xx_new(:));  % 更新成像值
        end
        % 小波变异策略 variation
        for ii = 1:size(TTPS, 2)
            group_product_variation = zeros(1, nn*mm);
            group_product_variation(xx_new_variation(:)) = 1; 
            for jj = 1:size(TTPS, 2)
                ntp = round((TTPS(xx_new_variation, ii) - TTPS(xx_new_variation, jj)) / dt) + k0;  % 计算传播时间差
                ntp = min(max(ntp, 1), nnt);  % 限制 ntp 范围
                fit_group_product_variation = fit_group_product_variation .* recW5(ntp, jj, ii)';  % 更新适应度值  
                group_product_variation(xx_new_variation(:)) = group_product_variation(xx_new_variation(:)) .* recW5(ntp, jj, ii)';% 更新成像值
            end
            fit_new_variation = fit_new_variation +fit_group_product_variation;
            recm1z_new(xx_new_variation(:)) = recm1z_new(xx_new_variation(:)) + group_product_variation(xx_new_variation(:));  % 更新成像值
        end
        
        % 贪婪操作：比较变异前后的适应度值并更新 fit_new 和 X
        for i = 1:SearchAgents_no
            if fit_new_variation(i) > fit_new(i)  % 如果变异后的适应度值更大
                fit_new(i) = fit_new_variation(i);  % 更新适应度值
                X(:, i) = X_new(:, i);  % 更新位置坐标
            end
        end
        
        % 更新全局最优
        [current_best, bestI] = max(fit_new);
        if current_best > y_global_best
            y_global_best = current_best;
            x_global_best = X(:, bestI);
        end

        % 记录当前迭代的全局最优适应度值
        Convergence_curve(t) = y_global_best;
    end
end
