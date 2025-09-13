function [y_global_best, x_global_best, Convergence_curve, recm1z] = ESOA(SearchAgents_no, Max_iter, lb, ub, dim, pop_init, nn, mm, recm1z, TTPS, Fs, k0, nnt, recW5)
    % ESOA 参数设置
    beta1 = 0.9;
    beta2 = 0.99;

    % 初始化种群和适应度值
    x = pop_init;  % 种群位置
    Convergence_curve = zeros(1, Max_iter);  % 记录每轮迭代的全局最优适应度值
    w = 2 * rand(dim, SearchAgents_no) - 1;  % 初始化权重
    m = zeros(dim, SearchAgents_no);  % 初始化动量
    v = zeros(dim, SearchAgents_no);  % 初始化二阶动量
    y = zeros(1,SearchAgents_no);  % 适应度值
    x_hist_best = x;  % 每个粒子的历史最优位置
    y_hist_best = y;  % 每个粒子的历史最优适应度值
    g_hist_best = zeros(dim, SearchAgents_no);  % 每个粒子的历史最优梯度

    % 计算初始适应度值 第一轮单独做
    dt = 1 / Fs;
    xx = sub2ind([nn, mm], round(x(1, :)), round(x(2, :)));  % 将位置转换为索引
%     %******原始成像逻辑*******
%     recm1z(xx(:)) = 0;  % 初始化成像值为零
%     for ii = 1:size(TTPS, 2)  % 遍历每个台站 ii
%         for jj = 1:size(TTPS, 2)  % 遍历每个台站 jj
%             ntp = round((TTPS(xx, ii) - TTPS(xx, jj)) / dt) + k0;  % 计算传播时间差
%             ntp = min(max(ntp, 1), nnt);  % 限制 ntp 范围
%             y = y + recW5(ntp, jj, ii)'.^2;  % 更新适应度值
%             recm1z(xx(:)) = recm1z(xx(:)) + recW5(ntp, jj, ii)'.^2;  % 更新成像值
%         end
%     end
    
    %******更改成像逻辑*******
    recm1z(xx(:)) = 0;  % recm1z 初始值为 0
    % 外层循环：组间（ii）
    for ii = 1:size(TTPS, 2)
        group_product = zeros(1, nn*mm);
        group_product(xx(:)) = 1; 
        for jj = 1:size(TTPS, 2)
            ntp = round((TTPS(xx, ii) - TTPS(xx, jj)) / dt) + k0;
            ntp = min(max(ntp, 1), nnt);  % 限制 ntp 范围
            group_product(xx(:)) = group_product(xx(:)) .* recW5(ntp, jj, ii)';
            y = y + recW5(ntp, jj, ii)';
        end
        recm1z(xx(:)) = recm1z(xx(:)) + group_product(xx(:));
    end

    % 初始化全局最优
    [y_global_best, bestI] = max(y);  % 找到全局最优适应度值
    x_global_best = x(:, bestI);  % 全局最优位置
    Convergence_curve(1) = y_global_best;  % 记录第一轮的最优适应度值

    % Main loop 剩余轮次
    for l = 2:Max_iter
        % 更新粒子位置
        for i = 1:SearchAgents_no
            % 计算梯度和其他更新逻辑
            p_y = sum(w(:, i) .* x(:, i));
            p = p_y - y(i);
            g_temp = p .* x(:, i);

            % 个体方向
            p_d = x_hist_best(:, i) - x(:, i);
            f_p_bias = y_hist_best(i) - y(i);
            p_d = p_d .* f_p_bias;
            p_d = p_d ./ ((sum(p_d) + eps) .* (sum(p_d) + eps));

            d_p = p_d + g_hist_best(:, i);

            % 群体方向
            c_d = x_global_best - x(:, i);
            f_c_bias = y_global_best - y(i);
            c_d = c_d .* f_c_bias;
            c_d = c_d ./ ((sum(c_d) + eps) .* (sum(c_d) + eps));

            d_g = c_d + g_hist_best(:, i);

            % 梯度估计
            r1 = rand(dim, 1);
            r2 = rand(dim, 1);
            g = (1 - r1 - r2) .* g_temp + r1 .* d_p + r2 .* d_g;
            g = g ./ (sum(g) + eps);

            % 更新动量
            m(:, i) = beta1 .* m(:, i) + (1 - beta1) .* g;
            v(:, i) = beta2 .* v(:, i) + (1 - beta2) .* g.^2;
            w(:, i) = w(:, i) - m(:, i) ./ (sqrt(v(:, i)) + eps);

            % 更新位置：Advice Forward
            x_o = x(:, i) + exp(-l / (0.1 * Max_iter)) * 0.1 .* (ub(1) - lb(1)) .* g;
            x_o = Bounds(x_o, lb, ub);  % 边界处理

            % 更新粒子位置
            x(:, i) = x_o;
        end

        % 计算新位置的适应度值
        xx = sub2ind([nn, mm], round(x(1, :)), round(x(2, :)));
       %******原始成像逻辑*******
%         recm1z(xx(:)) = 0;  % 初始化成像值为零
%         y = zeros(1,SearchAgents_no);  % 重置适应度值
%         for ii = 1:size(TTPS, 2)  % 遍历每个台站 ii
%             for jj = 1:size(TTPS, 2)  % 遍历每个台站 jj
%                 ntp = round((TTPS(xx, ii) - TTPS(xx, jj)) / dt) + k0;  % 计算传播时间差
%                 ntp = min(max(ntp, 1), nnt);  % 限制 ntp 范围
%                 y = y + recW5(ntp, jj, ii)'.^2;  % 更新适应度值
%                 recm1z(xx(:)) = recm1z(xx(:)) + recW5(ntp, jj, ii)'.^2;  % 更新成像值
%             end
%         end
        
            %******更改成像逻辑*******
        recm1z(xx(:)) = 0;  % recm1z 初始值为 0
        y = zeros(1,SearchAgents_no);
        % 外层循环：组间（ii）
        for ii = 1:size(TTPS, 2)
            group_product = zeros(1, nn*mm);
            group_product(xx(:)) = 1; 
            for jj = 1:size(TTPS, 2)
                ntp = round((TTPS(xx, ii) - TTPS(xx, jj)) / dt) + k0;
                ntp = min(max(ntp, 1), nnt);  % 限制 ntp 范围
                group_product(xx(:)) = group_product(xx(:)) .* recW5(ntp, jj, ii)';
                y = y + recW5(ntp, jj, ii)';
            end
            recm1z(xx(:)) = recm1z(xx(:)) + group_product(xx(:));
        end

        % 更新全局最优
        [y_current_best, bestI] = max(y);
        if y_current_best > y_global_best
            y_global_best = y_current_best;
            x_global_best = x(:, bestI);
        end

        % 记录当前迭代的全局最优适应度值
        Convergence_curve(l) = y_global_best;
    end
end

function s = Bounds(s, Lb, Ub)
    % 确保 s、Lb 和 Ub 维度一致
    s = s(:);  % 将 s 转为列向量
    Lb = Lb(:);  % 将 Lb 转为列向量
    Ub = Ub(:);  % 将 Ub 转为列向量

    % 生成随机数矩阵，大小与 s 相同
    randVec = rand(size(s));

    % 应用下界
    I = s < Lb;
    s(I) = Lb(I) + randVec(I) .* (Ub(I) - Lb(I));  % 随机生成下界内的值

    % 应用上界
    J = s > Ub;
    s(J) = Lb(J) + randVec(J) .* (Ub(J) - Lb(J));  % 随机生成上界内的值
end
