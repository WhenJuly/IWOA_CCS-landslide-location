function [y_global_best, x_global_best, Convergence_curve, recm1z] = SSA(SearchAgents_no, Max_iter, lb, ub, dim, pop_init, nn, mm, recm1z, TTPS, Fs, k0, nnt, recW5)
    % SSA 参数设置
    P_percent = 0.2;  % 生产者比例
    pNum = round(SearchAgents_no * P_percent);  % 生产者数量

    % 初始化种群和适应度值
    x = pop_init;  % 种群位置
    pmax = zeros(1, SearchAgents_no);  % 每个粒子的历史最优适应度值
    xpmax = x;  % 每个粒子的历史最优位置
    gmax = 0;  % 全局最优适应度值
    xgmax = zeros(dim, 1);  % 全局最优位置
    Convergence_curve = zeros(1, Max_iter);  % 记录每轮迭代的全局最优适应度值
    pmax = zeros(1, SearchAgents_no);  % 初始化每个粒子的历史最优适应度值
    ggmax = zeros(1, Max_iter);  % 初始化全局最优适应度记录
    dt = 1/Fs;
    %% 第一次迭代：计算初始适应度值和全局最优值
    xx = sub2ind([nn, mm], round(x(1, :)), round(x(2, :)));
%     %******原始成像逻辑*******
%     recm1z(xx(:)) = 0;
%     % 为所有粒子计算适应度值
%     for ii = 1:size(TTPS, 2)  % 遍历每个台站 ii
%         for jj = 1:size(TTPS, 2)  % 遍历每个台站 jj
%             ntp = round((TTPS(xx, ii) - TTPS(xx, jj)) / dt) + k0;  % 计算传播时间差并转换为采样点 xx是转换后的点位
%             ntp = min(max(ntp, 1), nnt);  % 限制 ntp 范围，确保不超出边界
%             pmax = pmax + recW5(ntp, jj, ii)'.^2;  % 累加每对台站的互相关结果的平方
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
            pmax = pmax + recW5(ntp, jj, ii)';
        end
        recm1z(xx(:)) = recm1z(xx(:)) + group_product(xx(:));
    end
    
    % 更新每个粒子的历史最优适应度值和位置
    gmax = max(pmax);  % 找到当前全局最优适应度值
    Convergence_curve(1) = gmax;  % 保存第一次迭代的最优适应度
    [gmax, maxpos] = max(recm1z(xx(:)));  % 找到最大成像值的位置
    xgmax = x(:, maxpos);  % 更新全局最优位置
    %% 从第二轮迭代开始
    for t = 2:Max_iter
        % 更新粒子位置
        [~, sortIndex] = sort(recm1z(xx(:)));  % 对适应度值排序
        worstIdx = sortIndex(end);  % 最差个体

        % 生产者更新
        for i = 1:pNum
            r1 = rand(1);
            x(:, sortIndex(i)) = x(:, sortIndex(i)) .* exp(-(i) / (r1 * Max_iter));
            x(:, sortIndex(i)) = Bounds(x(:, sortIndex(i)), lb, ub);  % 确保粒子位置在边界内
        end
%         disp(t)
        % 跟随者更新
        for i = (pNum + 1):SearchAgents_no
            A = floor(rand(1, dim) * 2) * 2 - 1;
            if i > SearchAgents_no / 2
                x(:, sortIndex(i)) = randn(dim, 1) .* exp((x(:, worstIdx) - x(:, sortIndex(i))) / (i)^2);
            else
                x(:, sortIndex(i)) = xgmax + (abs(x(:, sortIndex(i)) - xgmax)' * (A' * (A * A')^(-1)) * ones(1, dim))';
            end
%         disp(i)
        x(:, sortIndex(i)) = Bounds(x(:, sortIndex(i)), lb, ub);

        end

        % 计算当前种群的适应度值       
        xx = sub2ind([nn, mm], round(x(1, :)), round(x(2, :)));
%         %******原始成像逻辑*******
%         recm1z(xx(:)) = 0;  % 初始化成像值为零
%         % 为所有粒子计算适应度值
%         for ii = 1:size(TTPS, 2)  % 遍历每个台站 ii
%             for jj = 1:size(TTPS, 2)  % 遍历每个台站 jj
%                 ntp = round((TTPS(xx, ii) - TTPS(xx, jj)) / dt) + k0;  % 计算传播时间差并转换为采样点 xx是转换后的点位
%                 ntp = min(max(ntp, 1), nnt);  % 限制 ntp 范围，确保不超出边界
%                 recm1z(xx(:)) = recm1z(xx(:)) + recW5(ntp, jj, ii)'.^2;  % 更新成像值
%             end
%         end
%         
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
            end
            recm1z(xx(:)) = recm1z(xx(:)) + group_product(xx(:));
        end

        ind2 = recm1z(xx(:)) > pmax;  % 如果当前适应度优于历史最优，则更新
        pmax(ind2) = recm1z(ind2);  % 更新历史最优适应度 每个粒子对应起来的
        xpmax(:, ind2) = x(:, ind2);  % 更新历史最优位置

        [gmax, maxpos] = max(recm1z);  % 计算当前最优适应度
        [xg1, xg2] = ind2sub([nn, mm], maxpos);  % 获取最优适应度的位置
        xgmax = [xg1; xg2] .* ones(2, 1);  % 更新全局最优位置
        % 记录当前迭代的全局最优适应度值
    Convergence_curve(t) = gmax;
    end

    % 返回全局最优适应度值和位置
    y_global_best = gmax;
    x_global_best = xgmax;
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
