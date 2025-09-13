function new_position = wavelet_mutation(position, lb, ub, t, T)
    % position: 当前被囊个体的位置 (1xD 向量)
    % p_v: 变异概率
    % lb: 下界 (1xD 向量)
    % ub: 上界 (1xD 向量)
    % t: 当前迭代次数
    % T: 最大迭代次数
    
    D = length(position);  % 维度
    new_position = position;  % 初始化新位置
    N = fix(T/20);
    % 计算变异概率
    v = 7;  % 自由度
    gamma_function = @(z) gamma(z); % 伽马函数
    p_v = 4 * ((t/N).^(v/2 - 1)) .* exp(-t/(2*N)) / (2^(v/2) * gamma_function(v/2));
    
    % 随机数 c1
    c1 = rand();
    
    if c1 < p_v
        % 计算伸缩参数 a
        g = 10000;  % 上限
        xi = 5;  % 形状参数
        a = exp(-log(g) * (1 - t/T)^xi) + log(g);
        
        % 计算小波变异系数 sigma
        phi = -2.5*a + 5*a*rand();  % 随机数 phi
        sigma = (1/sqrt(a)) * exp(-(phi/a)^2 / 2) * cos(5 * phi/a);
        
        % 对每个维度进行变异
        for j = 1:D
            if sigma > 0
                new_position(j) = position(j) + sigma * (ub - position(j));
            else
                new_position(j) = position(j) + sigma * (position(j) - lb);
            end
                     % 边界检查
            if new_position(j) < lb || new_position(j) > ub
                new_position(j) = lb + (ub - lb) * rand();  % 在边界内随机生成一个新位置
            end
        end

%         % 边界处理
%         new_position = max(new_position, lb);
%         new_position = min(new_position, ub);
    end
end
