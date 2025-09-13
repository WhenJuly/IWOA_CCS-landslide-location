function pop_init = initialization(xmin, xmax, dim, SearchAgents_no, init_method)
% initialization - 优化初始化函数
%
% 输入参数：
%   xmin: 粒子位置的下界 (dim x SearchAgents_no)
%   xmax: 粒子位置的上界 (dim x SearchAgents_no)
%   dim: 搜索空间的维度
%   SearchAgents_no: 种群数量
%   init_method: 初始化方法 ('PWLCM', 'Sine', 'Tent')
%
% 输出参数：
%   pop_init: 初始化种群 (dim x SearchAgents_no)

% 检查 init_method 是否合法
if ~ismember(lower(init_method), {'pwclm', 'sine', 'tent'})
    error('未知的初始化方法。请选择 "PWLCM", "Sine", 或 "Tent"。');
end

% 根据 init_method 选择不同的初始化方法
switch lower(init_method)
    case 'pwclm'
        % PWLCM 初始化
        pop_init = zeros(dim, SearchAgents_no);
        P = 0.4;  % PWLCM 参数
        for d = 1:dim
            x = rand;  % 初始值
            for i = 2:SearchAgents_no
                if x >= 0 && x < P
                    x = x / P;
                elseif x >= P && x < 0.5
                    x = (x - P) / (0.5 - P);
                elseif x >= 0.5 && x < 1-P
                    x = (1 - P - x) / (0.5 - P);
                elseif x >= 1-P && x < 1
                    x = (1 - x) / P;
                end
                pop_init(d, i) = x;  % 存储混沌值
            end
        end
    case 'sine'
        % Sine 初始化
        pop_init = zeros(dim, SearchAgents_no);
        for d = 1:dim
            x = rand;  % 初始值
            for i = 2:SearchAgents_no
                x = 0.99 * sin(pi * x);  % Sine 映射
                pop_init(d, i) = x;  % 存储混沌值
            end
        end
    case 'tent'
        % Tent 初始化
        pop_init = zeros(dim, SearchAgents_no);
        Alpa = 0.499;  % Tent 参数
        for d = 1:dim
            x = rand;  % 初始值
            for i = 2:SearchAgents_no
                if x < Alpa
                    x = x / Alpa;
                else
                    x = (1 - x) / (1 - Alpa);
                end
                pop_init(d, i) = x;  % 存储混沌值
            end
        end
end

% 将混沌值映射到 [xmin, xmax] 范围内
pop_init = xmin + (xmax - xmin) .* pop_init;  % 线性映射
end
