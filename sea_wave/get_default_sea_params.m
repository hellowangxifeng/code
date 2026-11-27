function sea_params = get_default_sea_params()
    % =====================================================================
    % 动态海面默认参数生成器 (严格执行空间奈奎斯特频率绑定)
    % 核心逻辑：以物理最高频率 f_max 为主导，自适应反推并锁定网格密度
    % 设置频率范围、频率与方向分段数、空间范围，并由最高波浪频率计算网格密度。
    % =====================================================================
    sea_params.g = 9.81;
    sea_params.neumannConstant = 3.05;
    g = sea_params.g; % 重力加速度
    
    % --- 1. 空间域边界设定 ---
    sea_params.x_range      = [-5, 5];
    sea_params.y_range      = [-5, 5];
    
    % --- 2. 核心物理绑定：最高仿真频率 (主导参数) ---
    sea_params.f_max        = 4.0; % 最高仿真频率 (Hz)
    
    % --- 3. 空间离散化 (自适应计算网格密度) ---
    % 依据理论推导: lambda_min = 2*dx, 结合深水波色散关系 dx = g / (4*pi*f_max^2)
    % 确保生成的海面既不发生空间混叠，又满足奈奎斯特极限
    dx_required = g / (4 * pi * sea_params.f_max^2);
    
    % 取 X 和 Y 方向的最大跨度以保证全域安全
    L_x = sea_params.x_range(2) - sea_params.x_range(1);
    L_y = sea_params.y_range(2) - sea_params.y_range(1);
    L_max = max(L_x, L_y); 
    
    % 向上取整计算网格密度，加 1 以保证完美包含边界
    sea_params.grid_density = ceil(L_max / dx_required) + 1;
    
    % --- 4. 其他常规参数设定 ---
    sea_params.f_min        = 0.05;              % 频率下限 (Hz)，保留较大风速的长波能量
    sea_params.Fre_segment  = 50;                % 频率子带划分数
    sea_params.Dir_segment  = 36;                % 方向子带划分数
    sea_params.U            = 10;                % 风速 (m/s)
    sea_params.t0           = 0;                 % 初始时刻
    sea_params.seed         = 42;                % 随机种子 (保证复现性)
end
