function [XX, YY, Z_surf, plane_data, selected_plane_data, wave] = oceanSurfaceModel_phiwind_260315(f_min, f_max, Fre_segment, Dir_segment, U, phi_wind, x_range, y_range, grid_density, t0, seed, output_mode, plane_index)
% 基于双叠加法和平面拟合的海面建模 - 整合修正版本 (可见光水下通信专版)
% 按 Neumann 方向谱计算波浪振幅，加入随机相位和风向旋转，叠加出海面高度；再拟合网格小平面，输出法向等数据。
% 
% 输入参数:
%   f_min, f_max: 频率范围 (Hz)
%   Fre_segment, Dir_segment: 频率和方向分段数
%   U: 风速 (m/s)
%   phi_wind: 风向与光束的夹角(弧度)。0表示顺风(风沿Y轴)，pi/2表示侧风(风沿X轴) % [新增]
%   x_range, y_range: X和Y轴范围 [min, max] (建议设置能够容纳主波长的大小，如[-25, 25])
%   grid_density: 网格密度 (每个方向的网格点数)
%   t0: 时间参数
%   seed: 随机种子
%   output_mode: 输出模式 ('all' 或 'single')
%   plane_index: 当output_mode='single'时，指定要输出的平面索引 [i,j] 或线性索引
% 
% 输出:
%   XX, YY: 网格坐标
%   Z_surf: 波面高度
%   plane_data: 所有平面的参数结构体
%   selected_plane_data: 选中的单个平面数据（当output_mode='single'时）
%   wave: 可选第六输出，原海面同一组 omega/kx/ky/coefficient；复振幅参考全局原点及 t=0。

    % 参数默认值处理
    if nargin < 11
        seed = 1;
    end
    if nargin < 12
        output_mode = 'all'; % 默认输出所有平面
    end
    if nargin < 13
        plane_index = [1, 1]; % 默认选择第一个平面
    end
    
    % 验证输出模式参数
    if ~ismember(output_mode, {'all', 'single'})
        error('output_mode 必须是 "all" 或 "single"');
    end
   
    sea = get_default_sea_params();
    p = struct('g',sea.g,'fMinHz',f_min,'fMaxHz',f_max,'frequencyBins',Fre_segment, ...
        'directionBins',Dir_segment,'windSpeed',U,'windAngleDeg',rad2deg(phi_wind),'seed',seed);
    wave = RIS_wave_components(p);
    [XX,YY] = meshgrid(linspace(x_range(1),x_range(2),grid_density), ...
        linspace(y_range(1),y_range(2),grid_density));
    field = struct('wave',wave,'parameters',struct('startTime',t0));
    Z_surf = reshape(RIS_wave_field(field,[XX(:),YY(:)],0),size(XX));
    plane_data = calculatePlaneParameters(XX,YY,Z_surf);
    selected_plane_data = [];
    if strcmp(output_mode,'single')
        selected_plane_data = getSelectedPlaneData(plane_data,plane_index);
    end

    function plane_data = calculatePlaneParameters(XX, YY, Z_surf)
        % 计算每个小平面的ABC系数和法向量
        [rows, cols] = size(XX);
        
        % 初始化平面数据存储
        plane_data = struct();
        plane_data.A = zeros(rows-1, cols-1);
        plane_data.B = zeros(rows-1, cols-1);
        plane_data.C = zeros(rows-1, cols-1);
        plane_data.normal_x = zeros(rows-1, cols-1);
        plane_data.normal_y = zeros(rows-1, cols-1);
        plane_data.normal_z = zeros(rows-1, cols-1);
        plane_data.centroid_x = zeros(rows-1, cols-1);
        plane_data.centroid_y = zeros(rows-1, cols-1);
        plane_data.centroid_z = zeros(rows-1, cols-1);
        plane_data.fit_error = zeros(rows-1, cols-1);
        plane_data.vertices = cell(rows-1, cols-1); 
        
        % 计算每个小平面
        for idx_i = 1:rows-1
            for idx_j = 1:cols-1
                vertices = [
                    XX(idx_i,idx_j),     YY(idx_i,idx_j),     Z_surf(idx_i,idx_j);      
                    XX(idx_i,idx_j+1),   YY(idx_i,idx_j+1),   Z_surf(idx_i,idx_j+1);    
                    XX(idx_i+1,idx_j),   YY(idx_i+1,idx_j),   Z_surf(idx_i+1,idx_j);    
                    XX(idx_i+1,idx_j+1), YY(idx_i+1,idx_j+1), Z_surf(idx_i+1,idx_j+1)   
                ];
                
                [A, B, C, fitError, normal] = fitPlaneLS(vertices);
                
                plane_data.A(idx_i,idx_j) = A;
                plane_data.B(idx_i,idx_j) = B;
                plane_data.C(idx_i,idx_j) = C;
                plane_data.normal_x(idx_i,idx_j) = normal(1);
                plane_data.normal_y(idx_i,idx_j) = normal(2);
                plane_data.normal_z(idx_i,idx_j) = normal(3);
                plane_data.centroid_x(idx_i,idx_j) = mean(vertices(:,1));
                plane_data.centroid_y(idx_i,idx_j) = mean(vertices(:,2));
                plane_data.centroid_z(idx_i,idx_j) = mean(vertices(:,3));
                plane_data.fit_error(idx_i,idx_j) = fitError;
                plane_data.vertices{idx_i,idx_j} = vertices;
            end
        end
        
        function [A, B, C, fitError, normal] = fitPlaneLS(vertices)
            x = vertices(:,1); y = vertices(:,2); z = vertices(:,3);
            X_design = [x, y, ones(size(x))];
            coefficients = X_design \ z;
            
            A = coefficients(1); B = coefficients(2); C = coefficients(3);
            
            z_pred = X_design * coefficients;
            fitError = sqrt(mean((z - z_pred).^2));
            
            normal_vec = [A, B, -1];
            normal_length = norm(normal_vec);
            if normal_length > 0
                normal = normal_vec / normal_length;
            else
                normal = [0, 0, -1]; 
            end
            if normal(3) < 0
                normal = -normal;
            end
        end
    end

    function selected_data = getSelectedPlaneData(plane_data, plane_index)
        [rows, cols] = size(plane_data.A);
        if isscalar(plane_index)
            [i, j] = ind2sub([rows, cols], plane_index);
        else
            i = plane_index(1); j = plane_index(2);
        end
        
        selected_data = struct();
        selected_data.index = [i, j];
        selected_data.linear_index = sub2ind([rows, cols], i, j);
        selected_data.equation = struct('A', plane_data.A(i,j), 'B', plane_data.B(i,j), 'C', plane_data.C(i,j));
        selected_data.normal_vector = [plane_data.normal_x(i,j), plane_data.normal_y(i,j), plane_data.normal_z(i,j)];
        selected_data.centroid = [plane_data.centroid_x(i,j), plane_data.centroid_y(i,j), plane_data.centroid_z(i,j)];
        selected_data.fit_error = plane_data.fit_error(i,j);
        selected_data.vertices = plane_data.vertices{i,j};
        selected_data.equation_string = sprintf('z = %.6f*x + %.6f*y + %.6f', selected_data.equation.A, selected_data.equation.B, selected_data.equation.C);
    end

end
