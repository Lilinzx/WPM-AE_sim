function conf = get_sim_config()
% GET_SIM_CONFIG  加载仿真所需的全局物理参数与系统配置
%
% 该函数定义了房间声学环境、扬声器阵列几何、控制区域分布以及频域处理参数。
% 所有长度单位均为：米 (m)
% 所有时间单位均为：秒 (s)
%
% 输出:
%   conf: 包含所有仿真参数的结构体

    %% 1. 基础物理与房间声学参数 (Basic Physics & Room Acoustics)
    conf.c = 343;           % 声速 (m/s), 空气中约 343 m/s
    conf.fs = 16384;        % 采样率 (Hz), 决定了仿真的频宽
    
    % 房间尺寸 [长, 宽, 高]
    % 坐标系定义: 房间角落为 (0,0,0)，长宽高分别沿 x,y,z 轴
    conf.room_dim = [6, 5, 4]; 
    
    % 混响时间 RT60 (s)
    % 影响 ISM 镜像声源法中的反射系数 beta
    % 0.2s = 强吸声(录音棚); 0.5s = 普通会议室; >1.0s = 厅堂
    conf.rt60 = 0.2;        
    
    % RIR (房间脉冲响应) 的长度 (采样点数)
    % 2048 / 16384 ≈ 0.125秒，足以覆盖主要的反射声
    conf.nsamples_rir = 2048; 
    
    %% 2. 扬声器阵列设置 (Loudspeaker Array Setup)
    % 这里定义默认值，main 脚本中可能会根据用户设置进行覆盖
    conf.L = 16;            % 扬声器数量
    conf.R_array = 1.5;     % 阵列半径 (米)
    
    % 几何布局: 弧形阵列 (Arc)
    % 以房间中心为圆心，朝向听众区域
    d_theta_deg = 3;        % 扬声器间隔角度 (度)
    
    % 计算角度分布 (以 90度/正Y轴 为中心对称)
    span_deg = (conf.L - 1) * d_theta_deg;
    theta_deg = linspace(90 - span_deg/2, 90 + span_deg/2, conf.L);
    theta_rad = deg2rad(theta_deg);
    
    % 扬声器坐标计算 (默认高度 z=2m)
    % 注意: 这里先计算相对坐标，generate_rirs 可能会加上房间中心偏移
    % 但为了简化，我们在 ism_core 中通常假设输入已经是绝对坐标。
    % 这里我们定义相对于"仿真中心"的坐标。
    z_plane = 2.0; 
    conf.mic_height = z_plane;
    
    conf.ls_coords = zeros(conf.L, 3);
    conf.ls_coords(:, 1) = conf.R_array * cos(theta_rad); % x
    conf.ls_coords(:, 2) = conf.R_array * sin(theta_rad); % y
    conf.ls_coords(:, 3) = z_plane;                       % z
    
    %% 3. 控制区域设置 (Control Zones)
    % 定义亮区(Bright Zone)和暗区(Dark Zone)的控制点
    
    conf.Mb = 32;           % 亮区控制点数量 (采样点)
    conf.Md = 32;           % 暗区控制点数量
    conf.r_zone = 0.3;      % 区域半径 (米), 模拟人头大小范围
    conf.d_zones = 1.0;     % 两区域圆心间距 (米)
    
    % 区域中心坐标
    % 亮区在右 (0.5, 0), 暗区在左 (-0.5, 0)
    center_bright = [conf.d_zones/2, 0, z_plane];  
    center_dark   = [-conf.d_zones/2, 0, z_plane]; 
    
    % 生成亮区圆周上的采样点
    phi_b = linspace(0, 2*pi, conf.Mb+1); phi_b(end) = []; % 0~2pi
    conf.bright_coords = zeros(conf.Mb, 3);
    conf.bright_coords(:, 1) = center_bright(1) + conf.r_zone * cos(phi_b);
    conf.bright_coords(:, 2) = center_bright(2) + conf.r_zone * sin(phi_b);
    conf.bright_coords(:, 3) = z_plane;
    
    % 生成暗区圆周上的采样点
    phi_d = linspace(0, 2*pi, conf.Md+1); phi_d(end) = [];
    conf.dark_coords = zeros(conf.Md, 3);
    conf.dark_coords(:, 1) = center_dark(1) + conf.r_zone * cos(phi_d);
    conf.dark_coords(:, 2) = center_dark(2) + conf.r_zone * sin(phi_d);
    conf.dark_coords(:, 3) = z_plane;
    
    %% 4. 算法与频域参数 (Algorithm Parameters)
    conf.J = 2048;          % 控制滤波器长度 / FFT点数
    conf.target_ls_idx = 4; % 目标声源对应的扬声器索引 (用于定义期望声场)
    
    % 分析频率范围
    conf.f_range = [100, 4000]; % Hz (人耳敏感频段)
    
    % 频率轴计算
    % linspace(0, fs, J)
    conf.f_axis = linspace(0, conf.fs, conf.J); 
    
    % 标记有效频点 (Valid Bins)
    % 仅在有效频段内计算性能指标，避开 DC 和 Nyquist 附近的无效数据
    conf.valid_bins = find(conf.f_axis >= conf.f_range(1) & conf.f_axis <= conf.f_range(2));
    
end