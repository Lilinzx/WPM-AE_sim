%% 多区域声场控制算法综合仿真平台
% 包含算法：ACC, WPM-Fixed, WPM-LWE, WPM-AE (Proposed)
% 功能：生成RIR -> 运行算法 -> 绘制性能对比图 -> 绘制声场热力图
% 修复：解决了可视化与训练物理模型不一致的问题
clc; clear; close all;

%% === Part 1: 用户可调参数配置 (User Parameters) ===

% --- 1. 物理与几何配置 ---
rt60_val = 0.2;         % 混响时间 (秒) [推荐: 0.2 - 0.5]
mic_z    = 2.0;         % 麦克风/声源平面的高度 (米)
L_spk    = 16;          % 扬声器数量
dist_r   = 1.5;         % 阵列半径 (米)

% --- 2. 算法目标参数 ---
target_AC_dB  = 30;     % 目标声学对比度 (dB) [推荐: 20-40]
target_AE_dB  = 0;      % 阵列努力约束上限 (dB, 0dB = 单扬声器能量)
target_LWE_dB = 0;      % 白噪声增益约束上限 (dB)
delta_fixed   = 1e-3;   % WPM-Fixed 算法的固定正则化参数

% --- 3. 绘图配置 ---
plot_freq     = 4000;   % 画声场图时使用的代表性频率 (Hz)
vis_ism_order = 1;      % 可视化时的反射阶数 (1=仅一阶反射+直达声, 设高了会很慢!)

%% === Part 2: 系统初始化 ===
fprintf('=== 初始化系统 ===\n');

% 获取基础配置结构体
conf = get_sim_config(); 
% 覆盖用户自定义参数
conf.rt60 = rt60_val;
conf.L = L_spk;
conf.R_array = dist_r;
% 更新坐标 (因为L可能变了)
d_theta = 3 * pi/180; % 间隔弧度
span = (conf.L - 1) * d_theta;
angles = linspace(pi/2 - span/2, pi/2 + span/2, conf.L);
conf.ls_coords = [conf.R_array * cos(angles'), conf.R_array * sin(angles'), ones(conf.L,1)*mic_z];

% 生成或加载房间脉冲响应 (RIR)
fprintf('正在生成/计算 RIR (RT60=%.2fs)...\n', conf.rt60);
[H_bright, H_dark, d_bright] = generate_rirs(conf);

%% === Part 3: 运行四种对比算法 ===
fprintf('\n=== 开始算法仿真 ===\n');

% 1. ACC (Baseline: 最大化对比度，不顾及能量效率)
[q_acc, m_acc] = algo_acc(H_bright, H_dark, d_bright, conf);

% 2. WPM-Fixed (Baseline: 固定正则化参数)
[q_fix, m_fix] = algo_wpm_fixed(H_bright, H_dark, d_bright, conf, target_AC_dB, delta_fixed);

% 3. WPM-LWE (Constraint: 限制白噪声增益)
[q_lwe, m_lwe] = algo_wpm_lwe(H_bright, H_dark, d_bright, conf, target_AC_dB, target_LWE_dB);

% 4. WPM-AE (Proposed: 限制阵列努力)
[q_ae, m_ae]   = algo_wpm_ae(H_bright, H_dark, d_bright, conf, target_AC_dB, target_AE_dB);

%% === Part 4: 性能对比绘图 (Fig 1) ===
fprintf('\n=== 绘制性能指标对比图 ===\n');
f_axis = conf.f_axis(conf.valid_bins);
idx = conf.valid_bins;

figure('Name', 'Performance Comparison', 'Position', [100, 100, 1200, 800], 'Color', 'w');
t = tiledlayout(3, 1, 'Padding', 'compact', 'TileSpacing', 'compact');

% Subplot 1: Acoustic Contrast (AC)
nexttile;
plot(f_axis, m_acc.AC(idx), 'k-', 'LineWidth', 1.0, 'DisplayName', 'ACC'); hold on;
plot(f_axis, m_fix.AC(idx), 'g-.', 'LineWidth', 1.5, 'DisplayName', ['Fixed \delta=', num2str(delta_fixed)]);
plot(f_axis, m_lwe.AC(idx), 'r--', 'LineWidth', 1.5, 'DisplayName', 'WPM-LWE');
plot(f_axis, m_ae.AC(idx),  'b-',  'LineWidth', 2.0, 'DisplayName', 'WPM-AE (Prop)');
yline(target_AC_dB, 'm:', 'Target AC', 'LineWidth', 2);
ylabel('AC (dB)'); title('声学对比度 (Acoustic Contrast)'); 
legend('Location', 'best'); grid on; ylim([0, 40]); xlim([100, 4000]);

% Subplot 2: Array Effort (AE)
nexttile;
plot(f_axis, m_acc.AE(idx), 'k-', 'LineWidth', 1.0); hold on;
plot(f_axis, m_fix.AE(idx), 'g-.', 'LineWidth', 1.5);
plot(f_axis, m_lwe.AE(idx), 'r--', 'LineWidth', 1.5);
plot(f_axis, m_ae.AE(idx),  'b-',  'LineWidth', 2.0);
yline(target_AE_dB, 'k--', 'Constraint Limit');
ylabel('AE (dB)'); title('阵列努力 (Array Effort)'); 
grid on; ylim([-10, 20]); xlim([100, 4000]);

% Subplot 3: Spectral Distortion (SD) - 重建误差
nexttile;
% 注意: ACC通常不关注SD，所以可能很大。我们主要看WPM类的。
plot(f_axis, m_fix.SD(idx), 'g-.', 'LineWidth', 1.5); hold on;
plot(f_axis, m_lwe.SD(idx), 'r--', 'LineWidth', 1.5);
plot(f_axis, m_ae.SD(idx),  'b-',  'LineWidth', 2.0);
ylabel('SD (dB)'); xlabel('Frequency (Hz)'); 
title('频谱失真 (Spectral Distortion)'); 
legend({'WPM-Fixed', 'WPM-LWE', 'WPM-AE'}, 'Location', 'best');
grid on; xlim([100, 4000]); ylim([-30, 10]);

%% === Part 5: 声场可视化 (Fig 2) ===
fprintf('\n=== 计算与绘制 %.0f Hz 处的声场 (ISM Order=%d) ===\n', plot_freq, vis_ism_order);
% 找到最接近 plot_freq 的频点索引
[~, k_idx] = min(abs(conf.f_axis - plot_freq));
real_freq = conf.f_axis(k_idx);
fprintf('   实际绘制频率: %.1f Hz\n', real_freq);

% 1. 定义绘图网格
x_range = -1.0 : 0.05 : 1.0; 
y_range = -1.0 : 0.05 : 1.5; 
[X, Y] = meshgrid(x_range, y_range);
grid_pts = [X(:), Y(:), ones(numel(X), 1) * conf.mic_height];

% 2. 计算网格传递函数 (最耗时步骤)
% 必须使用 ism_core 确保与训练一致
H_grid_k = compute_grid_atf_single_freq(grid_pts, conf, real_freq, vis_ism_order);

% 3. 计算各算法声压分布
P_acc = reshape(H_grid_k * q_acc(k_idx, :).', size(X));
P_fix = reshape(H_grid_k * q_fix(k_idx, :).', size(X));
P_lwe = reshape(H_grid_k * q_lwe(k_idx, :).', size(X));
P_ae  = reshape(H_grid_k * q_ae(k_idx, :).',  size(X));

% 4. 绘图
figure('Name', 'Sound Field Comparison', 'Position', [150, 150, 1000, 800], 'Color', 'w');
alg_names = {'ACC', 'WPM-Fixed', 'WPM-LWE', 'WPM-AE (Prop)'};
maps = {P_acc, P_fix, P_lwe, P_ae};

% 计算归一化基准 (以亮区中心能量为 0dB)
center_dist = sqrt((X - 0.5).^2 + Y.^2);
mask_center = center_dist < 0.1;

for i = 1:4
    subplot(2, 2, i);
    
    E_mag = abs(maps{i}).^2;
    ref_val = mean(E_mag(mask_center));
    E_dB = 10*log10(E_mag / (ref_val + 1e-12));
    
    imagesc(x_range, y_range, E_dB);
    axis xy equal; colormap('jet'); caxis([-30, 5]);
    title(alg_names{i}); 
    if i==1, ylabel('y (m)'); end
    if i==3, xlabel('x (m)'); ylabel('y (m)'); end
    if i==4, xlabel('x (m)'); colorbar; end
    
    hold on;
    % 画区域圆圈
    viscircles([0.5, 0], conf.r_zone, 'Color', 'k', 'LineWidth', 1);
    viscircles([-0.5, 0], conf.r_zone, 'Color', 'w', 'LineWidth', 1, 'LineStyle', '--');
    % 画扬声器
    scatter(conf.ls_coords(:,1), conf.ls_coords(:,2), 30, 'k', 'filled', 's');
end

%% === 辅助函数: 单频点网格计算 ===
function H_out = compute_grid_atf_single_freq(grid_coords, conf, target_freq, order)
    % 专为绘图优化的单频点计算，调用 ism_core
    N = size(grid_coords, 1);
    L = conf.L;
    H_out = zeros(N, L);
    
    % 预计算常数
    omega = 2*pi*target_freq;
    V = prod(conf.room_dim);
    S = 2*(conf.room_dim(1)*conf.room_dim(2) + ...
           conf.room_dim(2)*conf.room_dim(3) + ...
           conf.room_dim(1)*conf.room_dim(3));
    alpha = 0.161 * V / (S * conf.rt60);
    beta = sqrt(1 - alpha);
    
    % 并行计算加速
    parfor n = 1:N
        rcv = grid_coords(n, :);
        row_h = zeros(1, L);
        for l = 1:L
            src = conf.ls_coords(l, :);
            % 调用独立文件 ism_core
            % 注意：ism_core 返回时域脉冲，我们需要做 FFT 提取特定频点
            % 为了速度，这里我们手动计算频域响应会更快？
            % 不，ism_core 比较复杂，还是通过 h -> fft 稳妥，或者修改 ism_core 支持频域。
            % 为了代码复用，我们这里调用 ism_core 得到 h，然后计算单点 DTFT。
            
            h_time = ism_core(src, rcv, conf, beta, order, false);
            
            % 计算单频点响应: sum(h[n] * exp(-j*w*n/fs))
            % 这是一个简单的离散傅里叶变换
            t_idx = (0:length(h_time)-1)';
            row_h(l) = sum(h_time .* exp(-1j * omega * t_idx / conf.fs));
        end
        H_out(n, :) = row_h;
    end
end