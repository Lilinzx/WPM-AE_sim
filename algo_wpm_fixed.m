function [q_fixed, metrics] = algo_wpm_fixed(H_bright, H_dark, d_bright, conf, AC_target_dB, delta)
% ALGO_WPM_FIXED  加权声压匹配 (WPM) - 固定正则化版
% 
% 这是一个鲁棒的基线算法 (Baseline)。
% 它不进行复杂的迭代搜索，而是直接求解一个带正则化的最小二乘问题。
%
% 目标函数 J = ||Hb*q - d||^2 + mu * ||Hd*q||^2 + delta * ||q||^2
%
% 输入:
%   delta: 正则化参数 (Regularization factor)，用于限制阵列能量，防止过载。
%          值越大，AE越低，但重构误差越大。
%   AC_target_dB: 在此函数中未直接使用，仅为了保持接口一致性而保留。

    fprintf('算法: 运行 WPM-Fixed (鲁棒基线, delta=%.1e)...\n', delta);
    
    [K, Mb, L] = size(H_bright);
    [~, Md, ~] = size(H_dark);
    
    q_fixed = zeros(K, L);
    
    % 初始化指标存储
    metrics.AC = zeros(K, 1);
    metrics.AE = zeros(K, 1);
    metrics.SD = zeros(K, 1);
    
    % 参考权重 (用于计算 AE 的归一化基准)
    q0 = zeros(L, 1); q0(conf.target_ls_idx) = 1;
    
    % 暗区抑制权重 (mu)
    % mu = 0   -> 纯声压匹配 (PM), 不管暗区
    % mu > 0   -> 加权声压匹配 (WPM), 压制暗区
    % 这里设为 1.0，表示暗区静音和亮区发声同等重要
    mu = 1.0; 
    
    %% 频域循环求解
    for k = 1:K
        % 1. 提取当前频点数据 [Mb x L]
        Hb = squeeze(H_bright(k, :, :)); 
        Hd = squeeze(H_dark(k, :, :));   
        db = reshape(d_bright(k, :, :), Mb, 1);   
        
        % 2. 构建相关矩阵
        % 为了公平对比，我们尽量保持与 WPM-AE 相似的矩阵缩放
        Rb = (Hb' * Hb) / Mb;
        Rd = (Hd' * Hd) / Md;
        r_db = (Hb' * db) / Mb; % 互相关向量
        
        % 3. 闭式求解 (Closed-form Solution)
        % 求解方程: (Rb + mu*Rd + delta*I) * q = r_db
        % 这是一个标准的 Tikhonov 正则化最小二乘解
        
        % 构建系统矩阵 A
        A_sys = Rb + mu * Rd + delta * eye(L);
        
        % 求解 q
        % 使用反斜杠算子，Matlab 会自动选择最优解法 (Cholesky 或 LU)
        q_sol = A_sys \ r_db;
        
        % 4. 存储结果
        q_fixed(k, :) = q_sol.';
        
        % 5. 计算性能指标 (Metrics)
        
        % (1) 声学对比度 AC
        Eb = real(q_sol' * Rb * q_sol);
        Ed = real(q_sol' * Rd * q_sol);
        metrics.AC(k) = 10*log10(Eb / (Ed + 1e-15));
        
        % (2) 阵列努力 AE
        % g0 是单扬声器在亮区产生的参考能量
        g0 = real(q0' * Rb * q0); 
        % AG (Array Gain) = 亮区能量 / 阵列总输入电功率
        AG = Eb / (real(q_sol' * q_sol) + 1e-15);
        metrics.AE(k) = 10*log10(g0 / (AG + 1e-15));
        
        % (3) 频谱失真 SD (Spectral Distortion)
        % 衡量实际生成的声场与期望声场 db 的误差
        % SD = 10 log10( ||Hb*q - d||^2 / ||d||^2 )
        err_energy = norm(Hb * q_sol - db)^2;
        target_energy = norm(db)^2;
        metrics.SD(k) = 10*log10(err_energy / (target_energy + 1e-15));
    end
    
    fprintf('    WPM-Fixed 完成. 平均 AC: %.2f dB\n', mean(metrics.AC));
end