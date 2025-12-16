function [q_ae, metrics] = algo_wpm_ae(H_bright, H_dark, d_bright, conf, AC_target_dB, AE_max_dB)
% ALGO_WPM_AE  WPM with Array Effort Constraint (WPM-AE)
% 
% 求解目标：最小化重构误差 (WPM)，同时满足 声学对比度(AC) 和 阵列努力(AE) 约束。
% 数学本质：求解带有两个二次约束的二次规划 (QCQP)。
%
% 本实现包含了论文 Eq.(20) 提到的 Trace-Normalization 加权，
% 确保 AE 约束在优化中具有比 AC 约束更高的优先级。
%
% 输入:
%   H_bright: 亮区传递函数矩阵
%   H_dark:   暗区传递函数矩阵
%   d_bright: 期望声场向量
%   conf:     系统配置
%   AC_target_dB: 目标 AC 值 (例如 30dB)
%   AE_max_dB:    允许的最大 AE 值 (通常为 0dB)

    fprintf('算法: 运行 WPM-AE (AC目标=%.0fdB, AE上限=%.0fdB, 含优先级加权)...\n', AC_target_dB, AE_max_dB);
    
    [K, Mb, L] = size(H_bright);
    [~, Md, ~] = size(H_dark);
    
    q_ae = zeros(K, L);
    metrics.AC = zeros(K, 1);
    metrics.AE = zeros(K, 1);
    metrics.SD = zeros(K, 1);
    
    % 线性域转换
    AC_lin = 10^(AC_target_dB/10);
    AE_limit_lin = 10^(AE_max_dB/10); 
    
    % 参考权重 (仅开启目标扬声器) 用于归一化 AE
    q0 = zeros(L, 1); q0(conf.target_ls_idx) = 1;
    
    % 优化器选项 (使用 SQP 算法)
    opts = optimoptions('fmincon', 'Display', 'off', ...
        'Algorithm', 'sqp', ...
        'SpecifyObjectiveGradient', false, ... 
        'OptimalityTolerance', 1e-6, ...
        'StepTolerance', 1e-8);

    % 进度条变量
    msg_len = 0;
    
    for k = 1:K
        % 简单的进度显示
        if mod(k, 200) == 0 || k == K
            fprintf(repmat('\b', 1, msg_len)); 
            msg = sprintf('    进度: %3d%%', round(k/K*100));
            fprintf(msg); msg_len = length(msg);
        end
        
        % 1. 提取当前频点数据
        Hb = squeeze(H_bright(k, :, :));
        Hd = squeeze(H_dark(k, :, :));
        db = reshape(d_bright(k, :, :), Mb, 1);
        
        % 空间相关矩阵
        Rb = (Hb' * Hb) / Mb;
        Rd = (Hd' * Hd) / Md;
        HbH_Hb = Hb' * Hb;
        HbH_db = Hb' * db;
        
        % 2. 构建对偶优化所需的矩阵 (转为实数域 2L x 2L)
        % 目标函数: Min ||Hb*q - db||^2  =>  q' H0 q - 2*y'*q
        H0 = [real(HbH_Hb), -imag(HbH_Hb); imag(HbH_Hb), real(HbH_Hb)];
        H0 = H0 + 1e-7 * eye(2*L); % 微量对角加载，防止完全奇异
        y_vec = [real(HbH_db); imag(HbH_db)];
        
        % 约束1 (AC): q'(AC*Rd - Rb)q <= 0
        S1 = AC_lin * Rd - Rb;
        H1 = [real(S1), -imag(S1); imag(S1), real(S1)];
        
        % 约束2 (AE): q'(eta*I - Rb)q <= 0
        % 推导: AE = g0 / AG <= Limit => AG >= g0/Limit
        % 整理得: q' ( (g0/Limit)*I - Rb ) q <= 0
        
        g0 = real(q0' * Rb * q0) / real(q0' * q0); % 系统的参考增益
        eta = g0 / (AE_limit_lin + 1e-12);
        
        S2_base = eta * eye(L) - Rb;
        
        % [论文复现修正] 优先级加权因子 (Priority Weighting)
        % 引用论文 Eq.(20): H2 = L/tr(Rb) * ...
        % 物理意义: 通过迹归一化，使得 AE 约束的梯度尺度与白噪声增益特性一致，
        % 从而在求解器中获得比 AC 约束更高的优先级。
        weight_factor = L / (trace(Rb) + 1e-15);
        S2_weighted = S2_base * weight_factor;
        
        H2 = [real(S2_weighted), -imag(S2_weighted); imag(S2_weighted), real(S2_weighted)];
        
        % 3. 对偶最大化 (Dual Maximization)
        % 寻找最佳的拉格朗日乘子 lambda = [lam1, lam2]
        objective = @(lam) dual_objective_safe(lam, H0, H1, H2, y_vec);
        
        % 初始猜测: 给一点小初值防止 H_sum 初始奇异
        x0 = [1e-5; 1e-5]; 
        
        try
            % 约束 lambda >= 0
            lambda_opt = fmincon(objective, x0, [], [], [], [], [0; 0], [], [], opts);
        catch
            lambda_opt = [0; 0];
        end
        
        % 4. 恢复原变量 (Primal Recovery)
        % x = (H0 + lam1*H1 + lam2*H2)^-1 * y
        H_final = H0 + lambda_opt(1)*H1 + lambda_opt(2)*H2;
        
        % 鲁棒求解线性方程
        if rcond(H_final) < 1e-12
             x_sol = pinv(H_final) * y_vec;
        else
             x_sol = H_final \ y_vec;
        end
        
        q_sol = x_sol(1:L) + 1j * x_sol(L+1:end);
        
        % 5. 亮区能量重缩放 (Rescaling)
        % 论文 Eq.(13): 确保实际输出能量与期望一致
        E_bright_curr = real(q_sol' * Rb * q_sol);
        E_bright_des  = real(q0' * Rb * q0);
        
        % 计算缩放因子 (防止除零)
        if E_bright_curr > 1e-15
            scale_factor = sqrt(E_bright_des / E_bright_curr);
        else
            scale_factor = 1.0; 
        end
        
        % 应用缩放 (注意：WPM-AE 必须保证缩放后仍满足 AE 约束，
        % 但由于我们是在优化 normalized AE，只要比例一致，AE 值本身是不变的)
        q_final = q_sol * scale_factor;
        q_ae(k, :) = q_final.';
        
        % 6. 计算指标
        E_bright = real(q_final' * Rb * q_final);
        E_dark   = real(q_final' * Rd * q_final);
        q_norm_sq = real(q_final' * q_final);
        
        metrics.AC(k) = 10*log10(E_bright / (E_dark + 1e-15));
        
        AG = E_bright / (q_norm_sq + 1e-15);
        metrics.AE(k) = 10*log10(g0 / (AG + 1e-15));
        
        % 计算 SD (Spectral Distortion)
        err_energy = norm(Hb * q_final - db)^2;
        target_energy = norm(db)^2;
        metrics.SD(k) = 10*log10(err_energy / (target_energy + 1e-15));
    end
    fprintf('\n    WPM-AE 完成. 平均 AC: %.2f dB\n', mean(metrics.AC));
end

function f = dual_objective_safe(lam, H0, H1, H2, y)
    % 鲁棒的对偶目标函数
    % 目标: 最小化 f(lam) = y' * H(lam)^-1 * y
    % (对应对偶问题的最大化 -y' H^-1 y)
    
    H_sum = H0 + lam(1)*H1 + lam(2)*H2;
    
    % [检查重点] 正定性检查 (Cholesky 分解)
    % 物理意义: 只有当 H_sum 正定，原问题才有下界，对偶函数才有意义。
    [R, flag] = chol(H_sum);
    
    if flag ~= 0
        % 矩阵非正定！返回惩罚值，迫使优化器调整 lambda
        f = 1e10 + sum(abs(lam)) * 1e5;
        return;
    end
    
    % 如果矩阵正定，使用 Cholesky 快速求解 x = H^-1 * y
    x = R \ (R' \ y);
    
    f = y' * x;
end