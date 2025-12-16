function [q_lwe, metrics] = algo_wpm_lwe(H_bright, H_dark, d_bright, conf, AC_target_dB, LWE_max_dB)
% ALGO_WPM_LWE  WPM with LWE Constraint
% 约束: 限制滤波器权重能量 (LWE, white noise gain) 不超过阈值。
% 类似于传统的鲁棒波束形成。

    fprintf('算法: 运行 WPM-LWE (AC目标=%.0fdB, LWE上限=%.0fdB)...\n', AC_target_dB, LWE_max_dB);
    
    [K, Mb, L] = size(H_bright);
    [~, Md, ~] = size(H_dark);
    
    q_lwe = zeros(K, L);
    metrics.AC = zeros(K, 1);
    metrics.AE = zeros(K, 1);
    metrics.SD = zeros(K, 1);
    
    AC_lin = 10^(AC_target_dB/10);
    LWE_limit_lin = 10^(LWE_max_dB/10); 
    
    q0 = zeros(L, 1); q0(conf.target_ls_idx) = 1;
    Limit_Val = LWE_limit_lin; % 这里假设参考能量为1
    
    opts = optimoptions('fmincon', 'Display', 'off', 'Algorithm', 'sqp');

    msg_len = 0;
    
    for k = 1:K
        if mod(k, 200) == 0 || k == K
            fprintf(repmat('\b', 1, msg_len)); 
            msg = sprintf('    进度: %3d%%', round(k/K*100));
            fprintf(msg); msg_len = length(msg);
        end
        
        Hb = squeeze(H_bright(k, :, :));
        Hd = squeeze(H_dark(k, :, :));
        db = reshape(d_bright(k, :, :), Mb, 1);
        
        Rb = (Hb' * Hb) / Mb;
        Rd = (Hd' * Hd) / Md;
        HbH_Hb = Hb' * Hb;
        HbH_db = Hb' * db;
        
        % 实数矩阵构建
        H0 = [real(HbH_Hb), -imag(HbH_Hb); imag(HbH_Hb), real(HbH_Hb)];
        H0 = H0 + 1e-7 * eye(2*L);
        y_vec = [real(HbH_db); imag(HbH_db)];
        
        % 约束1 AC
        S1 = AC_lin * Rd - Rb;
        H1 = [real(S1), -imag(S1); imag(S1), real(S1)];
        
        % 约束2 LWE: q'q <= Limit
        H2 = eye(2*L); % 单位阵
        c2 = -Limit_Val;
        
        % 目标函数
        objective = @(lam) dual_obj_safe_lwe(lam, H0, H1, H2, y_vec, c2);
        
        try
            lambda_opt = fmincon(objective, [0.1; 0.1], [], [], [], [], [0; 0], [], [], opts);
        catch
            lambda_opt = [0; 0];
        end
        
        H_final = H0 + lambda_opt(1)*H1 + lambda_opt(2)*H2;
        x_sol = H_final \ y_vec;
        
        if any(isnan(x_sol)), x_sol = zeros(2*L, 1); end
        
        q_sol = x_sol(1:L) + 1j * x_sol(L+1:end);
        q_lwe(k, :) = q_sol.';
        
        % 指标
        Eb = real(q_sol' * Rb * q_sol);
        Ed = real(q_sol' * Rd * q_sol);
        metrics.AC(k) = 10*log10(Eb / (Ed + 1e-15));
        
        g0 = real(q0' * Rb * q0);
        AG = Eb / (real(q_sol' * q_sol) + 1e-15);
        metrics.AE(k) = 10*log10(g0 / (AG + 1e-15));
        
        err_energy = norm(Hb * q_sol - db)^2;
        metrics.SD(k) = 10*log10(err_energy / (norm(db)^2 + 1e-15));
    end
    fprintf('\n    WPM-LWE 完成. 平均 AC: %.2f dB\n', mean(metrics.AC));
end

function f = dual_obj_safe_lwe(lam, H0, H1, H2, y, c2)
    H_sum = H0 + lam(1)*H1 + lam(2)*H2;
    [R, flag] = chol(H_sum);
    
    if flag ~= 0
        f = 1e10 + sum(abs(lam))*1e5; % Penalty
        return;
    end
    
    x = R \ (R' \ y);
    term1 = y' * x;
    term2 = lam(2) * c2;
    f = term1 - term2; % Minimized
end