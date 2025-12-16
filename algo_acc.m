function [q_acc_all, metrics] = algo_acc(H_bright, H_dark, d_bright, conf)
% ALGO_ACC 声学对比度控制 (Baseline)
% 通过广义特征值分解最大化亮暗能量比。

    fprintf('算法: 运行 ACC (Acoustic Contrast Control)...\n');
    
    [K, Mb, L] = size(H_bright);
    [~, Md, ~] = size(H_dark);
    
    q_acc_all = zeros(K, L);
    metrics.AC = zeros(K, 1);
    metrics.AE = zeros(K, 1);
    metrics.SD = zeros(K, 1);
    
    q0 = zeros(L, 1); q0(conf.target_ls_idx) = 1;
    
    for k = 1:K
        Hb = squeeze(H_bright(k, :, :));
        Hd = squeeze(H_dark(k, :, :));
        db = reshape(d_bright(k, :, :), Mb, 1);
        
        Rb = (Hb' * Hb) / Mb;
        Rd = (Hd' * Hd) / Md;
        
        % 广义特征值分解: Rb * q = lambda * Rd * q
        % 加上微量正则化防止 Rd 奇异
        [V, D] = eig(Rb, Rd + 1e-6*eye(L), 'vector');
        
        [~, idx] = max(real(D));
        q_tmp = V(:, idx);
        
        % 缩放: 匹配目标能量
        E_bright_des = real(q0' * Rb * q0);
        E_b_curr = real(q_tmp' * Rb * q_tmp);
        scale = sqrt(E_bright_des / (E_b_curr + 1e-15));
        q_final = q_tmp * scale;
        
        q_acc_all(k, :) = q_final.';
        
        % Metrics
        Eb = real(q_final' * Rb * q_final);
        Ed = real(q_final' * Rd * q_final);
        metrics.AC(k) = 10*log10(Eb / (Ed + 1e-15));
        
        g0 = E_bright_des / 1; 
        AG = Eb / (real(q_final' * q_final) + 1e-15);
        metrics.AE(k) = 10*log10(g0 / (AG + 1e-15));
        
        err_energy = norm(Hb * q_final - db)^2;
        metrics.SD(k) = 10*log10(err_energy / (norm(db)^2 + 1e-15));
    end
    fprintf('    ACC 完成.\n');
end