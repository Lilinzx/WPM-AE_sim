function [H_bright, H_dark, d_bright_direct] = generate_rirs(conf)
% GENERATE_RIRS 生成仿真所需的房间脉冲响应和传递函数
% 调用: ism_core.m

    fprintf('M2: 正在计算 RIR (使用独立物理引擎 ism_core)...\n');
    
    % ISM 参数
    V = prod(conf.room_dim);
    S = 2 * (conf.room_dim(1)*conf.room_dim(2) + ...
             conf.room_dim(2)*conf.room_dim(3) + ...
             conf.room_dim(1)*conf.room_dim(3));
    alpha = 0.161 * V / (S * conf.rt60);
    beta = sqrt(1 - alpha);
    
    order = 12; % 训练用的反射阶数 (高精度)
    
    % 预分配
    rir_b = zeros(conf.nsamples_rir, conf.Mb, conf.L);
    rir_d = zeros(conf.nsamples_rir, conf.Md, conf.L);
    rir_tgt = zeros(conf.nsamples_rir, conf.Mb, 1);
    
    % 1. 亮区 RIR
    for l = 1:conf.L
        src = conf.ls_coords(l, :);
        for m = 1:conf.Mb
            rcv = conf.bright_coords(m, :);
            rir_b(:, m, l) = ism_core(src, rcv, conf, beta, order, false);
        end
    end
    
    % 2. 暗区 RIR
    for l = 1:conf.L
        src = conf.ls_coords(l, :);
        for m = 1:conf.Md
            rcv = conf.dark_coords(m, :);
            rir_d(:, m, l) = ism_core(src, rcv, conf, beta, order, false);
        end
    end
    
    % 3. 目标声场 (仅直达声)
    src_tgt = conf.ls_coords(conf.target_ls_idx, :);
    for m = 1:conf.Mb
        rcv = conf.bright_coords(m, :);
        rir_tgt(:, m, 1) = ism_core(src_tgt, rcv, conf, beta, 0, true);
    end
    
    % 4. FFT 变换
    K = conf.J;
    H_bright = permute(fft(rir_b, K, 1), [1, 2, 3]);
    H_dark   = permute(fft(rir_d, K, 1), [1, 2, 3]);
    d_bright_direct = permute(fft(rir_tgt, K, 1), [1, 2, 3]);
    
    fprintf('    RIR 生成完毕.\n');
end