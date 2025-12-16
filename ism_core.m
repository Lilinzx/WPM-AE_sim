function h = ism_core(src, rcv, conf, beta, max_order, direct_only)
% ISM_CORE 镜像声源法核心计算引擎 (独立文件)
% 该函数用于计算从源点到接收点的房间脉冲响应(RIR)。
% 独立出来是为了保证"训练数据生成"和"可视化绘图"使用完全一致的物理模型。
%
% 输入:
%   src: 声源坐标 [x, y, z] (米)
%   rcv: 接收点坐标 [x, y, z] (米)
%   conf: 配置结构体 (包含房间尺寸 room_dim, 采样率 fs 等)
%   beta: 墙壁反射系数 (标量)
%   max_order: 最大反射阶数 (可视化时建议设小一点以提高速度)
%   direct_only: 布尔值，如果为 true，只计算直达声 (Order=0)

    Lx = conf.room_dim(1);
    Ly = conf.room_dim(2);
    Lz = conf.room_dim(3);
    
    h = zeros(conf.nsamples_rir, 1);
    c = conf.c;
    fs = conf.fs;
    
    % 坐标系原点修正：假设 conf 中定义的坐标是以房间中心(或者特定参考点)为基准
    % 这里我们假设 conf.ls_coords 等是绝对坐标，直接位于 0~Lx 范围内。
    % 如果你的原始逻辑是相对坐标，请确保输入前已转换。
    % (根据之前的逻辑，这里假设输入已经是房间内的绝对坐标)
    
    s_room = src;
    r_room = rcv;

    % 如果只算直达声，循环范围设为0
    if direct_only
        nx_range = 0; ny_range = 0; nz_range = 0;
    else
        nx_range = -max_order:max_order;
        ny_range = -max_order:max_order;
        nz_range = -max_order:max_order;
    end
    
    for nx = nx_range
        for ny = ny_range
            for nz = nz_range
                
                % 简单的阶数剪枝
                if abs(nx)+abs(ny)+abs(nz) > max_order
                    continue;
                end
                
                % 8个镜像组合 (x,y,z 正反向)
                for qx = 0:1
                for qy = 0:1
                for qz = 0:1
                    
                    % 镜像声源位置计算 (Standard Allen & Berkley Formula)
                    % Rp = [ x +/- xs + 2nxLx, ... ]
                    
                    % X分量
                    if qx==0
                        rx_val = r_room(1) - (s_room(1) + 2*nx*Lx);
                    else
                        rx_val = r_room(1) - (-s_room(1) + 2*nx*Lx);
                    end
                    
                    % Y分量
                    if qy==0
                        ry_val = r_room(2) - (s_room(2) + 2*ny*Ly);
                    else
                        ry_val = r_room(2) - (-s_room(2) + 2*ny*Ly);
                    end
                    
                    % Z分量
                    if qz==0
                        rz_val = r_room(3) - (s_room(3) + 2*nz*Lz);
                    else
                        rz_val = r_room(3) - (-s_room(3) + 2*nz*Lz);
                    end
                    
                    dist = sqrt(rx_val^2 + ry_val^2 + rz_val^2);
                    
                    % 计算延迟和幅值
                    delay = dist / c;
                    sample_idx = round(delay * fs) + 1;
                    
                    if sample_idx <= conf.nsamples_rir
                        % 估算反射次数 (用于衰减计算)
                        % 这种计算并不完美精确，但在ISM中是常用的近似
                        % 每次镜像操作(q=1)代表一次反射，nx代表穿墙次数(2次反射)
                        num_bounces = abs(nx)*2 + qx + abs(ny)*2 + qy + abs(nz)*2 + qz;
                        
                        % 幅值定律: 1/d * beta^bounces
                        amp = (beta^num_bounces) / (4*pi*dist);
                        
                        h(sample_idx) = h(sample_idx) + amp;
                    end
                    
                end
                end
                end
            end
        end
    end
end