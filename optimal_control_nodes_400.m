clear

%
categories = {'Anomic_average', 'Broca_average', 'Conduction_average', 'Global_average', 'None_average'};
% categories = {'Healthy'};
% 参考模型参数
Ar = load("/Users/dongkangli/Desktop/time_series/400/Healthy_average/Average_400_J/average_healthy_mri_400_J.txt"); % 参考模型矩阵 F_r(t)，这里假设为一个负单位阵的一半
% Ar = load('D:\average_healthy_mri_400.txt');
Qr = load("/Users/dongkangli/Desktop/time_series/400/Healthy_average/Average_400_J/average_healthy_mri_400_Sigma.txt"); % 参考模型状态权重矩阵 Q_r(t)
% Qr = load('D:\average_healthy_mri_400.txt');
basePath = "/Users/dongkangli/Desktop/time_series/400/Average";

% 遍历每个类别→每个受试者
for lei = 1:length(categories)

type = categories{lei};
folderPath = fullfile(basePath, type);
fileInfo = dir(folderPath);
filePaths = {};
for i = 1:length(fileInfo)
    % 排除目录和 .DS_Store 文件
    if ~fileInfo(i).isdir && ~strcmp(fileInfo(i).name, '.DS_Store')
        filePath = fullfile(folderPath, fileInfo(i).name);
        filePaths{end+1} = filePath;
    end
end
pattern = '([^/]+)(?=\.[^.]+$)'; 
num_sub = size(filePaths,2);
    for sub = 1:num_sub
        name = regexp(filePaths{sub}, pattern, 'match', 'once');
        fprintf('类别%s %s\n',type,name);
    
        J_folder = fullfile(basePath, [type,'_J']);
        A = load(fullfile(J_folder,name+'_J.txt'));
        Q_noise = load(fullfile(J_folder,name+'_Sigma.txt'));
        
        tic;
        iter = 20; t = 1000:-1:1; n = 400;%系统维度
        control_nodes = [1:n]; %控制的节点
        E_u = [];
        x_history = zeros(iter, length(t), n); xn_history = zeros(iter, length(t), n); xr_history = zeros(iter, length(t), n);
         
        parfor loop = 1:iter
%             loop
           
            % 定义控制矩阵 B
            B = zeros(n,n);
            for i = 1:length(control_nodes)
                B(control_nodes(i), i) = 1;
            end
            Q = eye(n); % 状态权重矩阵 Q(t)
            R = 1;%输入权重R
            
            % 终止条件
            tf = 1000;
            t0 = 1;
            dt = 1;
            P = zeros(n, n); % 终止时间的P矩阵
            Pr = zeros(n, n); 
            P_result = cell(tf,1);
            Pr_result = cell(tf,1);
            % 向后积分
            for t = tf:-dt:t0
                P_result{t} = P;
                Pr_result{t} = Pr;
                % 计算黎卡提方程的右侧
                dP_dt = P * A + A' * P - P * B * (1/R) * B' * P + Q;
                dPr_dt = Pr * Ar + Ar' * Pr - P * B * (1/R) * B' * Pr - Q;
                % 更新P矩阵
                P = P + dP_dt * dt;
                Pr = Pr + dPr_dt * dt;
            end
            
            % 初始状态
            x0 = zeros(n, 1); % 系统初始状态，随机初始化
            xr0 = zeros(n, 1); % 参考模型初始状态，随机初始化
            x = x0;
            xn = x0;
            xr = xr0;
            P = eye(n); % 初始协方差矩阵
            Pr = eye(n); % 参考模型初始协方差矩阵
            u = zeros(length(t), length(control_nodes));
            % 模拟系统
            for k = 1:dt:tf
                % 计算控制输入
                K = R\ (B' * P_result{k}); Kr = -R\ (B' * Pr_result{k});
                u(k, :) = -K * x + Kr * xr;
                xn = xn + dt * (A * xn + mvnrnd(zeros(n,1), Q_noise, 1)');
                x = x + dt * (A * x + B * u(k, :)');
                xr = xr + dt * (Ar * xr + mvnrnd(zeros(n,1), Qr, 1)');
            end
            % E_u
            for i = 1:n
                E_u(loop,i) = sum(u(:,i).^2);
            end
        end
        %         [~, idx] = sort(mean(E_u,1), 'descend'); 
        [~, idx] = sort(mean(E_u,1), 'ascend'); 

        new_iter = 20;
        pick_nodes = 20:20:n;
        parfor loop = 1:new_iter
             loop
            E_u_mean = []; 
            KL = []; 
            KL_num = [];
           kkk = 0;
           for node_num = pick_nodes
                kkk = kkk+1;
                control_nodes = idx(1:node_num);
                % 定义控制矩阵 B
                B_ = eye(n,n);
                B = B_(:,control_nodes);
                Q = eye(n); % 状态权重矩阵 Q(t)
                R = 1;%输入权重R
                
                % 终止条件
                tf = 1000;
                t0 = 1;
                dt = 1;
                P = zeros(n, n); % 终止时间的P矩阵
                Pr = zeros(n, n); 
                P_result = cell(tf,1);
                Pr_result = cell(tf,1);
                % 向后积分
                for t = tf:-dt:t0
                    P_result{t} = P;
                    Pr_result{t} = Pr;
                    % 计算黎卡提方程的右侧
                    dP_dt = P * A + A' * P - P * B * (1/R) * B' * P + Q;
                    dPr_dt = Pr * Ar + Ar' * Pr - P * B * (1/R) * B' * Pr - Q;
                    % 更新P矩阵
                    P = P + dP_dt * dt;
                    Pr = Pr + dPr_dt * dt;
                end
                
                % 初始状态
                x0 = zeros(n, 1); % 系统初始状态，随机初始化
                xr0 = zeros(n, 1); % 参考模型初始状态，随机初始化
                x = x0;
                xn = x0;
                xr = xr0;
                P = eye(n); % 初始协方差矩阵
                Pr = eye(n); % 参考模型初始协方差矩阵
                u = zeros(tf, length(control_nodes));
                x_history = zeros(tf, n); xn_history = zeros(tf, n); xr_history = zeros(tf, n);
                % 模拟系统
                for k = 1:dt:tf
                    % 计算控制输入
                    K = R\ (B' * P_result{k}); Kr = -R\ (B' * Pr_result{k});
                    u(k, :) = -K * x + Kr * xr;
                    xn = xn + dt * (A * xn + mvnrnd(zeros(n,1), Q_noise, 1)');
                    x = x + dt * (A * x + B * u(k, :)');
                    xr = xr + dt * (Ar * xr + mvnrnd(zeros(n,1), Qr, 1)');
                    % 存储状态变量
                    x_history(k, :) = x';
                    xn_history(k, :) = xn';
                    xr_history(k, :) = xr';
                end
                E_u = [];
                for i = 1:length(control_nodes)
                    E_u(i) = sum(u(:,i).^2);
                end
                E_u_mean(kkk) = sum(E_u);
                [K_x, K_num] = cal_KL(x_history,xr_history,xn_history);
                KL(kkk) = K_x; KL_num(kkk) = K_num;
           end
           E_u_all(lei,loop,:) = E_u_mean; 
           KL_all(lei,loop,:) = KL; 
           KL_num_all(lei,loop,:) = KL_num;
        end

        
        fprintf(' ... %.2f s\n',toc);
%         node_ave = ave_control(A);
%         node_modal= modal_control(A);
%         %保存
%         resultPath = fullfile(basePath, [type,'_result']);
%         if ~exist(resultPath, 'dir') 
%              mkdir(resultPath);    
%         end
%         %E_u ，x_history，xn_history，xr_history
%          save(fullfile(resultPath,name + '_optimal_result.mat'), "E_u_all", "KL_all", "KL_num_all");
    end

end


function [K_x, K_num] = cal_KL(x_history,xr_history,xn_history)
    % 假设 x_history, xr_history 和 xn_history 是 m x 400 的矩阵
    [~, n] = size(x_history);  % 假设所有矩阵的列数相同，且为 n（400）
    % 初始化 KL 散度的结果
    KL_x = zeros(1, n);
    KL_xn = zeros(1, n);  % 存储 x_history 和 xn_history 的 KL 散度
    % 平滑因子，防止零值
    epsilon = 100;
    % 计算每一列的 KL 散度
    for i = 1:n
        % 获取第 i 列
        p = xr_history(1:end, i);  % x_history 第 i 列
        q_x = x_history(1:end, i);  % xr_history 第 i 列
        q_xn = xn_history(1:end, i);  % xn_history 第 i 列
        % 添加平滑因子，防止零值
        p = p + epsilon;
        q_x = q_x + epsilon;
        q_xn = q_xn + epsilon;
        % 归一化
        p = p / sum(p);
        q_x = q_x / sum(q_x);
        q_xn = q_xn / sum(q_xn);
        % 计算 KL 散度
        % KL(p || q) = sum(p .* log(p ./ q))
        KL_x(i) = sum(p .* log(p ./ q_x), 'omitnan');  % 忽略 NaN 值
        KL_xn(i) = sum(p .* log(p ./ q_xn), 'omitnan');  % 忽略 NaN 值
    end
    % KL_x = log(KL_x);KL_xn = log(KL_xn);
    K_x = sum(KL_x - KL_xn);
    K_num = sum(KL_x < KL_xn);
end





function [values] = modal_control(A) 
A = A + eye(size(A));   % For functional normalized matrix
[U, T] = schur(A,'real');   % Schur stability
eigVals = diag(T);
N = size(A,1);
phi = zeros(N,1);
for i = 1 : N
    phi(i) = (U(i,:).^2) * (1 - eigVals.^2);
end
values = phi;
end




function [values] = ave_control(A)
% FUNCTION:
%         Returns values of AVERAGE CONTROLLABILITY for each node in a
%         network, given the adjacency matrix for that network. Average
%         controllability measures the ease by which input at that node can
%         steer the system into many easily-reachable states.
%
% INPUT:
%         A is the structural (NOT FUNCTIONAL) network adjacency matrix, 
% 	      such that the simple linear model of dynamics outlined in the 
%	      reference is an accurate estimate of brain state fluctuations. 
%	      Assumes all values in the matrix are positive, and that the 
%	      matrix is symmetric.
%
% OUTPUT:
%         Vector of average controllability values for each node
%
% Bassett Lab, University of Pennsylvania, 2016. 
% Reference: Gu, Pasqualetti, Cieslak, Telesford, Yu, Kahn, Medaglia,
%            Vettel, Miller, Grafton & Bassett, Nature Communications
%            6:8414, 2015.

% A = A./(1+svds(A,1));     % Matrix normalization 
A = A + eye(size(A));   % For functional normalized matrix
[U, T] = schur(A,'real'); % Schur stability
midMat = (U.^2)';
v = diag(T);
P = repmat(diag(1 - v*v'),1,size(A,1));
values = sum(midMat./P)';
end
