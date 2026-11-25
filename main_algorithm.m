% =========================================================================
% 整合版：迭代双边滤波（KD-tree加速）+ 脊柱候选点提取 + B样条拟合 + Cobb角
% =========================================================================

clc; clear; close all;

%% ==================== 0) 用户参数设置 ====================
% ---- 双边滤波参数 ----
numIterations = 3;      % 迭代平滑次数（可修改）
sigma_s = 3.0;          % 空间标准差 (mm)
sigma_r = 5.0;          % 值域标准差 (mm)
k = 20;                 % 邻域点数（K近邻）

% ---- 脊柱分析参数（原第二段代码中的参数） ----
% 全局字体大小
titleFontSize = 15;      
labelFontSize = 13;     
legendFontSize = 15;    
textFontSize = 15;      
tickFontSize = 14;      

% 切片范围选择
skip_start = 10;   % 忽略前几条切片
skip_end   = 10;   % 忽略后几条切片

% 切片提取参数
Delta     = 4;
dPlatform = 4;
dyBin     = 2;
minBins   = 25;
gapBreak  = 7;
surfaceMode = 'max';
qlenth = 140;           % 轮廓最小长度 (mm)
yRangeHalfWidth = 55;   % Y范围筛选半宽（mm）
mergeDist = 1;          % 合并最近点对的Y坐标距离阈值（mm）
cobb_dist = 15;         % Cobb角计算时两点间最小X距离阈值 (mm)，可自行调整

% ---- 文件路径（请根据实际情况修改） ----
inputFile = 'sample_data.txt';
% 输出文件名自动生成
[path, name, ext] = fileparts(inputFile);
outputFile = fullfile(path, [name, '_双边滤波迭代', num2str(numIterations), '次', ext]);

%% ==================== 1) 读取原始点云 ====================
P = readmatrix(inputFile);
P = P(:,1:3);
P = P(~any(isnan(P),2),:);
P_original = P;  % 备份原始点云用于对比

fprintf('原始点云点数：%d\n', size(P,1));
fprintf('开始迭代双边滤波，迭代次数：%d\n', numIterations);

%% ==================== 2) 迭代双边滤波（KD-tree加速） ====================
for iter = 1:numIterations
    fprintf('  第 %d 次迭代...\n', iter);
    
    % 构建KD-Tree索引（基于当前点云的XY坐标）
    idx_tree = createns(P(:,1:2), 'NSMethod', 'kdtree');
    
    % 批量搜索每个点的K近邻（K+1包含自身）
    [neighbor_indices_all, ~] = knnsearch(idx_tree, P(:,1:2), 'K', k+1);
    
    % 初始化新点云
    P_new = zeros(size(P));
    
    % 对每个点进行双边滤波
    for i = 1:size(P,1)
        % 获取K个邻居（排除自身）
        idx = neighbor_indices_all(i, 2:end);
        if length(idx) < 3
            P_new(i,:) = P(i,:);
            continue;
        end
        
        % 邻居点坐标
        neigh_xy = P(idx, 1:2);
        neigh_z  = P(idx, 3);
        
        % 空间距离权重
        dist_xy = sqrt((neigh_xy(:,1) - P(i,1)).^2 + (neigh_xy(:,2) - P(i,2)).^2);
        w_s = exp(-dist_xy.^2 / (2 * sigma_s^2));
        
        % 值域权重（高度差）
        dz = abs(neigh_z - P(i,3));
        w_r = exp(-dz.^2 / (2 * sigma_r^2));
        
        % 综合权重并归一化
        w = w_s .* w_r;
        w = w / sum(w);
        
        % 加权平均得到新Z值
        P_new(i,3) = sum(w .* neigh_z);
        P_new(i,1:2) = P(i,1:2);   % X,Y坐标保持不变
    end
    
    % 更新点云用于下一次迭代
    P = P_new;
end

P_smoothed = P;  % 最终平滑结果

%% ==================== 3) 可视化滤波前后对比 ====================
figure('Name', '点云平滑前后对比（迭代双边滤波）', 'Color','w');
subplot(1,2,1);
scatter3(P_original(:,1), P_original(:,2), P_original(:,3), 1, P_original(:,3), 'filled');
title('原始点云（彩色=Z值）'); axis equal; colorbar; colormap(jet); view(3);
xlabel('X (mm)'); ylabel('Y (mm)'); zlabel('Z (mm)');

subplot(1,2,2);
scatter3(P_smoothed(:,1), P_smoothed(:,2), P_smoothed(:,3), 1, P_smoothed(:,3), 'filled');
title(sprintf('迭代%d次双边滤波平滑后', numIterations)); axis equal; colorbar; colormap(jet); view(3);
xlabel('X (mm)'); ylabel('Y (mm)'); zlabel('Z (mm)');

%% ==================== 4) 保存滤波后的点云 ====================
writematrix(P_smoothed, outputFile, 'Delimiter', 'tab');
fprintf('平滑后点云已保存为：%s\n', outputFile);

%% ==================== 5) 脊柱分析（基于滤波后点云） ====================
% 使用滤波后的点云P_smoothed替换原有P，继续执行脊柱候选点提取等分析

fprintf('\n===== 开始脊柱候选点提取与分析（基于平滑后点云） =====\n');

% 重新赋值给后续使用的变量
P = P_smoothed;           
xAll = P(:,1);
yAll = P(:,2);
zAll = P(:,3);
P_original = P;   % 用于显示的点云改为平滑后点云（也可保留原始对比，这里用平滑后）

xlimit = [min(xAll), max(xAll)];

tic;

%% 5.1) 生成所有X轴切片列表
xcList = xlimit(1):dPlatform:xlimit(2);
totalSlices = numel(xcList);
fprintf('总切片数：%d\n', totalSlices);
fprintf('切片范围选择：忽略前 %d 条，忽略后 %d 条，有效切片序号范围 [%d, %d]\n', ...
    skip_start, skip_end, skip_start+1, totalSlices - skip_end);

%% 5.2) 初始化存储所有切片极小值点的变量
minima_by_slice = cell(totalSlices, 1);

%% 5.3) 遍历所有切片，提取极小值点
for sliceIdx = 1:totalSlices
    try
        xc = xcList(sliceIdx);
        idx = abs(xAll - xc) < Delta/2;
        ys = yAll(idx);
        zs = zAll(idx);

        if numel(ys) < 200
            fprintf('第%d条切片：点数不足（%d），跳过\n', sliceIdx, numel(ys));
            continue;
        end

        yMin = min(ys); yMax = max(ys);
        edges = yMin:dyBin:yMax;
        if numel(edges) < (minBins+1)
            fprintf('第%d条切片：bin数不足，跳过\n', sliceIdx);
            continue;
        end

        [~,~,bin] = histcounts(ys, edges);
        yb = edges(1:end-1) + dyBin/2;
        zb = nan(size(yb));

        for b = 1:numel(yb)
            idb = (bin == b);
            if ~any(idb), continue; end
            zvals = zs(idb);
            switch surfaceMode
                case 'max', zb(b) = max(zvals);
                case 'min', zb(b) = min(zvals);
                otherwise, zb(b) = median(zvals);
            end
        end

        ok = ~isnan(zb);
        yb = yb(ok);
        zb = zb(ok);

        if numel(yb) < minBins
            fprintf('第%d条切片：有效点不足，跳过\n', sliceIdx);
            continue;
        end

        w = max(5, round(numel(zb)*0.05));
        zb = smoothdata(zb, 'movmean', w);

        d = sqrt(diff(yb).^2 + diff(zb).^2);
        breakIdx = find(d > gapBreak);
        starts = [1, breakIdx+1];
        ends   = [breakIdx, numel(yb)];
        valid = (ends - starts + 1) >= 2;
        starts = starts(valid);
        ends   = ends(valid);
        if isempty(starts)
            fprintf('第%d条切片：无有效段，跳过\n', sliceIdx);
            continue;
        end

        best = 1; bestScore = -inf;
        for s = 1:numel(starts)
            i1 = starts(s); i2 = ends(s);
            score = abs(yb(i2) - yb(i1));
            if score > bestScore
                bestScore = score;
                best = s;
            end
        end
        i1 = starts(best); i2 = ends(best);
        yb_back = yb(i1:i2);
        zb_back = zb(i1:i2);

        y_span = max(yb_back) - min(yb_back);
        if y_span < qlenth
            fprintf('第%d条切片：背部轮廓Y方向跨度%.1f mm < %d mm，放弃该切片\n', sliceIdx, y_span, qlenth);
            continue;
        end

        y_fine = linspace(min(yb_back), max(yb_back), 1000);
        z_fine = [];
        method_used = '';

        try
            if length(yb_back) >= 4
                pp = spapi(4, yb_back, zb_back);
                z_fine = fnval(pp, y_fine);
                method_used = 'B-spline';
            else
                error('点数不足4');
            end
        catch
            z_fine = interp1(yb_back, zb_back, y_fine, 'pchip');
            method_used = 'pchip';
        end

        if isempty(z_fine)
            fprintf('第%d条切片：插值失败，跳过\n', sliceIdx);
            continue;
        end

        y_mean_fit = mean(y_fine);
        if y_mean_fit >= min(y_fine) && y_mean_fit <= max(y_fine)
            z_at_mean = interp1(y_fine, z_fine, y_mean_fit, 'pchip');
        else
            [~, idx_closest] = min(abs(y_fine - y_mean_fit));
            z_at_mean = z_fine(idx_closest);
            fprintf('第%d条切片：y_mean_fit超出范围，使用最近点\n', sliceIdx);
        end

        y_min_limit = y_mean_fit - yRangeHalfWidth;
        y_max_limit = y_mean_fit + yRangeHalfWidth;

        dz_dy = gradient(z_fine, y_fine);
        minima_2d = [];
        for i = 2:length(dz_dy)-1
            if dz_dy(i-1) < -1e-6 && dz_dy(i+1) > 1e-6
                minima_2d = [minima_2d; y_fine(i), z_fine(i)];
            end
        end

        if ~isempty(minima_2d)
            y_in_range = (minima_2d(:,1) >= y_min_limit) & (minima_2d(:,1) <= y_max_limit);
            minima_2d = minima_2d(y_in_range, :);
        end

        if ~isempty(minima_2d) && size(minima_2d,1) >= 2
            merged = true;
            while merged
                merged = false;
                n = size(minima_2d,1);
                minDist = inf;
                bestPair = [0,0];
                for i = 1:n-1
                    for j = i+1:n
                        dist = abs(minima_2d(i,1) - minima_2d(j,1));
                        if dist < minDist
                            minDist = dist;
                            bestPair = [i, j];
                        end
                    end
                end
                if minDist <= mergeDist
                    newY = (minima_2d(bestPair(1),1) + minima_2d(bestPair(2),1)) / 2;
                    newZ = (minima_2d(bestPair(1),2) + minima_2d(bestPair(2),2)) / 2;
                    keepIdx = setdiff(1:n, bestPair);
                    minima_2d = [minima_2d(keepIdx,:); newY, newZ];
                    merged = true;
                    fprintf('  第%d条切片：合并距离=%.3f mm的两个极小值点\n', sliceIdx, minDist);
                end
            end
        end

        if ~isempty(minima_2d) && size(minima_2d,1) > 2
            dist_to_mean = abs(minima_2d(:,1) - y_mean_fit);
            [~, idx_sorted] = sort(dist_to_mean);
            minima_2d = minima_2d(idx_sorted(1:2), :);
            fprintf('  第%d条切片：点数超过2，保留距离拟合曲线Y均值最近的两个点\n', sliceIdx);
        end

        if ~isempty(minima_2d)
            minima_3d = [xc * ones(size(minima_2d,1),1), minima_2d];
            minima_by_slice{sliceIdx} = minima_3d;
            fprintf('第%d条切片（%s）：初步保留%d个极小值点（轮廓跨度%.1f mm，拟合曲线Y均值%.1f）\n', ...
                sliceIdx, method_used, size(minima_2d,1), y_span, y_mean_fit);
        else
            minima_by_slice{sliceIdx} = [];
            fprintf('第%d条切片（%s）：无极小值点\n', sliceIdx, method_used);
        end

    catch ME
        fprintf('第%d条切片处理出错：%s\n', sliceIdx, ME.message);
        minima_by_slice{sliceIdx} = [];
        continue;
    end
end

%% 5.4) 相邻切片距离比较
for sliceIdx = 1:totalSlices
    points = minima_by_slice{sliceIdx};
    if isempty(points) || size(points,1) ~= 2
        continue;
    end
    
    prev_points = [];
    next_points = [];
    if sliceIdx > 1
        prev_points = minima_by_slice{sliceIdx-1};
    end
    if sliceIdx < totalSlices
        next_points = minima_by_slice{sliceIdx+1};
    end
    
    if isempty(prev_points) && isempty(next_points)
        y_mean_cur = mean(points(:,2));
        dist_to_mean = abs(points(:,2) - y_mean_cur);
        [~, idx] = min(dist_to_mean);
        minima_by_slice{sliceIdx} = points(idx, :);
        fprintf('第%d条切片：两个点无相邻切片参考，保留离当前Y均值最近的点\n', sliceIdx);
        continue;
    end
    
    best_idx = 1;
    best_diff = inf;
    best_sum = inf;
    for p = 1:2
        cur_point = points(p, :);
        dist_prev = inf;
        if ~isempty(prev_points)
            dists_prev = sqrt(sum((prev_points - cur_point).^2, 2));
            dist_prev = min(dists_prev);
        end
        dist_next = inf;
        if ~isempty(next_points)
            dists_next = sqrt(sum((next_points - cur_point).^2, 2));
            dist_next = min(dists_next);
        end
        if isempty(prev_points) && ~isempty(next_points)
            diff_val = dist_next;
            sum_val = dist_next;
        elseif ~isempty(prev_points) && isempty(next_points)
            diff_val = dist_prev;
            sum_val = dist_prev;
        elseif isempty(prev_points) && isempty(next_points)
            diff_val = 0;
            sum_val = 0;
        else
            diff_val = abs(dist_prev - dist_next);
            sum_val = dist_prev + dist_next;
        end
        if diff_val < best_diff - 1e-6
            best_diff = diff_val;
            best_sum = sum_val;
            best_idx = p;
        elseif abs(diff_val - best_diff) < 1e-6
            if sum_val < best_sum - 1e-6
                best_sum = sum_val;
                best_idx = p;
            end
        end
    end
    minima_by_slice{sliceIdx} = points(best_idx, :);
    fprintf('第%d条切片：两个点中保留与相邻切片距离更均衡的点\n', sliceIdx);
end

%% 5.5) 合并所有切片结果，应用切片范围选择
all_minima_3d = [];
for sliceIdx = 1:totalSlices
    if sliceIdx <= skip_start || sliceIdx > totalSlices - skip_end
        continue;
    end
    if ~isempty(minima_by_slice{sliceIdx})
        all_minima_3d = [all_minima_3d; minima_by_slice{sliceIdx}];
    end
end

fprintf('有效切片范围 [%d, %d] 内共筛选出 %d 个极值点。\n', ...
    skip_start+1, totalSlices - skip_end, size(all_minima_3d,1));

%% 5.6) 绘制原始点云 + 所有极小值点（三维）
figure('Name', '所有切片极小值点（三维）', 'Color','w');
hold on; grid on; view(3);
scatter3(P_original(:,1), P_original(:,2), P_original(:,3), 1, [0.8 0.8 0.8], 'filled');
if ~isempty(all_minima_3d)
    scatter3(all_minima_3d(:,1), all_minima_3d(:,2), all_minima_3d(:,3), ...
             50, 'r', 'filled', 'MarkerEdgeColor', 'k');
end
xlabel('X (mm)', 'FontSize', labelFontSize);
ylabel('Y (mm)', 'FontSize', labelFontSize);
zlabel('Z (mm)', 'FontSize', labelFontSize);
title(sprintf('Smoothed point cloud + effective slice minimum points（total %d）', size(all_minima_3d,1)), ...
      'FontSize', titleFontSize);
axis equal;
legend('Smoothed point cloud', 'minimum points', 'Location', 'best', 'FontSize', legendFontSize);
set(gca, 'FontSize', tickFontSize);

%% 5.7) 基于最终极值点的B样条拟合曲线，展示三个视角并计算Cobb角
if size(all_minima_3d, 1) >= 4
    [x_sort, idx_sort] = sort(all_minima_3d(:,1));
    y_sort = all_minima_3d(idx_sort, 2);
    z_sort = all_minima_3d(idx_sort, 3);
    
    x_fine = linspace(min(x_sort), max(x_sort), 500);
    
    try
        sp_y = spapi(4, x_sort, y_sort);
        y_fine = fnval(sp_y, x_fine);
        sp_z = spapi(4, x_sort, z_sort);
        z_fine = fnval(sp_z, x_fine);
        method_curve = 'B-spline';
    catch
        warning('B样条拟合脊柱曲线失败，使用pchip插值');
        y_fine = interp1(x_sort, y_sort, x_fine, 'pchip');
        z_fine = interp1(x_sort, z_sort, x_fine, 'pchip');
        method_curve = 'pchip';
    end
    
    % 窗口1：三维曲线
    figure('Name', 'B-Spline fitting of spinal midline(3D)', 'Color','w');
    hold on; grid on; view(3);
    scatter3(P_original(:,1), P_original(:,2), P_original(:,3), 1, [0.8 0.8 0.8], 'filled');
    scatter3(all_minima_3d(:,1), all_minima_3d(:,2), all_minima_3d(:,3), 50, 'b', 'filled', 'MarkerEdgeColor','k');
    plot3(x_fine, y_fine, z_fine, 'r-', 'LineWidth', 2);
    xlabel('X (mm)', 'FontSize', labelFontSize);
    ylabel('Y (mm)', 'FontSize', labelFontSize);
    zlabel('Z (mm)', 'FontSize', labelFontSize);
    title(sprintf('Spinal midline fitting（%s）', method_curve), 'FontSize', titleFontSize);
    legend('Smoothed cloud', 'Discrete extremum', 'fitting curve', 'Location', 'best', 'FontSize', legendFontSize);
    axis equal;
    set(gca, 'FontSize', tickFontSize);
    
    % 窗口2：XZ视角（侧视图）
    figure('Name', 'Spinal midline XZ projection', 'Color','w');
    hold on; grid on;
    plot(x_fine, z_fine, 'r-', 'LineWidth', 2);
    scatter(x_sort, z_sort, 30, 'b', 'filled');
    xlabel('X (mm)', 'FontSize', labelFontSize);
    ylabel('Z (mm)', 'FontSize', labelFontSize);
    title(sprintf('Sagittal projection（%s）', method_curve), 'FontSize', titleFontSize);
    legend('fitting curve', 'discrete points', 'Location', 'best', 'FontSize', legendFontSize);
    set(gca, 'FontSize', tickFontSize);
    
    % 窗口3：XY视角（冠状面俯视图）+ Cobb角计算
    figXY = figure('Name', 'Spinal midline XY projection（including Cobb angles）', 'Color','w');
    hold on; grid on; axis equal;
    plot(x_fine, y_fine, 'r-', 'LineWidth', 2);
    scatter(x_sort, y_sort, 30, 'b', 'filled');
    xlabel('X (mm)', 'FontSize', labelFontSize);
    ylabel('Y (mm)', 'FontSize', labelFontSize);
    title(sprintf('Coronal projection（%s）', method_curve), 'FontSize', titleFontSize);
    legend('fitting curve', 'discrete points', 'Location', 'best', 'FontSize', legendFontSize);
    set(gca, 'FontSize', tickFontSize);
    
    % ========== Cobb角计算 ==========
    dx = diff(x_fine);
    dy = diff(y_fine);
    slope = dy ./ dx;
    theta_deg = atan(slope) * 180 / pi;
    x_mid = (x_fine(1:end-1) + x_fine(2:end)) / 2;
    
    local_max = false(size(theta_deg));
    local_min = false(size(theta_deg));
    for i = 2:length(theta_deg)-1
        if theta_deg(i) > theta_deg(i-1) && theta_deg(i) > theta_deg(i+1)
            local_max(i) = true;
        elseif theta_deg(i) < theta_deg(i-1) && theta_deg(i) < theta_deg(i+1)
            local_min(i) = true;
        end
    end
    
    idx_max = find(local_max);
    idx_min = find(local_min);
    all_extrema = sort([idx_max, idx_min]);
    
    if length(all_extrema) >= 2
        cobb_candidates = [];
        for k = 1:length(all_extrema)-1
            idx1 = all_extrema(k);
            idx2 = all_extrema(k+1);
            delta = abs(theta_deg(idx2) - theta_deg(idx1));
            x_dist = x_mid(idx2) - x_mid(idx1);
            if delta > 5 && x_dist > cobb_dist
                cobb_candidates = [cobb_candidates; idx1, idx2, delta, x_mid(idx1), x_mid(idx2)];
            end
        end
        
        if ~isempty(cobb_candidates)
            [~, sortIdx] = sort(cobb_candidates(:,3), 'descend');
            top4 = cobb_candidates(sortIdx, :);
            n_top = min(4, size(top4,1));
            top4 = top4(1:n_top, :);
            
            fprintf('\n========== Cobb角计算结果（XY投影） ==========\n');
            for k = 1:n_top
                idx1 = top4(k,1);
                idx2 = top4(k,2);
                cobb_val = top4(k,3);
                x_start = top4(k,4);
                x_end = top4(k,5);
                y_start = y_fine(idx1);
                y_end = y_fine(idx2);
                fprintf('第%d大Cobb角：%.2f° （X范围 [%.1f, %.1f] mm，距离=%.1f mm）\n', ...
                    k, cobb_val, x_start, x_end, x_end - x_start);
                
                figure(figXY);
                plot([x_start, x_end], [y_start, y_end], 'g-', 'LineWidth', 2);
                mid_x = (x_start + x_end) / 2;
                mid_y = (y_start + y_end) / 2;
                text(mid_x, mid_y + 5, sprintf('%.1f°', cobb_val), ...
                    'Color', 'red', 'FontSize', textFontSize, 'FontWeight', 'bold', ...
                    'HorizontalAlignment', 'center', 'BackgroundColor', 'w');
            end
            legend('fitting curve', 'discrete points', 'Cobb corner line segment', ...
                   'Location', 'best', 'FontSize', legendFontSize);
            fprintf('=============================================\n');
        else
            fprintf('未检测到有效Cobb角（角度变化<5°或X距离≤%d mm）。\n', cobb_dist);
        end
    else
        fprintf('切线角度极值点不足，无法计算Cobb角。\n');
    end
else
    fprintf('有效范围内极值点数量不足4个，无法进行B样条曲线拟合。\n');
end

toc;

% ==================== 结束 ====================
disp('全部处理完成。');