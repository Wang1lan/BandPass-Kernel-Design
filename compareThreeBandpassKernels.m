% compareThreeBandpassKernels.m
% 依赖：仓库 code_bandPassFilter/bandpassFilter.m，以及 inspectBandpassKernel.m
% 实验起点参数，不是已标定的最终频带；请以运行后输出的实际频响为准。
clear; clc; close all;

names = {'低频大尺度', '中频纹理'};
% 列顺序：fRadius, sigmaHighFreq, sigmaLowFreq, scaleLowFreq
params = [5, 0.55, 0.2, 0.85; 3, 1, 0.5, 0.95]; detailGains = [6, 3];
detailSign = 1;

periods = [2 3 4 5 6 8 10 12 16 20 24 32];
anglesDeg = [0 45 90];

results = cell(1, 2);
kernels = cell(1, 2);
for k = 1:2
    p = params(k,:);
    h = bandpassFilter(p(1),p(2),p(3),p(4));
    kernels{k} = h;
    fprintf('\n========== %s (%dx%d) ==========%s', ...
        names{k},2*p(1)+1,2*p(1)+1,newline);
    results{k} = inspectBandpassKernel(h,periods,anglesDeg);
    % 分析函数保留绝对增益；仅将条纹对比图改为按各核二维峰值归一化。
    stripePeak = max(results{k}.magnitude2D(:));
    close(results{k}.figures(2));
    results{k}.figures(2) = plotNormalizedStripes(results{k},stripePeak,names{k});
    fprintf('%s: 条纹归一化基准（二维幅频峰值）= %.6g\n',names{k},stripePeak);
    fprintf('半功率(-3dB)主频带：r=%.3f ~ %.3f，条纹周期=%.1f ~ %.1f px\n', ...
        results{k}.mainBand,results{k}.periodRange);
    if any(results{k}.bandTouchesLimit)
        warning('%s 的 -3dB 主频带触及 r=0 或 r=1，需检查是否偏离带通目标。',names{k});
    end
end

% 构造合成核：与 demo.m 一致，包含原图直通项与两个细节核的带符号增益叠加。
maxSize = max(cellfun(@(x) size(x,1), kernels));
compositeDetailKernel = zeros(maxSize, maxSize);
for k = 1:2
    compositeDetailKernel = compositeDetailKernel + ...
        detailSign * detailGains(k) * centerPadKernel(kernels{k}, maxSize);
end
compositeKernel = compositeDetailKernel;
centerIdx = floor(maxSize/2) + 1;
compositeKernel(centerIdx, centerIdx) = compositeKernel(centerIdx, centerIdx) + 1;

compResp = analyzeKernelResponse(compositeKernel, 1024);

% 合成核对不同周期与方向条纹的滤波增益。
% inspectBandpassKernel 保留实际增益；仅对条纹可视化及归一化增益表做峰值归一化。
% 合成核和下方合成频响仍是 1 + gain1*H1 + gain2*H2，不改变实际增强逻辑。
fprintf('\n========== 合成核条纹滤波增益（含原图直通） ==========\n');
compositeStripeResult = inspectBandpassKernel(compositeKernel,periods,anglesDeg);
set(compositeStripeResult.figures(1),'Name','合成核频响与条纹增益（实际值）');
compositeStripePeak = max(compositeStripeResult.magnitude2D(:));
close(compositeStripeResult.figures(2));
compositeStripeResult.figures(2) = plotNormalizedStripes( ...
    compositeStripeResult,compositeStripePeak,'合成核');

% 同时输出原始增益和归一化增益：原始增益反映增强强度，归一化增益反映频带形状。
fprintf('\n合成核二维幅频峰值（条纹归一化基准）= %.6g\n',compositeStripePeak);
fprintf('合成核条纹增益：理论绝对/实测绝对 | 理论归一化/实测归一化\n');
fprintf('%8s %8s', 'T(px)', 'r=2/T');
for a = 1:numel(anglesDeg)
    fprintf(' | %g°绝对(理/实) %g°归一(理/实)',anglesDeg(a),anglesDeg(a));
end
fprintf('\n');
for p = 1:numel(periods)
    fprintf('%8g %8.3f', periods(p),2/periods(p));
    for a = 1:numel(anglesDeg)
        theoretical = compositeStripeResult.theoreticalGain(p,a);
        measured = compositeStripeResult.measuredGain(p,a);
        fprintf(' | %6.3f/%6.3f %6.3f/%6.3f', ...
            theoretical,measured,theoretical/compositeStripePeak,measured/compositeStripePeak);
    end
    fprintf('\n');
end

% 两核径向频响同图比较：看主峰分离、带宽以及相邻频段重叠。
figure('Name','两个带通核及合成响应比较','Color','w');
hold on;
colors = lines(3);
for k = 1:2
    r = results{k};
    curve = r.radialRMS / max(r.radialRMS);
    plot(r.radialFrequency,curve,'LineWidth',2,'Color',colors(k,:), ...
         'DisplayName',sprintf('%s (%dx%d), 峰值 %.3f', ...
         names{k},size(r.kernel,1),size(r.kernel,2),r.peakFrequency));
end
compCurve = compResp.radialRMS / max(compResp.radialRMS);
plot(compResp.radialFrequency, compCurve, 'LineWidth', 2.5, 'Color', colors(3,:), ...
    'DisplayName', sprintf('合成响应 (1 + %gH_1 + %gH_2)', detailGains(1), detailGains(2)));
yline(1/sqrt(2),'k--','-3dB');
xlim([0 1]); ylim([0 1.05]);
xlabel('归一化径向频率 r（条纹周期 T=2/r 像素）');
ylabel('径向 RMS / 各响应自身峰值（仅比较频段形状）');
title('两个带通核与合成响应的归一化频响对比');
legend('Location','best'); grid on;

% 新增：合成核的二维频率响应热力图
figure('Name','合成核二维频率响应热力图','Color','w');
imagesc(compResp.freqAxis, compResp.freqAxis, compResp.responseNorm);
axis image;
set(gca,'YDir','normal');
colormap(jet);
colorbar;
xlabel('归一化频率 f_x / (Nyquist)');
ylabel('归一化频率 f_y / (Nyquist)');
title(sprintf('合成核二维频率响应热力图 |1 + %gH_1 + %gH_2|（归一化）', ...
    detailGains(1), detailGains(2)));

fprintf('\n======= 两核汇总（请依据实际响应决定是否修改参数） =======\n');
for k = 1:2
    r = results{k};
    fprintf('%s: %dx%d; peak=%.3f; -3dB=[%.3f, %.3f]; T=[%.1f, %.1f]px; peakGain=%.3f; whiteNoiseGain=%.3f\n', ...
        names{k},size(r.kernel,1),size(r.kernel,2),r.peakFrequency, ...
        r.mainBand,r.periodRange,max(r.radialRMS),norm(r.kernel(:),2));
end
fprintf('合成响应: %dx%d; peak=%.3f; peakGain=%.3f; spatialL2=%.3f\n', ...
    size(compositeKernel,1), size(compositeKernel,2), compResp.peakFrequency, ...
    max(compResp.radialRMS), norm(compositeKernel(:),2));

function fig = plotNormalizedStripes(result,peakGain,figLabel)
% 条纹归一化只用于显示：以该核二维幅频峰值为 1，不改变滤波核或实测绝对增益。
% 去掉直流偏置后叠加 0.5，保证不同周期共用相同灰度色标。
assert(isfinite(peakGain) && peakGain > 0,'条纹归一化基准必须为正数。');
periods = result.periods;
nPeriod = numel(periods);
fig = figure('Name',[figLabel '：归一化条纹增益'],'Color','w', ...
    'Position',[60 120 max(1100,140*nPeriod) 430]);
layout = tiledlayout(2,nPeriod,'TileSpacing','compact','Padding','compact');
showRows = 1:min(64,size(result.inputStripes{1,1},1));
showCols = 1:min(64,size(result.inputStripes{1,1},2));
for p = 1:nPeriod
    inputStripe = result.inputStripes{1,p}(showRows,showCols);
    outputFull = result.detailStripes{1,p};
    outputStripe = outputFull(showRows,showCols);
    outputStripe = 0.5+(outputStripe-mean(outputFull(:)))/peakGain;
    nexttile(p);
    imagesc(inputStripe,[0.25 0.75]);
    axis image off;
    title(sprintf('T=%g px',periods(p)));
    nexttile(nPeriod+p);
    imagesc(outputStripe,[0.25 0.75]);
    axis image off;
    title(sprintf('归一化增益 %.3f',result.measuredGain(p,1)/peakGain));
end
colormap(fig,gray(256));
title(layout,sprintf([ ...
    '%s：%g°条纹；上：相同幅值输入；下：去直流、按峰值归一化后的输出\n' ...
    '峰值基准 |H|_{max}=%.4g；所有条纹统一灰度范围 [0.25, 0.75]'], ...
    figLabel,result.anglesDeg(1),peakGain));
end

function out = centerPadKernel(kernel, targetSize)
% 将较小奇数方核按中心对齐，零填充到 targetSize×targetSize。
currentSize = size(kernel,1);
assert(size(kernel,1) == size(kernel,2), 'kernel 必须为方阵。');
assert(mod(currentSize,2) == 1 && mod(targetSize,2) == 1, '尺寸必须为奇数。');
assert(currentSize <= targetSize, 'targetSize 必须不小于当前核尺寸。');
out = zeros(targetSize, targetSize);
offset = (targetSize - currentSize)/2;
idx = offset + (1:currentSize);
out(idx, idx) = kernel;
end

function resp = analyzeKernelResponse(kernel, fftSize)
% 计算二维频率响应热力图，以及径向 RMS 频响曲线。
if nargin < 2
    fftSize = 1024;
end

H = fftshift(fft2(kernel, fftSize, fftSize));
mag = abs(H);

% FFT 频率网格：-1 ~ 1-2/N，避免 linspace(-1,1,N) 引入频率偏差。
axisVals = (-fftSize/2:fftSize/2-1) * (2/fftSize);
[fx, fy] = meshgrid(axisVals, axisVals);
r = hypot(fx, fy);
radialFrequency = linspace(0, 1, 400);
radialRMS = zeros(size(radialFrequency));

for i = 1:numel(radialFrequency)
    if i == 1
        mask = r <= radialFrequency(2) / 2;
    else
        leftEdge = (radialFrequency(i-1) + radialFrequency(i)) / 2;
        if i == numel(radialFrequency)
            % 末端仅统计单位圆内的 Nyquist 径向频率，不能延伸到 sqrt(2)。
            % 原来使用 inf 会混入单位圆外角落的低响应点，造成 r=1 突降。
            rightEdge = 1;
        else
            rightEdge = (radialFrequency(i) + radialFrequency(i+1)) / 2;
        end
        mask = (r > leftEdge) & (r <= rightEdge);
    end
    vals = mag(mask);
    if isempty(vals)
        radialRMS(i) = NaN;
    else
        radialRMS(i) = sqrt(mean(vals.^2));
    end
end

[~, peakIdx] = max(radialRMS);
resp.kernel = kernel;
resp.response2D = mag;
resp.responseNorm = mag / max(mag(:));
resp.freqAxis = axisVals;
resp.radialFrequency = radialFrequency;
resp.radialRMS = radialRMS;
resp.peakFrequency = radialFrequency(peakIdx);
end
