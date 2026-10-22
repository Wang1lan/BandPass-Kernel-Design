function result = inspectBandpassKernel(h, periods, anglesDeg)
% inspectBandpassKernel 单核的“频率—像素尺度—细节输出”调试。
% 用法：h=bandpassFilter(3,1.0,0.3,0.8); result=inspectBandpassKernel(h);
% 或：load filter_fx.mat f; result=inspectBandpassKernel(double(f)/8192);
% h：最终实数奇数尺寸核；函数不改变核增益，不执行额外归一化。
% periods：条纹周期(px)，默认[2 3 4 6 8 12 16 24 32]。
% anglesDeg：亮度变化方向(度)，默认[0 45 90]；与条纹走向垂直。
% 输出result：频响、径向RMS、主响应区间、理论/实测增益及条纹图。
% 图1：二维幅度热图与频响曲线；图2：第一个方向的输入/细节条纹。
% 归一化径向频率r=2/T；r<=1用于全方向比较，四角高频看二维热图。
% 基础MATLAB即可运行；绘图使用tiledlayout（R2019b及以上）。

if nargin < 2 || isempty(periods)
    periods = [2 3 4 6 8 12 16 24 32];
end
if nargin < 3 || isempty(anglesDeg)
    anglesDeg = [0 45 90];
end
validateattributes(h, {'numeric'}, {'2d','real','finite','nonempty'});
assert(all(mod(size(h),2)==1) && any(h(:)~=0), ...
    'h必须是非零的奇数尺寸核。');
validateattributes(periods, {'numeric'}, {'vector','real','finite','>=',2});
validateattributes(anglesDeg, {'numeric'}, {'vector','real','finite','nonempty'});
h = double(h);
periods = double(periods(:).');
anglesDeg = double(anglesDeg(:).');

% 热图保留实际增益；FFT中的核位置只影响相位，不影响幅度。
nFFT = max(512, 2^nextpow2(max(size(h))));
f = (-nFFT/2:nFFT/2-1)*2/nFFT;
mag2D = abs(fftshift(fft2(h,nFFT,nFFT)));
r = linspace(0,1,501).';
theta = 0:2:358;
Hradial = kernelResponse(h, r*cosd(theta), r*sind(theta));
radialRMS = sqrt(mean(abs(Hradial).^2,2));
Hdirection = kernelResponse(h, r*cosd(anglesDeg), r*sind(anglesDeg));

% -3dB只取含最大峰值的连续区间；触及端点时不把端点当作截止。
[peakGain, peakIndex] = max(radialRMS);
threshold = peakGain/sqrt(2);
lo = peakIndex;
hi = peakIndex;
while lo > 1 && radialRMS(lo-1) >= threshold
    lo = lo-1;
end
while hi < numel(r) && radialRMS(hi+1) >= threshold
    hi = hi+1;
end
mainBand = [r(lo),r(hi)];

% 每组条纹均为0.5+0.25*cos，单独滤波，避免分区接缝产生伪响应。
padR = (size(h,1)-1)/2;
padC = (size(h,2)-1)/2;
testSize = 256+2*max(padR,padC);
[x,y] = meshgrid(0:testSize-1);
rows = padR+1:testSize-padR;
cols = padC+1:testSize-padC;
padRows = min(max((1-padR):(testSize+padR),1),testSize);
padCols = min(max((1-padC):(testSize+padC),1),testSize);
nAngle = numel(anglesDeg);
nPeriod = numel(periods);
inputStripes = cell(nAngle,nPeriod);
detailStripes = cell(nAngle,nPeriod);
measuredGain = zeros(nPeriod,nAngle);
theoreticalGain = zeros(nPeriod,nAngle);
amplitude = 0.25;
for a = 1:nAngle
    for p = 1:nPeriod
        phase = 2*pi*(x*cosd(anglesDeg(a))+y*sind(anglesDeg(a)))/periods(p);
        input = 0.5+amplitude*cos(phase);
        padded = input(padRows,padCols);
        % conv2翻转核一次后执行相关，与当前imfilter默认方向一致。
        detail = conv2(padded,rot90(h,2),'valid');
        phaseCore = phase(rows,cols);
        detailCore = detail(rows,cols);
        basis = [cos(phaseCore(:)),sin(phaseCore(:)),ones(numel(phaseCore),1)];
        % 拟合振幅，避开边界；pinv也适用于T=2时退化的正弦列。
        fit = pinv(basis)*detailCore(:);
        measuredGain(p,a) = hypot(fit(1),fit(2))/amplitude;
        freq = 2/periods(p);
        theoreticalGain(p,a) = abs(kernelResponse(h, ...
            freq*cosd(anglesDeg(a)),freq*sind(anglesDeg(a))));
        inputStripes{a,p} = input(rows,cols);
        detailStripes{a,p} = detailCore;
    end
end

figResponse = figure('Name','单核频率响应','Color','w','Position',[60 80 1300 580]);
layout = tiledlayout(1,2,'TileSpacing','compact','Padding','compact');
nexttile;
imagesc(f,f,mag2D);
axis image xy;
colormap(gca,parula(256));
colorbar;
xlabel('归一化水平频率 f_x');
ylabel('归一化垂直频率 f_y');
title('最终核的二维幅度响应 |H|');

nexttile;
hold on;
topGain = max([mag2D(:);radialRMS(:);abs(Hdirection(:));measuredGain(:)])*1.1;
patch([mainBand fliplr(mainBand)],[0 0 topGain topGain], ...
    [0.8 0.9 0.85],'FaceAlpha',0.3,'EdgeColor','none','HandleVisibility','off');
plot(r,radialRMS,'k','LineWidth',2,'DisplayName','径向RMS');
colors = lines(nAngle);

lineStyles = {'-', '--', ':'};
markers = {'o', 's', '^'};

for a = 1:nAngle
    plot(r, abs(Hdirection(:,a)), ...
        'Color', colors(a,:), ...
        'LineStyle', lineStyles{a}, ...
        'LineWidth', 1.8, ...
        'DisplayName', sprintf('%g° 理论', anglesDeg(a)));

    plot(2./periods, measuredGain(:,a), ...
        markers{a}, ...
        'Color', colors(a,:), ...
        'MarkerSize', 6, ...
        'LineStyle', 'none', ...
        'DisplayName', sprintf('%g° 条纹实测', anglesDeg(a)));
end

plot([0 1],[threshold threshold],':','Color',[0.5 0.5 0.5],'HandleVisibility','off');
xlim([0 1]);
ylim([0 topGain]);
xticks([0 0.125 0.25 0.5 0.75 1]);
xticklabels({'0 / Inf','0.125 / 16','0.25 / 8','0.5 / 4','0.75 / 2.67','1 / 2'});
xlabel('归一化径向频率 r / 条纹周期 T (px)');
ylabel('幅度增益（未归一化到峰值）');
title(sprintf('径向RMS峰值 r=%.3f；主响应区间 %.3f–%.3f', ...
    r(peakIndex),mainBand));
grid on;
legend('Location','best','FontSize',8);
title(layout,'实际频率响应与已知周期条纹验证；阴影为主峰的-3dB区间');

% 所有输出共用对称色标，零细节为中灰；只展示第一个方向。
limits = cellfun(@(d) max(abs(d(:))),detailStripes);
detailLimit = max(limits(:));
detailLimit = max(detailLimit,eps);
figStripes = figure('Name','像素尺度与细节输出','Color','w', ...
    'Position',[60 120 max(1100,140*nPeriod) 430]);
layout = tiledlayout(2,nPeriod,'TileSpacing','compact','Padding','compact');
% 显示内部64×64小块，避免整幅条纹缩小后掩盖最细周期。
showRows = 1:min(64,size(inputStripes{1,1},1));
showCols = 1:min(64,size(inputStripes{1,1},2));
for p = 1:nPeriod
    nexttile(p);
    imagesc(inputStripes{1,p}(showRows,showCols),[0 1]);
    axis image off;
    title(sprintf('T=%g px',periods(p)));
    nexttile(nPeriod+p);
    imagesc(detailStripes{1,p}(showRows,showCols),[-detailLimit detailLimit]);
    axis image off;
    title(sprintf('增益 %.3f',measuredGain(p,1)));
end
colormap(figStripes,gray(256));
title(layout,sprintf(['方向 %g°；上：同振幅输入；下：有符号细节，统一色标±%.3g\n'],...
                        anglesDeg(1),detailLimit));

result.kernel = h;
result.kernelSum = sum(h(:));
result.frequency2D = f;
result.magnitude2D = mag2D;
result.radialFrequency = r;
result.radialRMS = radialRMS;
result.directionResponse = Hdirection;
result.peakFrequency = r(peakIndex);
result.mainBand = mainBand; % 采样网格上的主峰-3dB区间，不是硬截止。
result.bandTouchesLimit = [lo==1,hi==numel(r)];
result.periodRange = 2./fliplr(mainBand);
result.periods = periods;
result.anglesDeg = anglesDeg;
result.theoreticalGain = theoreticalGain; % 行：periods；列：anglesDeg。
result.measuredGain = measuredGain;
result.inputStripes = inputStripes; % 单元格行：方向；列：周期；已裁去边界。
result.detailStripes = detailStripes;
result.detailDisplayLimit = detailLimit;
result.figures = [figResponse figStripes];
fprintf('核总和：%.6g；径向RMS峰值：r=%.3f，增益=%.4f\n', ...
    result.kernelSum,result.peakFrequency,peakGain);
fprintf('主响应区间：[%.3f, %.3f]；对应周期：[%.2f, %.2f] px\n', ...
    mainBand,result.periodRange);
end

function H = kernelResponse(h,fx,fy)
% 以核中心为原点计算相关响应，保留有符号/复数响应。
rowOffset = (1:size(h,1))-(size(h,1)+1)/2;
colOffset = (1:size(h,2))-(size(h,2)+1)/2;
H = zeros(size(fx));
for i = 1:size(h,1)
    for j = 1:size(h,2)
        H = H+h(i,j)*exp(1i*pi*(fx*colOffset(j)+fy*rowOffset(i)));
    end
end
end
