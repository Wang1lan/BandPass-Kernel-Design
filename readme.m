% BANDPASS-KERNEL-DESIGN 带通核分析与双带通图像增强
%
% 仓库：https://github.com/Wang1lan/BandPass-Kernel-Design
% 本文按 main 分支源码 020906f8e1e0884b5ad3388ae5a47db506d83aeb 整理。
% 参数说明以该版本为基准，后续调参请以脚本顶部的实际设置为准。
% 本文件全部为 MATLAB 注释，可在编辑器阅读，也可在仓库根目录执行 help readme。
%
% 1. 项目内容
% 本项目包含两类 MATLAB 浮点实验：
%   1) 分析自定义带通核的二维频响、径向频带与不同周期/方向的条纹响应。
%   2) 在 RGB 或 YCbCr-Y 上进行双带通细节增强，并支持增强前去噪与高光保护。
% 去噪支持双边滤波和自引导滤波；BOTH 模式同时显示 RGB 与 Y 增强结果。
% 核参数和增益是可调实验配置，不代表已经标定的最终频带或图像质量标准。
%
% 2. 文件与入口
%   demo_doubleK.m                   双带通增强、去噪、显示及实图评价脚本。
%   compareThreeBandpassKernels.m     当前实际比较两个细节核及其合成核。
%   inspectBandpassKernel.m          单核频响及正弦条纹分析函数。
%   Functions/bandpassFilter.m       仓库自定义带通核生成函数。
%   Functions/filterShiftSum.m       整幅图错位累加卷积，支持二维及 RGB double。
%   Functions/bilateralFilterShiftSum.m  二维/联合 RGB 双边去噪。
%   Functions/rgbToYCbCrFull.m       double 全范围 RGB 转 YCbCr。
%   Functions/yCbCrFullToRgb.m       与正变换配套的 YCbCr 转 RGB。
%   Functions/buildHighlightMask.m   从原图构造 hardMask、softMask。
%   Functions/evaluateTEGGlobal.m    非高光区域的梯度均值比。
%   Functions/evaluateNoisePCA.m     共享 ROI 上的 PCA 噪声估计。
%   testImg/olymRaw.jpg              随仓库提供的默认测试图片。
%
% bandpassFilter 是本仓库自定义函数，不是 MATLAB 内置的 bandpass API。
% compareThreeBandpassKernels 的文件名虽保留 Three，当前配置只包含两个细节核。
%
% 3. 环境与快速运行
% 建议使用 MATLAB R2022a，并安装 Image Processing Toolbox。
% 代码使用 fsamp2、padarray、bwdist、imfilter；引导去噪使用 imguidedfilter。
% 交互选择 PCA ROI 还需要 drawrectangle。绘图使用 tiledlayout。
%
% 仓库提供默认测试图片 testImg/olymRaw.jpg，可直接用于 demo 的 imgPath。
% 使用其他图像时，请修改 demo_doubleK.m 顶部的 imgPath。
% 输入须为 H×W×3 的 RGB uint8；当前仓库未附带结果图或自动化验证脚本。
%
% 将 MATLAB 当前文件夹切换到仓库根目录后，执行：
%   demo_doubleK
%
% demo 会根据自身路径定位资源并加入 Functions，无须依赖启动时的当前文件夹。
% demo 和 compare 脚本均执行 clear、clc、close all，会清理工作区及已有图窗。
% 参数在脚本顶部设置；仅在命令行预先给同名变量赋值，不能替代脚本内配置。
%
% 4. demo 当前配置与调参位置
% 当前主流程默认值：
%   enhanceMode            = 'RGB'       增强方式：'RGB'、'Y'、'BOTH'。
%   denoise.enable         = true        去噪开启。
%   denoise.method         = 'guided'    方法：'bilateral' 或 'guided'。
%   denoise.mode           = 'Y'         去噪通道：'RGB' 或 'Y'。
%   denoise.windowSize     = 5           实际窗口边长，5 表示 5×5，不是半径。
%   denoise.sigmaSpatial   = 1.5         仅双边使用，单位为像素。
%   denoise.sigmaRange     = 10          仅双边使用，单位为 8 位灰度级。
%   denoise.guidedEps      = 1           仅引导使用，单位为灰度级平方。
%   detailSign             = 1           两路细节共用的符号。
%   enableNoise            = false       仅控制 PCA 评价，与去噪开关独立。
%   noiseROI               = []          开启 PCA 时，在原图交互框选。
%   detailDisplayLimit     = 30          细节热图共用范围 [-30,30]。
%
% 两路核的参数列为 [fRadius, sigmaHighFreq, sigmaLowFreq, scaleLowFreq]：
%   第一路：[5, 0.6, 0.2, 0.8]，核尺寸 11×11，detailGain_1 = 8。
%   第二路：[3, 1.0, 0.5, 0.95]，核尺寸 7×7，detailGain_2 = 4。
% fRadius 决定核边长 2*fRadius+1；两个 sigma 是频域高斯参数，不能直接当作
% 硬通带边界或条纹周期。scaleLowFreq 是低频项系数，实际频带应通过分析工具查看。
%
% 常用设置示例：修改脚本顶部对应字段后重新运行。
%   关闭去噪：denoise.enable = false;
%   双边 RGB 去噪：denoise.enable = true; denoise.method = 'bilateral';
%                  denoise.mode = 'RGB';
%   引导 Y 去噪并对比增强：denoise.method = 'guided'; denoise.mode = 'Y';
%                          denoise.enable = true; enhanceMode = 'BOTH';
%
% 增强方式、去噪方法、去噪通道三者独立。去噪关闭时，方法及专用参数均旁路。
% 去噪开启时窗口须为有限正奇整数；窗口为 1 合法。双边仅检查两个 sigma 为正，
% 引导仅检查 guidedEps 为有限正标量，并要求窗口不超过图像的高和宽。
% guidedEps 不与 sigmaRange 自动绑定；不同方法使用相同窗口也不代表平滑强度相同。
%
% 5. 图像处理流程与数值约定
%   原始 RGB uint8 -> 原图高光掩模 -> 可选一次去噪得到 Iwork
%                 -> RGB 和/或 Y 双带通增强 -> 各分支限幅、舍入、输出 uint8。
%
% img 始终保留原图；Iwork 为共享的 RGB double 输入，底图与两路细节均从它计算。
% 两路核并行处理同一个输入，不串联。BOTH 的两个增强分支也互不串联、不融合。
% RGB 分支的重建为：
%   Iwork + detailGain_1.*(1-softMask).*detailRGB_1
%         + detailGain_2.*(1-softMask).*detailRGB_2
% 其中 detailRGB_i = detailSign .* filterShiftSum(Iwork,kernel_i,'replicate')。
% Y 分支采用相同的两路公式，但只更新 Y，保留增强分支输入的 Cb、Cr。
%
% YCbCr 使用本仓库配套的全范围转换：Y=0.299R+0.587G+0.114B，Cb/Cr 中性值128。
% RGB 和 Y 均保持 0～255 灰度级单位的 double；转换本身不裁剪、不舍入。
% 色度保持应在逆转换前检查，最终 RGB 限幅后重新计算的色度可能发生变化。
% 最终输出统一按 uint8(round(min(max(raw,0),255))) 生成。
%
% 双边 RGB 去噪采用通道均方颜色距离，三个通道共享一张权重图，边界为 replicate。
% 引导去噪使用 imguidedfilter(src,src,...)，明确传入窗口及 DegreeOfSmoothing。
% RGB 自引导逐通道分别处理，Y 自引导仅处理亮度；引导沿用函数默认边界。
% 当前输入不归一化为 0～1，故 guidedEps 也使用灰度级平方单位。
%
% 高光由原始 img 的最大/最小通道值检测，经膨胀和平滑形成 softMask。
% softMask=1 表示禁止叠加增强增量，不会撤销之前的去噪效果。
% 原始检测 hardMask 还用于排除 TEG 评价中的高光邻域。
%
% 6. 显示、输出变量与评价
% 整图窗口并排显示原图及实际运行的增强分支；BOTH 为三幅图。
% 细节热图显示掩模及增益之前的两路带符号亮度细节；BOTH 为 2×2 布局。
% 所有热图共用蓝-白-红色图与 detailDisplayLimit，标题标明分支、核及去噪方式。
%
% RGB 分支变量：imgRawRGB（裁剪前 double）、imgEnhRGB（RGB uint8）、tegRGB。
% Y 分支变量：imgRawY（转回 RGB 后的 double）、imgEnhY（RGB uint8）、tegY。
% 开启 PCA 后另有 noiseRGB/noiseY；未运行的分支不生成对应结果。
%
% TEGGlobal = 有效区域增强图平均 Sobel 梯度 / 原图平均 Sobel 梯度。
% 有效区域由原图 hardMask、highlightDistance 和 borderWidth 共同确定。
% 无有效区域或原图平均梯度接近零时返回 NaN。比值增大表示梯度增强，
% 不能单凭该数值认定细节质量提高，因为噪声也可能贡献梯度。
%
% enableNoise=true 时，noiseROI=[] 会在原图框选并双击确认；也可填写一基整数
% [x y width height]。BOTH 只选择一次 ROI，两份结果复用相同 patch 中心。
% 当前 patchSize=5、stride=3、tailNum=5；ROI 须容纳至少 patchSize^2+1 个完整 patch。
% PCA 返回 sigmaOrig、sigmaEnh、noiseGain、noiseGainDB、roiPos、validPatchNum：
%   noiseGain = sigmaEnh/sigmaOrig；noiseGainDB = 20*log10(noiseGain)。
% 原图 sigma 接近零时 gain 为 NaN；原图 sigma 为正且输出 sigma 为零时，dB 可为 -Inf。
% ROI 应选择平坦、无饱和高光的区域；估计值依赖 ROI 和 patch 配置。
% 两个评价函数内部沿用亮度近似 (77R+150G+29B)/256，始终以原始 img 为基准。
%
% 7. 带通核与合成响应分析
% 在仓库根目录执行以下命令，分析不需要实际图像：
%   addpath(fullfile(pwd,'Functions'));
%   compareThreeBandpassKernels
%
% compare 脚本本身不加入 Functions 路径，运行前不能省略上述 addpath。
% 其顶部旧注释中的 code_bandPassFilter 不是当前仓库目录；核函数位于 Functions。
% compare 当前独立设置为：params=[5 .55 .2 .85; 3 1 .5 .95]，detailGains=[6 3]。
% 它不会自动读取 demo 的参数；如需分析相同实验，应手动同步核参数、符号及增益。
% 合成核包含原图直通项，响应为 1 + detailSign*(gain1*H1+gain2*H2)。
% 该线性合成响应不包含实图高光掩模和非线性去噪的影响。
% 归一化曲线/条纹仅用于比较频段形状；实际幅值增益仍会单独输出。
%
% 单独分析一个核：
%   addpath(fullfile(pwd,'Functions'));
%   h = bandpassFilter(5,0.6,0.2,0.8);
%   r = inspectBandpassKernel(h,[2 3 4 6 8 12 16 24 32],[0 45 90]);
%
% h 应为有限实数、非零、奇数尺寸二维核。当前绘图请最多使用三个方向。
% r 包含 magnitude2D、radialRMS、peakFrequency、mainBand、periodRange、
% theoreticalGain、measuredGain、inputStripes、detailStripes 和两个 figures 句柄。
% 归一化径向频率 rFreq 与条纹周期 T 的关系为 rFreq=2/T，T 单位为像素。
% mainBand 是采样网格上主峰的连续 -3 dB 区间，不是硬截止频率。
% 理论/实测幅值增益不保留响应的符号，不能单凭幅值图判断 detailSign 极性。
% compare 中的 whiteNoiseGain=norm(h(:),2) 是线性核的理论白噪声增益，
% 与实图 ROI 上估计的 PCA noiseGain 含义不同。
%
% 8. 保存结果
% 默认只显示及评价，不自动写出增强图像。demo 的保存示例当前位于脚本前部，
% 并被块注释关闭；应在增强完成后执行该段，或移至处理流程末尾再启用。
% 不能只在原位置取消注释，否则 runRGB/runY 及增强结果尚未生成。
% 使用原保存示例前还需确认 savePath 存在；示例本身不会自动创建目录。
% 原示例文件名区分增强分支、去噪方法/通道及核参数、增益。
%
% 也可先运行完 demo，再在命令行按当前分支保存（以下为独立示例）：
%   saveDir = fullfile(projectDir,'result');
%   if ~isfolder(saveDir), mkdir(saveDir); end
%   tag = datestr(now,'yyyymmdd_HHMMSS');
%   if runRGB, imwrite(imgEnhRGB,fullfile(saveDir,['RGB_' tag '.png'])); end
%   if runY, imwrite(imgEnhY,fullfile(saveDir,['Y_' tag '.png'])); end
%
% PNG 为无损格式，可用于与工作区 uint8 数组逐像素核对。
%
% 9. 使用说明的证据范围
% 本文根据仓库源码说明接口与行为，不把本地历史实验结果当作本提交的验收结论。
% 随附样例已于 2026-10-23 在 MATLAB R2022a 按上述默认配置运行，
% 输出为 1080×1350×3 RGB uint8，并完成尺寸、类型和浮点有限性检查。
% 仓库提供一张默认测试图；其他图像、参数组合和 ROI 的验证需按实际用途进行。
% 调整图像、核、增益、去噪参数或 ROI 后，应重新检查输出、细纹理及评价值。
%
