# 视觉认知脑电响应分析

本程序只使用四个 MAT 文件的原始 Fz、F3、F4 通道（第 1–3 行）。第 8 行解析视觉提示，第 9 行解析目标和点击，第 10 行检查时间戳。设备自带的 Decon 通道不参与计算。

## 运行

在 MATLAB R2024a 中，先切换到本仓库目录，再执行：

```matlab
addpath(pwd);
sourceDir = fullfile(pwd, "data");
run_q1(sourceDir, fullfile(pwd, "results"));
```

依赖 MATLAB Signal Processing Toolbox（`butter`、`filtfilt`）和 Statistics and Machine Learning Toolbox（`prctile`）。

## 处理方法

1. 按非零提示的起点定位每次试验，并在下一次提示前匹配目标和点击；不把提示方向当成目标方向。
2. 现行主分析将连续信号经 0.5–30 Hz 零相位滤波，之后按提示起点截取 −0.2 至 0.8 秒，使用提示前数据校正基线。新增的高通专项在同一批 350 个保留试次上比较 0.10、0.20、0.25、0.50 Hz，其他预处理及 Huber 拟合参数固定；现行主分析设置尚未更改。
3. 原始试次中任一通道达到绝对值 999 的标记为疑似采集截断。其余试次根据峰峰值与相邻样本最大跳变的记录内中位数绝对偏差标记强伪影。
4. 对通过质量检查的左右试次分别计算普通平均和 Huber 稳健平均，并对稳健平均施加二阶差分平滑，得到可画图的拟合曲线。
5. 输出 100–250 ms 和 250–500 ms 时间窗的均值、bootstrap 区间与左右差值；计算左右差波的分半重复性。分类验证包括同记录按试次时间顺序划分的五折测试（以提示序列循环移位构造零假设）和同任务 A/B 两份记录互测。
6. 同时比较原三电极特征与仅保留两种空间对比的特征：`Fz-(F3+F4)/2` 和 `F3-F4`。这是分析用空间重参考，并不等于证明被减掉的共同分量全是伪影。

滤波频段、伪影阈值、平滑强度均集中在 `run_q1.m` 的 `default_config` 中，可在正式论文的敏感性分析中调整。解码仅是检测视觉信息是否保留的辅助检验，不代替脑电波形与统计分析。

## 输出

- `events_and_quality.csv`：逐试次事件、提示/目标/点击方向、时间和质量标记。
- `quality_summary.csv`：各文件左右试次数、截断和异常试次数。
- `erp_curves.csv`：三电极左右提示的普通、稳健和拟合响应曲线。
- `spatial_mode_curves.csv`：拟合曲线的共同分量、中线对比和左右对比。
- `window_statistics.csv`：两段预定时间窗的响应和左右差值及 bootstrap 区间。
- `split_half_reliability.csv`：同一记录内两半试次的波形相关性。
- `contrast_split_half.csv`：左右差波的分半相关性，比单纯共模波形相关更贴近“保留形状特征”的要求。
- `filter_sensitivity.csv`：既有 0.1/0.5/1.0 Hz 诊断；使用 −0.8 至 0.8 秒扩展试次，并额外排除扩展区间内的饱和试次，因此不用于本轮固定试次的频率对照。
- `cross_record_decoding.csv`：同任务 A/B 互测的平衡准确率、AUC 和置换检验结果。
- `within_record_decoding.csv`、`within_record_folds.csv`：按时间留出的同记录五折验证；提示序列循环移位保持其连续性，用于构造零假设。
- `late_wave_summary.csv`、`late_wave_common_mode.csv`：0.55–0.80 秒晚期波动的幅值、峰时及电极波形相关性；电极相关性是平均波形沿时间轴的相关，不是独立试次的跨脑区连接强度。
- `figures/*.png`：各记录的三电极左右响应波形。
- `results_highpass/q1_trial_balance.csv`：四组记录中左右条件的原始与保留试次数。
- `results_highpass/q1_highpass_sensitivity_real.csv`、`q1_highpass_waveforms.csv`：固定 350 试次的 12 个“记录×Task×电极”组合及左右、差波的四频率对照。
- `results_highpass/q1_filter_simulation_recovery.csv`：两类半模拟参考模板、12 次重复的滤波恢复指标；同目录保留拆分、校准参数和示例波形。
- `results_highpass/figures/*.png`：滤波响应、真实波形、完整矩阵和半模拟恢复图。方法解释与工作设置建议见 `q1_filter_implementation_audit.md`、`q1_highpass_decision_report.md`。

高通专项可在 MATLAB R2024a 中从仓库目录运行：

```matlab
addpath(pwd);
q1_highpass_real(fullfile(pwd,"data"), fullfile(pwd,"results_highpass"), fullfile(pwd,"results","events_and_quality.csv"));
q1_filter_simulation(fullfile(pwd,"data"), fullfile(pwd,"results_highpass"), fullfile(pwd,"results","events_and_quality.csv"));
q1_future_event_leakage(fullfile(pwd,"data"), fullfile(pwd,"results_highpass"), fullfile(pwd,"results","events_and_quality.csv"));
```

随后执行 `python q1_highpass_figures.py results_highpass` 生成 PNG 图。该脚本需要 NumPy、SciPy、pandas 和 Matplotlib。

## 当前数据诊断

四份记录各有 100 次提示，主分析分别保留 80、93、87、90 次，共 350/400 次。既有扩展窗诊断中，A 的 Task 2、F3 的右减左 250–500 ms 均值从 0.1 Hz 下约 45.6 变为 0.5 Hz 下约 8.6 个源数据单位；该诊断因额外筛除 17 次而实际使用 333 次，不能与 350 次主分析直接比较。严格固定 350 次后，四种高通下该组合依次为 61.33、44.50、37.72、18.54；记录 B 的对应值为 19.55、20.51、22.77、3.29，未复现 A 的幅度变化。0.20–0.25 Hz 是当前值得考虑的折中范围，0.20 Hz 只是待确认的工作设置；现行 0.5 Hz 主分析仍保留。约 0.55 秒开始的同步大幅波动不能仅凭本数据断言它是 P300、视觉形状信息或眼动伪影。

A/B 互测的左右提示平衡准确率约 0.41–0.52；加入空间对比特征后约 0.42–0.47，均未显著高于随机水平。同记录按时间留出的验证也没有通过循环移位检验，八项经 Holm 校正的 p 值均为 1。项目二 0.55–0.80 秒的大幅波动在四个“记录×左右提示”组合中都存在，三电极平均曲线峰值约在 0.66–0.74 秒；它更像是共同任务相关信号，而非已验证的左右形状特征，但现有三电极数据不能进一步确定其来源。当前结果可以用于描述和拟合观测响应，不足以证明可跨记录稳定区分左右三角形。Fz/F3/F4 都在额区，也不宜直接把 250–500 ms 正波称为典型中央—顶区 P300。原数据未给物理标定单位，图和表中的幅值保留为“源数据单位”，不能写成 μV。

简明结果和论文表述见 `results_summary.md`。

程序不会写入或修改 `data/` 中的原始 MAT 文件。`cross_record_decoding.csv` 是跨记录验证；只有确认 A/B 确为不同受试者后，才能称作跨受试者验证。

`results/` 是当前 0.5 Hz 主分析输出；`results_checked/` 保留 0.1 Hz 参数下的探索性运行；`results_highpass/` 是固定试次的高通专项输出。项目描述与数据校验信息见 `DATA_PROVENANCE.md`；题目原文见 `problem/`。
