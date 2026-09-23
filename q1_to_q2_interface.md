# Q1 → Q2 标准接口

## 1. 接口目标

Q2 的观测目标不是某个窗口的单一均值，而是 Q1 最终估计器输出的完整、带不确定性的 ERP waveform。为避免把当前 Q1 草稿路径写死，Q2 只读取一个逻辑接口：

```text
Q1_FINAL_EXPORT/
├── erp_waveform.csv
├── metadata.json
└── provenance.json
```

目录名可以在运行配置中改变，但字段和语义必须保持不变。只有满足本接口 schema、状态和溯源校验的导出目录才能作为 `Q1_FINAL_EXPORT`。

## 2. 标准 waveform 表

`erp_waveform.csv` 每行表示一个 `record × task × condition × channel × time` 单元，至少包含：

| 字段 | 类型 | 约束与含义 |
|---|---|---|
| `record` | string | 记录标识；不得把它解释为受试者，除非有独立受试者映射 |
| `task` | string | 例如 `Task1`、`Task2`；使用稳定词表 |
| `condition` | string | 例如 `Left`、`Right`；保留原始条件语义 |
| `channel` | string | 只能为 `F3`、`Fz`、`F4`，顺序在 metadata 中声明 |
| `time_s` | float | 相对刺激 onset 的秒数；单调递增且在同一组内等间隔或显式记录不等间隔 |
| `erp_value` | float | Q1 最终 estimator 的点估计，单位必须在 metadata 中给出（通常为 µV） |
| `ci_low` | float | bootstrap 或其他预先声明方法得到的下界 |
| `ci_high` | float | 与 `ci_low` 同定义的上界，且 `ci_low ≤ erp_value ≤ ci_high` 应在允许数值容差内成立 |
| `n_trials` | integer | 该单元实际纳入的 trial 数，必须为非负整数；零 trial 不应伪装成零波形 |

允许追加字段（如 `se`、`estimator`、`filter_id`），但不得改变上述字段含义。

## 3. `metadata.json` 最小内容

```json
{
  "schema_version": "Q1_FINAL_EXPORT.v1",
  "status": "draft_or_candidate_or_final",
  "time_reference": "stimulus_onset",
  "time_unit": "s",
  "erp_unit": "uV",
  "channels": ["F3", "Fz", "F4"],
  "channel_reference": "describe_average_or_reference",
  "preprocessing": {
    "filter_type": "declare",
    "highpass_hz": null,
    "lowpass_hz": null,
    "notch_hz": null,
    "baseline_window_s": null,
    "artifact_rule": "declare"
  },
  "estimator": {
    "name": "declare_final_method",
    "version": "declare",
    "bootstrap": "declare_resampling_definition"
  },
  "trial_set_definition": "declare_inclusion_exclusion_rule",
  "source_data_hash": "sha256_or_equivalent",
  "created_at": "ISO-8601"
}
```

`channel_reference` 必须写清平均参考、双乳突参考或其他参考方式；不能只写“EEG”。`preprocessing` 需要记录滤波器类型、截止频率、阶数/相位方式（如适用）、基线区间和伪迹规则。这样 Q2 才能区分真实波形差异与预处理差异。

## 4. `provenance.json` 与版本边界

建议至少记录：

- 原始 MAT 文件的哈希、文件名和读取脚本版本；
- Q1 生成命令、MATLAB/Python 版本；
- estimator 代码 commit；
- 生成时间和数据表行数；
- trial-set 哈希或纳入 trial ID 的摘要哈希；
- 是否为训练、验证或最终汇总数据。

如果 Q1 方法比较尚未完成，导出文件必须标记 `status: draft`，Q2 不得将其作为最终观测目标。

## 5. 验证规则

在 Q2 读取前执行以下只读检查：

1. 必需列全部存在，列名大小写一致；
2. 每个 `record/task/condition/channel` 组的 `time_s` 唯一、递增，时间网格与其他组兼容；
3. `erp_value`、CI 和 `n_trials` 无非法缺失；CI 顺序和单位有效；
4. 每组 `n_trials` 与 Q1 trial 统计一致，不能用窗口内有效点数替代 trial 数；
5. `Left/Right`、`Task1/Task2` 的词表不发生隐式重命名；
6. 记录 `channel_reference`、预处理、估计器和 bootstrap 定义；
7. 重新计算导出文件哈希并与 `provenance.json` 对比；
8. 发现任何 schema/单位/哈希错误时停止模型拟合并报告错误。

## 6. Q2 使用方式

Q2 先把表重塑为

\[
Y_{r,\tau,q}\in\mathbb{R}^{T\times3},
\]

其中每个时间点含 `[F3,Fz,F4]`。模型预测 `\hat Y(θ)` 与 `Y` 在同一时间网格上比较，完整波形优先：

\[
J_{wave}=\sum_{r,\tau,q,t}
\|\hat Y_{r,\tau,q}(t)-Y_{r,\tau,q}(t)\|_2^2,
\qquad q\equiv\texttt{condition}\in\{Left,Right\}.
\]

若使用 Q1 CI，应在训练阶段依据 CI 或估计协方差进行标准化，并在报告中说明；不能只挑 250–500 ms 窗口或只挑一个通道拟合。窗口幅度、峰潜伏期、左右差异波形和空间模式作为预注册诊断项。

## 7. 训练/验证边界

- 预处理参数、标准化、模型结构选择和拟合只使用训练记录/训练条件。
- 验证集用于比较 Model A/B/C 和复杂度，不用于再次调参。
- 最终汇总只在模型和参数化冻结后生成。
- 不把四个 record 自动称为四个受试者，也不把 F3/Fz/F4 自动称为三个脑区。

## 8. Q1 输出与 Q2 状态约定

| 状态 | 允许用途 |
|---|---|
| `draft` | schema 检查、绘图和接口联调 |
| `candidate` | 候选 estimator 比较；可做预注册的模型模拟，不作最终 Q2 结论 |
| `final` | Q1 方法、trial-set、预处理和不确定性已冻结，可作为 Q2 waveform target |

Q2 当前阶段只完成接口设计；在 `status=final` 且通过验证规则前，不开始最终参数拟合。
