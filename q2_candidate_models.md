# Q2 候选模型与四层结构

## 1. 统一记号和边界

条件索引为 `r`（record）、`τ`（task）和 `q∈{L,R}`（condition：Left/Right）。三通道输出固定为

\[
y_{r,\tau,q}(t)= [F3(t),Fz(t),F4(t)]^\top .
\]

本文定义候选结构、参数角色和后续拟合协议；最终观测目标由通过校验的 `Q1_FINAL_EXPORT` 提供。

整体计算链为

\[
u_q(t)\rightarrow x_L(t)\rightarrow x_C(t)\rightarrow J(t)\rightarrow y(t).
\]

`Left/Right` 首先解释为形状/方向类别。除非实验设计额外证明空间视野左右位置，否则不把标签写成左、右视野。

## 2. 四层最小模型

### Layer 1：visual stimulus encoding

以刺激时间包络 `p(t)` 和低维特征向量 `φ_q` 表示输入：

\[
u_q(t)=p(t)\,\phi_q,
\qquad \phi_L\neq\phi_R.
\]

第一版可采用二维编码并固定范数，例如 `φ_L=[1,0]^T`、`φ_R=[0,1]^T`，或由图形特征构成的二维/三维向量。该选择是 D 级计算抽象，不声称大脑只有两个神经元。为了允许左右差异很弱，可令 `φ_R=φ_L+Δφ`，并在模型比较中检验 `Δφ=0` 的受限模型。

Stringer 的结果提醒我们：真实视觉群体编码可以是高维和分布式的，因此低维 `φ_q` 只承担“本题两个刺激类别的有效驱动方向”，不能被解释为完整 V1 表征。

### Layer 2：LGN dynamics

使用一阶、带输入延迟的低维滤波器：

\[
\tau_L\dot{x}_L(t)=-x_L(t)+W_Lu_q(t-\delta_L)+b_L+\eta_L(t),
\]

其中 `x_L∈R^{d_L}`，`W_L` 为特征到 LGN 的增益，`τ_L>0` 为有效时间常数，`δ_L≥0` 为输入延迟。主模型优先 `d_L=1`；只有在留出波形确实增加解释力时才考虑 `d_L=2`。

Lien–Scanziani 支持“丘脑输入可以具有相位/空间结构并被皮层放大”，因此保留 `W_L` 和延迟通道是机制上合理的；他们没有提供可直接用于本题人类 EEG 的 `τ_L` 或 `δ_L`。

### Layer 3：cortical population dynamics

#### Model A：经验线性基线

\[
\dot{x}(t)=A x(t)+B u_q(t)+b+\eta(t),
\qquad
J(t)=H x(t),
\qquad
y(t)=LJ(t)+d+\varepsilon(t).
\]

`x∈R^{d}` 是未标记的经验状态，`A` 需稳定或受到稳定性正则约束。Model A 不把 `x` 命名为 LGN 或 E/I，只回答“一个低维线性动态系统是否已足够解释 ERP 波形”。

- 参数：`A,B,b,H,L,d` 与噪声协方差。
- 物理含义：最小复杂度的时间滤波和三通道混合。
- 文献支持：Horrocks 的 FA 说明低维潜在轨迹可作描述性工具（B/C）；Breakspear 支持状态空间与观测模型的一般框架（C）。
- 优点：参数少、可线性化、最适合作为泛化基线。
- 风险：不能区分 LGN 驱动与皮层内部动力学，`x` 的生理解释弱。
- Left/Right：通过 `B(φ_L-φ_R)` 产生波形差；若差异不稳定，模型可以估计接近零的驱动差。
- 与 Q3：仅提供可复用的潜在状态，不含记忆/决策环节。

#### Model B：LGN + linear cortex + low-dimensional source

\[
\begin{aligned}
\tau_L\dot{x}_L &= -x_L+W_Lu_q(t-\delta_L)+b_L+\eta_L,\\
\tau_C\dot{x}_C &= A_Cx_C+B_Cx_L+b_C+\eta_C,\\
J(t)&=H_Cx_C(t),\\
y(t)&=LJ(t)+d+\varepsilon(t).
\end{aligned}
\]

主设定为 `d_L=1`、`d_C=1`、`K=1` 或 `2` 个宏观源；`J∈R^K`，`L∈R^{3×K}`。`A_C` 必须稳定，`H_C` 负责把皮层状态映射到源活动。

- 参数：`τ_L,δ_L,W_L,b_L,τ_C,A_C,B_C,b_C,H_C,L,d`。
- 物理含义：可分离的丘脑输入、皮层滤波/耦合和头皮混合。
- 文献支持：Lien–Scanziani 对丘脑输入—皮层放大提供 B 级机制；Breakspear 对隐藏群体状态和 `y=LJ+ε` 提供 C 级方法依据。
- 优点：比 A 有可解释的上游输入，仍能通过线性系统进行稳定比较。
- 风险：`W_L`、`B_C`、`H_C` 和 `L` 存在乘积补偿；三通道不足以识别多个自由源。
- Left/Right：差异可来自 `φ_q`，也可在受限版本中只允许 `W_L` 的条件增益变化；不能同时放开所有层级的左右特异参数。
- 与 Q3：可把 `x_L,x_C,J` 作为上游状态摘要，但不引入海马、工作记忆或决策状态。

#### Model C：LGN + Wilson–Cowan E/I + low-dimensional source

LGN 层同 Model B；皮层层使用一个局部兴奋/抑制群体：

\[
\begin{aligned}
\tau_E\dot E &= -E+S_E(w_{EE}E-w_{EI}I+g_Ex_L+b_E)+\eta_E,\\
\tau_I\dot I &= -I+S_I(w_{IE}E-w_{II}I+g_Ix_L+b_I)+\eta_I,\\
J_k(t)&=\alpha_{E,k}E(t)-\alpha_{I,k}I(t),\\
y(t)&=LJ(t)+d+\varepsilon(t).
\end{aligned}
\]

`S_E,S_I` 是平滑 sigmoid。主版本只使用一个 E/I 群体和 `K=1` 或 `2` 个源；多群体网络不作为首选。`J_k=α_EE-α_II` 是把群体活动当作宏观电流代理的建模假设，并非四篇论文给出的头皮源方程。

- 参数：`τ_E,τ_I,w_{EE},w_{EI},w_{IE},w_{II},g_E,g_I,b_E,b_I,α_E,α_I,L,d`。
- 物理含义：E/I 相互作用、丘脑驱动、非线性增益和宏观源代理。
- 文献支持：Breakspear 明确介绍 neural mass、sigmoid、Wilson–Cowan/Jansen–Rit、隐藏状态和 forward model（C）；Lien–Scanziani 支持丘脑驱动被皮层回路放大的方向（B）。
- 优点：可解释增益、抑制、振荡与状态调制，适合检验非线性是否改善波形和左右差异解释。
- 风险：参数多、易产生多稳态/补偿，三通道 ERP 不能唯一确定全部 E/I 参数；不允许用手工调参追求某个记录的拟合。
- Left/Right：允许输入方向差异传播到 E/I 状态；`φ_L=φ_R` 的消融应给出相同驱动下的受限比较。
- 与 Q3：只输出共享的动态状态和可解释增益；Q3 若需要记忆/决策必须另建模型，不在本模型中偷加。

### Layer 4：brain source → scalp EEG observation

三个候选模型统一使用

\[
y(t)=LJ(t)+d+\varepsilon(t),
\qquad L\in\mathbb{R}^{3\times K}.
\]

`J(t)` 是低维宏观源代理，`L` 是相对 lead-field/mixing matrix，`d` 是通道基线，`ε` 是观测噪声。由于没有个体 MRI、头模型、精确源坐标和高密度 EEG，该层只能解释为低维参数化混合，不能解释为精确 source localization。

主比较从 `K=1` 开始，再比较 `K=2`；`K=3` 只作为复杂度上限。每列固定 `||L_:k||₂=1`，并固定符号约定，从而把无法观测的绝对尺度留给 source amplitude。F3、Fz、F4 只作为三个混合观测通道，不对应三个独立脑区。

## 3. Task1/Task2：共享骨架 + 小调制

两个 task 共用 `φ_q` 的编码形式、LGN/皮层状态维度、源数候选和 `L` 的拓扑。允许变化的最小集合为：

| 参数组 | 共享/调制策略 | 理由 |
|---|---|---|
| `τ_L, τ_E, τ_I` | 默认共享；仅在留出数据支持时作小幅 task modulation | 它们是有效动力学时间尺度，独立估计会制造不可辨识自由度 |
| `L` 的列方向与符号规范 | 共享 | 头皮几何不应随 task 改变 |
| 输入增益 `g_{in,τ}`、baseline `b_τ` | 可调制 | 对应注意/警觉状态的简化表示；Horrocks 支持状态改变群体时间轨迹，但不提供本题数值 |
| 少量耦合调制 `Δw_τ` 或阻尼项 | 可选、强正则 | 只在波形形状/振荡差异稳定时启用 |
| `φ_L-φ_R` | 跨 task 共享方向，幅度可弱调制 | 避免为每个 task 重新定义刺激编码 |

## 4. 未来目标函数（只定义，不执行）

在 Q1 提供完整 waveform 后，候选模型的预测为 `ŷ_{r,τ,q}(t;θ)`。建议先定义未加权的多通道波形误差：

\[
J_{wave}(θ)=\sum_{r,τ,q}\sum_t
\|\hat y_{r,τ,q}(t;θ)-y^{Q1}_{r,τ,q}(t)\|_2^2.
\]

可报告但暂不拍定权重的诊断项包括：250–500 ms 窗口幅度、峰潜伏期、左右差异波形能量和 F3/Fz/F4 空间模式。若使用 Q1 bootstrap CI，可将残差按 CI 或估计协方差标准化，而不是人为放大某一时间窗。

## 5. 候选特征

| 特征 | 类型 | 解释 |
|---|---|---|
| `A_F3, A_Fz, A_F4`（250–500 ms） | 直接 EEG | 三通道窗口平均/峰值；不宣称是特定 ERP 成分 |
| `D_LR=y_R-y_L` 的波形能量、峰值和峰潜伏期 | 直接 EEG | 检验左右差异是否跨记录稳定 |
| `D_M=A_Fz-(A_F3+A_F4)/2` | 直接 EEG | 描述中线相对两侧的空间模式 |
| `max x_L`, `max E/I`, `max J_k` | 模型派生 | 只能作为潜在状态特征，不能当作直接测量 |
| latent trajectory speed/tangling | 模型派生/描述性 | 借鉴 Horrocks 的群体轨迹思想；需在三通道低维状态上谨慎解释 |

## 6. 消融和选模规则

至少预注册以下消融：

1. 去掉 LGN 动力学，直接由 `u_q` 驱动皮层；
2. 用线性皮层替代 E/I neural mass（C → B）；
3. 强制 `φ_L=φ_R`；
4. 强制 Left/Right 共享 source spatial pattern；
5. 比较 `K=1`、`K=2`（必要时 `K=3`）并加入复杂度惩罚；
6. 去掉 task-specific gain/baseline。

选模依据为留出记录/条件上的完整波形误差、校准后的不确定性、复杂度惩罚和左右特征的可复现性。若复杂模型只改善训练记录而不改善留出记录，保留更简单模型；若 B 与 C 接近，优先 B。
