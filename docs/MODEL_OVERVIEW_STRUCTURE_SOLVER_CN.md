# 模型求解 README：经济结构与求解器

Last updated: 2026-05-21.

本文档说明 Bridge 项目当前结构模型的经济设定、数据接口、求解流程和校准矩构造。它面向论文写作和代码复现，主要对应 `code/03_model/functions.jl` 与 `code/03_model/main.jl`。

参考文档：`docs/MODEL_NOTES.md`、`docs/USER_CORRECTIONS.md`。本文档只描述 Bridge 项目当前实现，不引入其他项目或旧草稿中的模型设定。

## 目录

1. [用途与范围](#1-用途与范围)
2. [经济环境](#2-经济环境)
   - [2.1 主体与索引](#21-主体与索引)
   - [2.2 输入、状态变量与模型对象](#22-输入状态变量与模型对象)
   - [2.3 参数](#23-参数)
3. [模型结构](#3-模型结构)
   - [3.1 工人选择](#31-工人选择)
   - [3.2 劳动供给弹性与 markdown](#32-劳动供给弹性与-markdown)
   - [3.3 企业问题](#33-企业问题)
   - [3.4 均衡条件](#34-均衡条件)
   - [3.5 桥梁/隧道反事实](#35-桥梁隧道反事实)
4. [求解器接口](#4-求解器接口)
   - [4.1 数据预处理与归一化](#41-数据预处理与归一化)
   - [4.2 反推 firm amenity 与生产率](#42-反推-firm-amenity-与生产率)
   - [4.3 基准均衡](#43-基准均衡)
   - [4.4 反事实均衡](#44-反事实均衡)
5. [模型矩与参数校准](#5-模型矩与参数校准)
   - [5.1 模型矩](#51-模型矩)
   - [5.2 目标函数](#52-目标函数)
   - [5.3 参数边界与起点](#53-参数边界与起点)
6. [算法流程](#6-算法流程)
7. [实现约定与维护事项](#7-实现约定与维护事项)

## 1. 用途与范围

本 README 记录当前 Julia 求解器所实现的一层 logit 结构模型。模型用于研究桥梁/隧道带来的通勤时间变化如何影响企业面对的劳动供给弹性、工资、就业和 markdown。

论文文字中的术语约定如下：

- $z$ 表示地区，即街道/乡镇；在代码中对应通勤矩阵和人口权重中的 worker origin。
- $j$ 表示企业。
- $a_j$ 称为 firm amenity。
- $\nu_j = 1 + 1 / \varepsilon_j$ 称为 markdown，且 $\nu_j > 1$。
- 最终模型使用当前一层 logit 结构，不使用旧草稿中的嵌套 logit 参数 $\sigma$。
- 模型侧 $bigMA$ 阈值固定为 0.5 分钟。

$\beta_{\text{target}}$ 来自实证结果。当前本机实证代码可能尚未更新到直接生成同一组目标值，后续实证代码更新后应同步核对。

## 2. 经济环境

### 2.1 主体与索引

| 符号 | 维度 | 代码名 | 含义 |
|---|---:|---|---|
| $z$ | $Z$ | 行维度 | 地区，即街道/乡镇 |
| $j$ | $J$ | 列维度 | 企业 |
| $l_z$ | $Z$ | `l` | 地区人口权重，归一化为总和 1 |
| $d_{zj}$ | $Z \times J$ | `d` | 无桥/隧道情形下，从地区到企业的通勤时间，单位为分钟 |
| $d'_{zj}$ | $Z \times J$ | `d′` | 有桥/隧道情形下，从地区到企业的通勤时间，单位为分钟 |

### 2.2 输入、状态变量与模型对象

| 对象 | 代码名 | 含义 |
|---|---|---|
| $w_j^{data}$ | `wⱼ_data` | 基准企业工资，归一化为基准总工资账单为 1 |
| $l_j^{data}$ | `lⱼ_data` | 基准企业就业份额，归一化为总和 1 |
| $\pi_{zj}$ | `π_zj` | 地区 $z$ 工人选择企业 $j$ 的概率 |
| $l_j$ | `lⱼ` | 企业就业 |
| $\varepsilon_{zj}$ | `ε_zj` | 地区-企业层面的劳动供给弹性 |
| $\varepsilon_j$ | `εⱼ` | 企业层面的劳动供给弹性 |
| $w_j$ | `wⱼ` | 均衡工资 |
| $a_j$ | `aⱼ` | firm amenity，均值归一化为 0 |
| $z_j$ | `zⱼ` | 企业生产率，均值归一化为 1 |
| $\nu_j$ | `νⱼ` | markdown，等于 $1 + 1 / \varepsilon_j$ |

### 2.3 参数

| 参数 | 代码名 | 含义 | 当前设置 |
|---|---|---|---|
| $\alpha$ | `α` | 企业生产函数的劳动收益递减参数 | `0.4` |
| $\eta$ | `η` | 通勤成本参数 | 校准参数 |
| $\theta$ | `θ` | 偏好离散度/选择响应强度 | 校准参数 |

$\eta$ 越大，通勤时间对工人效用的惩罚越大。$\theta$ 越大，工人对工资、firm amenity 和通勤时间差异越敏感。

## 3. 模型结构

### 3.1 工人选择

当前代码中，通勤时间以分钟计量，但进入效用的是通勤时间的对数。工人效用为：

$$
\begin{aligned}
U_{izj}
&= \log w_j + a_j - \eta \log d_{zj}
+ \frac{1}{\theta}\varepsilon_{ij}.
\end{aligned}
$$

地区 $z$ 的工人选择企业 $j$ 的概率为：

$$
\begin{aligned}
\pi_{zj}
&=
\frac{
\exp\{\theta(\log w_j + a_j - \eta \log d_{zj})\}
}{
\sum_k \exp\{\theta(\log w_k + a_k - \eta \log d_{zk})\}
}.
\end{aligned}
$$

等价地：

$$
\begin{aligned}
\pi_{zj}
&=
\frac{
(w_j e^{a_j} d_{zj}^{-\eta})^\theta
}{
\sum_k (w_k e^{a_k} d_{zk}^{-\eta})^\theta
}.
\end{aligned}
$$

企业就业由各地区选择概率加总得到：

$$
\begin{aligned}
l_j &= \sum_z \pi_{zj} l_z.
\end{aligned}
$$

代码接口 `WorkerChoice(wⱼ, l, d, aⱼ, params)` 要求显式传入 $a_j$。求解器不再为缺失 firm amenity 自动使用零向量。

### 3.2 劳动供给弹性与 markdown

一层 logit 结构下，地区-企业层面的工资弹性为：

$$
\begin{aligned}
\varepsilon_{zj} &= \theta(1 - \pi_{zj}).
\end{aligned}
$$

企业 $j$ 的工人来源权重为：

$$
\begin{aligned}
\gamma_{zj} &= \frac{\pi_{zj} l_z}{l_j}.
\end{aligned}
$$

企业层面的劳动供给弹性为：

$$
\begin{aligned}
\varepsilon_j &= \sum_z \gamma_{zj}\varepsilon_{zj}.
\end{aligned}
$$

结构模型中的 markdown 定义为：

$$
\begin{aligned}
\nu_j &= 1 + \frac{1}{\varepsilon_j}, \qquad \nu_j > 1.
\end{aligned}
$$

有效劳动供给弹性越低，markdown 越高。桥梁/隧道通过改变选择概率和工人来源结构，进而改变企业劳动供给弹性与 markdown。

### 3.3 企业问题

企业生产函数为：

$$
\begin{aligned}
Y_j(l_j) &= z_j l_j^\alpha.
\end{aligned}
$$

工资一阶条件为：

$$
\begin{aligned}
w_j
&=
\alpha z_j l_j^{\alpha - 1}
\frac{\varepsilon_j}{1 + \varepsilon_j}.
\end{aligned}
$$

Firm amenity 不直接进入企业一阶条件。它通过工人选择概率、企业就业和劳动供给弹性影响均衡工资和就业。

### 3.4 均衡条件

给定参数、地区人口、通勤矩阵、企业生产率和 firm amenity，均衡满足：

1. 工人选择概率由一层 logit 给出。
2. 企业就业等于各地区选择概率的加权和。
3. 企业劳动供给弹性由地区-企业弹性按工人来源权重加总。
4. 企业工资满足企业一阶条件。
5. 求解器归一化工资，使总工资账单为 1。

`SolveModel` 通过工资固定点迭代求解上述条件。

### 3.5 桥梁/隧道反事实

数据准备代码将无桥/隧道通勤时间记为 `dzj`，将有桥/隧道通勤时间记为 `dzj_prime`。模型代码中相应数组为 `d` 和 `d′`。

市场可达性改善定义为人口加权节省分钟数：

$$
\begin{aligned}
dMA_j &= \sum_z l_z(d_{zj} - d'_{zj}).
\end{aligned}
$$

模型矩中的处理变量为：

$$
\begin{aligned}
bigMA_j &= \mathbf{1}\{dMA_j \geq 0.5\}.
\end{aligned}
$$

反事实只改变通勤时间矩阵，从 $d$ 改为 $d'$；地区人口、生产率、firm amenity 和参数保持不变。

## 4. 求解器接口

### 4.1 数据预处理与归一化

`main.jl` 读取模型样本后执行以下归一化：

$$
\begin{aligned}
l_z &= \frac{pop_z}{\sum_{z'} pop_{z'}}, \\
l_j^{data} &= \frac{employ_j}{\sum_k employ_k}, \\
w_j^{data} &= \frac{wage_j}{\sum_k wage_k l_k^{data}}.
\end{aligned}
$$

归一化后，地区人口权重总和为 1，企业就业份额总和为 1，基准总工资账单为 1。通勤时间单位为分钟；零值在取对数前替换为 `1e-2`。

### 4.2 反推 firm amenity 与生产率

`SolveAmenitiesFromEmployment` 给定基准工资、地区人口、基准通勤矩阵和基准企业就业，求解 firm amenity，使模型基准就业匹配数据就业：

$$
\begin{aligned}
l_j^{data}
&=
\sum_z
\pi_{zj}(w^{data}, a, d; \eta,\theta)l_z.
\end{aligned}
$$

Firm amenity 的归一化为：

$$
\begin{aligned}
\frac{1}{J}\sum_j a_j &= 0.
\end{aligned}
$$

数值实现中，代码在 $\log q_j = \theta a_j$ 上迭代：

$$
\begin{aligned}
\log q_j^{new}
&=
\log q_j^{old}
+ \log l_j^{data}
- \log l_j^{model}.
\end{aligned}
$$

每次迭代后减去 $\log q_j$ 的均值，以维持归一化。

`SolveFirmPrimitivesFromData` 在得到 firm amenity 后，用企业一阶条件反推生产率：

$$
\begin{aligned}
z_j
&=
\frac{
w_j^{data}(1+\varepsilon_j)
}{
\alpha l_j^{\alpha - 1}\varepsilon_j
}.
\end{aligned}
$$

随后将生产率均值归一化为 1。

### 4.3 基准均衡

在 `ComputeModelMoments` 中，基准状态直接采用观测工资和数据就业：

$$
\begin{aligned}
w_j &= w_j^{data}, \\
l_j &= l_j^{data}.
\end{aligned}
$$

这样处理是因为 firm amenity 已经被构造为匹配基准就业，生产率也被构造为匹配基准工资一阶条件。校准矩计算不需要重复求解基准固定点。

在非校准模拟部分，`main.jl` 仍可调用 `SolveModel` 求解基准均衡，并使用 $w_j^{data}$ 作为 warm start。

### 4.4 反事实均衡

反事实保持 $l_z$、$z_j$、$a_j$ 和参数不变，只将通勤矩阵从 $d$ 改为 $d'$。为提高数值稳定性，`ComputeModelMoments` 使用 continuation 路径：

$$
\begin{aligned}
d(s) &= (1-s)d + s d',
\qquad
s \in \left\{\frac{1}{S}, \frac{2}{S}, \ldots, 1\right\}.
\end{aligned}
$$

每一步以前一步工资为初始值调用 `SolveModel`。工资更新规则为：

$$
\begin{aligned}
w_j^{new}
&=
\alpha z_j l_j^{\alpha - 1}
\frac{\varepsilon_j}{1 + \varepsilon_j}.
\end{aligned}
$$

随后执行工资账单归一化：

$$
\begin{aligned}
w_j^{new}
&\leftarrow
\frac{w_j^{new}}{\sum_k w_k^{new} l_k}.
\end{aligned}
$$

## 5. 模型矩与参数校准

### 5.1 模型矩

`ComputeModelMoments` 使用基准和反事实企业就业构造两期 DID 型模型矩。定义：

$$
\begin{aligned}
\Delta \log l_j &= \log l'_j - \log l_j, \\
w\_diff_j &= \log w_j - \overline{\log w}.
\end{aligned}
$$

模型矩对应的两期面板回归为：

$$
\begin{aligned}
\log l_{jt}
&=
\beta_1(bigMA_j \times post_t)
+ \beta_2(w\_diff_j \times bigMA_j \times post_t) \\
&\quad
+ \lambda_j
+ \tau_t
+ \rho_t w\_diff_j
+ u_{jt}.
\end{aligned}
$$

在两期差分形式下，代码通过 Frisch-Waugh-Lovell 残差化估计：

$$
\begin{aligned}
\Delta \log l_j
&=
\beta_1 bigMA_j
+ \beta_2(w\_diff_j \times bigMA_j)
+ c
+ \rho w\_diff_j
+ e_j.
\end{aligned}
$$

模型返回两个矩：

$$
\begin{aligned}
\beta_{\text{model}}
&=
\begin{bmatrix}
\beta_1 \\
\beta_2
\end{bmatrix}.
\end{aligned}
$$

`wage_center` 可设为 `:all` 或 `:treated`。`moment_firm_mask` 可将模型矩限制在匹配的经验回归样本企业上。

### 5.2 目标函数

当前目标矩为：

$$
\begin{aligned}
\beta_{\text{target}}
&=
\begin{bmatrix}
-0.092242 \\
0.0939565
\end{bmatrix}.
\end{aligned}
$$

目标函数为：

$$
\begin{aligned}
Objective(\eta,\theta)
&=
\sum_m
\left[
\beta_{\text{model},m}(\eta,\theta)
- \beta_{\text{target},m}
\right]^2.
\end{aligned}
$$

当前不加入外部劳动供给弹性矩或其他额外矩来约束 $\theta$。

### 5.3 参数边界与起点

当前 `main.jl` 使用网格搜索选择局部优化起点，并用 `Optim.NelderMead()` 在参数边界内优化。

| 对象 | 当前值 |
|---|---|
| `eta_bounds` | `[0.25, 15.0]` |
| `eta_grid` | `[1.5, 3.0, 4.5, 6.0, 7.5]` |
| `theta_bounds` | `[0.25, 8.0]` |
| `theta_grid` | `[0.75, 1.5, 2.5, 3.5, 4.5]` |
| `x0` | `[8.0, 0.75]` |

由于 firm amenity 会随每个候选参数重新反推，仅使用当前两个就业矩时，$\theta$ 的约束可能较弱。若局部优化把 $\theta$ 推向边界，应将其解释为边界解或识别不足，而不是精确的内点估计。

## 6. 算法流程

给定 $l$、$d$、$d'$、$w^{data}$、$l^{data}$、$\alpha$、$\beta_{\text{target}}$、参数网格和参数边界，校准流程如下：

1. 对每个候选 $\eta$ 和 $\theta$ 设定参数组 $(\eta,\theta,\alpha)$。
2. 用基准工资、基准通勤矩阵和基准企业就业反推 firm amenity。
3. 计算劳动供给弹性，并由企业一阶条件反推生产率。
4. 将观测工资和数据就业作为基准状态。
5. 沿 continuation 路径从 $d$ 过渡到 $d'$，逐步求解反事实均衡。
6. 由基准就业和反事实就业计算企业就业变化。
7. 用人口加权通勤节省计算 $dMA$，并设定 $BIG = \mathbf{1}\{dMA > 0.5\}$。
8. 按 `calculate_calibration_4_moments.do` 构造去均值后的基准工资异质性和 `lndma`。
9. 用 FWL 回归得到四个模型矩。
10. 计算模型矩和目标矩之间的平方距离。
11. 选取网格搜索中目标函数较小的点作为 starts。
12. 用 Nelder-Mead 在参数边界内最小化目标函数。

最终输出包括校准参数、反推的 firm amenity、反推的生产率、基准与反事实的工资、就业、劳动供给弹性、markdown，以及模型矩。

## 7. 实现约定与维护事项

- `SolveModel`、`SolveZfromW` 和相关函数要求 `vars` 中显式包含 $a_j$。
- $a_j$ 是求解器原语，不是事后报告变量。
- 通勤时间以分钟为单位存储；工人效用使用 $\log d_{zj}$。
- 最终模型为一层 logit，不实现旧草稿中的嵌套参数 $\sigma$。
- $z$ 在论文文字中统一解释为地区，即街道/乡镇。
- $\nu_j$ 在论文中统一称为 markdown，且为大于 1 的对象。
- 模型矩中的 $bigMA$ 阈值固定为 0.5 分钟。
- 后续需要与经验部分同步的是 treatment 的文字定义、样本限制，以及本机实证代码对 $\beta_{\text{target}}$ 来源和数值的更新。
