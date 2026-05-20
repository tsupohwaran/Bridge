# 模型概览、经济结构与求解器

Last updated: 2026-05-20.

> 主要依据：`docs/MODEL_NOTES.md`、`code/03_model/functions.jl`、`code/03_model/main.jl`。  
> 参考格式：`/Users/pohwaran/Doctorate/Paper/Electrification/manuscript/algorithms/solver.md`。  
> 说明：本文档只整理 Bridge 项目当前 Julia 实现。参考文件只用于章节组织，不引入其贸易模型设定。

## 目录

1. [模型概览](#1-模型概览)
2. [经济结构](#2-经济结构)
   - 2.1 [主体与维度](#21-主体与维度)
   - 2.2 [变量](#22-变量)
   - 2.3 [参数](#23-参数)
   - 2.4 [工人选择](#24-工人选择)
   - 2.5 [劳动供给弹性与市场势力](#25-劳动供给弹性与市场势力)
   - 2.6 [企业问题](#26-企业问题)
   - 2.7 [均衡条件](#27-均衡条件)
   - 2.8 [桥梁/隧道反事实](#28-桥梁隧道反事实)
3. [求解器](#3-求解器)
   - 3.1 [求解思路](#31-求解思路)
   - 3.2 [数据预处理与归一化](#32-数据预处理与归一化)
   - 3.3 [反推企业 amenity 与生产率](#33-反推企业-amenity-与生产率)
   - 3.4 [基准均衡](#34-基准均衡)
   - 3.5 [反事实均衡](#35-反事实均衡)
   - 3.6 [模型矩构造](#36-模型矩构造)
   - 3.7 [参数校准](#37-参数校准)
   - 3.8 [收敛与数值稳定](#38-收敛与数值稳定)
4. [伪代码](#4-伪代码)
5. [实现说明与待确认事项](#5-实现说明与待确认事项)

## 1. 模型概览

已确认：

- 当前结构模型在 `code/03_model/functions.jl` 中实现，由 `code/03_model/main.jl` 调用。
- 模型主体包括居住地/乡镇 `z` 和企业 `j`。
- 工人根据工资、企业非货币 amenity 和通勤时间选择企业。
- 企业面对向上倾斜的劳动供给曲线，工资由边际产出和劳动供给弹性共同决定。
- 当前通勤成本直接使用分钟数 `d_zj`，而不是 `ln(d_zj)`。
- 当前效用函数包含企业层面的非货币 amenity `a_j`，并在基准数据中由企业就业份额反推。

推断：

- 模型的核心机制是：桥梁/隧道降低部分居住地到企业的通勤时间，改变工人选择概率和企业劳动力来源结构，从而改变企业劳动供给弹性、工资和就业。
- 校准目标试图让模型生成的就业变化在市场可达性改善和初始工资异质性上的模式接近经验回归矩。

未确定：

- 论文最终应将 `a_j` 称为 firm amenity、残差企业吸引力、补偿性差异，还是其他术语。
- 最终模型是否继续使用当前一层 logit 结构，还是回到旧 draft 中的嵌套结构。

## 2. 经济结构

### 2.1 主体与维度

已确认：

| 符号 | 维度 | 含义 |
|---|---:|---|
| `z` | `Z` | 工人居住地/乡镇/原点位置 |
| `j` | `J` | 企业 |
| `l_z` | `Z` | 居住地人口权重，当前在 `main.jl` 中归一化为总和 1 |
| `d_zj` | `Z x J` | 从居住地 `z` 到企业 `j` 的通勤时间，单位为分钟 |
| `d'_zj` | `Z x J` | 桥梁/隧道情形下的通勤时间，单位为分钟 |

未确定：

- `z` 在论文文字中应严格解释为居住地、乡镇政府点，还是市场可达性计算中的人口加权原点。

### 2.2 变量

已确认：

**外生或由数据给定：**

| 符号 | 代码名 | 含义 |
|---|---|---|
| `l_z` | `l` | 居住地人口权重 |
| `d_zj` | `d` | 基准通勤时间，分钟 |
| `d'_zj` | `d′` | 反事实通勤时间，分钟 |
| `w_j^data` | `wⱼ_data` | 基准企业工资，当前归一化为总工资账单为 1 |
| `l_j^data` | `lⱼ_data` | 基准企业就业份额，当前归一化为总和 1 |

**内生或由模型反推：**

| 符号 | 代码名 | 含义 |
|---|---|---|
| `pi_zj` | `π_zj` | 居住地 `z` 工人选择企业 `j` 的概率 |
| `l_j` | `lⱼ` | 企业就业 |
| `epsilon_zj` | `ε_zj` | 居住地-企业层面的劳动供给弹性 |
| `epsilon_j` | `εⱼ` | 企业层面的劳动供给弹性 |
| `w_j` | `wⱼ` | 均衡工资 |
| `a_j` | `aⱼ` | 企业非货币 amenity，均值归一化为 0 |
| `z_j` | `zⱼ` | 企业生产率，均值归一化为 1 |

### 2.3 参数

已确认：

| 参数 | 含义 | 当前设置/说明 |
|---|---|---|
| `alpha` / `α` | 企业生产函数的劳动收益递减参数 | `main.jl` 中设为 `0.4` |
| `eta` / `η` | 每分钟通勤时间的效用成本 | 校准参数 |
| `theta` / `θ` | 偏好离散度/工人选择响应强度 | 校准参数 |

推断：

- `η` 越大，额外通勤分钟对工人效用的惩罚越大。
- `θ` 越大，工人对工资、amenity 和通勤时间差异越敏感。

### 2.4 工人选择

已确认：

当前 Julia 实现对应的工人效用为：

```text
U_izj = ln(w_j) + a_j - η d_zj + (1 / θ) ε_ij
```

因此，居住地 `z` 的工人选择企业 `j` 的概率为：

```text
π_zj = exp(θ [ln(w_j) + a_j - η d_zj])
       / sum_k exp(θ [ln(w_k) + a_k - η d_zk])
```

等价写法为：

```text
π_zj = (w_j exp(a_j - η d_zj))^θ
       / sum_k (w_k exp(a_k - η d_zk))^θ
```

企业就业由各居住地选择概率加总得到：

```text
l_j = sum_z π_zj l_z
```

已确认：

- `WorkerChoice(wⱼ, l, d, aⱼ, params)` 显式要求传入 `aⱼ`。
- 当前代码不再为缺失 amenity 自动使用零向量。

### 2.5 劳动供给弹性与市场势力

已确认：

在当前一层 logit 结构下，居住地-企业层面的工资弹性为：

```text
ε_zj = θ (1 - π_zj)
```

企业 `j` 的工人来源权重为：

```text
γ_zj = π_zj l_z / l_j
```

企业层面的劳动供给弹性为：

```text
ε_j = sum_z γ_zj ε_zj
```

代码中还计算或使用与市场势力相关的对象：

```text
ν_j = 1 + 1 / ε_j
```

推断：

- 如果企业面对的有效劳动供给弹性较低，则企业工资相对边际产出的折价更大，劳动市场势力更强。
- 桥梁/隧道通过改变 `π_zj` 和 `γ_zj`，可以改变 `ε_j`，即使企业生产率和 amenity 本身保持不变。

未确定：

- 论文中应将 `ν_j` 叫作 markdown、inverse markdown、monopsony wedge 还是 labor-market-power wedge。

### 2.6 企业问题

已确认：

当前模型使用的企业生产函数为：

```text
Y_j(l_j) = z_j l_j^α
```

工资一阶条件在代码中写为：

```text
w_j = α z_j l_j^(α - 1) ε_j / (1 + ε_j)
```

推断：

- 这是一条带有买方劳动市场势力的工资条件：工资等于劳动边际产出乘以弹性调整项。
- `a_j` 不直接进入企业一阶条件；它通过工人选择概率、企业就业和劳动供给弹性影响均衡。

### 2.7 均衡条件

已确认：

给定参数 `(η, θ, α)`、人口 `l_z`、通勤矩阵 `d_zj`、企业生产率 `z_j` 和企业 amenity `a_j`，均衡满足：

1. 工人选择概率由 logit 公式给出。
2. 企业就业满足 `l_j = sum_z π_zj l_z`。
3. 企业劳动供给弹性由 `ε_zj` 和 `γ_zj` 加权得到。
4. 企业工资满足 `w_j = α z_j l_j^(α - 1) ε_j / (1 + ε_j)`。
5. 工资在求解器中归一化，使总工资账单为 1。

当前 `SolveModel` 通过工资固定点迭代求解上述条件。

### 2.8 桥梁/隧道反事实

已确认：

- 数据准备代码将无桥/隧道通勤时间记为 `dzj`，将有桥/隧道通勤时间记为 `dzj_prime`。
- 当前模型代码中相应数组为 `d` 和 `d′`。
- 市场可达性改善按分钟计算：

```text
dMA_j = sum_z l_z (d_zj - d'_zj)
```

- 当前模型矩中使用：

```text
bigMA_j = 1[dMA_j >= 0.5]
```

推断：

- 正的 `dMA_j` 表示桥梁/隧道降低了到企业 `j` 的人口加权通勤时间。
- 反事实中通常保持 `z_j` 和 `a_j` 固定，只改变通勤时间矩阵。

未确定：

- `bigMA` 的最终阈值是否应固定为 `0.5` 分钟，还是应与经验回归中的处理变量定义完全一致。

## 3. 求解器

### 3.1 求解思路

已确认：

- 当前 Bridge 项目使用 level fixed-point solver，不是 Electrification 参考文档中的 exact hat-algebra solver。
- 校准时，对每一个候选参数 `(η, θ)`，先用基准就业和工资反推企业 amenity 与生产率，再求解通勤时间变化后的反事实均衡。
- `ComputeModelMoments` 在校准路径上直接使用反推得到的基准对象，不额外求解一次基准固定点。

求解流程为：

1. 读取并归一化基准人口、就业、工资和通勤时间。
2. 给定候选 `(η, θ)`。
3. 由 `l_j^data` 和 `w_j^data` 反推 `a_j`。
4. 由工资一阶条件反推 `z_j`。
5. 使用 `d` 作为基准通勤矩阵，使用 `d′` 作为反事实通勤矩阵。
6. 求解反事实均衡工资、就业和劳动供给弹性。
7. 用基准和反事实企业就业构造两期 DID 型模型矩。
8. 选择使模型矩接近经验目标矩的 `(η, θ)`。

### 3.2 数据预处理与归一化

已确认：

`main.jl` 当前使用：

```text
l_z       = pop / sum_z pop
l_j^data  = employ_j / sum_j employ_j
w_j^data  = wage_j / sum_j wage_j l_j^data
```

其中：

- `l_z` 总和为 1。
- `l_j^data` 总和为 1。
- `w_j^data` 被归一化，使基准总工资账单为 1。
- `d_zj` 和 `d'_zj` 的单位是分钟，零值被替换为 `1e-2`。

### 3.3 反推企业 amenity 与生产率

已确认：

`SolveAmenitiesFromEmployment` 给定 `w_j^data`、`l_z`、`d_zj` 和 `l_j^data`，寻找 `a_j` 使模型就业匹配数据就业：

```text
l_j^data = sum_z π_zj(w^data, a, d) l_z
```

amenity 归一化为：

```text
mean_j(a_j) = 0
```

数值上，代码在 `logq_j = θ a_j` 上迭代：

```text
logq_j(new) = logq_j(old) + log(l_j^data) - log(l_j^model)
```

每次迭代后减去均值以保留归一化。

已确认：

`SolveFirmPrimitivesFromData` 在得到 `a_j` 后，使用工资一阶条件反推生产率：

```text
z_j = w_j^data (1 + ε_j) / [α l_j^(α - 1) ε_j]
```

并将 `z_j` 均值归一化为 1。

推断：

- `a_j` 吸收了在工资和通勤时间之外使企业就业更高或更低的企业层面吸引力。
- `z_j` 是使观测工资与企业一阶条件一致的生产率对象，尺度由归一化决定。

### 3.4 基准均衡

已确认：

在 `ComputeModelMoments` 中，基准状态直接使用：

```text
w_j = w_j^data
l_j = l_j^data  # 实际来自 amenity inversion 后的 primitives.l_j
```

原因是 `a_j` 被构造为匹配基准就业，`z_j` 被构造为匹配基准工资一阶条件。因此，校准矩计算中不需要重复求解基准固定点。

已确认：

在非校准模拟部分，`main.jl` 仍可调用 `SolveModel` 求解基准均衡，并使用 `w_j^data` 作为 warm start。

### 3.5 反事实均衡

已确认：

反事实保持 `l_z`、`z_j`、`a_j` 和参数不变，只将通勤时间从 `d` 改为 `d′`。

`ComputeModelMoments` 为了数值稳定，采用通勤时间路径 continuation：

```text
d_path(s) = (1 - s) d + s d′,  s = 1/S, 2/S, ..., 1
```

每一步都以前一步的工资作为初始值，调用：

```text
SolveModel((; l, d = d_path, z_j, a_j), params)
```

已确认：

`SolveModel` 的工资更新为：

```text
w_j(new) = α z_j l_j^(α - 1) ε_j / (1 + ε_j)
```

随后归一化：

```text
w_j(new) = w_j(new) / sum_j w_j(new) l_j
```

### 3.6 模型矩构造

已确认：

当前 `ComputeModelMoments` 使用基准和反事实企业就业构造两期 DID 型矩。令：

```text
Δlnl_j = ln(l'_j) - ln(l_j)
w_diff_j = ln(w_j) - mean(ln(w_j))
```

模型矩对应两期面板回归：

```text
lnl_jt ~ bigMA_j × post_t
       + w_diff_j × bigMA_j × post_t
       + firm fixed effects
       + year fixed effects
       + w_diff_j × year_t slopes
```

在两期差分形式下，代码通过 Frisch-Waugh-Lovell 残差化实现：

```text
Δlnl_j = β1 bigMA_j + β2 (w_diff_j × bigMA_j)
         + absorbed constant and w_diff_j + residual
```

模型返回：

```text
β_model = [β1, β2]
```

已确认：

- `wage_center` 可设为 `:all` 或 `:treated`。
- `moment_firm_mask` 可将模型矩限制在匹配的经验回归样本企业上。

### 3.7 参数校准

已确认：

当前校准目标为：

```text
β_target = [-0.092242, 0.0939565]
```

目标函数为：

```text
Objective(η, θ) = sum_m (β_model,m(η, θ) - β_target,m)^2
```

当前 `main.jl` 使用网格搜索寻找局部优化起点，并用 `Optim.NelderMead()` 在参数边界内优化。

当前参数边界和网格为：

```text
η_bounds = [0.005, 0.25]
η_grid   = [0.02, 0.04, 0.05, 0.06, 0.08, 0.10, 0.12, 0.16, 0.20]

θ_bounds = [0.25, 5.0]
θ_grid   = [0.75, 1.5, 2.5, 3.5, 4.5]
```

已确认：

- `CalibrateEtaTheta` 可以接受 `x0 = nothing`。
- 如果提供 `starts`，局部优化从这些起点开始。
- 优化器内部通过 logit 变换处理有界参数。

推断：

- 因为 `a_j` 会随每个候选 `(η, θ)` 重新反推，当前两个就业矩对 `θ` 的约束较弱。
- 如果局部优化持续把 `θ` 推向上界，应将其解释为边界解或识别不足，而不是精确的内点估计。

未确定：

- `β_target` 的最终来源和是否为论文最终目标矩仍需确认。
- 是否应加入外部劳动供给弹性矩或其他矩来约束 `θ`。

### 3.8 收敛与数值稳定

已确认：

- `WorkerChoice` 使用 log-sum-exp 技巧计算选择概率，避免 `exp(large number)` 溢出。
- `SolveModel` 使用 damping 和 warm start。
- `ComputeModelMoments` 对反事实通勤变化使用 continuation，而不是一次性从 `d` 跳到 `d′`。
- 如果 amenity inversion、反事实固定点或模型矩出现 `NaN`、`Inf` 或过大绝对值，目标函数返回惩罚值。

## 4. 伪代码

```text
输入：
    l, d, d′, w_data, l_data, α, β_target
    η_grid, θ_grid, η_bounds, θ_bounds

网格搜索：
    对每个 η in η_grid:
        对每个 θ in θ_grid:
            params = (η, θ, α)

            # 反推基准原语
            求 a，使 sum_z π_zj(w_data, a, d; params) l_z = l_data_j
            计算 ε_j
            由工资 FOC 反推 z_j

            # 基准状态
            w_base = w_data
            l_base = l_data

            # 反事实状态
            w_start = w_base
            对 s = 1/S, ..., 1:
                d_path = (1 - s) d + s d′
                解工资固定点：
                    π_zj = WorkerChoice(w, l, d_path, a; params)
                    l_j = sum_z π_zj l_z
                    ε_j = sum_z (π_zj l_z / l_j) θ(1 - π_zj)
                    w_j = α z_j l_j^(α - 1) ε_j / (1 + ε_j)
                    归一化 w_j
                w_start = w_j

            # 模型矩
            Δlnl = ln(l_cf) - ln(l_base)
            dMA = sum_z l_z (d_zj - d′_zj)
            bigMA = 1[dMA >= 0.5]
            w_diff = ln(w_base) - center(ln(w_base))
            β_model = FWL_regression(Δlnl on bigMA and bigMA*w_diff,
                                     absorbing constant and w_diff)

            objective = sum((β_model - β_target)^2)

局部优化：
    选取网格搜索中 objective 最小的若干点作为 starts
    用 Nelder-Mead 在 η_bounds 和 θ_bounds 内最小化 objective

输出：
    η_est, θ_est
    反推的 a_j, z_j
    基准与反事实的 w_j, l_j, ε_j
    模型矩 β_model
```

## 5. 实现说明与待确认事项

已确认：

- `SolveModel`、`SolveZfromW` 和相关函数现在要求 `vars` 中显式包含 `a_j`。
- `a_j` 是求解器需要的模型原语之一，不是事后报告变量。
- 当前校准和反事实全部使用分钟单位的 `d`。
- 当前模型不实现旧草稿中的嵌套参数 `σ`。

推断：

- 新增 `a_j` 后，基准就业可以被更精确地匹配，但也使 `θ` 的识别更依赖反事实矩或额外外部矩。
- 使用分钟单位后，`η` 的数量级应比旧的 `ln(d)` 设定小很多。

未确定：

- 最终 treatment 定义、样本限制、`β_target` 来源和 `bigMA` 阈值仍需与经验部分统一。
- 论文正文是否应报告 `a_j` 的分布或只把它作为校准残差。
- 是否应为 `θ` 增加外部目标，例如平均劳动供给弹性、markdown 水平，或文献中的弹性范围。
