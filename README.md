# Bridge Effect on Labor Reallocation

**胶州湾大桥对劳动力重配置效应的量化空间经济学研究**

## 项目概述

本项目研究胶州湾大桥（2011年通车）对青岛地区企业区位选择、劳动力重配置及土地价格的影响。采用双重差分法(DID)进行实证分析，并构建量化空间经济学模型进行结构估计。

## 目录结构

```
Bridge/
├── code/                    # 所有代码
│   ├── 00_setup/            # 环境配置脚本
│   ├── 01_data_prep/        # 数据准备 (Stata/Python)
│   ├── 02_empirical/        # 实证分析 (Stata)
│   ├── 03_model/            # 结构模型 (Julia)
│   ├── utils/               # 公共函数
│   └── sandbox/             # 草稿/实验代码
├── data/                    # 所有数据
│   ├── raw/                 # 原始数据 (只读)
│   │   ├── cied/            # 工企数据 1998-2014
│   │   ├── ctsd/            # 税调数据 2007-2020
│   │   └── geo/             # 地理数据
│   ├── processed/           # 处理后数据
│   └── temp/                # 临时数据 (已忽略)
├── output/                  # 输出结果
│   ├── figures/             # 图表
│   │   ├── empirical/       # 实证图 (land/firm)
│   │   ├── model/           # 模型图
│   │   └── markdown/        # 其他图表
│   ├── tables/              # 回归表格
│   └── logs/                # 运行日志
├── manuscript/              # 论文手稿
│   ├── paper/               # 正文 (LyX/TeX)
│   ├── archive/             # 历史版本
│   └── notes/               # 研究笔记
├── gis/                     # GIS项目文件
├── Project.toml             # Julia依赖
└── Manifest.toml            # Julia锁文件
```

## 环境配置

### Julia (结构模型)

```bash
cd Bridge
julia --project=. -e 'using Pkg; Pkg.instantiate()'
```

### Stata (实证分析)

```stata
do code/00_setup/PkgInstall.do
```

### Python (Market Access计算)

```bash
pip install -r code/00_setup/requirements.txt
```

## 运行流程

### 1. 数据准备

```stata
// 按顺序执行
do code/01_data_prep/CIEDPrep.do      // 工企数据清洗
do code/01_data_prep/CTSDPrep.do      // 税调数据清洗
do code/01_data_prep/MAPrep.do        // Market Access计算
do code/01_data_prep/MdEst.do         // Markdown估计
```

### 2. 实证分析

```stata
do code/02_empirical/reg_land.do      // 土地价格回归
do code/02_empirical/reg_firm.do      // 企业回归
```

### 3. 结构模型

```julia
using Pkg; Pkg.activate(".")
include("code/03_model/Main.jl")
```

## 数据说明

| 数据源 | 时间范围 | 说明 |
|--------|----------|------|
| CIED (工企数据库) | 1998-2014 | 规模以上工业企业 |
| CTSD (税调数据) | 2007-2020 | 全国税收调查企业 |
| 土地交易数据 | 2007-2019 | 青岛市土地出让 |
| 路网数据 | 2011前后 | OpenStreetMap + 官方数据 |

## 关键变量

- **Market Access (MA)**: 基于通勤时间计算的市场可达性指标
- **Bridge Effect**: 大桥通车前后MA变化量
- **Markdown**: 企业成本加成率

## 作者

[姓名], [机构]

## 许可证

本项目仅供学术研究使用，数据使用需遵守相关保密协议。

---

*Last updated: 2026-03-30*
