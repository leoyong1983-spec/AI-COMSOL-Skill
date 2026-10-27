# AI-COMSOL-Skill

[![License: MIT](https://img.shields.io/badge/License-MIT-blue.svg)](LICENSE)
[![Validation](https://github.com/leoyong1983-spec/AI-COMSOL-Skill/actions/workflows/validate.yml/badge.svg)](https://github.com/leoyong1983-spec/AI-COMSOL-Skill/actions/workflows/validate.yml)

让 AI 通过 **COMSOL 原生 Java API + batch** 建立、求解并核对乙二醇电解研究模型。提供可加载的 `SKILL.md`、可运行控制器、反应示例、催化剂/动力学来源和明确标注的参数假设。

Native COMSOL automation for ethylene-glycol electrolysis research: 0D balances, 2D/3D planar-channel scenarios, source provenance and saved-model readback. Windows-first; Python standard library; separately licensed COMSOL installation required.

## 已实现的范围

| 能力 | 实现与边界 |
|---|---|
| 0D反应 | 共用电解液或双CSTR；给定支路FE或前向Tafel；元素、电荷和两侧电流闭合；独立RK4核对 |
| 2D/3D反应 | 无膜直流道；充分发展层流速度PDE + EG/产物输运 + 二次电流分布；3D包含有限宽度侧壁 |
| 参数准备 | 文献催化剂/反应器/表观动力学，来源分级与假设区间；4/6/8电子反应支路 |
| 原生验证 | Java生成、求解、表格/场数据导出、保存mph、重开回读、解析流动与网格比较 |
| 可选DWSIM | 中性水/EG的参考物性与逐属性有效性检查 |

当前空间模型适用于声明假设下的研究情景。一般歧管/入口CFD、离子迁移、OH⁻局部耗竭、多孔电极内部、气泡、传热、脉冲失活再生和实验校准需要另建模型。`navier-stokes`选项仍为实验模式。

## 快速运行

环境：Windows、Python **3.12**、用户自有COMSOL安装和对应物理接口许可。当前验证版本及数值证据见[验证说明](docs/verification.md)。无需安装MPh、MATLAB或MCP；基础测试不需要COMSOL。

```powershell
git clone https://github.com/leoyong1983-spec/AI-COMSOL-Skill.git
Set-Location AI-COMSOL-Skill
py -3.12 -B tools/validate_repo.py
```

0D参考计算可直接运行，无需仿真器：

```powershell
$s = '.\skills\ai-comsol-skill'
py -3.12 -B "$s\scripts\electrosynthesis.py" reference --config "$s\assets\literature_eg_equivalent.json" --out '.\AGENTCODING\07_handoff\eg-reference-01'
```

原生0D、2D、3D：先将下例路径改成实际安装目录，仅设置当前进程环境。

```powershell
$env:COMSOL_ROOT = 'C:\Program Files\COMSOL\COMSOL64\Multiphysics'
py -3.12 -B "$s\scripts\electrosynthesis.py" native --config "$s\assets\literature_eg_equivalent.json" --out '.\AGENTCODING\07_handoff\eg-0d-01'
py -3.12 -B "$s\scripts\spatial_reactor.py" native --config "$s\assets\spatial_pdnimoo4.json" --dimension 2 --out '.\AGENTCODING\07_handoff\eg-2d-01'
py -3.12 -B "$s\scripts\spatial_reactor.py" native --config "$s\assets\spatial_pdnimoo4.json" --dimension 3 --out '.\AGENTCODING\07_handoff\eg-3d-01'
```

每次使用新运行目录，程序拒绝覆盖。模型、输入快照和结果在`AGENTCODING/07_handoff`，日志在`05_qa`，临时Java/配置在`04_work`。这些本地文件由Git忽略。

技能可放在支持`SKILL.md`的AI工具的技能目录；完整入口是[skills/ai-comsol-skill/SKILL.md](skills/ai-comsol-skill/SKILL.md)。安装位置与计算项目分开：在用户项目根目录运行控制器，或显式传`--project-root`。不要把输出指向原始资料或正式设计文件。

示例请求：

> 使用 $ai-comsol-skill，以公开文献为依据建立乙二醇制乙醇酸盐的二维流道情景，逐项标注假设，计算浓度与电流分布，并做守恒、保存重开和网格检查。

## 示例与参数

| 文件 | 用途 |
|---|---|
| `assets/literature_eg_equivalent.json` | 文献锚定0D等效电荷情景；不是脉冲动力学复现 |
| `assets/demo_eg_her.json` | EG氧化/HER合成动力学测试 |
| `assets/demo_kinetic.json` | 含可选糠醛阴极支路的合成测试 |
| `assets/spatial_pdnimoo4.json` | Pd/NiMoO₄/NF表观阳极代理与2D/3D平板流道假设 |
| `assets/project_input_template.json` | 故意留空的项目输入模板；未补齐时拒绝求解 |

以上路径相对`skills/ai-comsol-skill`。[0D反应说明](skills/ai-comsol-skill/references/reaction-model.md)、[空间模型与网格检查](skills/ai-comsol-skill/references/spatial-models.md)、[催化剂与参数来源](skills/ai-comsol-skill/references/catalyst-reactor-data-v2.json)、[数据准备与DWSIM边界](skills/ai-comsol-skill/references/data-and-validation.md)提供具体口径。

另保留[33项COMSOL自动化来源与案例线索](skills/ai-comsol-skill/references/automation-sources.json)，为2026-10-26的检索快照。第三方工具的声明与本仓库实际运行证据分开，未将其自动安装为依赖。

空间示例的通道间隙、HER参数、产品FE分配、部分物性和扩散系数为假设或相邻证据先验。原论文未补偿iR的阳极参数属于表观代理，不能称为本征动力学。数值收敛不等于实验验证，盐形态产物也不等于已分离的游离乙醇酸。

## 贡献与许可

修改反应式、关键参数、几何或物理接口时保留依据，并重新运行相应检查。提交前运行`py -3.12 -B tools/validate_repo.py`。GitHub Actions仅执行无需商业软件的检查；原生COMSOL验证须在自有许可环境执行。

原创代码和说明采用[MIT许可证](LICENSE)。软件和论文的权利边界见[第三方说明](THIRD_PARTY_NOTICES.md)。仓库不包含私有项目资料、论文全文、COMSOL/DWSIM二进制或许可文件。
