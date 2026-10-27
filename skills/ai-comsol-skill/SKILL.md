---
name: ai-comsol-skill
description: Build and run auditable COMSOL 6.4 ethylene-glycol electrolysis research models through native Java and batch. Use for 0D reaction balances, 2D/3D planar-channel reaction scenarios, literature parameter provenance, optional DWSIM reference properties, conservation checks and saved-model readback. General manifold CFD, ionic migration and experimental calibration require additional models and evidence.
---

# AI-COMSOL-Skill

面向乙二醇电解研究的原生 COMSOL 自动化技能。先确认问题、参数依据和物理边界，再生成并运行模型。已有实现使用 Java API 与 batch，不要求 MPh、MATLAB 或 MCP。

## 按问题选择入口

- **0D物料/电荷衡算**：读[反应模型](references/reaction-model.md)，使用 `scripts/electrosynthesis.py`。支持共用电解液或双CSTR、给定支路FE或前向Tafel竞争支路。独立RK4与电流代数根用作核对。
- **2D/3D直流道情景**：读[空间模型](references/spatial-models.md)、[参数依据](references/catalyst-reactor-data-v2.json)与[校核契约](references/spatial-model-contract.json)，使用 `scripts/spatial_reactor.py`。默认求充分发展Stokes–Poiseuille速度PDE，再耦合有机物输运与二次电流分布；3D含有限宽度侧壁。
- **物性参考或拟合准备**：读[数据与验证](references/data-and-validation.md)和[0D参数表](references/parameter-sources.json)。`scripts/run_dwsim_properties.py`可选调用用户安装的DWSIM计算中性水/EG参考物性。
- **新几何、三级电流、多孔电极、传热或气泡**：读[扩展路线](references/spatial-model-routing.md)。先核接口许可与输入，再建立新模型和独立验证。仓库中的 `navier-stokes` 选项为实验模式，未形成已验证的一般流道CFD能力。

## 项目与运行

技能安装目录与用户项目根目录分开。每次明确当前项目；`--project-root`省略时取当前工作目录。公开示例没有私有项目绑定。自有配置若声明`project_root`，必须与当前选择相同，否则拒绝。

用户新生成的模型/配置快照/结果进入所选项目的`AGENTCODING/07_handoff`；Java、类、临时配置进入`04_work`；日志进入`05_qa`。使用唯一运行目录；拒绝覆盖既有输出。原始资料和既有模型只读，修改时另建副本。不得关闭COMSOL安全策略来绕过导出限制。

从仓库根目录运行：

```powershell
$s = '.\skills\ai-comsol-skill'
$env:COMSOL_ROOT = 'C:\Program Files\COMSOL\COMSOL64\Multiphysics' # 改成实际安装路径，仅本进程
py -3.12 -B "$s\scripts\electrosynthesis.py" validate --config "$s\assets\literature_eg_equivalent.json"
py -3.12 -B "$s\scripts\electrosynthesis.py" native --config "$s\assets\literature_eg_equivalent.json" --out '.\AGENTCODING\07_handoff\eg-0d-01'
py -3.12 -B "$s\scripts\spatial_reactor.py" native --config "$s\assets\spatial_pdnimoo4.json" --dimension 2 --out '.\AGENTCODING\07_handoff\eg-2d-01'
py -3.12 -B "$s\scripts\spatial_reactor.py" native --config "$s\assets\spatial_pdnimoo4.json" --dimension 3 --out '.\AGENTCODING\07_handoff\eg-3d-01'
```

技能装在别处时，`$s`改为技能的完整路径，在当前用户项目根目录运行，或显式传`--project-root`及项目内输出路径。程序使用Python标准库；原生运行需要Windows与匹配的COMSOL安装/模块许可。`COMSOL_ROOT`也可用`--comsol-root`指定。

## 参数与电化学约束

分开标识文献测量、软件计算、工程假设、合成测试值和待补输入。缺少实验数据时，可执行用户授权范围内的文献情景与有依据假设；每项假设保留单位、依据、范围、影响和替换条件。不同催化剂的最佳FE、最低电位和最高电流不能拼成一组测量数据。

- 阳极EG→乙醇酸根、两个甲酸根、草酸根分别按4/6/8电子计量；先做元素与电荷检查。
- 阴阳极共用总电流，FE分别在各侧归一；不能把两侧FE加起来要求≤1。未测支路的分配必须显式列为假设。
- 净新鲜进料与内部循环分开。参比电极电位、过电位和整槽电压分开。
- 原论文未补偿iR的极化数据只能作为表观代理；当前空间阳极边界可能重复计入参考欧姆损失，不宣称本征动力学辨识。
- 空间动力学参考EG浓度固定于文献锚点，改变入口浓度不得跟随改变参考浓度。温度须与已声明动力学参考匹配；Re和入口发展长度检查不通过时不得沿用充分发展模型。
- 空间模型以恒定支持电解质电导率近似，GLY/FMA/OX按被动产物标量输运；没有空间OH-/Na+守恒或离子迁移。高局部产物浓度下必须重评该近似。
- DWSIM中性参考不提供电解液电导率、离子活度、反应速率或FE。拒绝非有限、零或负物性；某属性失败不伪造替代值。
- 盐形态出口不是已分离游离乙醇酸。槽压、热负荷、纯度和工业放大缺少相应基准时保持未知。

## 验证与交付

返回码0与生成mph均不足以证明成功。检查原生输出、有限值、浓度极值、物料/元素/电流闭合及保存重开一致。0D还要和独立积分比较；空间模型另做解析流动核对和同物理输入的网格加密。少量数值欠冲保留原值并按显式容限判断，不裁剪。

保存回读按每个数值自身量级比较，空间模型包括完整导出场。网格比较重查有效状态及文件哈希，必须加密实际使用的方向。解析或网格检查不通过的退出码为2；历史结果不能自动获得新契约的通过状态。支持电解质近似警示必须随数值通过结果一并交付。

输出含`config_used.json`、`result_card.json`、原生表格与mph；0D另有历史与结果说明，空间另有场数据和流动模型。完整模型可重开，不依赖临时Java源文件。结果卡区分软件检查、实验校准与工程验收；仅前者可能由程序确认。

软件检查：`py -3.12 -B -m unittest discover -s skills/ai-comsol-skill/scripts -p 'test_*.py' -v`。发布版本实际验证范围见仓库`docs/verification.md`。不得把历史运行或其他催化剂论文当作新配置的运行证据。
