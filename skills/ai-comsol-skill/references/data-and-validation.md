# 数据准备与结果验证

参数保留来源、原始单位、SI值、催化剂/电解液/电位参考、适用范围和性质。文献测量、软件计算与工程假设不能互换；图上读取要注明方法和读数精度。关系量、区间及文献性能汇总不能冒充逐时实验点。

## 0D示例

`literature_eg_equivalent.json`以[Han等2025论文](https://www.nature.com/articles/s41467-025-58813-3)为等效电荷情景依据，使用0.10M EG、1M KOH及装置库存量级。由0.475mmol/(cm²·h)、1cm²、名义FE0.92推得平均等效电流约0.05535A。它不是脉冲电流轨迹或失活再生模型。500mL是库存；8mL/min是内部循环，不能设置为连续新鲜进料。

剩余阳极电流分配给OER、阴极FE=1均为情景假设。`demo_eg_her.json`和`demo_kinetic.json`使用合成动力学，只用于软件验证。`project_input_template.json`故意未填完整，必须被输入检查拒绝。

## DWSIM参考

可选原生适配器使用`CalculatorInterface.Calculator.CalcProp`与Raoult's Law，计算中性水/EG的参考物性。固定状态为298.15/313.15/333.15K、101325Pa绝压、EG摩尔分数0/0.0018/0.02/1。

密度、黏度和导热系数仅为中性参考，不能直接代表KOH/NaOH+EG电解液。在已测接口中液体Cp返回0；适配器逐属性拒绝该值并保留null。不得用未验证的纯组分Cp代替混合电解液Cp。

```powershell
py -3.12 -B .\skills\ai-comsol-skill\scripts\run_dwsim_properties.py --dwsim-root '<实际DWSIM安装目录>' --out '.\AGENTCODING\07_handoff\dwsim-reference-01'
```

## 验证层次

1. 无许可依赖的测试：计量/电子数、FE、Faraday产量、解析CSTR极限、耗尽拒绝、异常输入、输出边界与安装位置独立性。
2. 0D原生求解：与独立RK4和代数电流根比较，保留结果表及保存重开差异。显式设置实际Time求解器容差，不能只设置Study显示值。
3. 空间求解：流量、两极电流、物种总通量与碳闭合，浓度极值、解析流动和网格加密，详见`spatial-model-contract.json`。
4. 实验校准：需要同一电极/电解液的极化、产品谱、转化与动态响应。当前公开示例未完成这一层。

累计转化以初始库存加累计进料为基准，不是CSTR稳态单程转化；净生成碳选择性避免重复累计中间体，存在其他含碳进料时不自动给出EG基准值。

公开版本实际运行记录的脱敏摘要见仓库`docs/verification.md`。完整本地模型/日志留在用户项目的AGENTCODING，不随开源包分发。

## 官方方法参考

- [Time求解器容差](https://doc.comsol.com/6.4/doc/com.comsol.help.comsol/comsol_api_solver.51.51.html)
- [Windows运行参数](https://doc.comsol.com/6.4/doc/com.comsol.help.comsol/comsol_ref_running.38.31.html)
- [临时/配置目录](https://www.comsol.com/support/knowledgebase/1083)
