# 二维/三维反应模型：范围、方程与运行

本模式在本项目内生成原生COMSOL几何、网格、流动基准、物种场与电流分布，并保存、重开回读。它是有来源参数及假设下的研究情景，没有完成项目实验校准。参数性质逐项见[catalyst-reactor-data-v2.json](catalyst-reactor-data-v2.json)。

## 可运行几何与流动

- 主几何是无膜平板直流道。2D采用实际面外宽度W；3D包含有限宽度的两侧壁，解的浓度/电势可以沿x、y、z变化。
- 经过验证的默认模式是fully-developed-poiseuille：在原生Coefficient Form PDE内求充分发展、单向层流的横截面速度。以截面积分归一到给定体积流量；由黏度与速度解计算压降。2D可用解析抛物线、3D可用矩形管Fourier级数独立核对。
- 该模式明确假设有足够长的上游直流道。当前基准假设20mm，按Re与水力直径估算入口发展长度约8mm。它不解歧管、入口回流、复杂流场和电极孔道。
- navier-stokes保留为实验模式；早期模型的守恒检查未通过，不能把该模式列为已验证一般流道CFD。需要实际流道时，以受控几何另建模型并重新验证。
- 对沿流向和间隙方向的两条相对边采用一致分布，避免仅一侧加密后网格横向变化。3D使用映射面及扫掠六面体；侧壁方向也加密。对流保持SUPG，当前2D保留交风稳定化，3D在流向对齐网格上关闭非线性交风项；各例均另查极值和守恒。

## 反应与传输

速度先求解并冻结，再耦合Transport of Diluted Species与Secondary Current Distribution。EG为中性底物；GLY、FMA、OX是用于产品衡算的被动标量，不在本模式中声称已经描述其电迁移。1M NaOH的电导率和活度按支持电解质近似固定。

稳态底物/产物方程为u·grad(c_i)=div(D_i grad(c_i))；壁面通量由显式电子计量给出。阳极电流采用：

    j_a = j_anchor (c_EG/c_ref)^m × 10^[(E_a - phi_l - E_anchor)/b_a]

该式是表观局部边界假设。Pd/NiMoO4/NF基准取100mA/cm²@0.79VvsRHE及207.7mV/dec；原论文未补偿iR，不能把这些参数称为已去除传质/欧姆影响的本征动力学。m=0.3来自另一Pd电极的反应级数先验。阴极采用另行标注的HER经验假设，不能因此宣称真实槽压得到验证。

`kinetics.anode_EG_reference_mol_m3`固定为该锚点的1000 mol/m³，独立于`operation.feed_EG_mol_m3`。改变进料浓度时不得同步改动参考浓度，除非换用了有依据的新动力学锚点。程序逐点核对导出的电流密度与该公式。

GLY、两个FMA、OX分别按4、6、8电子计量生成。FE分配必须显式且和为1；当前0.8/0.15/0.05是工程产品谱假设。未在低于析氧热力学电位的基准中把所有剩余电流分配给OER。

没有空间解析OH-/Na+耗竭、盐活度、催化吸附覆盖度、失活、脉冲再生、孔内传输、气泡、析盐或热量；这些用途需新增相应物理与数据。温度基准25°C为室温近似，非任意温度动力学拟合。

输入检查拒绝超出动力学参考温度的直接外推，并以Re<2000、假设上游长度≥0.05ReDh筛查充分发展层流适用性。入口长度关联来自圆管层流的近似工程估算，用水力直径移用于矩形通道，不构成实际入口流场验证。结果另外报告`max(c_GLY+c_FMA+2c_OX)/c_support`；超过0.1时提示重评支持电解质近似。0.1为工程审查提示值，不是经过验证的误差上限。

## 运行

在已选项目根目录执行，输出目录每次唯一。先运行输入校验：

    $s = '.\skills\ai-comsol-skill'
    py -3.12 -B "$s\scripts\spatial_reactor.py" validate --config "$s\assets\spatial_pdnimoo4.json"
    py -3.12 -B "$s\scripts\spatial_reactor.py" native --config "$s\assets\spatial_pdnimoo4.json" --dimension 2 --out '.\AGENTCODING\07_handoff\COMSOL空间运行\eg-2d-new'
    py -3.12 -B "$s\scripts\spatial_reactor.py" native --config "$s\assets\spatial_pdnimoo4.json" --dimension 3 --mesh-factor 0.75 --out '.\AGENTCODING\07_handoff\COMSOL空间运行\eg-3d-new'

--mesh-factor改变各方向单元数；--flow-ml-min改变模型局部通流量，不代表整个循环系统的新鲜进料。需改单个参数时在AGENTCODING内建立带来源/假设的新配置，保留每次config_used.json。

结果包括原生reactor_spatial.mph、flow_solved.mph、原生场数据、积分结果、重开回读表、配置及result_card.json。打开最终mph不依赖临时Java或外部过程文件；重算需COMSOL安装及本包控制器。所有临时Java、类、软件配置进入04_work，日志进入05_qa。

## 检查

结果卡须满足有限值、EG正值/最大值、流量及两极电流闭合、4种物种总通量及碳衡算，并保存重开一致。少量产品数值负向欠冲保留原值与极值记录，按明确数值容限检查，不裁剪成零。阈值见[模型契约](spatial-model-contract.json)。只有result_card.json明确通过的运行才进入有效成果索引。

当前`spatial-2`结果还检查流量是否等于指令值、各支路是否满足法拉第计量，并逐列核验保存前后积分值和完整导出场：绝对容差1e-14（对应列SI单位）加该数值自身量级的1e-8相对容差。不得用统一的1作为微量产物流率的相对比较分母。

    py -3.12 -B "$s\scripts\spatial_review.py" analytic --case '<运行目录>' --out '<07_handoff内新JSON>'
    py -3.12 -B "$s\scripts\spatial_review.py" mesh --coarse '<粗网格目录>' --fine '<细网格目录>' --out '<07_handoff内新JSON>'

网格比较必须保持几何、操作、物性、动力学及维数相同；报告总电流、最小表面EG、最大局部电流、乙醇酸根生成及压降的差异，不能只看场图。相同模型网格收敛不替代实验验证。

比较前重查通过状态、输入/模型/场文件哈希及原始导出值，求解器设置和网格分布比保持一致。细网格必须至少加密一个实际使用的方向且不粗化其他方向；2D只改未使用的跨度网格无效。历史`spatial-1`结果需重跑才能按新契约复核。解析流动检查设定速度最大误差/均速≤1%、压降误差≤0.5%；网格关键指标差异≤2%。任一未满足时返回非零退出码2，不能被批处理当作通过。

## 官方方法参照

- [COMSOL6.4流动电池示例](https://doc.comsol.com/6.4/doc/com.comsol.help.models.battery.pb_flow_battery/pb_flow_battery.html)：支持先求稳态流场、再耦合电化学传输的工作方式；其化学体系为相邻范例，不能当作本项目反应参数。
- [电化学接口](https://doc.comsol.com/6.4/doc/com.comsol.help.echem/echem_introduction.02.02.html)：本模式是支持电解质近似，不把它改称三级电流分布或浓电解液机理验证。
- [入口和出口边界的放置](https://www.comsol.com/blogs/how-to-place-inlet-and-outlet-boundary-conditions-in-cfd-simulations)：层流入口长度0.05ReDh的来源；矩形流道应用属于近似筛查。
