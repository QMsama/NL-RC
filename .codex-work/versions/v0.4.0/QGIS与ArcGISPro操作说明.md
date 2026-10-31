# 2024年TN关键源区识别与GIS操作

## 本版输入与方法

采用v0.3.1的water_2024_point21全季RM/RCC栅格，计算期6月1日至9月23日。按两模式实际种植掩膜分别识别田间系统出口TN流失关键源区候选。

文献框架：Zhuang et al. (2016), Identification of critical source areas for nonpoint source pollution in the Danjiangkou Reservoir Basin, China, DOI [10.1080/10402381.2016.1204396](https://www.tandfonline.com/doi/full/10.1080/10402381.2016.1204396)：对负荷—面积曲线拟合，依据导数1及2定位分级界限。本版使用端点固定为(0,0)/(1,1)、导数非负且递减的三次B样条最小二乘拟合；样条函数是本版的实现选择。归一化面积、负荷都用0—1或同时用0—100%，斜率不变。

先按kg/ha降序排序，横轴为累计面积比例，纵轴为累计质量比例。等面积像元的面积比例等于像元数比例，质量由kg/ha×像元面积计算，并与kg/像元文件核对。有效零值参与面积，NoData排除。主结果内部均匀结点数24，另加端部结点；具体结点和系数见thresholds_and_fit.json。在2001个等距面积样点及端部补充样点上拟合，单调凹约束使导数阈值对应连续的高负荷排名段。

解f'(x)=2和f'(x)=1，按累计面积比例换算像元排名，再提取该排名的负荷值作为GIS阈值。排名取最近整数，阈值同值像元整体纳入较高等级；因此实际面积和拟合面积可以有细微差异。实际负荷贡献重新由分类像元汇总，不能直接把拟合曲线值当作统计结果。

## 阈值与结果

|处理|1/2级界限 kg/ha|2/3级界限 kg/ha|关键区面积占比|关键区负荷占比|
|---|---:|---:|---:|---:|
|RM|4.893727779388428|9.774025917053223|6.8180%|15.4907%|
|RCC|5.467080116271973|10.915754318237305|4.0483%|8.6672%|

类别1为其他有效种植区域（L<T1），2为次关键源区（T1≤L<T2），3为关键源区候选（L≥T2）。0为分类文件NoData。两个模式独立归一化，类别表示模式内相对贡献。

RM关键区面积38825.77 ha、负荷430.98 t；RCC关键区面积13708.57 ha、负荷160.18 t。完整三级统计见class_statistics.csv。阈值来自本次拟合求导，不直接用2倍平均负荷代替。

## QGIS：直接查看

1. 加载正式成果目录“模型/2026.9.3/面上推广/2024/point21_关键源区”中的两张TN_source_class GeoTIFF。
2. 图层属性→符号系统→调色板/唯一值，读取类别。1浅灰蓝、2橙色、3红色；0透明。GeoTIFF已内嵌颜色表。
3. 加载县界作为透明填充、细线边界参照。展示地图时可以动态重投影；分类统计保持原始网格。
4. 图例明确写出“全季TN关键源区候选，模式内分级”。展示分类数据采用最近邻，保持类别整数。

## QGIS：从全季负荷栅格重现分类

1. 从“point21_三期累计TN流失”加载TN_RM_00_season_0601_0923_kg_ha.tif，将图层名改为RM；RCC对应图改为RCC。
2. 打开“栅格→栅格计算器”。RM表达式：

```text
1 + ("RM@1" >= 4.893727779388428) + ("RM@1" >= 9.774025917053223)
```

RCC表达式：

```text
1 + ("RCC@1" >= 5.467080116271973) + ("RCC@1" >= 10.915754318237305)
```

3. 输出范围、CRS、行列数都使用对应输入图层；保持原像元对齐。保存为新的GeoTIFF。
4. 设置唯一值颜色及图例。使用“栅格图层唯一值报告”检查有效分类数量：RM的1/2/3分别78774/37737/8525，RCC分别42262/29081/3010。原始负荷NoData不参加计算，不能改成0。

这两个表达式已在QGIS 3.44.12的栅格计算器引擎实际运行，分类和NoData掩膜与交付结果逐像元一致。参考：[QGIS栅格计算器](https://doc.qgis.org/3.44/en/docs/user_manual/working_with_raster/raster_analysis.html)。

## ArcGIS Pro：直接查看和重现分类

1. Map→Add Data加载两张分类GeoTIFF；Symbology选择Unique Values，按1/2/3设置图例和颜色。
2. 如需重新分类，加载全季kg/ha输入，分别将图层名改为RM和RCC。
3. 在Geoprocessing→Raster Calculator中计算（需相应Spatial Analyst许可）：

```text
Con("RM" >= 9.774025917053223, 3, Con("RM" >= 4.893727779388428, 2, 1))
```

```text
Con("RCC" >= 10.915754318237305, 3, Con("RCC" >= 5.467080116271973, 2, 1))
```

4. Environments设置：Snap Raster、Cell Size、Processing Extent、Output Coordinate System均以对应全季输入为准。条件栅格NoData应继续为NoData。
5. 对整数输出执行Build Raster Attribute Table，核对COUNT与上述QGIS像元数量一致。也可用Reclassify工具，但须明确边界等值归属，不能直接采用自动自然断点。

本次未启动ArcGIS Pro执行表达式，步骤依据Esri工具文档：[Reclassify](https://doc.esri.com/en/arcgis-pro/latest/tool-reference/spatial-analyst/reclassify.html)。直接加载已生成GeoTIFF即可查看结果。

## 如何汇总各等级负荷

分类图作为分区图，匹配的全季kg/像元图作为值栅格。ArcGIS Pro使用Zonal Statistics as Table（分区字段Value，统计SUM），得到各级TN质量kg，除1000为t。COUNT×4.554342663494959为面积ha。QGIS可用分区栅格统计工具；若先转面再做矢量分区统计，保留精确像元边界并核对统计差异。

不可将kg/ha像元直接求和当作kg；也不可将分类值1/2/3求和当作负荷。需要县域归属时另指定县界和跨界像元分摊规则，本版统计为模式内全区三级统计。

## 如何重跑拟合和重新得到阈值

GIS重分类工具负责把阈值映射到地图，样条拟合与求导使用随版脚本。在QGIS自带Python环境执行，需NumPy、SciPy、GDAL、Matplotlib；本机已实跑：

```powershell
& 'C:/Program Files/QGIS 3.44.12/bin/python-qgis-ltr.bat' -B '.codex-work/exports/identify_source_areas.py' --input '完整稻季栅格目录' --output '新的分类栅格目录' --records '新的核验目录'
```

从NL-RC仓库根目录执行。每次换年份、管理规则、研究区或输入掩膜，重新运行拟合，采用新thresholds_and_fit.json中的阈值；本表阈值只适用于当前point21全季结果。使用新输出目录保存各次版本。

## 核验与解释

有效像元、投影、行列、像元对齐、分区面积和质量闭合、GeoTIFF读回均通过。拟合预测与实际空间累计负荷差最大0.002632个百分点。16/24/32均匀内部结点配置下，主结果与对照分类一致率99.8480%—99.9360%（详细各项见JSON）；曲线拟合R²接近1仅表示对已有累计曲线的逼近。

识别依据仍为2024气象、2021土地分类及2018施氮代理的情景输出。空间出口负荷候选区与实际入河贡献的对应，需要另外核验河网连通和输移。
