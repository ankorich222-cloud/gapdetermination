
# RAW体数据空隙分析功能
当前SurfaceDeterminationProject目录下游以下代码
VolumeData.h 负责体数据加载，
SurfaceDetermination.h 负责灰度阈值辅助估计，
MarchingCubes.h 做初始等值面提取，
main.cpp 实现简单的人机交互流程。

下面是整技术流程，再按模块逐步补齐。
总体业务流程：

1. 体数据读取与预处理
目标：把输入CT/工业扫描体数据整理成后续算法可直接使用的统一格式。
具体步骤：
读取 .raw 数据，获得 dims、spacing、origin、voxels。
统计灰度范围 minVal/maxVal，用于阈值映射、直方图显示和后续参数约束。
可选预处理（后续补充开发）：
降噪：高斯、均值或中值滤波，降低散斑噪声。
灰度归一化：仅用于显示，不建议破坏原始灰度。
ROI 裁剪：如果用户只关心局部区域，可减少计算量。
建立体素坐标系与物理坐标系转换函数：
体素坐标 (i,j,k)
物理坐标 (x,y,z) = origin + index * spacing
后续所有距离类参数都必须支持两种单位：
voxel
mm
这是后面“法向搜索距离”和“空隙体积/等效直径”能做准确计算的前提。
2. 表面测定与初始三维表面模型提取
目标：先得到目标物体的闭合初始表面。
具体步骤：
计算全局灰度直方图。
自动估计背景峰、材料峰，给出初始 isoValue。
在二维切片中叠加阈值轮廓，允许用户人工微调：
background
material
isoValue
使用 Marching Cubes 提取初始等值面，得到体数据后续需要根据得到的处理

3. 高级模式：法向搜索距离优化边界
目标：让初始等值面更贴近真实边界，而不是仅依赖单一全局阈值。
3.1 输入参数
高级模式增加参数：
normalSearchDistance
distanceUnit：voxel 或 mm
可选：
searchStep
gradientThreshold
maxVertexShift
normalSmoothIterations
3.2 核心思想
初始表面出来后，对网格上每个顶点沿法向方向取样，寻找更合理的“真实边缘位置”。
3.3 算法流程
对每个边缘线的点执行：
取顶点位置 p 和单位法向 n。
将搜索距离统一换算到物理长度或体素长度。
沿 [-d, +d] 范围采样灰度剖面：
p + t * n
t 按固定步长递增
采样时使用三线性插值读取体数据灰度，不能直接最近邻取值。
在剖面上寻找最优边界点，方法：
先找灰度穿越，再在穿越附近找最大梯度点
若找到有效边界点，则把顶点沿法向移动到新位置。
若没有找到可信点，则保留原位置。
对全部顶点更新后，再做一次轻微几何平滑，避免局部尖刺。
3.4得到高级模式后的边缘表面，但其实后面空隙操作是对表面包含的体数据进行空隙分析
4. 根据表面模型限定“内部分析区域”
目标：后续空隙分割只在材料内部进行，避免把外部空气误识别为空隙。
这是整个“空隙提取”流程里很关键的一步。
做法：
Mesh → Voxel Rasterization

把三角面转为体掩膜

VTK直接支持：

vtkPolyDataToImageStencil
vtkImageStencil

流程：

surface (polydata)
      ↓
vtkPolyDataToImageStencil
      ↓
vtkImageStencil
      ↓
insideMask

这一步本质是：

Mesh → Binary Volume
从优化后的闭合表面生成一个“材料内部掩膜”。

5. 空隙候选区域生成
目标：从物体内部挑出可能属于空隙的低灰度体素。
5.1 输入参数
voidGrayMin
voidGrayMax
可选：
seedStrategy
minVoidVolume
connectivity：6/18/26
5.2 基本逻辑
对所有体素扫描：
如果体素在 insideObject 内部；
且灰度落在用户设定的空隙灰度范围内；
则标记为“空隙候选体素”。
得到：
candidateVoidMask
5.3 可选增强
为了提高稳定性，可以在候选掩膜上增加预处理：
去除孤立噪声点
小区域开运算
小孔闭运算
边界附近置信度过滤
第一版先不做太多形态学，避免引入额外误差。

6. 种子点选择
目标：从候选空隙区域中得到可靠的起始点。
这里可以设计成两种模式。
6.1 自动种子模式
系统自动在 candidateVoidMask 中寻找连通块内部点作为种子。
流程：
找所有候选体素。
对连续低灰度区域做初步聚类。
在每个区域内部选择一个或多个中心体素作为种子。
6.2 手动/半自动种子模式
用户在二维切片上点击一个或多个点作为种子，系统检查：
是否在空隙灰度范围内
是否靠近边界
6.3 建议
第一版先做“自动种子”，因为最终目标是“全部空隙提取”，不是只提单个空隙。
实际上，如果 candidateVoidMask 已经定义得足够好，后面直接做三维连通域分析就够了，种子点更多是为了交互和调试。
7. 三维连通域分析与空隙分割
目标：把所有候选空隙自动拆分成一个个独立缺陷。
7.1 输入
candidateVoidMask
insideObjectMask
connectivity
7.2 流程
在 candidateVoidMask 上执行 3D 连通域标记。
得到每个独立连通块的标签 labelId。
对每个连通块统计基础信息：
体素数
包围盒
中心位置
删除不满足条件的小连通块：
体素数太少
体积太小
明显为噪声
保留剩余连通块作为有效空隙对象。
7.3 连通性建议
第一版用 26-连通，更适合空隙这种不规则结构。
以后可以开放成参数。
7.4 特别注意
如果某个低灰度区域和物体外部通过表面漏洞相连，它不应当被算作内部空隙。
所以表面闭合性检查非常重要。
8. 空隙几何参数计算
目标：对每个识别出的空隙做定量分析。
建议至少输出以下参数。
8.1 位置类
质心坐标 centroid（可以是小数形式）
包围盒中心
在切片中的投影位置
相对物体坐标位置
距离外表面的最近距离
8.2 尺寸类
体积 V = voxelCount * spacingX * spacingY * spacingZ
等效直径
定义为与空隙体积等体积球体的直径：
Deq = (6V/pi)^(1/3)
三向尺寸：
dx
dy
dz
最大包围盒对角线
8.3 形状类
表面积
球形度
长宽高比
PCA 主轴长度
扁平度/细长度
紧致度
8.4 拓扑类
连通域编号
体素数
是否触碰边界
是否靠近外表面
8.5 推荐第一版必须实现
空隙编号
质心位置
体积
等效直径
包围盒尺寸
体素数
这部分最先有业务价值。
9. 二维与三维高亮显示
目标：把识别结果稳定、清晰地显示出来。
9.1 二维切片显示
在 XY/XZ/YZ 切片上叠加：
空隙区域 mask
空隙轮廓
空隙编号
当前选中空隙的参数信息
建议颜色区分：
表面轮廓：黄色
空隙区域：红色/橙色半透明
当前选中空隙：亮绿色或青色
9.2 三维显示
三维视图中同时显示：
物体外表面：半透明灰色
空隙表面：红色实体或高亮色
可选显示空隙质心点
可选显示编号标签
9.3 三维空隙表面生成
对每个连通域可单独再跑一次 Marching Cubes：
输入是单个空隙的二值 mask
输出是空隙网格
这样三维高亮会比直接画体素块更好看。
10. 结果导出
目标：支持后续分析和复现。
建议导出：
surface.stl：外表面模型
void_xxx.stl：每个空隙表面
voids.csv：空隙参数表
void_mask.raw/mhd 或 label_volume.raw/mhd：标签体
配置文件 json：保存本次阈值、搜索距离、连接方式等参数
11. 建议的数据结构设计
为了后面代码不乱，建议尽快补几个结构体。
struct AdvancedSurfaceParams {
    bool enabled = false;
    float normalSearchDistance = 0.0f;
    bool useMillimeter = false;
    float searchStep = 0.5f;
    float maxVertexShift = 2.0f;
    float gradientThreshold = 0.0f;
};

struct VoidDetectionParams {
    float grayMin = 0.0f;
    float grayMax = 0.0f;
    int connectivity = 26;
    double minVolumeMM3 = 0.0;
};

struct VoidRegion {
    int id = 0;
    std::vector<size_t> voxelIndices;
    size_t voxelCount = 0;
    double volumeMM3 = 0.0;
    double equivalentDiameterMM = 0.0;
    std::array<double, 3> centroidMM = {0, 0, 0};
    std::array<int, 6> bbox = {0, 0, 0, 0, 0, 0};
};
