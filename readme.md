# B样条曲面拟合球

## 一、编程环境

### 1.1 开发环境

* Windows 10

* Microsoft Visual Studio Community 2019 版本 16.8.4

* C++ 11

  

### 1.2 第三方库

* glm：c++数学库，用于存储向量和矩阵、计算矩阵变换，并将其传输至glsl
* eigen：c++数学库，用于求解线性方程组
* glad：用于访问OpenGL规范接口，加载驱动中的OpenGL实现函数
* glfw：一个轻量级的跨平台OpenGL窗口框架
* ImGUI：即时渲染GUI库，用于交互



## 二、使用方法

<img src="C:\Users\supuo\AppData\Roaming\Typora\typora-user-images\image-20210124153059792.png" alt="image-20210124153059792" style="zoom: 50%;" />

本程序可以直接使用GUI进行交互，可以选择b样条曲线拟合圆或b样条曲面拟合球，提供`uniform`、`chordal`、`centripetal`三种参数化方式，`uniform`、`average`两种节点向量生成方式，`open`、`clamped`、`closed`三种b样条类型，可以自选插入的数据点和b样条阶数。

- 按`z`切换显示拟合出的几何图形

- 按`x`切换显示数据点

- 按`c`切换显示控制点和控制多边形

- 按`v`键后可以用`wasd`、`空格`、`左shift`进行前后左右上下移动，可移动鼠标实现视角改变，滚动滑轮进行缩放。

- 按`esc`退出

  


## 三、数据结构

### 3.1 b样条曲线

```c++
struct BSpline {
	enum class BSplineType {
		Open, Clamped, Closed
	};

	BSpline(int _p, const std::vector<double>& us, BSplineType type = BSplineType::Clamped);
	std::vector<double> computeCoefficients(double u) const;
	glm::vec3 operator()(const std::vector<glm::vec3>& controlPoints, double u) const;

	int n = 0, m = 0, p = 0;
	std::vector<double> us;
	BSplineType type;
};
```



### 3.2 b样条曲面

```c++
struct BSplineSurface {
	BSplineSurface(int _p,
	               int _q,
	               const std::vector<double>& us,
	               const std::vector<double>& vs,
	               BSpline::BSplineType utype = BSpline::BSplineType::Clamped,
	               BSpline::BSplineType vtype = BSpline::BSplineType::Clamped);
	glm::vec3 operator()(const std::vector<std::vector<glm::vec3>>& controlPoints, double u, double v) const;

	BSpline ubsp, vbsp;
};
```



### 3.3 拟合器

```c++
class Fitter {
public:
	enum class ParametrizationMethod {
		Uniform, Chordal, Centripetal
	};

	enum class KnotGenerationMethod {
		Uniform, Average
	};

	Fitter(ParametrizationMethod pType = ParametrizationMethod::Uniform,
	       KnotGenerationMethod kType = KnotGenerationMethod::Uniform):
		parametrizationType(pType),
		knotGenerationType(kType) {}

	// ParametrizationMethod
	static std::vector<double> uniformParametrization(const std::vector<glm::vec3>& dataPoints);
	static std::vector<double> chordalParametrization(const std::vector<glm::vec3>& dataPoints);
	static std::vector<double> centripetalParametrization(const std::vector<glm::vec3>& dataPoints);

	// KnotGenerationMethod
	std::vector<double> generateKnots(int p, std::vector<double>& ts, BSpline::BSplineType bspType) const;
	std::vector<double> generateUniformKnots(int p, std::vector<double>& ts, BSpline::BSplineType bspType) const;
	std::vector<double> generateAverageKnots(int p, std::vector<double>& ts, BSpline::BSplineType bspType) const;

	// Parametrization
	std::vector<double> curveParametrization(const std::vector<glm::vec3>& dataPoints) const;
	std::vector<std::vector<double>>
	surfaceParametrization(const std::vector<std::vector<glm::vec3>>& dataPoints) const;

	// Interpolation
	static std::vector<glm::vec3> computeControlPoints(const BSpline& bsp,
	                                                   const std::vector<glm::vec3>& dataPoints,
	                                                   const std::vector<double>& ts);
	BSpline interpolateCurve(int p,
	                         const std::vector<glm::vec3>& dataPoints,
	                         std::vector<glm::vec3>& controlPoints,
	                         BSpline::BSplineType type) const;
	BSplineSurface interpolateSurface(int p,
	                                  int q,
	                                  const std::vector<std::vector<glm::vec3>>& dataPoints,
	                                  std::vector<std::vector<glm::vec3>>& controlPoints,
	                                  BSpline::BSplineType utype,
	                                  BSpline::BSplineType vtype) const;

	ParametrizationMethod parametrizationType;
	KnotGenerationMethod knotGenerationType;
};

```



## 四、算法描述

### 4.1 参数化方法

#### 4.1.1 均匀参数化

有$n+1$个数据点参数均匀分布在$[0,1]$区间内，将其划分为$n$个子区间，每段子区间长度位$1/n$，划分点位份别为$0,1/n,2/n,3/n,...,(n-1)/n,1$。



#### 4.1.2 弦长参数化

为了使数据点参数的分布情况与数据点将曲线划分的情况相联系，利用每两个数据点之间的弦长除以总弦长以确定数据点参数。



#### 4.1.3 向心参数化

为了减小长弦的影响，增大短弦的影响，减少拟合曲线的凸出，将弦长求幂（一般取0.5）。



### 4.2 曲面参数化

对每列分别在u（行）方向进行曲线参数化，将所得的uNet在每行进行平均得到的是曲面在u方向上的参数；对每行分别在v（列）方向进行曲线参数化，将所得的vNet在每lie进行平均得到的是曲面在v方向上的参数。

![img](https://pages.mtu.edu/~shene/COURSES/cs3621/NOTES/INT-APP/FIG-SURF-INT-parameter-1.jpg)

```cpp
for (int j = 0; j <= n; ++j) {
		vector<vec3> columnData(m + 1);
		for (int i = 0; i <= m; ++i) {
			columnData[i] = dataPoints[i][j];
		}
		uNet.emplace_back(curveParametrization(columnData));
	}
for (int i = 0; i <= n; ++i) {
    for (int j = 0; j <= m; ++j) {
        parameter[0][j] += uNet[i][j] / (n + 1);  // u方向参数
    }
}
for (int i = 0; i <= m; ++i) {
    vector<vec3> rowData(n + 1);
    for (int j = 0; j <= n; ++j) {
        rowData[j] = dataPoints[i][j];
    }
    vNet.emplace_back(curveParametrization(rowData));
}
for (int i = 0; i <= m; ++i) {
    for (int j = 0; j <= n; ++j) {
        parameter[1][j] += vNet[i][j] / (m + 1); // v方向参数
    }
}
```



### 4.3 节点向量生成

#### 4.3.1 均匀方法

使参数定义域所在的节点区间均匀分布在$[0,1]$区间内。



#### 4.3.2 平均方法

利用将参数平均的方法求节点的位置，可以解决均匀方法在弦长参数化插值时出现的线性方程系统奇异的问题，但仅适用于clamped B样条。

<img src="https://pages.mtu.edu/~shene/COURSES/cs3621/NOTES/INT-APP/EQN-knot-generation-2.jpg" alt="img" style="zoom: 80%;" />



### 4.4 计算混合函数（Ni,p）

<img src="https://pages.mtu.edu/~shene/COURSES/cs3621/NOTES/spline/B-spline/bs-basis.jpg" alt="img" style="zoom:80%;" />

设置长度位$p+1$的数组，在计算时需找到$u$所在区间$i$，令`us[i]= 1`，之后利用动态规划的方法和三角形递推形式得到$N_{i,p}(u)$

<img src="https://pages.mtu.edu/~shene/COURSES/cs3621/NOTES/spline/B-spline/bspline-triangle.jpg" alt="img" style="zoom: 67%;" />



### 4.5 曲线拟合

#### 4.5.1 插值

利用![img](https://pages.mtu.edu/~shene/COURSES/cs3621/NOTES/INT-APP/bspline-eqn.jpg)公式求出控制点$P_i$，可列矩阵方程![img](https://pages.mtu.edu/~shene/COURSES/cs3621/NOTES/INT-APP/GLOBAL-INT-matrix-3.jpg)。其中

<img src="https://pages.mtu.edu/~shene/COURSES/cs3621/NOTES/INT-APP/GLOBAL-INT-matrix-1.jpg" alt="img" style="zoom:80%;" />，<img src="https://pages.mtu.edu/~shene/COURSES/cs3621/NOTES/INT-APP/GLOBAL-INT-matrix-2.jpg" alt="img" style="zoom:80%;" />

对于open（不经过首尾控制点）、clamped（经过首尾控制点）和closed（曲线处处$C^{k-1}$连续）三种样条曲线有两种情况

- clamped情况：

  对参数t进行参数化时得到$t$的区间为$[0,1]$，由于参数$t$定义域要求在$[u_p,u_{n+1]}]$内，而clamped曲线的节点向量为
  $$
  \left\{u_0,u_1,...u_p,u_{p+1},u_{p+2},...,u_{n},u_{n+1},...,u_{m-1},u_m\right\}=\left\{0,0,...0,u_{p+1},u_{p+2},...u_{n},1,...1,1\right\}
  $$
  因此可以直接利用上述矩阵方程求解矩阵$P$即为控制点坐标矩阵

- open和closed情况：

  由于参数化后$t$的区间为$[0,1]$，且要求$t$的定义域在$[u_p,u_{n+1]}]$内，因此需要将t的区间转变为$[u_p,u_{n+1]}]$。首先将$t$重新参数化

  ```cpp
  inline void reparameterize(std::vector<double>& ts, const std::vector<double>& us) {
  	int n = static_cast<int>(ts.size()) - 1, m = static_cast<int>(us.size()) - 1, p = m - n - 1;
  	for (auto& t : ts) {
  		t = reparameterize(t, us[n + 1] - us[p], us[p]);
  	}
  }
  inline double reparameterize(const double t, const double factor, const double offset) {
  	return offset + factor * t;
  }
  ```

  随后对于open情况直接求解上述矩阵方程即可。而对于closed情况，由于需要处处满足$C^{k-1}$连续，因此还需要折叠控制点，使
  $$
  p_0 = p_{n-p+1}, p_1 = p_{n-p+2},...,p_{p-2}=p_{n-1},p_{p-1}=p_n
  $$
  因此矩阵方程$D=N\cdot P=N\cdot M\cdot P'$，其中
  
  
  $$
  M=\begin{bmatrix}
  1&0&\ldots&0\\
  0&1&\ldots&0\\
  \vdots&\vdots&&\vdots\\
  0&0&\ldots&1\\
  1&0&\ldots&0\\
  0&1&\ldots&0\\
  \ldots&\ldots&\ldots&\ldots
  \end{bmatrix},
  P'=\begin{bmatrix}
  p_0\\
  p_1\\
  \vdots\\
  p_{n-p}\\
  \end{bmatrix}
  $$
  由由于多项式函数不能完美逼近圆，因此对圆插值的求解过程中会出现约束条件过多矩阵方程无解情况，此时Eigen库会自动解出尽量满足矩阵方程的解。
  
  

### 4.6 曲面拟合

先对曲面的每列进行曲线拟合，得到每列（d）的控制节点$Q_{i,d}$ ,将其作为中间结果数据点组成网格，利用$Q$的每行（c）生成行控制点$P_{c,j}$ ，此为b样条曲面的控制点。

```cpp
	for (int d = 0; d <= n; ++d) {
		vector<vec3> columnDataPoints(m + 1);
		for (int i = 0; i <= m; ++i) {
			columnDataPoints[i] = dataPoints[i][d];
		}
		vector<vec3> intermediate = computeControlPoints(bs.ubsp, columnDataPoints, ts[0]);
		for (int i = 0; i <= m; ++i) {
			Q[i][d] = intermediate[i];
		}
	}
	controlPoints.resize(m + 1);
	for (int c = 0; c <= m; ++c) {
		controlPoints[c] = computeControlPoints(bs.vbsp, Q[c], ts[1]);
	}
```



### 4.7 计算b样条曲线位置点(利用定义)

![img](https://pages.mtu.edu/~shene/COURSES/cs3621/NOTES/spline/B-spline/bspline-curve-eqn.jpg)



### 4.8 计算b样条曲面位置点

将曲面的每行作为b样条曲线，利用曲面每行的控制点和knotsv求出参数v所在的位置，将每行的v位置点作为新的控制点，利用knotsu求出参数u所在的位置即为b样条曲面在参数（u，v）的位置

```cpp
int m = static_cast<int>(controlPoints.size()) - 1, n = static_cast<int>(controlPoints[0].size()) - 1;
vector<vec3> intermediate(m + 1);
for (int i = 0; i <= m; ++i) {
    vector<vec3> rowControlPoints(n + 1);
    for (int j = 0; j <= n; ++j) {
        rowControlPoints[j] = controlPoints[i][j];
    }
    intermediate[i] = vbsp(rowControlPoints, v);
}
return ubsp(intermediate, u);
```



## 五、实验结果

- 由于对球面利用极坐标等角度间隔采样,所以利用三种参数化形式所得误差几乎没有区别.

- 其他条件相同情况下，插值点数越多，所得几何图形与原图形误差越小

- 在一定范围内，阶数越高，所得几何图形与原图形误差越小，结束过高会产生大幅波动（龙格现象）。由于绘制点间距为0.01，因此可能部分拟合失效情况不能直观显示，但控制面板中有对每维度均匀随机采样100个数据点的误差可显示该情况。

- 对于不同B样条类型，open在阶数变高情况下最早失效，clamped情况次之，closed类型b样条在每维度100数据点，双向20阶情况下仍能保证非常完美的拟合，效果最好。

  <img src="C:\Users\supuo\AppData\Roaming\Typora\typora-user-images\image-20210124153835130.png" alt="image-20210124153835130" style="zoom: 50%;" />

  <img src="C:\Users\supuo\AppData\Roaming\Typora\typora-user-images\image-20210124153900163.png" alt="image-20210124153900163" style="zoom: 50%;" />

  <img src="C:\Users\supuo\AppData\Roaming\Typora\typora-user-images\image-20210124153924398.png" alt="image-20210124153924398" style="zoom: 50%;" />
  
  



## 六、声明

### 6.1 参考资料

- 理论和公式：https://pages.mtu.edu/~shene/COURSES/cs3621/NOTES/
- 程序组织结构：https://github.com/shihaoL/bspline



### 6.2 创新点

- 完成了Open、Clamped、Closed的三种形式的B样条
- 使用GUI交互