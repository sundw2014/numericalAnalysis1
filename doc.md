# 数值分析大作业一
> 孙大伟 2014010782

## 需求分析
 所需功能：输入一张图片，对其进行变换，最后输出变换后的图像
 插值算法：可以选择最近邻、双线性、双三次插值方法中的一种
 变换：可以选择扭曲、畸变、TPS变形中的一种

## 方案设计
#### 符号说明
 $I_i$表示输入图像，$I_O$表示输出图像，$f(x,y)$表示变换所对应的映射，将输出图像中的坐标映射到输入图像上。
#### 基本思路
 算法的目的是求输出图像中各像素点的值$I_O(x,y)$，其中$x,y$是两个整数值，即像素的索引。
 对于坐标为$(x,y)$的像素处理如下：
  1. 得到$I_I$中对应的坐标$(x_I,y_I) = f(x,y)$
  2. 通过插值函数和原图像得到$\tilde I_I(x_I,y_I)$
  3. 令$I_O(x,y)=\tilde I_I(x_I,y_I)$
## 方案具体原理
### 插值算法
插值算法输入的是坐标$(x_I,y_I)$，需要根据原图像$I_I$中的某些像素点得到插值结果$\tilde I_I(x_I,y_I)$
#### 最近邻插值
 对于坐标$(x_I,y_I)$，与其最为接近的整数点为$(round(x_I),round(y_I))$,其中$round(x)$是对$x$求近似（四舍五入），所以得到 $$\tilde I_I(x_I,y_I) = I_I(round(x_I),round(y_I))$$
#### 双线性插值
 双线性插值算法可以用矩阵运算来描述
 $$ \tilde I_I(i+u,j+v) =
 \left[ \begin{matrix}
 1-u & u
 \end{matrix}\right]
 \left[ \begin{matrix}
 I_I(i,j) & I_I(i,j+1)\\
 I_I(i+1,j) & I_I(i+1,j+1)
 \end{matrix}\right]
 \left[\begin{matrix}
 1-v \\
 v
 \end{matrix}\right]
  $$
#### 双三次插值
 双三次插值算法可以用矩阵运算来描述
 $$\tilde I_I(i+u,j+v) = ABC^T$$
 其中
 $$A = \left[ \begin{matrix}
 S(u+1) & S(u) & S(u-1) & S(u-2)
  \end{matrix}\right]\\
  C = \left[ \begin{matrix}
 S(v+1) & S(v) & S(v-1) & S(v-2)
  \end{matrix}\right]\\
 B=I_I(i-1:i+2,j-1:j+2)$$
 $S(x)$为三次插值的核函数，可以按照如下公式近似
 $$S(x) =
 \begin{cases}
	 (a+2)\left|x\right|^3 - (a+3)\left|x\right|^2 + 1 & \left|x\right|\leq1\\
	 a\left|x\right|^3 - 5a\left|x\right|^2 + 8a\left|x\right| - 4a & 1<\left|x\right|<2\\
	 0 & \mbox{otherwise}
 \end{cases}$$
 其中，根据经验$a$取$-0.5\sim-0.7$，本次试验中取$-0.6$
### 变换算法
###### 符号说明
 $(x_O,y_O)$是$I_O$中的直角坐标，$(x_I,y_I)$是$I_I$中的直角坐标；$(r_O,\alpha_O)$是$I_O$中的极坐标，$(r_I,\alpha_I)$是$I_I$中的极坐标；
#### 扭曲
 扭曲变换的公式如下：
 $$x_O=r_I*{cos\left({\alpha_I+\theta*\frac{row-r_I}{row}}\right)}, y_O=r_I*{sin\left({\alpha_I+\theta*\frac{row-r_I}{row}}\right)}$$
其中$row$为最大旋转半径，$\theta$为最大旋转角度。
 从上述公式可知$\alpha_O = {\alpha_I+\theta*\frac{row-r_I}{row}}，r_O=r_I$。从而易得反变换的公式为
 $$\alpha_I =  \alpha_O - \theta*\frac{row-r_O}{row} ， r_I = r_O\\
 x_I = r_I * cos\alpha_I , y_I = r_I * sin\alpha_I$$
#### 畸变
 在本次试验中畸变的模型选为三个参数的径向畸变，描述如下
 $$r_I = r_O(1+k_1r_O^2+k_2r_O^4+k_3r_O^6),\alpha_i=\alpha_O$$
 所以得
 $$x_I=r_icos(\theta_I),y_i=r_Isin(\theta_I)$$
#### TPS
 用户通过点击鼠标输入两个点集，$\left\{(x_O^1,y_O^1)(x_O^2,y_O^2)(x_O^3,y_O^3)...(x_O^N,y_O^N)\right\}$为输出图像中的坐标点，$\left\{(x_I^1,y_I^1),(x_I^2,y_I^2),(x_I^3,y_I^3)...(x_I^N,y_I^N)\right\}$为输入图像中的坐标点。通过算法可以得到一个函数，从而将输出图像中的任意点$(x_O,y_O)$映射到输入图像中的点$(x_I,y_I)$，为了得到这样的映射，我们分别考虑$x$和$y$坐标。
 先考虑$x$坐标，我们要找到一个函数$f(x,y)$，满足:
 $$f(x_O^i,y_O^i) = x_I^i - x_O^i ，i=1,2,3...N$$
 为了使映射光滑，假设函数$f(x,y)$具有以下形式：
 $$f(x,y) = a_1 + a_xx + a_yy + \sum\limits_{i=1}^{N}\omega_iU(\left\|(x_O^i,y_O^i)-(x,y)\right\|)$$
 其中$U(x)=x^2log(x)$，为了使映射形成的曲面能量最小添加一些限制条件，最终得到：
  $$\begin{cases}
  f(x_O^i,y_O^i) = x_I^i - x_O^i ，i=1,2,3...N\\
  \sum\limits_{i=0}^{N}\omega_i = 0\\
  \sum\limits_{i=0}^{N} \omega_ix_i = 0\\
  \sum\limits_{i=0}^{N} \omega_iy_i = 0
  \end{cases}$$
  共$N+3$个方程，$N+3$个未知数，可以解出$f(x,y)$，假设得到函数$f_X(x,y)$，同理也可以通过相同方法得到$f_Y(x,y)$。
  从而可以根据以下关系得到输出图像中点在输入图像中的对应点：
  $$x_I = f_X(x_I,y_I) + x_O，y_I = f_Y(x_O,y_O) + y_O$$
## 误差分析

## 结果分析
 
