# 射线与三角形求交：Möller–Trumbore 算法详解与证明

## 前言

这篇文章是我和 AI 一来一回聊出来的。

起因只是一个很小的问题：射线与三角形求交，除了「先求射线与平面的交点、再判断点是否在三角形内」这个老办法，Möller–Trumbore 那个一口气解出 $t,u,v$ 的做法究竟是怎么来的？结果问题越追越多——为什么这里冒出一个克莱姆法则？克莱姆法则本身为什么成立？混合积为什么能写成行列式？由于很多数学知识遗忘了，于是就这么一路问下去，把原本只会照抄公式的地方，逐个变成了能自己推出结论的地方。

收获比预期大。以前看这类实现，眼里只是一串叉乘和点乘；现在知道每一步在算什么、为什么可以那样重排，以及符号最容易在哪一步翻车。怕自己忘，也怕以后只记得结论不记得原因，就把这次讨论整理成了这篇。

全文按「问题 → 数学工具 → 推导 → 实现 → 工程细节」的顺序展开。想补推导的话，§4（克莱姆法则）和 §5（混合积）是关键；只想拿代码，直接跳到 §8 即可。

---

## 1. 问题定义与记号

三维空间里最常被问的一个问题：**「这条射线打中了哪个三角形？」**

教科书上的标准答案是分两步：先求射线与三角形所在平面的交点，再判断这个交点是否落在三角形内部。Möller 与 Trumbore 在 1997 年给出了另一条路径——只用一个 3×3 线性方程组，一次解出射线参数与重心坐标，把叉乘次数减半。下面先统一记号。

**射线**：

$$
\mathbf{R}(t) = \mathbf{O} + t\,\mathbf{D}, \qquad t \in [0, +\infty)
$$

- $\mathbf{O}$：起点
- $\mathbf{D}$：方向向量。**不要求是单位向量**，这一点后面 §9.3 会用到

**三角形**：三个顶点 $\mathbf{P}_0, \mathbf{P}_1, \mathbf{P}_2$，两条边向量

$$
\mathbf{E}_1 = \mathbf{P}_1 - \mathbf{P}_0, \qquad \mathbf{E}_2 = \mathbf{P}_2 - \mathbf{P}_0
$$

从射线起点指向 $\mathbf{P}_0$ 的向量：

$$
\mathbf{S} = \mathbf{O} - \mathbf{P}_0
$$

### 1.1 重心坐标

三角形所在平面上的任意一点都可以写成三个顶点的仿射组合：

$$
\mathbf{P}(u,v) = (1 - u - v)\,\mathbf{P}_0 + u\,\mathbf{P}_1 + v\,\mathbf{P}_2
$$

$(1-u-v,\; u,\; v)$ 就是**重心坐标**（barycentric coordinates），三者恒和为 1。

几何上它们是**面积比**：

$$
u = \frac{S_{\triangle \mathbf{P}_0 Q \mathbf{P}_2}}{S_{\triangle \mathbf{P}_0 \mathbf{P}_1 \mathbf{P}_2}},
\qquad
v = \frac{S_{\triangle \mathbf{P}_0 \mathbf{P}_1 Q}}{S_{\triangle \mathbf{P}_0 \mathbf{P}_1 \mathbf{P}_2}},
\qquad
1-u-v = \frac{S_{\triangle Q \mathbf{P}_1 \mathbf{P}_2}}{S_{\triangle \mathbf{P}_0 \mathbf{P}_1 \mathbf{P}_2}}
$$

这也解释了为什么它天然适合做**属性插值**：给三个顶点各挂一个法线 / UV / 颜色，按 $(1-u-v, u, v)$ 加权就能得到交点处的插值结果。


### 1.2 判定条件

相交 $\iff$ 存在 $(t,u,v)$ 使

$$
\mathbf{O} + t\mathbf{D} = (1-u-v)\mathbf{P}_0 + u\mathbf{P}_1 + v\mathbf{P}_2
$$

且

$$
\boxed{\;t > 0,\qquad u \ge 0,\qquad v \ge 0,\qquad u + v \le 1\;}
$$

前一个条件说交点在射线正方向，后三个说交点在三角形内部（含边界）。

---

## 2. 朴素做法：两步法

### 2.1 第一步，与平面求交

平面法向

$$
\mathbf{N} = \mathbf{E}_1 \times \mathbf{E}_2
$$

过 $\mathbf{P}_0$ 的平面方程 $\mathbf{N}\cdot(\mathbf{X}-\mathbf{P}_0)=0$。代入 $\mathbf{X}=\mathbf{O}+t\mathbf{D}$：

$$
\mathbf{N}\cdot(\mathbf{O} + t\mathbf{D} - \mathbf{P}_0) = 0
\;\Longrightarrow\;
\mathbf{N}\cdot\mathbf{S} + t\,\mathbf{N}\cdot\mathbf{D} = 0
$$

$$
\boxed{\;t = \frac{\mathbf{N}\cdot(\mathbf{P}_0-\mathbf{O})}{\mathbf{N}\cdot\mathbf{D}} = -\frac{\mathbf{N}\cdot\mathbf{S}}{\mathbf{N}\cdot\mathbf{D}}\;}
$$

$\mathbf{N}\cdot\mathbf{D}=0$ 表示射线平行于平面：无交点，或整条射线躺在平面内。

### 2.2 第二步，判断是否在三角形内

得到交点 $\mathbf{Q}=\mathbf{O}+t\mathbf{D}$ 后，两种常见写法。

**边函数（三次叉乘同向）**：

$$
\mathbf{c}_0 = (\mathbf{P}_1-\mathbf{P}_0)\times(\mathbf{Q}-\mathbf{P}_0),\quad
\mathbf{c}_1 = (\mathbf{P}_2-\mathbf{P}_1)\times(\mathbf{Q}-\mathbf{P}_1),\quad
\mathbf{c}_2 = (\mathbf{P}_0-\mathbf{P}_2)\times(\mathbf{Q}-\mathbf{P}_2)
$$

若 $\mathbf{c}_0\cdot\mathbf{N},\ \mathbf{c}_1\cdot\mathbf{N},\ \mathbf{c}_2\cdot\mathbf{N}$ 同号，则 $\mathbf{Q}$ 在内部。

**面积比**：直接按 §1.1 的定义算 $u,v$，代价是三次叉乘 + 三次取模 + 除法。

### 2.3 代价

| 项目 | 两步法 | Möller–Trumbore |
|---|---|---|
| 向量减法 | 6 | 3 |
| 叉乘 | **4** | **2** |
| 点乘 | 5 | 4 |
| 除法 | 1 | 1 |
| 是否构造交点 $\mathbf{Q}$ | 需要 | 不需要 |

### 2.4 三个隐性成本

1. **中间量 $\mathbf{Q}$ 没有复用价值**。它被算出来只为了做内测，而 $t$ 本身已经是最终答案。
2. **内测要重新算一遍边**。$\mathbf{Q}-\mathbf{P}_i$ 是三次减法加三次叉乘，和第一步毫无共享。
3. **信息没利用干净**。$\mathbf{S},\mathbf{D},\mathbf{E}_1,\mathbf{E}_2$ 全都在场，$t$ 里其实已经隐含了重心坐标的信息，但要先求出 $\mathbf{Q}$ 才能把它「读」出来——多绕了一圈。

反过来想：既然 $(t,u,v)$ 是三个未知数，而等式在三维空间中本来是三个标量方程，那它本来就该被**一次性解出来**。

---

## 3. 换个思路：一个 3×3 线性方程组

从核心等式出发，把未知数归到左边。先展开右端：

$$
(1-u-v)\mathbf{P}_0 + u\mathbf{P}_1 + v\mathbf{P}_2
= \mathbf{P}_0 + u(\mathbf{P}_1-\mathbf{P}_0) + v(\mathbf{P}_2-\mathbf{P}_0)
= \mathbf{P}_0 + u\mathbf{E}_1 + v\mathbf{E}_2
$$

代回并移项：

$$
t\mathbf{D} - u\mathbf{E}_1 - v\mathbf{E}_2 = \mathbf{P}_0 - \mathbf{O} = -\mathbf{S}
$$

两边同乘 $-1$：

$$
\boxed{\;-t\,\mathbf{D} + u\,\mathbf{E}_1 + v\,\mathbf{E}_2 = \mathbf{S}\;}
\tag{3.1}
$$

写成矩阵形式（**三列依次是 $-\mathbf{D},\ \mathbf{E}_1,\ \mathbf{E}_2$**）：

$$
\underbrace{
\begin{bmatrix}
\vert & \vert & \vert\\
-\mathbf{D} & \mathbf{E}_1 & \mathbf{E}_2\\
\vert & \vert & \vert
\end{bmatrix}
}_{\mathbf{M}}
\begin{bmatrix} t \\ u \\ v \end{bmatrix}
=
\begin{bmatrix} S_x \\ S_y \\ S_z \end{bmatrix}
\tag{3.2}
$$

$$
\mathbf{M} = \big[\,-\mathbf{D}\;\big|\;\mathbf{E}_1\;\big|\;\mathbf{E}_2\,\big]
$$

于是问题变成：**解一个 3×3 线性方程组**。

式 (3.1) 左端的负号来自移项，随后被矩阵的第一列 $-\mathbf{D}$ 吸收，所以式 (3.2) 的右端是 $+\mathbf{S}$ 而不是 $-\mathbf{S}$。顺着移项走一遍，符号不会错。

---

## 4. 数学工具一：克莱姆法则

### 4.1 定理

设 $\mathbf{A}$ 是 $n\times n$ 方阵，$\det(\mathbf{A})\neq 0$。若 $\mathbf{x}$ 是

$$
\mathbf{A}\mathbf{x} = \mathbf{b}
$$

的解，把 $\mathbf{A}$ 的**第 $i$ 列**替换为 $\mathbf{b}$ 得到 $\mathbf{A}_i$，则

$$
\boxed{\;x_i = \frac{\det(\mathbf{A}_i)}{\det(\mathbf{A})}\;}
$$

**记忆**：分母永远是原矩阵的行列式；分子是把「该未知数对应的那一列」换成常数列。

**它的价值**：不需要消元（无需处理主元为 0 的情形），每个未知数都有一个**独立、闭式**的表达式。当这些行列式能借助向量运算快速求出时，它就是最快的路径。

### 4.2 证明一：换列 + 行列式的列线性

以 $2\times2$ 为例：

$$
\mathbf{A} = \begin{bmatrix} a_{11} & a_{12} \\ a_{21} & a_{22}\end{bmatrix},
\qquad
\mathbf{b} = \begin{bmatrix} b_1 \\ b_2\end{bmatrix},
\qquad
\mathbf{A}_1 = \begin{bmatrix} b_1 & a_{12} \\ b_2 & a_{22}\end{bmatrix}
$$

由 $\mathbf{b}=\mathbf{A}\mathbf{x}$，即 $b_1=a_{11}x_1+a_{12}x_2$、$b_2=a_{21}x_1+a_{22}x_2$，代入：

$$
\mathbf{A}_1 = \begin{bmatrix} a_{11}x_1 + a_{12}x_2 & a_{12} \\ a_{21}x_1 + a_{22}x_2 & a_{22}\end{bmatrix}
$$

**关键性质**：行列式对**每一列/行**是线性的（固定其余列，对当前列线性）。于是把第一列拆开：

$$
\det(\mathbf{A}_1)
= \det\begin{bmatrix} a_{11}x_1 & a_{12}\\ a_{21}x_1 & a_{22}\end{bmatrix}
+ \det\begin{bmatrix} a_{12}x_2 & a_{12}\\ a_{22}x_2 & a_{22}\end{bmatrix}
$$

第一项提出 $x_1$ 就是 $x_1\det(\mathbf{A})$；第二项两列相同，行列式为 $0$。故

$$
\det(\mathbf{A}_1) = x_1\det(\mathbf{A})
\;\Longrightarrow\;
x_1 = \frac{\det(\mathbf{A}_1)}{\det(\mathbf{A})}
\qquad\blacksquare
$$

对 $x_2$ 完全对称。推广到 $n$ 维只是「拆成 $n$ 项、其余 $n-1$ 项因列重复而消失」，论证相同。

### 4.3 证明二：伴随矩阵

$\mathbf{A}$ 可逆时

$$
\mathbf{x} = \mathbf{A}^{-1}\mathbf{b} = \frac{1}{\det(\mathbf{A})}\,\operatorname{adj}(\mathbf{A})\,\mathbf{b},
\qquad
(\operatorname{adj}\mathbf{A})_{ij} = C_{ji}
$$

其中 $C_{ji}$ 是代数余子式。取第 $i$ 个分量：

$$
x_i = \frac{1}{\det(\mathbf{A})}\sum_{j=1}^{n} C_{ji}\,b_j
$$

而右端这个求和正是**把第 $i$ 列换成 $\mathbf{b}$ 后按第 $i$ 列展开**的公式：

$$
\det(\mathbf{A}_i) = \sum_{j=1}^{n} b_j\,C_{ji}
\qquad\Longrightarrow\qquad
x_i = \frac{\det(\mathbf{A}_i)}{\det(\mathbf{A})}
\qquad\blacksquare
$$

> 两个证明的分工：证明一让你**看见**为什么是「换列」；证明二让你**理解**它在矩阵语言里的位置。写代码时两者都不需要，但想清楚「分母为什么总是同一个」，靠的是它们。


## 5. 数学工具二：混合积

### 5.1 定义

对三个三维向量 $\mathbf{a},\mathbf{b},\mathbf{c}$：

$$
[\mathbf{a},\mathbf{b},\mathbf{c}] \;\equiv\; \mathbf{a}\cdot(\mathbf{b}\times\mathbf{c})
$$

先叉乘得向量，再点乘得标量，所以叫**标量三重积**（scalar triple product）；三个向量「混」在一起，也叫**混合积**。

### 5.2 定理：以三个向量为列的行列式 = 混合积

$$
\mathbf{A} = \big[\,\mathbf{a}\;\big|\;\mathbf{b}\;\big|\;\mathbf{c}\,\big] =
\begin{bmatrix}
a_x & b_x & c_x\\
a_y & b_y & c_y\\
a_z & b_z & c_z
\end{bmatrix}
$$

按第一列展开：

$$
\det(\mathbf{A})
= a_x\begin{vmatrix} b_y & c_y\\ b_z & c_z\end{vmatrix}
- a_y\begin{vmatrix} b_x & c_x\\ b_z & c_z\end{vmatrix}
+ a_z\begin{vmatrix} b_x & c_x\\ b_y & c_y\end{vmatrix}
$$

而叉乘的分量是

$$
\mathbf{b}\times\mathbf{c} =
\begin{bmatrix}
b_y c_z - b_z c_y\\
b_z c_x - b_x c_z\\
b_x c_y - b_y c_x
\end{bmatrix}
$$

逐项对照：

$$
\begin{vmatrix} b_y & c_y\\ b_z & c_z\end{vmatrix} = b_yc_z - b_zc_y = (\mathbf{b}\times\mathbf{c})_x,
\qquad
\begin{vmatrix} b_x & c_x\\ b_z & c_z\end{vmatrix} = b_xc_z - b_zc_x = -(\mathbf{b}\times\mathbf{c})_y,
\qquad
\begin{vmatrix} b_x & c_x\\ b_y & c_y\end{vmatrix} = b_xc_y - b_yc_x = (\mathbf{b}\times\mathbf{c})_z
$$

中间那个负号正好与展开式里的 $-\,$ 抵消：

$$
\det(\mathbf{A}) = a_x(\mathbf{b}\times\mathbf{c})_x + a_y(\mathbf{b}\times\mathbf{c})_y + a_z(\mathbf{b}\times\mathbf{c})_z
= \mathbf{a}\cdot(\mathbf{b}\times\mathbf{c})
\qquad\blacksquare
$$

这里之所以选**列**来证，是因为后面要用的都是列：§6 的 $\mathbf{M}=[-\mathbf{D}\,|\,\mathbf{E}_1\,|\,\mathbf{E}_2]$ 就是列拼，直接套用即可。若改以三个向量为**行**拼成矩阵，结论同样成立——$\det(\mathbf{A}^\top)=\det(\mathbf{A})$，且按第一行展开与按第一列展开逐项对应。

### 5.3 性质

| 性质 | 表达式 |
|---|---|
| 轮换不变 | $\mathbf{a}\cdot(\mathbf{b}\times\mathbf{c}) = \mathbf{b}\cdot(\mathbf{c}\times\mathbf{a}) = \mathbf{c}\cdot(\mathbf{a}\times\mathbf{b})$ |
| 交换变号 | $\mathbf{a}\cdot(\mathbf{b}\times\mathbf{c}) = -\,\mathbf{a}\cdot(\mathbf{c}\times\mathbf{b})$ |
| 有重复为零 | 任意两个向量相同或共线 $\Rightarrow 0$ |
| 共面判据 | 三向量共面 $\iff [\mathbf{a},\mathbf{b},\mathbf{c}]=0$ |
| 体积 | $\bigl\vert[\mathbf{a},\mathbf{b},\mathbf{c}]\bigr\vert =$ 平行六面体体积 |
| 点叉可换 | $\mathbf{a}\cdot(\mathbf{b}\times\mathbf{c}) = (\mathbf{a}\times\mathbf{b})\cdot\mathbf{c}$ |

**轮换性的证明**：把 $[\mathbf{a},\mathbf{b},\mathbf{c}]$ 看成以三向量为列的行列式 $\det[\mathbf{a}\mid\mathbf{b}\mid\mathbf{c}]$。

行列式有一条基本性质：**交换任意两行（或两列），行列式变号。** 于是把 $(\mathbf{a},\mathbf{b},\mathbf{c})$ 轮换成 $(\mathbf{b},\mathbf{c},\mathbf{a})$ 只需两步——

$$
(\mathbf{a},\mathbf{b},\mathbf{c}) \;\longrightarrow\; (\mathbf{b},\mathbf{a},\mathbf{c}) \;\longrightarrow\; (\mathbf{b},\mathbf{c},\mathbf{a})
$$

两次交换，变了两次号，负负得正，行列式不变。同理，交换任意两个向量只需一次交换，符号翻转。

两种情形合起来，就是「6 种排列里 3 种相等、3 种互为相反数」的全部理由：**交换次数为偶则保号，为奇则变号**。$\blacksquare$

**与点叉互换的关系**：轮换不变还有一条等价说法，即「叉号可以往前挪一位」：

$$
\mathbf{u}\cdot(\mathbf{v}\times\mathbf{w}) = (\mathbf{u}\times\mathbf{v})\cdot\mathbf{w}
$$

这一条很容易被误当成轮换性的「原因」，其实它是轮换性的**另一种写法**：由轮换式直接代入可得它；反过来用它配合点积的交换律，也能推出轮换式。所以它不提供任何新的信息，只是一条更方便的记账规则。

一行证明：左边是 $\det[\mathbf{u}\mid\mathbf{v}\mid\mathbf{w}]$，右边是 $\mathbf{w}\cdot(\mathbf{u}\times\mathbf{v})=\det[\mathbf{w}\mid\mathbf{u}\mid\mathbf{v}]$，而后者正是前者的三列轮换。

**这条性质属于谁**：既不属于点积，也不属于叉积。

- 点积是对称的（$\mathbf{u}\cdot\mathbf{v}=\mathbf{v}\cdot\mathbf{u}$），它只负责「投影求值」，不产生符号；
- 叉积是反对称的（$\mathbf{u}\times\mathbf{v}=-(\mathbf{v}\times\mathbf{u})$），符号从它这里来；
- 但**符号随排列如何变化**，由行列式的「对换变号」规定，与两个算子各自的对称性都无关。

一个反例说明它不能归给叉积：三重叉积 $\mathbf{a}\times(\mathbf{b}\times\mathbf{c})$ 轮换之后并不相等。取 $\mathbf{a}=(1,0,0)$、$\mathbf{b}=(1,1,0)$、$\mathbf{c}=(0,1,1)$：

$$
\mathbf{a}\times(\mathbf{b}\times\mathbf{c})=(0,-1,-1),\qquad
\mathbf{b}\times(\mathbf{c}\times\mathbf{a})=(-1,1,1),\qquad
\mathbf{c}\times(\mathbf{a}\times\mathbf{b})=(1,0,0)
$$

三者互不相等——**只有「一叉一点」的嵌套才有轮换不变性**。

**工程上的用处**：这条性质让你可以自由调整顺序、挑最省力的写法。§6.1 里那一步

$$
\det\big[-\mathbf{D}\mid\mathbf{E}_1\mid\mathbf{E}_2\big]
= -\mathbf{D}\cdot(\mathbf{E}_1\times\mathbf{E}_2)
= (\mathbf{D}\times\mathbf{E}_2)\cdot\mathbf{E}_1
= \mathbf{S}_1\cdot\mathbf{E}_1
$$

中间的等号就是轮换性。改写成最后的形式，是为了复用 $\mathbf{S}_1=\mathbf{D}\times\mathbf{E}_2$（算 $u$ 时还要用），从而省掉一次叉乘——这正是「叉乘从 4 次降到 2 次」的来源。

### 5.4 几何意义：为什么是体积

1. $\mathbf{b}\times\mathbf{c}$ 的模长 = $\mathbf{b},\mathbf{c}$ 张成的平行四边形**底面积**，方向垂直于该平行四边形；
2. $\mathbf{a}\cdot(\mathbf{b}\times\mathbf{c}) = |\mathbf{a}|\,|\mathbf{b}\times\mathbf{c}|\cos\theta$；
3. 而 $|\mathbf{a}|\cos\theta$ 是 $\mathbf{a}$ 在法向上的**投影长度**，即平行六面体的**高** $h$；
4. 于是 $\big|[\mathbf{a},\mathbf{b},\mathbf{c}]\big| = S_{\text{底}}\cdot h = V$。$\blacksquare$

**轮换的几何解释**：三种写法只是**换了底面**。

| 写法 | 底面 | 高由谁提供 |
|---|---|---|
| $\mathbf{a}\cdot(\mathbf{b}\times\mathbf{c})$ | $\mathbf{b},\mathbf{c}$ 张成 | $\mathbf{a}$ |
| $\mathbf{b}\cdot(\mathbf{c}\times\mathbf{a})$ | $\mathbf{c},\mathbf{a}$ 张成 | $\mathbf{b}$ |
| $\mathbf{c}\cdot(\mathbf{a}\times\mathbf{b})$ | $\mathbf{a},\mathbf{b}$ 张成 | $\mathbf{c}$ |

盒子是同一个盒子。底面换了，高也随之改变，但「底面积 × 高」始终相同——这就是 §5.3 的轮换性在几何上的说法。

> 注意：这个几何解释只直接说明**绝对值**不变。体积一样大，不等于定向一样；符号不变还得靠 §5.3 的「变两次号」。只用「体积相等」去论证轮换不变，其实只证了一半。

**符号的意义**：正 ⟹ $\{\mathbf{a},\mathbf{b},\mathbf{c}\}$ 是右手系；负 ⟹ 左手系；零 ⟹ 共面（体积为 0）。

在图形学里，法线朝向、三角形绕序（winding order）、背面剔除、光线打的是正面还是背面——判据全是混合积的**符号**。

---

## 6. 合成：解出 $t, u, v$

### 6.1 分母

$$
\det(\mathbf{M}) = \det\big[-\mathbf{D}\,\big|\,\mathbf{E}_1\,\big|\,\mathbf{E}_2\big]
= (-\mathbf{D})\cdot(\mathbf{E}_1\times\mathbf{E}_2)
= -\,\mathbf{D}\cdot(\mathbf{E}_1\times\mathbf{E}_2)
$$

用轮换性整理成最好算的形式。定义

$$
\boxed{\;\mathbf{S}_1 = \mathbf{D}\times\mathbf{E}_2\;}
$$

则

$$
\mathbf{S}_1\cdot\mathbf{E}_1 = (\mathbf{D}\times\mathbf{E}_2)\cdot\mathbf{E}_1
= \mathbf{D}\cdot(\mathbf{E}_2\times\mathbf{E}_1)
= -\,\mathbf{D}\cdot(\mathbf{E}_1\times\mathbf{E}_2)
= \det(\mathbf{M})
$$

$$
\boxed{\;\det = \mathbf{S}_1\cdot\mathbf{E}_1\;}
$$

### 6.2 分子

| 未知数 | 替换的列 | 分子行列式 | 化简结果 |
|---|---|---|---|
| $t$ | 第 1 列 | $\det[\mathbf{S}\,\vert\,\mathbf{E}_1\,\vert\,\mathbf{E}_2] = \mathbf{S}\cdot(\mathbf{E}_1\times\mathbf{E}_2)$ | $= \mathbf{S}_2\cdot\mathbf{E}_2$ |
| $u$ | 第 2 列 | $\det[-\mathbf{D}\,\vert\,\mathbf{S}\,\vert\,\mathbf{E}_2] = (-\mathbf{D})\cdot(\mathbf{S}\times\mathbf{E}_2)$ | $= \mathbf{S}_1\cdot\mathbf{S}$ |
| $v$ | 第 3 列 | $\det[-\mathbf{D}\,\vert\,\mathbf{E}_1\,\vert\,\mathbf{S}] = (-\mathbf{D})\cdot(\mathbf{E}_1\times\mathbf{S})$ | $= \mathbf{S}_2\cdot\mathbf{D}$ |

其中

$$
\boxed{\;\mathbf{S}_2 = \mathbf{S}\times\mathbf{E}_1\;}
$$

### 6.3 结果

$$
\boxed{
\begin{aligned}
\mathbf{S}_1 &= \mathbf{D}\times\mathbf{E}_2, &\quad
\mathbf{S}_2 &= \mathbf{S}\times\mathbf{E}_1\\[4pt]
\det &= \mathbf{S}_1\cdot\mathbf{E}_1, &\quad
\mathrm{inv} &= \frac{1}{\det}\\[4pt]
t &= (\mathbf{S}_2\cdot\mathbf{E}_2)\,\mathrm{inv}, &\quad
u &= (\mathbf{S}_1\cdot\mathbf{S})\,\mathrm{inv}, &\quad
v &= (\mathbf{S}_2\cdot\mathbf{D})\,\mathrm{inv}
\end{aligned}}
\tag{6.1}
$$

即

$$
\begin{bmatrix} t\\ u\\ v\end{bmatrix}
= \frac{1}{\mathbf{S}_1\cdot\mathbf{E}_1}
\begin{bmatrix}
\mathbf{S}_2\cdot\mathbf{E}_2\\
\mathbf{S}_1\cdot\mathbf{S}\\
\mathbf{S}_2\cdot\mathbf{D}
\end{bmatrix}
$$

### 6.4 正确性证明（直接验证）

已知 $\det\neq0$，即 $\mathbf{D},\mathbf{E}_1,\mathbf{E}_2$ 线性无关。要验证式 (3.1)：

$$
-t\mathbf{D} + u\mathbf{E}_1 + v\mathbf{E}_2 = \mathbf{S}
$$

技巧：**用叉乘去消项**。拿 $\mathbf{S}_1=\mathbf{D}\times\mathbf{E}_2$ 点乘两端，含 $\mathbf{D}$ 和 $\mathbf{E}_2$ 的项自动消失。

**验证 $u$**：

$$
(-t\mathbf{D} + u\mathbf{E}_1 + v\mathbf{E}_2)\cdot\mathbf{S}_1
= -t\underbrace{\mathbf{D}\cdot(\mathbf{D}\times\mathbf{E}_2)}_{=0}
+ u\underbrace{\mathbf{E}_1\cdot(\mathbf{D}\times\mathbf{E}_2)}_{=\det}
+ v\underbrace{\mathbf{E}_2\cdot(\mathbf{D}\times\mathbf{E}_2)}_{=0}
= u\det
$$

右端 $\mathbf{S}\cdot\mathbf{S}_1 = \mathbf{S}_1\cdot\mathbf{S}$，故 $u = (\mathbf{S}_1\cdot\mathbf{S})/\det$。$\checkmark$

**验证 $v$**：两端点乘 $\mathbf{E}_1\times\mathbf{D}$（垂直于 $\mathbf{E}_1$ 与 $\mathbf{D}$）：

$$
(-t\mathbf{D} + u\mathbf{E}_1 + v\mathbf{E}_2)\cdot(\mathbf{E}_1\times\mathbf{D})
= v\,\mathbf{E}_2\cdot(\mathbf{E}_1\times\mathbf{D})
= v\big[-\mathbf{D}\cdot(\mathbf{E}_1\times\mathbf{E}_2)\big]
= v\det
$$

右端 $\mathbf{S}\cdot(\mathbf{E}_1\times\mathbf{D}) = \mathbf{D}\cdot(\mathbf{S}\times\mathbf{E}_1) = \mathbf{D}\cdot\mathbf{S}_2$，故 $v=(\mathbf{S}_2\cdot\mathbf{D})/\det$。$\checkmark$

**验证 $t$**：两端点乘 $\mathbf{E}_1\times\mathbf{E}_2$（垂直于 $\mathbf{E}_1$ 与 $\mathbf{E}_2$）：

$$
(-t\mathbf{D})\cdot(\mathbf{E}_1\times\mathbf{E}_2) = -t\big[-\det\big] = t\det
$$

右端 $\mathbf{S}\cdot(\mathbf{E}_1\times\mathbf{E}_2) = \mathbf{E}_2\cdot(\mathbf{S}\times\mathbf{E}_1) = \mathbf{S}_2\cdot\mathbf{E}_2$，故 $t = (\mathbf{S}_2\cdot\mathbf{E}_2)/\det$。$\checkmark$

三式全部验证通过；$\det\neq0$ 时解唯一，故这就是唯一解。$\blacksquare$

> 这个证明比「由克莱姆法则，显然」更有价值：它顺便解释了**为什么偏偏是这三个叉乘方向**——因为用 $(\mathbf{E}_i\times\mathbf{E}_j)$ 点乘恰好能消掉另外两个未知数的项。克莱姆法则是「怎么想到的」，这个证明是「为什么必然对」。

---

## 7. 它到底快在哪

先把话说明白：**混合积本身并不比行列式快**。展开一个 3×3 行列式是 9 乘 5 加；用一次叉乘加一次点乘算同一个混合积是 9 乘 + 3 加 + 3 加，并没有更省。

真正的收益在下面四点。

### 7.1 复用

四个行列式（分母 + 三个分子）共享同一批 $\mathbf{S}_1,\mathbf{S}_2$。按定义展开是 $4\times9=36$ 次乘法；复用后是 2 次叉乘 + 4 次点乘 $\approx 24$ 次。**叉乘从 4 次降到 2 次**是主要收益，而叉乘是这套运算里最贵的单目（每次 6 乘 3 减）。

### 7.2 只有一次除法

$\mathrm{inv}=1/\det$ 只算一次，三个未知数各乘一次——而不是三次除法。

### 7.3 提前退出（early-out）

$u$ 只依赖 $\mathbf{S}_1$，所以它可以在 $\mathbf{S}_2$ 之前就判定。$u$ 落在 $[0,1]$ 外时直接返回，此时**只花了 1 次叉乘 + 2 次点乘**。在光线追踪里，绝大多数射线是「未命中」，走的就是这条最便宜的路径。

### 7.4 延迟计算 $t$

$t$ 放到最后。未命中的射线根本不需要 $t$，省下来的不只是那一次点乘，还有一次结果写回。

于是推荐的计算顺序是：

```
S1   = cross(D, E2)
det  = dot(S1, E1)
if |det| < eps:            return miss      // 平行或退化三角形
inv  = 1 / det
u    = dot(S1, S) * inv
if u < 0 || u > 1:         return miss      // 拒绝率最高的关卡
S2   = cross(S, E1)
v    = dot(S2, D) * inv
if v < 0 || u + v > 1:     return miss
t    = dot(S2, E2) * inv
if t <= eps:               return miss      // 交点在射线背后
return hit(t, u, v)
```

顺序不是随意的：$u$ 在最前是因为它最便宜且拒绝率最高；$t$ 在最后是因为它最「不重要」——只在命中时才需要。

---

## 8. 参考实现

一份自包含的 C++ 版本，只依赖一个三维向量和点乘、叉乘：

```cpp
#include <cmath>

struct Vec3 {
    float x, y, z;
    Vec3 operator-(const Vec3& o) const { return {x - o.x, y - o.y, z - o.z}; }
};

inline Vec3 cross(const Vec3& a, const Vec3& b) {
    return {a.y * b.z - a.z * b.y,
            a.z * b.x - a.x * b.z,
            a.x * b.y - a.y * b.x};
}

inline float dot(const Vec3& a, const Vec3& b) {
    return a.x * b.x + a.y * b.y + a.z * b.z;
}

// 命中返回 true，并写出 t（射线参数）、u/v（重心坐标，对应 P1/P2 的权重）
// 未命中时不修改输出参数
bool rayTriangle(const Vec3& p0, const Vec3& p1, const Vec3& p2,
                 const Vec3& orig, const Vec3& dir,
                 float& t, float& u, float& v,
                 bool cullBackface = false)
{
    const Vec3  e1 = p1 - p0;
    const Vec3  e2 = p2 - p0;
    const Vec3  s  = orig - p0;              // 注意方向：orig - p0

    const Vec3  s1  = cross(dir, e2);
    const float det = dot(s1, e1);

    if (cullBackface) { if (det <  1e-9f) return false; }
    else              { if (std::fabs(det) < 1e-9f) return false; }

    const float inv = 1.0f / det;

    const float uu = dot(s1, s) * inv;
    if (uu < 0.0f || uu > 1.0f) return false;

    const Vec3  s2 = cross(s, e1);
    const float vv = dot(s2, dir) * inv;
    if (vv < 0.0f || uu + vv > 1.0f) return false;

    const float tt = dot(s2, e2) * inv;
    if (tt <= 1e-6f) return false;

    t = tt; u = uu; v = vv;
    return true;
}
```

几点说明：

- **背面剔除**用 $\det$ 的符号实现。$\det = -\mathbf{D}\cdot\mathbf{N}$，所以「射线从正面打来」$\iff \mathbf{D}\cdot\mathbf{N}<0 \iff \det>0$。哪一侧算正面取决于三角形绕序，改约定时只改这一处即可。
- 想做成**单面/双面可切换**，就是把 $\det$ 那一行换成上面那个分支。
- **最近命中**（光线追踪的常规需求）不要在这个函数里维护全局最小 $t$，而应把 `t` 返回给调用方比较——保持函数无状态，便于并行。

---

## 9. 边界情况与数值细节

### 9.1 $\det\approx0$ 的物理含义

$\det=0 \iff \mathbf{D},\mathbf{E}_1,\mathbf{E}_2$ 线性相关，涵盖三种情况：

- 射线**平行**于三角形所在平面；
- 三角形**退化**（三点共线，或两顶点重合，$\mathbf{E}_1\parallel\mathbf{E}_2$）；
- 射线方向退化（$\mathbf{D}=\mathbf{0}$）。

三种都该判「不命中」，代码里一句 `|det| < eps` 就能全部覆盖。这也是它相对两步法更省心的地方：两步法要分别处理 $\mathbf{N}=\mathbf{0}$（退化）和 $\mathbf{N}\cdot\mathbf{D}=0$（平行）。

### 9.2 $\epsilon$ 怎么取

- **判 $\det$**：绝对值阈值在三角形尺度接近 1 时够用；更稳的做法是相对量级，例如与 $|\mathbf{D}|\,|\mathbf{E}_1|\,|\mathbf{E}_2|$ 相比。
- **判 $u,v$**：一般用 `u < 0`、`u > 1`；若要容忍共享边被两个三角形重复命中，可加微小容差。注意**不要**给 $u$ 和 $u+v$ 用互相矛盾的容差，否则会出现「判定通过但插值权重为负」的诡异结果。
- **判 $t$**：常用 $10^{-4}\sim10^{-6}$ 量级，避免起点贴在三角形上时自相交。

### 9.3 $\mathbf{D}$ 不必归一化

式 (6.1) 对任意 $\mathbf{D}$ 成立。若 $\mathbf{D}$ 非单位向量，$t$ 的单位是「$\mathbf{D}$ 的倍数」，**不等于世界空间距离**：

$$
\text{真实距离} = t\,|\mathbf{D}|
$$

只比较 $t$ 大小来选最近命中时完全不受影响；但要拿 $t$ 做距离判断、介质吸收、或与别的系统交换数据时，必须记得这一点。

一个直观验证：把 $\mathbf{D}$ 乘以 2，$t$ 会减半，交点不变。

### 9.4 `float` 够不够

常规场景（三角形尺度 $10^{-1}\sim10^{3}$）用 `float` 足够。出现噪点状漏光时，优先怀疑两件事：

1. 射线起点与表面重合（阴影射线的经典问题），应沿法线偏移一个 $\epsilon$；
2. 世界坐标量级过大（$10^{6}$ 以上），`det` 的有效位被吃掉。

后者的对症做法是先把顶点平移到以射线起点为原点——本算法已经这么做了（$\mathbf{S}=\mathbf{O}-\mathbf{P}_0$ 正是这个平移），能显著改善数值表现。

---

## 10. 变体与演进

### 10.1 两个方向的等价表述

§2 的两步法与 §3 的一步法在数学上等价，只是解同一个方程组的不同路径。工程里还有第三种流行写法：**先把射线变换到三角形的局部坐标系**（用平面法向和两条边构造基），在二维里做边函数判定——它与 Möller–Trumbore 的差别只在于「在哪一步除以 $\det$」。

### 10.2 水密求交（watertight）

经典实现在三角形**共享边**上会遇到浮点误差问题：同一条边被相邻两个三角形用不同的顶点组合计算，舍入结果不一致，于是出现**亚像素级裂缝**，光线从缝里穿过去，画面呈现为稀疏噪点（漏光）。

Woop、Benthin、Wald 在 2013 年给出的解法，核心是让共享边两侧的判定**使用完全相同的浮点运算**：

- 先按 $\mathbf{D}$ 的**最大分量**选择主轴 $k_z$，并交换轴，使 $z$ 方向最陡（这样 $1/D_{k_z}$ 不会爆炸）；
- 把三个顶点平移到以射线起点为原点，再用 $1/D_{k_z}$ 缩放，投影到二维；
- 于是求交退化为二维的**边函数符号判定**，共享边两侧算出的是同一组数，符号必然一致。

代价是多几次轴交换和一次除法，收益是**彻底消除裂缝**。渲染器通常只在主光线或大型网格上开启水密模式。

### 10.3 SIMD 与批量求交

在 BVH 遍历中，一次要拿一条射线测 4 条或 8 条三角形。此时标量实现会被改写为：

- 结构体数组（SoA）布局，把多个三角形的 $\mathbf{E}_1,\mathbf{E}_2,\mathbf{P}_0$ 分量分别存成连续数组；
- 用 $\det$ 的符号做**掩码**，把 `if` 全换成 `blend` / `select`；
- 「提前退出」变成「掩码更新」，因为 SIMD 里不能真的提前退出。

这违背了 §7.3 的单射线最优策略——所以**标量场景与 SIMD 场景的最优计算顺序恰好不同**，是两套写法。

---

## 11. 常见错误清单

| # | 错误 | 症状 | 如何定位 |
|---|---|---|---|
| 1 | $\mathbf{S}$ 写成 $\mathbf{P}_0-\mathbf{O}$ | 全部未命中，或只有背面命中 | 手算 §12 的例子，看 $t$ 是否为 $1$ |
| 2 | 移项后左右两端的负号不匹配 | $t<0$、$u,v$ 反号，全部未命中 | 按 §3 的方程组逐项核对 |
| 3 | 漏判 $u+v\le1$ | 三角形外侧（同一平面）也被命中，轮廓「变胖」 | 用两个三角形拼一个正方形测试 |
| 4 | 漏判 $t>0$ | 射线背后的三角形也被命中 | 把相机放进物体内部 |
| 5 | 把 $u,v$ 当纹理坐标用反 | 几何正确但贴图错位 | 用带格子纹理的四边形 |
| 6 | $\mathbf{E}_1,\mathbf{E}_2$ 顺序反了，却没同步改 $\mathbf{S}_1,\mathbf{S}_2$ | $u,v$ 互换，插值错但形状对 | 打印 $u,v$ 与手算值对比 |
| 7 | 用 `det == 0` 判断 | 偶发 NaN，或 $t$ 爆炸到 $10^{30}$ | 改成 `fabs(det) < eps` |
| 8 | 先算 $t$ 再判 $u,v$ | 结果对但慢 2–3 倍 | 按 §7.4 重排 |

第 6 条最隐蔽：**形状看起来完全正常**，只有贴图或法线插值时才露馅——交换 $u,v$ 不改变「是否命中」，只改变「命中点的属性从哪来」。

---

## 12. 一个可以手算的例子

取

$$
\mathbf{O}=(0,0,0),\quad \mathbf{D}=(0,0,1),\quad
\mathbf{P}_0=(0,0,1),\quad \mathbf{P}_1=(1,0,1),\quad \mathbf{P}_2=(0,1,1)
$$

则 $\mathbf{E}_1=(1,0,0)$，$\mathbf{E}_2=(0,1,0)$，$\mathbf{S}=\mathbf{O}-\mathbf{P}_0=(0,0,-1)$。

$$
\mathbf{S}_1 = \mathbf{D}\times\mathbf{E}_2 = (0,0,1)\times(0,1,0) = (-1,0,0)
$$

$$
\det = \mathbf{S}_1\cdot\mathbf{E}_1 = (-1,0,0)\cdot(1,0,0) = -1 \neq 0 \quad\checkmark
$$

$$
u = \frac{\mathbf{S}_1\cdot\mathbf{S}}{\det} = \frac{0}{-1} = 0
$$

$$
\mathbf{S}_2 = \mathbf{S}\times\mathbf{E}_1 = (0,0,-1)\times(1,0,0) = (0,-1,0)
$$

$$
v = \frac{\mathbf{S}_2\cdot\mathbf{D}}{\det} = \frac{0}{-1} = 0,
\qquad
t = \frac{\mathbf{S}_2\cdot\mathbf{E}_2}{\det} = \frac{-1}{-1} = 1 \quad\checkmark
$$

结果 $t=1,\ u=0,\ v=0$：交点正是 $\mathbf{P}_0=(0,0,1)$。

### 调试用测试用例

| # | $\mathbf{O}$ | $\mathbf{D}$ | 期望结果 | 覆盖的行为 |
|---|---|---|---|---|
| 1 | $(0,0,0)$ | $(0,0,1)$ | 命中 $t=1,\ u=0,\ v=0$ | 顶点命中 |
| 2 | 同上，但 $\mathbf{P}_1,\mathbf{P}_2$ 互换 | $(0,0,1)$ | 命中 $t=1,\ u=0,\ v=0$ | 绕序翻转仍命中 |
| 3 | $(0.2,0.2,0)$ | $(0,0,1)$ | 命中 $t=1,\ u=0.2,\ v=0.2$ | 内部点，权重正确 |
| 4 | $(0,0,0)$ | $(0,0,-1)$ | 未命中 | $t<0$ |
| 5 | $(0,0,0)$ | $(1,0,0)$ | 未命中 | 平行，$\det=0$ |
| 6 | 三点共线 | 任意 | 未命中 | 退化三角形 |
| 7 | $(0,0,2)$ | $(0,0,1)$ | 未命中 | 射线背离平面 |

第 3 个用例值得手算：$u=0.2$ 正好对应面积比。

---

## 13. 公式速查

$$
\begin{aligned}
\mathbf{E}_1 &= \mathbf{P}_1-\mathbf{P}_0, \quad \mathbf{E}_2 = \mathbf{P}_2-\mathbf{P}_0, \quad \mathbf{S} = \mathbf{O}-\mathbf{P}_0\\[4pt]
\mathbf{S}_1 &= \mathbf{D}\times\mathbf{E}_2, \quad \mathbf{S}_2 = \mathbf{S}\times\mathbf{E}_1\\[4pt]
\det &= \mathbf{S}_1\cdot\mathbf{E}_1\\[4pt]
t &= \frac{\mathbf{S}_2\cdot\mathbf{E}_2}{\det},\quad
u = \frac{\mathbf{S}_1\cdot\mathbf{S}}{\det},\quad
v = \frac{\mathbf{S}_2\cdot\mathbf{D}}{\det}\\[4pt]
\text{命中} &\iff |\det| > \epsilon \ \wedge\ t>\epsilon \ \wedge\ u\ge0 \ \wedge\ v\ge0 \ \wedge\ u+v\le1
\end{aligned}
$$

一句话记忆：**一叉（$\mathbf{S}_1$）定分母与 $u$，二叉（$\mathbf{S}_2$）定 $v$ 与 $t$。**

---

## 参考

1. Tomas Möller, Ben Trumbore. *Fast, Minimum Storage Ray/Triangle Intersection*. Journal of Graphics Tools, 2(1), 1997.
2. Sven Woop, Carsten Benthin, Ingo Wald. *Watertight Ray/Triangle Intersection*. Journal of Computer Graphics Techniques, 2(1), 2013.
3. Ingo Wald, Sven Woop, Carsten Benthin, et al. *Embree: A Kernel Framework for Efficient CPU Ray Tracing*. ACM Transactions on Graphics, 33(4), 2014.
4. Matt Pharr, Wenzel Jakob, Greg Humphreys. *Physically Based Rendering: From Theory to Implementation*（第 3 版 §3.4 对光线-三角形求交的讨论）。
