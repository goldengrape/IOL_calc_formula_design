#!/usr/bin/env python
# coding: utf-8

# # 光线传输矩阵
# 
# 光线传输矩阵又称做[ABCD矩阵](https://en.wikipedia.org/wiki/Ray_transfer_matrix_analysis) 用来描述近轴光线的传播过程. 
# 
# ![](https://upload.wikimedia.org/wikipedia/commons/thumb/5/5e/RayTransferMatrixDefinitions.svg/600px-RayTransferMatrixDefinitions.svg.png)
# 
# 光线穿过一个光学系统的过程可以用一个2x2的矩阵来表示, 如果: 
# * 入射光线距离光轴是x1, 角度是 θ1,  
# * 出射光线距离光轴是x2, 角度是 θ2,
# 
# 那么: 
# $$
# \begin{bmatrix}x_2\\ \theta_2 \end{bmatrix} = \begin{bmatrix}A & B \\C & D \end{bmatrix} \begin{bmatrix}x_1\\ \theta_1 \end{bmatrix}
# $$
# 
# 对于不同的光学元件, ABCD的取值不同. 

# # ABCD Matrix
# 
# 内置上常用的参数字典, 其实最常用的也就是在自由空间内传播, 和在球面上折射, 而这两者可以合并成在球面上折射然后再传播一段距离. 
# 
# * 在自由空间传播距离为d的矩阵是: 
# $$
# free\space space\space matrix=\begin{bmatrix} 1 & d \\ 0 & 1 \end{bmatrix}
# $$
# * 在球面上折射的矩阵, 如果从n1折射率介质经过半径为R的球面, 折射到n2折射率介质: 
# $$
# curved\space interface\space matrix=\begin{bmatrix} 1 & 0 \\ \frac{(n_1-n_2)}{R n_2} & \frac{n_1}{n_2} \end{bmatrix}
# $$
# * 两者合并是相乘的关系, 注意乘法的顺序, 从右至左应当先经过曲面, 再自由传播: 
# 
# $$
# pass\space curved\space interface\space  matrix=
# \begin{bmatrix} 1 & d \\ 0 & 1 \end{bmatrix}
# \begin{bmatrix} 1 & 0 \\ \frac{(n_1-n_2)}{R n_2} & \frac{n_1}{n_2} \end{bmatrix}
# $$
# 
# 

# ## 薄透镜推导
# 
# 从$n_1$介质, 经过半径为$R_1$的球面进入$n_2$介质, 再经过忽略不计的厚度后, 从$n_2$介质经过半径为$R_2$的球面, 进入到$n_1$介质 
# 
# $$
# \begin{bmatrix} 1 & 0 \\ \frac{(n_2-n_1)}{R_2 n_1} & \frac{n_2}{n_1} \end{bmatrix}
# \begin{bmatrix} 1 & 0 \\ \frac{(n_1-n_2)}{R_1 n_2} & \frac{n_1}{n_2} \end{bmatrix}
# =
# \begin{bmatrix}  {1} & {0} \\ {-\frac{n_{2}-n_{1}}{n_{1}}\left(\frac{1}{R_{1}}-\frac{1}{R_{2}}\right)} & {1}\end{bmatrix}
# $$
# 
# 定义焦距$f$, 使$f$满足
# $$
# \frac{1}{f}=
# {\frac{n_{2}-n_{1}}{n_{1}}\left(\frac{1}{R_{1}}-\frac{1}{R_{2}}\right)}
# $$
# 
# 定义屈光度$P=\frac{1}{f}$ [此处存疑], 可能应当是$P=\frac{n}{f}$, n为两侧介质
# 

# $$
# thin\space lens\space matrix=
# \begin{bmatrix}  {1} & {0} \\ {-\frac{1}{f}} & {1}\end{bmatrix}
# =
# \begin{bmatrix}  {1} & {0} \\ {-P} & {1}\end{bmatrix}
# $$ 

# # 眼球的光线传输矩阵
# 
# 光线依次通过:
# 
# * 8 眼镜(薄透镜): 眼镜度数=R_t=目标设定值
# $$
# \begin{bmatrix}  {1} & {0} \\ {-R_t} & {1}\end{bmatrix} % 眼镜(薄透镜)
# $$
# * 7 眼镜到角膜前表面的自由空间(自由空间): 眼镜顶点距=$V_d=12mm$
# $$
# \begin{bmatrix} 1 & V_d \\ 0 & 1 \end{bmatrix} % 眼镜到角膜前表面的自由空间(自由空间)
# $$
# * 6 角膜前表面(球面上的折射): 空气折射率$n_{air}=1.00$, 角膜折射率=$n_{c}=1.376$, 角膜前表面曲率半径=$R_a$=测量值
# $$
# \begin{bmatrix} 1 & 0 \\ \frac{(n_{air}-n_{c})}{R_{a} n_{c}} & \frac{n_{air}}{n_{c}} \end{bmatrix} %角膜前表面(球面上的折射)
# $$
# * 5 角膜厚度(自由空间传播): 角膜厚度$T_c=0.5mm$
# $$
# \begin{bmatrix} 1 & T_c \\ 0 & 1 \end{bmatrix} %角膜厚度(自由空间传播)
# $$
# * 4 角膜后表面(球面上的折射): 角膜折射率=$n_{c}=1.376$, 房水折射率=$n_{h}$, 角膜后表面曲率半径=$R_p$=测量值
# $$
# \begin{bmatrix} 1 & 0 \\ \frac{(n_{c}-n_{h})}{R_{p} n_{h}} & \frac{n_{c}}{n_{h}} \end{bmatrix} % 角膜后表面(球面上的折射)
# $$
# * 3 前房(自由空间传播): 前房深度定义为角膜后表面到人工晶体前表面=$ELP$=需要估计的值
# $$
# \begin{bmatrix} 1 & ELP \\ 0 & 1 \end{bmatrix} %前房(自由空间传播)
# $$ 
# * 2 人工晶体(薄透镜)
# $$
# \begin{bmatrix}  {1} & {0} \\ {-IOL} & {1}\end{bmatrix} %人工晶体(薄透镜)
# $$
# * 1 人工晶体后到视网膜(自由空间传播): 眼轴长=AL
# $$
# \begin{bmatrix} 1 & AL-T_c-ELP \\ 0 & 1 \end{bmatrix} %人工晶体后到视网膜(自由空间传播)
# $$
# 
# 注意相乘时是从右至左连乘

# $$
# eye\space matrix=
# \begin{bmatrix} 1 & AL-T_c-ELP \\ 0 & 1 \end{bmatrix}                                              %人工晶体后到视网膜(自由空间传播)
# \begin{bmatrix}  {1} & {0} \\ {-IOL} & {1}\end{bmatrix}                                            %人工晶体(薄透镜)
# \begin{bmatrix} 1 & ELP \\ 0 & 1 \end{bmatrix}                                                     %前房(自由空间传播)
# \begin{bmatrix} 1 & 0 \\ \frac{(n_{c}-n_{h})}{R_{p} n_{h}} & \frac{n_{c}}{n_{h}} \end{bmatrix}     %角膜后表面(球面上的折射)
# \begin{bmatrix} 1 & T_c \\ 0 & 1 \end{bmatrix}                                                     %角膜厚度(自由空间传播)
# \begin{bmatrix} 1 & 0 \\ \frac{(n_{air}-n_{c})}{R_{a} n_{c}} & \frac{n_{air}}{n_{c}} \end{bmatrix} %角膜前表面(球面上的折射)
# \begin{bmatrix} 1 & V_d \\ 0 & 1 \end{bmatrix}                                                     %眼镜到角膜前表面的自由空间(自由空间)
# \begin{bmatrix}  {1} & {0} \\ {-R_t} & {1}\end{bmatrix}                                            %眼镜(薄透镜)
# $$

# ## 角膜K值
# 
# 由于测量设备原理的限制, 很多设备是无法测量角膜后表面曲率半径R_p的, 因此将角膜前后表面合并
# 
# * 6 角膜前表面(球面上的折射): 空气折射率$n_{air}=1.00$, 角膜折射率=$n_{c}=1.376$, 角膜前表面曲率半径=$R_a$=测量值
# $$
# \begin{bmatrix} 1 & 0 \\ \frac{(n_{air}-n_{c})}{R_{a} n_{c}} & \frac{n_{air}}{n_{c}} \end{bmatrix}
# $$
# * 5 角膜厚度(自由空间传播): 角膜厚度$T_c=0.5mm$
# $$
# \begin{bmatrix} 1 & T_c \\ 0 & 1 \end{bmatrix}
# $$
# * 4 角膜后表面(球面上的折射): 角膜折射率=$n_{c}=1.376$, 房水折射率=$n_{h}$, 角膜后表面曲率半径=$R_p$=测量值
# $$
# \begin{bmatrix} 1 & 0 \\ \frac{(n_{c}-n_{h})}{R_{p} n_{h}} & \frac{n_{c}}{n_{h}} \end{bmatrix}
# $$
# 
# 合并后是一个单一的球面, 其曲率半径仍是角膜前表面的曲率半径$R_a$, 折射率$n_k=1.3375$
# $$
# \begin{bmatrix} 1 & 0 \\ \frac{(n_{air}-n_{k})}{R_{a} n_{k}} & \frac{n_{air}}{n_{k}} \end{bmatrix}
# $$
# 其中, 又把$\frac{n_k-n_{air}}{R_a}=K$
# 于是有
# $$
# \begin{bmatrix} 1 & 0 \\ - \frac{K}{n_{k}} & \frac{n_{air}}{n_{k}} \end{bmatrix}                   %K值表示的角膜
# $$
# 

# 那么使用角膜K值来定义的, 
# 
# $$
# eye\space matrix\space by\space K=
# \begin{bmatrix} 1 & AL-T_c-ELP \\ 0 & 1 \end{bmatrix}                                              %人工晶体后到视网膜(自由空间传播)
# \begin{bmatrix}  {1} & {0} \\ {-IOL} & {1}\end{bmatrix}                                            %人工晶体(薄透镜)
# \begin{bmatrix} 1 & ELP \\ 0 & 1 \end{bmatrix}                                                     %前房(自由空间传播)
# \begin{bmatrix} 1 & 0 \\ - \frac{K}{n_{k}} & \frac{n_{air}}{n_{k}} \end{bmatrix}                   %K值表示的角膜
# \begin{bmatrix} 1 & V_d \\ 0 & 1 \end{bmatrix}                                                     %眼镜到角膜前表面的自由空间(自由空间)
# \begin{bmatrix}  {1} & {0} \\ {-R_t} & {1}\end{bmatrix}                                            %眼镜(薄透镜)
# $$
# 
# 其中需要令T_c=0

# ## 平行光聚焦于视网膜
# 
# 就是从高度为任意x, 角度为0入射到眼镜上的光线, 最终到达视网膜时, 高度=0, 角度无所谓.
# $$
# \begin{bmatrix} 0 \\ \theta \end{bmatrix} = 
# eye\space matrix 
# \begin{bmatrix}x \\ 0 \end{bmatrix}
# $$
# 
# 如果eye matrix也写成$2 \times 2$的矩阵
# $$
# eye\space matrix = \begin{bmatrix} eye_{0,0} & eye_{0,1} \\ eye_{1,0} & eye_{1,1} \end{bmatrix}                                              %人工晶体后到视网膜(自由空间传播)
# $$
# 
# 那么对于上述方程, 当$eye_{0,0}==0$时 可以求解出平行光聚焦到视网膜上的结果, 
# 
# * 当已知目标度数, 则IOL为未知数进行求解. 这是通常计算人工晶体度数的过程.
# * 当已知植入的IOL, 则使用R_t作为未知数进行求解. 这是对公式中的常数参数进行估计的过程. 
# 
# 注意, 在手术之前, 对于人工晶体的位置ELP, 实际上是不知道的, 只能通过各种已知的参数进行估计. 

# 
# 代入前面eye matrix
# 
# $$
# \begin{bmatrix} 0 \\ \theta \end{bmatrix} = 
# \begin{bmatrix} 1 & AL-T_c-ELP \\ 0 & 1 \end{bmatrix}
# \begin{bmatrix}  {1} & {0} \\ {-IOL} & {1}\end{bmatrix}
# \begin{bmatrix} 1 & ELP \\ 0 & 1 \end{bmatrix}
# \begin{bmatrix} 1 & 0 \\ \frac{(n_{c}-n_{h})}{R_{p} n_{h}} & \frac{n_{c}}{n_{h}} \end{bmatrix}
# \begin{bmatrix} 1 & T_c \\ 0 & 1 \end{bmatrix}
# \begin{bmatrix} 1 & 0 \\ \frac{(n_{air}-n_{c})}{R_{a} n_{c}} & \frac{n_{air}}{n_{c}} \end{bmatrix}
# \begin{bmatrix} 1 & V_d \\ 0 & 1 \end{bmatrix}
# \begin{bmatrix}  {1} & {0} \\ {-R_t} & {1}\end{bmatrix}
# \begin{bmatrix}x \\ 0 \end{bmatrix}
# $$
# 
# 或者使用角膜K值
# 
# $$
# \begin{bmatrix} 0 \\ \theta \end{bmatrix} = 
# \begin{bmatrix} 1 & AL-T_c-ELP \\ 0 & 1 \end{bmatrix}                                              %人工晶体后到视网膜(自由空间传播)
# \begin{bmatrix}  {1} & {0} \\ {-IOL} & {1}\end{bmatrix}                                            %人工晶体(薄透镜)
# \begin{bmatrix} 1 & ELP \\ 0 & 1 \end{bmatrix}                                                     %前房(自由空间传播)
# \begin{bmatrix} 1 & 0 \\ - \frac{K}{n_{k}} & \frac{n_{air}}{n_{k}} \end{bmatrix}                   %K值表示的角膜
# \begin{bmatrix} 1 & V_d \\ 0 & 1 \end{bmatrix}                                                     %眼镜到角膜前表面的自由空间(自由空间)
# \begin{bmatrix}  {1} & {0} \\ {-R_t} & {1}\end{bmatrix}                                            %眼镜(薄透镜)
# \begin{bmatrix}x \\ 0 \end{bmatrix}
# \\
# T_c=0
# $$
# 

# # 求解
# 
# 这么长的矩阵乘法, 虽然写出来容易, 但计算的时候当然要用程序

# 导入SymPy符号计算库, 定义所需要的符号变量

# In[1]:


from sympy import *
x,θ=symbols('x, θ')
n_air, n_c,n_h,n_k =symbols('n_air,n_c,n_h,n_k')
AL, T_c, ELP, IOL, R_a, R_p, V_d, R_t, K=symbols('AL, T_c, ELP, IOL, R_a, R_p, V_d, R_t, K')


# 分别列出每一步光线传播的过程矩阵, 顺便检查

# In[2]:


glass=Matrix([[1, 0], [-R_t,1]])
glass


# In[3]:


glass_to_cornea=Matrix([[1, V_d], [0,1]])
glass_to_cornea


# In[4]:


cornea_anterior=Matrix([[1,0],[(n_air-n_c)/(R_a*n_c), n_air/n_c]])
cornea_anterior


# In[5]:


cornea_inside=Matrix([[1, T_c], [0,1]])
cornea_inside


# In[6]:


cornea_posterior=Matrix([[1,0],[(n_c-n_h)/(R_p*n_h), n_c/n_h]])
cornea_posterior


# In[7]:


cornea_by_K=Matrix([[1,0],[-K/n_k, n_air/n_k]])
cornea_by_K


# In[8]:


anterior_chamber=Matrix([[1,ELP],[0,1]])
anterior_chamber


# In[9]:


pass_IOL=Matrix([[1, 0], [-IOL/n_h,1]])
pass_IOL


# In[10]:


vitreous=Matrix([[1, AL-T_c-ELP], [0,1]])
vitreous


# In[11]:


eye_matrix=vitreous * pass_IOL * anterior_chamber * cornea_posterior * cornea_inside * cornea_anterior * glass_to_cornea
eye_matrix[0,0]


# In[12]:


eye_matrix_by_K=(vitreous * pass_IOL * anterior_chamber * cornea_by_K * glass_to_cornea).subs(T_c,0)
eye_matrix_by_K[0,0]


# In[13]:


IOL_power_by_K= solve(eye_matrix_by_K[0,0], IOL)[0]
IOL_power_by_R= solve(eye_matrix[0,0], IOL)[0]


# In[14]:


IOL_power_by_K.subs({n_k:1.3375, AL:23/1000, K:40, ELP:4/1000, n_h:1.336})


# In[ ]:




