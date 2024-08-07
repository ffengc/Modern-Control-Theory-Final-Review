## 常用拉普拉斯变换

$$
\begin{array}{ccc}
\hline \text { 序号 } & \text { 拉氏变换 } \mathrm{F}(\mathrm{s}) & \text { 时间函数 } \mathrm{f}(\mathrm{t}) \\
\hline 1 & 1 & \delta(\mathrm{t}) \\
 2 & \frac{1}{1-e^{-5}} & \delta_{T}(t)=\sum_{n=0}^{\infty} \delta(t-n T) \\
 3 & \frac{1}{s} & 1(t) \\
 4 & \frac{1}{s^{2}} & \mathrm{t} \\
 5 & \frac{1}{s^{3}} & \frac{t^{2}}{2} \\
 6 & \frac{1}{s^{n+1}} & \frac{t^{n}}{n!} \\
 7 & \frac{1}{s+a} & e^{-a t} \\
 8 & \frac{1}{(s+a)^{2}} & t e^{-a t} \\
 9 & \frac{a}{s(s+a)} & 1-e^{-a t} \\
 10 & \frac{b-a}{(s+a)(s+b)} & e^{-a t}-e^{-b t} \\
 11 & \frac{\omega}{s^{2}+\omega^{2}} & \sin \omega t \\
 12 & \frac{s}{s^{2}+\omega^{2}} & \cos \omega t \\
 13 & \frac{\omega}{(s+a)^{2}+\omega^{2}} & e^{-a t} \sin \omega t \\
 14 & \frac{s+a}{(s+a)^{2}+\omega^{2}} & e^{-a t} \cos \omega t \\
 15 & \frac{1}{s-(1 / T) \ln a} & a^{t / T} \\
\hline
\end{array}
$$


齐次性 $L[a f(t)]=a F(s)$ 

叠加性 $L\left[f_{1}(t) \pm f_{2}(t)\right]=F_{1}(s) \pm F_{2}(s)$

微分定理一般形式:

$$
L\left[\frac{d f(t)}{d t}\right]=s F(s)-f(0) \\
L\left[\frac{d^{2} f(t)}{d t^{2}}\right]=s^{2} F(s)-s f(0)-f^{\prime}(0) \\
\vdots \\
L\left[\frac{d^{n} f(t)}{d t^{n}}\right]=s^{n} F(s)-\sum_{k=1}^{n} s^{n-k} f^{(k-1)}(0) \\
f^{(k-1)}(t)=\frac{d^{k-1} f(t)}{d t^{k-1}} \\
$$

微分定理（当初始条件为0的时候）:

$$
L\left[\frac{d^{n} f(t)}{d t^{n}}\right]=s^{n} F(s) \\
$$

积分定理一般形式: 

$$
L\left[\int f(t) d t\right]=\frac{F(s)}{s}+\frac{\left[\int f(t) d t\right]_{t-0}}{s} \\

L\left[\iint f(t)(d t)^{2}\right]=\frac{F(s)}{s^{2}}+\frac{\left[\int f(t) d t\right]_{t_{-} 0}}{s^{2}}+\frac{\left[\iint f(t)(d t)^{2}\right]_{t_{-} 0}}{s} \\

L\left[\int \cdots \int f(t)(d t)^{n}\right]=\frac{F(s)}{s^{n}}+\sum_{k=1}^{n} \frac{1}{s^{n-k+1}}\left[\int \cdots \int f(t)(d t)^{n}\right]_{t-0}
$$

积分定理（当初始条件为0的时候）:

$$
L\left[\int_{\cdots} \int f(t)(d t)^{n}\right]=\frac{F(s)}{s^{n}}$$
$$
## 电感电容

$$
U_L = L \frac{\mathrm{d} i}{\mathrm{~d} t} \\
C \frac{\mathrm{d} u_{C}}{\mathrm{~d} t}=i
$$

## 常用z变换

![image-20240705164132783](./assets/image-20240705164132783.png)
