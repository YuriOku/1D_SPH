In this section, we will discuss two methods of calculating the spatial derivative of the kernel function 
$\nabla_i W$.

## Standard gradient
One method is to compute the derivative

$$
\nabla_i W\left(r_{ij}, h \right) 
= \frac{\partial r_{ij}}{\partial \mathbf{r}_i}  
\frac{\partial W(r_{ij}, h)}{\partial r_{ij}} 
= \mathbf{e}_{ij} \frac{\partial W(r_{ij}, h)}{\partial r_{ij}},
$$

where

$$
r_{ij} = \left| \mathbf{r}_i - \mathbf{r}_j \right|,\ \mathbf{e}_{ij} = \frac{\mathbf{r}_i - \mathbf{r}_j}{\left|\mathbf{r}_i - \mathbf{r}_j \right|}.
$$

For 
$\frac{\partial W}{\partial r}$
, see [Kernel Functions](https://github.com/YuriOku/1D_SPH/wiki/Kernel-function).

## Integral Approach
The Integral Approach (IA), proposed by García-Senz et al. (2012), can handle Rayleigh-Taylor instability and Kelvin-Helmholtz instability with better accuracy than the conventional SPH method.
First, we have the integral

$$
I(\mathbf{r}) = \int \left(f(\mathbf{r}') - f(\mathbf{r}) \right) \left(\mathbf{r}' - \mathbf{r} \right) W\left(\left|\mathbf{r}' - \mathbf{r} \right|, h \right) d{r'}^d, \cdot\cdot\cdot\cdot\cdot (1)
$$

where 
$d$
 is the dimension of space.
Taylor expansion of 
$f(\mathbf{r}')$
around a point
$\mathbf{r}$
gives

$$
f(\mathbf{r}') - f(\mathbf{r}) = \nabla f (\mathbf{r}) \cdot (\mathbf{r}' - \mathbf{r}) + \mathcal{O}({\mathbf{r}'}^2)
\cdot\cdot\cdot\cdot\cdot (2)
$$

We ignore the higher-order terms, substitute Eq. (2) into Eq. (1), and consider solving for 
$\nabla f(\mathbf{r})$
.

### 1D
In one dimension, equation (1) becomes

$$
I(r) = \frac{\partial f(r)}{\partial r} \int  \left(r' - r \right)^2 W\left(\left|r' - r \right|, h \right) dr'.
\cdot\cdot\cdot\cdot\cdot (3)
$$

Solving equation (3) for
$\frac{\partial f(r)}{\partial r}$
and substituting equation (1) for 
$I(r)$
, we get

$$
\frac{\partial f(r)}{\partial r} =  \frac{\int \left(f(r') - f(r) \right) \left(r' - r \right) W\left(\left|r' - r \right|, h \right) dr'}{\int  \left(r' - r \right)^2 W\left(\left|r' - r \right|, h \right) dr'}.
\cdot\cdot\cdot\cdot\cdot (4)
$$

By using the kernel approximation, we can replace this integral with the kernel sum

$$
\frac{\partial f(r_i)}{\partial r} = 
\frac{\sum_j dV_j \left(f(r_j) - f(r_i) \right) \left(r_j - r_i \right) W\left(\left|r_j - r_i \right|, h \right)}{\sum_j  dV_j \left(r_j - r_i \right)^2 W\left(\left|r_j - r_i \right|, h \right)}.
\cdot\cdot\cdot\cdot\cdot (5)
$$

This derivative is correct when the field is linear. For example, the derivative of the density field 

$$
\rho_j = \rho_i + a(r_j - r_i)
$$

 is

$$
\begin{align*}
\frac{\partial \rho_i}{\partial r} 
&= \frac{\sum_j dV_j \left(\rho_j - \rho_i \right) \left(r_j - r_i \right) W\left(\left|r_j - r_i \right|, h \right)}{\sum_j  dV_j \left(r_j - r_i \right)^2 W\left(\left|r_j - r_i \right|, h \right)}\\
&= \frac{\sum_j  dV_j a \left(r_j - r_i \right)^2 W\left(\left|r_j - r_i \right|, h \right)}{\sum_j  dV_j \left(r_j - r_i \right)^2 W\left(\left|r_j - r_i \right|, h \right)}\\
&= a.
\end{align*}
$$

This method is equivalent to the linear-exact gradient (Price 2004; Rosswog 2015).
However, this method does not conserve momentum because the equation of motion is not antisymmetric for the exchange of 
$i,\,j$.
To be antisymmetric, the equation of motion must be of the form

$$
\begin{align*}
&\frac{\partial f_i}{\partial r} = \sum_j f_j G_{ij}\\
&(G_{ij} = -G_{ji}).
\end{align*}
$$

In the Integral Approach, we assume

$$
\sum_j dV_j \left(\mathbf{r}_j - \mathbf{r}_i \right) W\left(\left|\mathbf{r}_j - \mathbf{r}_i \right|, h \right) = 0.
$$

This is the assumption that the distribution of particles is unbiased and uniform. Using this assumption, equation (1) can be expressed as

$$
I(\mathbf{r}) \sim \sum_j dV_j f(\mathbf{r}_j) \left(\mathbf{r}_j - \mathbf{r}_i \right) W\left(\left|\mathbf{r}_j - \mathbf{r}_i \right|, h \right).
\cdot\cdot\cdot\cdot\cdot (7)
$$

As a result, the derivative (5) becomes

$$
\frac{\partial f(r_i)}{\partial r} = 
\frac{\sum_j dV_j f(r_j) (r_j - r_i ) W\left(\left|r_j - r_i \right|, h \right)}{\sum_j  dV_j (r_j - r_i )^2 W\left(\left|r_j - r_i \right|, h \right)}.
$$

By comparing with the usual derivative in SPH

$$
\nabla_i f(r_i) = \sum_j f(r_j) dV_j \nabla_i W(| r_i - r_j|, h),
$$

the spatial derivative of the kernel function in the Integral Approach

$$
\nabla_i W\left(\left| r_i - r_j \right|, h\right) = \frac{(r_j - r_i ) W\left(\left|r_j - r_i \right|, h \right)}{\sum_j  dV_j (r_j - r_i )^2 W\left(\left|r_j - r_i \right|, h \right)}
$$

can be obtained.

### 3-D
Equation (1) in three dimensions is

$$
\left[ \begin{array}{c}
I_1 (\mathbf{r})\\
I_2 (\mathbf{r})\\
I_3 (\mathbf{r})
\end{array} \right] = \int \left( \left[ \begin{array}{c}
\partial f(\mathbf{r})/\partial x_1\\
\partial f(\mathbf{r})/\partial x_2\\
\partial f(\mathbf{r})/\partial x_3
\end{array} \right] \cdot
\left[ \begin{array}{c}
x_1' - x_1\\
x_2' - x_2\\
x_3' - x_3
\end{array} \right] \right)
\left[ \begin{array}{c}
x_1' - x_1\\
x_2' - x_2\\
x_3' - x_3
\end{array} \right]
W\left(\left| \mathbf{r}' - \mathbf{r} \right|, h \right) d{r'}^3,
\cdot\cdot\cdot\cdot\cdot (11)
$$

where 

$$
\mathbf{r} = x_1 \mathbf{i} + x_2 \mathbf{j} + x_3 \mathbf{k}
$$

and

$$
I(\mathbf{r}) = I_1 (\mathbf{r})\mathbf{i} + I_2 (\mathbf{r})\mathbf{j} + I_3 (\mathbf{r})\mathbf{k}.
$$

Solving this for 
$\nabla f(\mathbf{r})$
, we obtain

$$
\left[ \begin{array}{c}
\partial f(\mathbf{r})/\partial x_1\\
\partial f(\mathbf{r})/\partial x_2\\
\partial f(\mathbf{r})/\partial x_3
\end{array} \right] =
\left[ \begin{array}{ccc}
\tau_{11} & \tau_{12} & \tau_{13}\\
\tau_{21} & \tau_{22} & \tau_{23}\\
\tau_{31} & \tau_{32} & \tau_{33}
\end{array} \right]^{-1}
\left[ \begin{array}{c}
I_1 (\mathbf{r})\\
I_2 (\mathbf{r})\\
I_3 (\mathbf{r})
\end{array} \right].
\cdot\cdot\cdot\cdot\cdot (12)
$$

where 

$$
\tau_{ij} = \int (x_i'  - x_i)(x_j' - x_j) W\left(\left| \mathbf{r}' - \mathbf{r}\right|, h\right) d{r'}^3
$$

and expressing it as a kernel sum, we get

$$
\tau_{ij} = \sum_b dV_b (x_{i,b}  - x_{i,a})(x_{j,b} - x_{j,a}) W\left(\left| \mathbf{r}_b - \mathbf{r}_a\right|, h\right).
$$

The inverse matrix in the first term on the right-hand side of equation (12) is

$$
C = \frac{1}{\tau_{11}\tau_{22}\tau_{33} - \tau_{11}\tau_{23}\tau_{32} + \tau_{12}\tau_{21}\tau_{33} -\tau_{12}\tau_{23}\tau_{31} + \tau_{13}\tau_{21}\tau_{32} - \tau_{13}\tau_{22}\tau_{31}}
\left[ \begin{array}{ccc}
\tau_{22}\tau_{33} - \tau_{23}\tau_{32} & \tau_{13}\tau_{32} - \tau_{12}\tau_{33} & \tau_{12}\tau_{23} - \tau_{13}\tau_{22}\\
\tau_{23}\tau_{31} - \tau_{21}\tau_{33} & \tau_{11}\tau_{33} - \tau_{13}\tau_{31} & \tau_{13}\tau_{21} - \tau_{11}\tau_{23}\\
\tau_{21}\tau_{32} - \tau_{22}\tau_{31} & \tau_{12}\tau_{31} - \tau_{11}\tau_{32} & \tau_{11}\tau_{22} - \tau_{12}\tau_{21}
\end{array} \right].
$$

The i-component of the integral (1)

$$
I_i(\mathbf{r}) = \int \left(f(\mathbf{r}') - f(\mathbf{r}) \right) \left(x_i' - x_i \right) W\left(\left|\mathbf{r}' - \mathbf{r} \right|, h \right) d{r'}^3
$$

 can be approximated by eliminating 
 $f(\mathbf{r})$ 
 from the kernel

$$
I_i(\mathbf{r}_a) = \sum_b dV_b f(\mathbf{r}_b) \left(x_{i,b} - x_{i,a} \right) W\left(\left|\mathbf{r}_b - \mathbf{r}_a \right|, h \right).
$$

The derivative of the function 
$f(\mathbf{r})$ 
 can be obtained from equation (12) as

$$
\frac{\partial f(\mathbf{r}_a)}{\partial x_i} = \sum_{l = 1}^3 C^{il}I_l(\mathbf{r}_a).
$$

Substituting Eq. (16) into Eq. (17) and comparing it with the derivative of the usual SPH, the derivative of the kernel function becomes

$$
\frac{\partial W\left(\left|\mathbf{r}_b - \mathbf{r}_a\right|, h\right)}{\partial x_i} = \sum_{l = 1}^3 C^{il} (x_{l,b} - x_{l,a}) W\left(\left|\mathbf{r}_b - \mathbf{r}_a\right|, h\right).
$$


## References
- García-Senz, D., Cabezón, R. M., and Escartín, J. A., "Improving smoothed particle hydrodynamics with an integral approach to calculating gradients", <i>Astronomy and Astrophysics</i>, vol. 538, 2012. doi:10.1051/0004-6361/201117939.
- Price, D. J., "Magnetic fields in Astrophysics", PhDT, 2004.
- Rosswog, S., "Boosting the accuracy of SPH techniques: Newtonian and special-relativistic tests", <i>Monthly Notices of the Royal Astronomical Society</i>, vol. 448, no. 4, pp. 3628-3664, 2015. doi:10.1093/mnras/stv225.
