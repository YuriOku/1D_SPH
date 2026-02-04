There are several ways to formulate the time evolution equations in the SPH method.
The following two systems of equations are implemented in 1D_SPH.
These methods are derived and explained in detail in Rosswog (2009). Please refer to Rosswog (2009) if you are interested.

## Vanilla ice
The standard method to obtain the time evolution equation for SPH particles is to discretize the system of fluid equations expressed in Lagrange form by a kernel approximation.
In the SPH method, the motion of the SPH particle represents the fluid flow, so the conservation of mass is automatically satisfied, and the continuity equation is not needed.
The equation of motion and energy equation is 

$$
\begin{align*}
\frac{dv}{dt} &= -\frac{\nabla P}{\rho}\\
\frac{du}{dt} &= -\frac{P}{\rho}\nabla \cdot v
\end{align*}
$$

,where
$v$
is velocity,
$P$
is pressure, 
$\rho$
is density and
$u$
is internal energy per unit mass. 
The spatial derivative is calculated using the kernel approximation.
Under the kernel approximation, the physical quantity
$f$
is expressed as

$$
f(x) = \int f(x') W(|x - x'|, h) dx'
\cdot\cdot\cdot\cdot\cdot(1)
$$

,where 
$W$
is the [kernel function](https://github.com/YuriOku/1D_SPH/wiki/kernel-function) for the particle spread.
When this is discretized using the volume element 
$dV$
of each SPH particle, it becomes

$$
f_i = \sum_j f_j dV_j W(| x_i - x_j|, h_i).
\cdot\cdot\cdot\cdot\cdot(2)
$$

The volume element is obtained from the physical quantity 
$Z$
 of the SPH particle and the corresponding continuous quantity 
$y$
  as
$dV = \frac{Z}{y}$
, where 
$Z = m,\,y = \rho$
 is used as standard. For more information, please see the [volume element](https://github.com/YuriOku/1D_SPH/wiki/volume-element) page.
The spatial derivative of 
$f$
 is expressed as 

$$
\nabla_i f_i = \sum_j f_j dV_j \nabla_i W(| x_i - x_j|, h_i).
\cdot\cdot\cdot\cdot\cdot(3)
$$

For the calculation of 
$\nabla_i W$
, see [Kernel Interpolation Method](https://github.com/YuriOku/1D_SPH/wiki/Kernel-interpolation).
Before applying eq. (3) to the right-hand side of the equation of motion, we use the relation 

$$
\begin{align*}
\frac{\nabla P}{\rho} = \frac{P}{\rho^\lambda}\nabla\left( \frac{1}{\rho^{1-\lambda}} \right) + \frac{1}{\rho^{2 - \lambda}}\nabla \left(\frac{P}{\rho^{\lambda - 1}} \right)
\end{align*}
$$

(Monaghan 1992) so that the equation is antisymmetric for the exchange of i, j where 
$\lambda$
 is a constant.
Traditionally, the form for the 
$\lambda = 2$
 case is used, which is naturally obtained from the derivative of 
$P/\rho$ 
, but
1D_SPH uses 
 $\lambda = 1$
 to suppress the noise that is generated when there are SPH particles of different masses (Ritchie & Thomas 2001).
The form with 
$\lambda = 1$ 
 is used in GASOLINE2 (Wadsley et al. 2017) and MAGMA2 (Rosswog 2020).
The SPH system of equations then becomes

$$
\begin{align*}
\frac{d\vec{v_i}}{dt} &= -\frac{1}{m_i}\sum_j Z_i Z_j \left( \frac{P_i + P_j}{y_i y_j}\right) \nabla_i \tilde{W}_{ij}\\
\frac{du_i}{dt} &= \frac{P_i}{m_i} \sum_j \frac{Z_i}{y_i}\frac{Z_j}{y_j}\vec{v_{ij}}\cdot\nabla_i \tilde{W}_{ij}
\end{align*}
$$

Here,
$\vec{v_{ij}}=\vec{v_i} - \vec{v_j}$
,
$\tilde{W_{ij}}= [W(r_{ij},h_i)+W(r_{ij}, h_j)]/2$
.

## Lagrangian
Another way to derive the SPH equation is to start with the Lagrangian and use the variational principle.
The Lagrangian of the SPH particle system is 

$$
L = \sum_i m_i(\frac{v_i^2}{2} - u_i).
\cdot\cdot\cdot\cdot\cdot(4)
$$

This is substituted into the Euler-Lagrange equation 

$$
\frac{d}{dt}\left(\frac{\partial L}{\partial \vec{v}_i} \right) - \frac{\partial L}{\partial \vec{r}_i} = 0.
\cdot\cdot\cdot\cdot\cdot(5)
$$

The first term is the left-hand side of the equation of motion.
Substituting the Lagrangian into the second term, we obtain 

$$
\frac{\partial L}{\partial \vec{r}_i} = -\sum_j m_j\frac{\partial u_j}{\partial \rho_j}
\frac{\partial \rho_j}{\partial \vec{r}_i}.
\cdot\cdot\cdot\cdot\cdot(6)
$$

 Now, using the first law of thermodynamics for isentropic flow

$$
du = - PdV = \frac{P}{\rho^2}d\rho
\cdot\cdot\cdot\cdot\cdot(7)
$$

, the Euler-Lagrange equation becomes 

$$
m_i\frac{d\vec{v}_i}{dt} = -\sum_j m_j\frac{P_j}{\rho_j^2}
\frac{\partial \rho_j}{\partial \vec{r}_i}.
\cdot\cdot\cdot\cdot\cdot(8)
$$

When we formulate the SPH equations from the density (the case with $Z = m,\,y = \rho$), substitute

$$
\rho_j = \sum_k m_k W(r_{jk}, h_j)
\cdot\cdot\cdot\cdot\cdot(9)
$$

into eq. (8). Then the derivative of the density becomes 

$$
\frac{\partial \rho_j}{\partial \vec{r}_i} = \sum_k m_k \left(\nabla_i W(r_{jk}, h_j)
+\frac{\partial W(r_{jk}, h_j)}{\partial h_j} \frac{\partial h_j}{\partial \rho_k}\frac{\partial \rho_k}{\partial \vec{r}_i} \right) = f_j \sum_k m_k \nabla_i W(r_{jk}, h_j)
\cdot\cdot\cdot\cdot\cdot(10)
$$

, where 


$$
f_i = \left(1 - \frac{\partial h_i}{\partial \rho_i} \sum_j m_j \frac{\partial}{\partial h_i} W(r_{ij}, h_i) \right)^{-1}
$$

 is a term to take into account changes in the smoothing length 
$h$
Called the "grad-h-term".
The equation of motion 

$$
\frac{d \vec{v}_i}{dt} = - \sum_j m_j \left(
f_i \frac{P_i}{\rho_i^2} \nabla_i W(r_{ij}, h_i) + 
f_j \frac{P_j}{\rho_j^2} \nabla_i W(r_{ij}, h_j)
\right)
$$

 can be obtained by transforming the equation, paying attention to the subscripts.

The time derivative of the energy is 

$$
\frac{du_i}{dt} = \frac{P_i}{\rho_i^2}\frac{d \rho_i}{dt}
= f_i \frac{P_i}{\rho_i^2} \sum_j m_j \vec{v_{ij}} \cdot \nabla_i W(r_{ij}, h_i)
$$

From the first law of thermodynamics.

For the case with a general volume element, we substitute 

$$
y_j = \sum_k Z_k W(r_{jk}, h_j)
\cdot\cdot\cdot\cdot\cdot(11)
$$

into Euler-Lagrange eq. (8) and obtain

$$
\begin{align*}
\frac{d \vec{v}_i}{dt} &= - \frac{1}{m_i}\sum_j Z_i Z_j \left(
f_i \frac{P_i}{y_i^2} \nabla_i W(r_{ij}, h_i) + 
f_j \frac{P_j}{y_j^2} \nabla_i W(r_{ij}, h_j)
\right)\\
\frac{du_i}{dt} &= f_i \frac{P_i}{m_i y_i^2} \sum_j Z_i Z_j \vec{v_{ij}} \cdot \nabla_i W(r_{ij}, h_i).
\end{align*}
$$

## References
- Monaghan, J. J., “Smoothed particle hydrodynamics.”, <i>Annual Review of Astronomy and Astrophysics</i>, vol. 30, pp. 543–574, 1992. doi:10.1146/annurev.aa.30.090192.002551.
- Ritchie, B. W. and Thomas, P. A., “Multiphase smoothed-particle hydrodynamics”, <i>Monthly Notices of the Royal Astronomical Society</i>, vol. 323, no. 3, pp. 743–756, 2001. doi:10.1046/j.1365-8711.2001.04268.x.
- Rosswog, S., ``Astrophysical smooth particle hydrodynamics'', <i>New Astronomy Reviews</i>, vol. 53, no. 4–6, pp. 78–104, 2009. doi:10.1016/j.newar.2009.08.007.
- Rosswog, S., “The Lagrangian hydrodynamics code MAGMA2”, <i>Monthly Notices of the Royal Astronomical Society</i>, vol. 498, no. 3, pp. 4230–4255, 2020. doi:10.1093/mnras/staa2591.
- Wadsley, J. W., Keller, B. W., and Quinn, T. R., “Gasoline2: a modern smoothed particle hydrodynamics code”, <i>Monthly Notices of the Royal Astronomical Society</i>, vol. 471, no. 2, pp. 2357–2369, 2017. doi:10.1093/mnras/stx1643.
