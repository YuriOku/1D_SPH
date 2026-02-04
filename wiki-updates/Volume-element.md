In the SPH method, the fluid is discretized into SPH particles and handled. To do this, we need the volume 
$dV$
 of each particle.
The volume is obtained from the physical quantity 
$Z$
 of the SPH particle and the corresponding continuous field 
$y$
, with 

$dV = \frac{Z}{y}$
.

## $Z = m,\,y = \rho$ 

Calculating the volume element using density as the fundamental quantity is the standard method in the SPH method, which is also called the **Standard SPH** (SSPH).
In this case, the density can be expressed as the kernel sum of the masses as follows:

$$
\begin{align*}
\rho_i &= \sum_j dV_j \rho_j W(r_{ij}, h_i)\\
 &= \sum_j \frac{m_j}{\rho_j} \rho_j W(r_{ij}, h_i) \\
 &= \sum_j m_j W(r_{ij}, h_i) 
\end{align*}
$$

In the standard SPH method, the density calculated in the above equation is used to determine other physical quantities.
For example, pressure can be obtained from the equation of state, which in 1D_SPH is the ideal gas equation of state, 

$$
P = (\gamma - 1)\rho u
$$


,where 
$\gamma$
 is the specific heat ratio and 
$u$
 is the internal energy per unit mass.

## $$
Z = mu,\,y = P/(\gamma - 1)
$$ 

In the standard SPH method, other physical quantities are calculated based on the density obtained as the kernel sum. Kernel summation leads to a field distortion, making it impossible to treat discontinuously varying densities such as contact discontinuities with high accuracy (Agertz et al., 2007). Saitoh & Makino (2013) proposed a method to formulate the SPH method using pressure as the basic quantity to solve this problem. This method is called density-independent SPH (DISPH). The basic quantity, pressure, can be obtained as 

$$
P_i= (\gamma - 1)\sum_j m_j u_j W(r_{ij}, h)
$$

and use this pressure as the basis for calculating other physical quantities.
Since the pressure is constant at the contact discontinuity, DISPH can handle the physical quantities at the contact discontinuity with high accuracy.

## References
- Agertz, O. et al., “Fundamental differences between SPH and grid methods”, <i>Monthly Notices of the Royal Astronomical Society</i>, vol. 380, no. 3, pp. 963–978, 2007. doi:10.1111/j.1365-2966.2007.12183.x.
- Saitoh, T. R. and Makino, J., “A Density-independent Formulation of Smoothed Particle Hydrodynamics”, <i>The Astrophysical Journal</i>, vol. 768, no. 1, 2013. doi:10.1088/0004-637X/768/1/44.