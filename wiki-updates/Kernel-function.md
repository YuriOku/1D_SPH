In the SPH method, the field is represented by a kernel approximation.
Under the kernel approximation, the physical quantity 
 $f$
 is represented as 

$$
f(x) = \int f(x') W(|x - x'|, h) dx'
$$

 where 
$W$
 is the kernel function that represents the spread of the particle, and 
$h$
 is the smoothing length that represents the spread of the kernel function.

The kernel function must satisfy the following three properties:

- It is normalized.

$$
\int W(\vec{r}, h) dV = 1
$$

- In the limit of 
 $h \rightarrow 0$,
it reduces to
$\delta$ function.

- Twice differentiable

![kernel approximation](https://user-images.githubusercontent.com/62641316/113504614-3140b000-9574-11eb-8b45-4bd6cbd5fede.png)

The length corresponding to the radius of the SPH particle (the distance at which the value of the kernel function becomes zero) is called the kernel support radius, and in general, the kernel support radius is 
$H = 2h$
.
However, the SPH simulation code GADGET (Springel et al. 2001; Springel 2005; Springel et al. 2020), which is commonly used in the field of astrophysics, uses 
$H = h$
. 1D_SPH also uses this definition.
Usually, in 3D simulations, the smoothing length is determined so that the number of particles contained within the kernel support radius is a constant value 
$N_H$.
In the 1D code, 1D_SPH, the smoothing length is 

$h = \eta \left(\frac{m}{\rho}\right)$

. This is obtained by an iterative method. Here, 
$\eta$
 is a parameter, and the standard value is 
$\eta = 2.4$
.

## Cubic spline function
The most standard kernel function is the cubic spline function (Monaghan & Lattanzio 1985) 

$$
W(r, h) = C
\begin{cases}
1 - 6\left(\frac{r}{h}\right)^2 + 6\left(\frac{r}{h}\right)^3 & \left( 0 < \frac{r}{h} \leq \frac{1}{2} \right)\\
2\left( 1- \frac{r}{h} \right)^3 &\left( \frac{1}{2} < \frac{r}{h} \leq 1\right)\\
0 &\left( 1 < \frac{r}{h} \right)
\end{cases}
$$

. where 
$C$
 is a normalization constant, and in one dimension it is 
$C = \frac{4}{3h}$
.
The radial derivative of the kernel function is 

$$
\frac{\partial W(r, h)}{\partial r} = \frac{C}{h}
\begin{cases}
-12\left(\frac{r}{h}\right) + 18\left(\frac{r}{h}\right)^2 & \left( 0 < \frac{r}{h} \leq \frac{1}{2} \right)\\
6\left( 1- \frac{r}{h} \right)^2 &\left( \frac{1}{2} < \frac{r}{h} \leq 1\right)\\
0 &\left( 1 < \frac{r}{h} \right).
\end{cases}
$$


The spline series also includes the higher-order Quartic spline (4th order), Quintic spline (5th order),... functions.

## Wendland C2 functions
The Wendland C2 function (Wendland 1995) 

$$
W(r, h) = C
\begin{cases}
\left(1 - \frac{r}{h} \right)^3 \left( 1 + 3\frac{r}{h} \right) & \left( \frac{r}{h} \leq 1 \right) \\
0 &  \left( \frac{r}{h} > 1 \right) 
\end{cases}
$$

 is another commonly used kernel function, where the normalization constant is 
$C = \frac{5}{4h}$
 in one dimension.
The Wendland kernel avoids a numerical instability called pairing instability, which can be seen when the number of neighbours is increased in the spline kernel (Dehnen & Aly 2012).
The derivative of the Wendland C2 function is 

$$
\frac{\partial W(r, h)}{\partial r} = \frac{C}{h}
\begin{cases}
-\frac{4r}{h}\left(1 - \frac{r}{h} \right)^2  & \left( \frac{r}{h} \leq 1 \right) \\
0 &  \left( \frac{r}{h} > 1 \right).
\end{cases}
$$

Higher-order C4 and C6 functions are also available in the Wendland series.

## References
- Dehnen, W. and Aly, H., “Improving convergence in smoothed particle hydrodynamics simulations without pairing instability”, <i>Monthly Notices of the Royal Astronomical Society</i>, vol. 425, no. 2, pp. 1068–1082, 2012. doi:10.1111/j.1365-2966.2012.21439.x.
- Monaghan, J. J. and Lattanzio, J. C., “A refined particle method for astrophysical problems”, <i>Astronomy and Astrophysics</i>, vol. 149, no. 1, pp. 135–143, 1985.
- Springel, V., Yoshida, N., and White, S. D. M., “GADGET: a code for collisionless and gasdynamical cosmological simulations”, <i>New Astronomy</i>, vol. 6, no. 2, pp. 79–117, 2001. doi:10.1016/S1384-1076(01)00042-2.
- Springel, V., “The cosmological simulation code GADGET-2”, <i>Monthly Notices of the Royal Astronomical Society</i>, vol. 364, no. 4, pp. 1105–1134, 2005. doi:10.1111/j.1365-2966.2005.09655.x.
- Springel, V., Pakmor, R., Zier, O., and Reinecke, M., “Simulating cosmic structure formation with the GADGET-4 code”, <i>arXiv e-prints</i>, 2020.
- Wendland, H. "Piecewise polynomial, positive definite and compactly supported radial functions of minimal degree." Advances in computational Mathematics 4.1, 389-396, 1995
