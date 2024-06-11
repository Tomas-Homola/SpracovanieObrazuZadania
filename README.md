# Image filtration app

Implementation of several image filtering methods such as convolution, linear/non-linear diffusion and curvature based models (MCF, GMCF) using **C++**, **Qt** (for UI design) and **OpenMP** library (for parallel computations). The filters are based on **partial differential equations** that are solved numerically using the finite volume method with explicit, implicit or semi-implicit schemes. The implicit and semi-implicit schemes utilize either SOR or BiCGStab linear solvers. All models are coupled with the homogenous Neumann boundary condition $\nabla u \cdot \mathbf{n} = 0$. The linear diffusion is formulated by the following equation:

$$
  \frac{\partial\ u}{\partial\ t} = \nabla^{2}\ u
$$

and it has a blurring effect, as can be seen bellow:

<img src="https://github.com/Tomas-Homola/SpracovanieObrazuZadania/assets/61438447/7d0368e3-4616-40a9-87ad-70a719f07f04" width="450"/> <img src="https://github.com/Tomas-Homola/SpracovanieObrazuZadania/assets/61438447/9d3a8f49-9f16-4d76-8f74-447f1e6e6c79" width="450"/> 

As for the non-linear diffusion, we have the regularized Perona-Malik model

$$
  \frac{\partial\ u}{\partial\ t} = \nabla \cdot \big( g( |\nabla G_{\sigma} * u| ) \nabla u \big),
$$

which is useful when we want to preserve edges:

<img src="https://github.com/Tomas-Homola/SpracovanieObrazuZadania/assets/61438447/7d0368e3-4616-40a9-87ad-70a719f07f04" width="450"/> <img src="https://github.com/Tomas-Homola/SpracovanieObrazuZadania/assets/61438447/367aea5e-b589-4c22-9f76-bf33fc165f3b" width="450"/> 

The curvature based filters, e.g. Mean Curvature Flow, are derived from the level-set formulation and it leads to the following equation:

$$
  \frac{1}{|\nabla u|} \frac{\partial\ u}{\partial\ t} = \nabla \cdot \left( \frac{\nabla u}{| \nabla u |} \right)
$$

and this filter is great at removing salt&pepper noise.

<img src="https://github.com/Tomas-Homola/SpracovanieObrazuZadania/assets/61438447/b2d6b10a-fe43-4ff7-8e9d-c4a706f63b68" width="450"/> <img src="https://github.com/Tomas-Homola/SpracovanieObrazuZadania/assets/61438447/426dc0b1-e859-40d8-bc07-1c10f1982714" width="450"/> 
