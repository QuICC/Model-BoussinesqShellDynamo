# Boussinesq thermal dynamo convection in a spherical shell

The model equations are

$$
(\partial_t - Pm \nabla^2){\bf u}  =  {\bf u}\times(\nabla\times{\bf u}) -\nabla \Pi  +\frac{Pm^2}{E}\widetilde{Ra} \Theta \frac{\bf r}{r_o}  - \frac{Pm}{E}\hat{\bf z}\times{\bf u}  + \frac{Pm}{E}(\nabla\times{\bf B})\times {\bf B},
$$

$$
    (\partial_t - \frac{Pm}{Pr}\nabla^2)\Theta =- ({\bf u}\cdot \nabla) \Theta -u_r \frac{d T_c}{dr},
$$

$$
    (\partial_t - \nabla^2) {\bf B} = \nabla\times({\bf u} \times {\bf B }.)
$$

with the parameters defined as

$$
Pr = \frac{\nu}{\kappa};\quad Pr = \frac{\nu}{\eta};\quad E=\frac{\nu}{2\Omega d^2};\quad Ra=\frac{ \alpha g_o \mathcal{T} d}{2\nu\Omega}
$$

where $\mathcal{T}$ is the temperature nondimensionalisation value and depends on the temperature boundary conditions and heating mode.

The system is solved using a Toroidal/Poloidal decomposition of the velocity

$$
{\bf u} = \nabla\times T_u {\bf r} + \nabla\times\nabla\times P_u {\bf r} .
$$

and the magnetic field:

$$
{\bf B} = \nabla\times T_B {\bf r} + \nabla\times\nabla\times P_B {\bf r} .
$$