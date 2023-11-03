function [T1,T2] = bm_geotherm(y)
u = 1e3 / (3600*24*365*1e6); %erosion rate
L = 100e3; %lithospheric thickness
k  = 2.5; %thermal conductivity [W m^{-1} K^{-1}]
rho = 2700; %density [kg m^-3]
Cp = 1000; %heat capacity [J kg^{-1} K^{-1}]
kappa = k/rho/Cp; %thermal diffusivity, [m^2 s^{-1}]
As = 3e-6; %volumetric heat production at the surface [W m^{-3}]
Ts = 10; %surface temperature (oC)
Tm = 1300; %mantle temperature (oC)
eta = u*L/kappa;
% steady-state geotherm with exponential heat production, no erosion, no topography
h = 10e3; %heat production exponent
a = As * h^2 / k;
x = linspace(0,50e3,100); %distance (m)
% y = linspace(0,10e3,100); %depth (m)
% [x2d,y2d]=meshgrid(x,y);
ph = 4.5e-3;%atmospheric adiabatic gradient [oC m^{-1}]
H0 = 0.75e3; %half amplitude topography [m]
lam = 30e3; %wavelength of topography [m]
yc = 30e3; %crustal thickness [m]
qm = 60e-3; %mantle heat flux [mW m^{-2}]

T0 = Ts + y/L * (Tm - Ts);

T1 = Ts + a * (1 - exp(- y/h) ) + y/L ...
    * ( (Tm - Ts) - a * ( 1 - exp(- L/h) ) );

b = As/rho/Cp * h^2 /(kappa - u *h );

ga = ((Tm-Ts) - b * (1 - exp(-L/h)) ) ...
    / (1 - exp(-u/kappa * L));

T2 = Ts + b * (1 - exp(- y/h) ) + ...
    + ga * ( 1 - exp(- u/kappa * y) ) ;



% m2 = 1/2 * (-u/kappa - sqrt( (u/kappa)^2 + (4*pi/lam)^2 ) );

% D = (u/kappa*ga + b/h - ph )*H0;

% T3 = repmat(T2',1,length(x)) + D * cos(2*pi/lam * x2d) .* exp(m2 *y2d);
% 
% T4 = Ts + qm*y/k + (As/k)*y.*(yc - y/2); 
% T4 = Ts + qm*y/k + h^2*(As/k)*(1 - exp(-y/h) - y/h*exp(-yc/h)); 