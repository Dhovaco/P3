
%Choose # Elements
N = 100

%Radius Definition
R =1

%Angle Definition  %ASSUME Ntheta = Nphi

N_theta = (1/2)*N
N_phi = (1/2)*N

%Define Linspace for j discretization

j_phi = linspace(1,N_phi,N_phi)
j_theta = linspace(1,N_theta,N_theta)

%Define Angle Matrix

phi = ((2*pi)/N_phi)*(j_phi - 0.5)
theta = ((pi)/N_theta)*(j_theta - 0.5)

%Position Vector Matrix

x = (sin(theta).').*(cos(phi).')
y = (sin(theta).').*(sin(phi).')
z = (cos(theta).')

vecR=[(x),(y),(z)]

normR = (R^-1).*vecR

%Area for each element

AreaR = (R^2)*((2*pi)/N_phi)*(cos(theta - (pi/(2*N_theta)))-cos(theta + (pi/(2*N_theta))))

