clear all; close all
%clear all will be removed once the files can work together

%Variables which will be defined elsewhere
%Number of elements
N=10
%The potential given by phi = E_0 z
E_0 = 1


%Define the implied potential - Scaled by electric field - Not it is
%transposed

phi_0  = (E_0.*ones(N,1))

%Define 0 matrix with N elements


%Define Drude Fuction - e1 e2 =  e_w , e(omega)
%Symbolic Variable Definition
    
e_w = sym('e_w')
e_m = sym('e_m')


%A question arises here, as for the final matrix to solve includes the
%fraction ew/e(omega). omega is however not a continuous function, and
%depends on the incident wavelength. How is the matrix then evaluated when
%at which value(s)?

% Define Parameters 
%All Units are SI
%h_bar = 1.054571817E-34;
%omega_p = (8.95*1.602E-19)/ h_bar;
%epsilon_water = 1.33;
%epsilon_0 = 8.85418782E-12;
%e_gamma = (65.8E-3*1.602E-19)/h_bar;
%c = 2.9979E8;
%mu_0 = 1.257E-6;
%r = 15E-9;
%lambda_0 = 450E-9 ;
%omega_0 = (2*pi*c)/lambda_0;
%omega_barp = (2.96*1.602E-19)/h_bar;
%gamma = (0.59*1.602E-19)/h_bar;
%gamma_lille = (0.59*1.602E-19)/h_bar;



%i = sqrt(-1);
%Nsteps=400;
%lambdav=linspace(10e-9,1000e-9,Nsteps);
%for j = 1:Nsteps;
 % lambda(j) = lambdav(j);
  
 % epsilon_drude(j) = 9 - (((omega_p)^2) / (((2*pi*c)/lambda(j))*(((2*pi*c)/lambda(j))+i*e_gamma))) - ((omega_barp)^2 / ((((2*pi*c)/lambda(j))^2)-((omega_0)^2)+(i*((2*pi*c)/lambda(j))*gamma)));   
%end 

%Define Input variable potential for N-element Matrix 
%A symbolic matrix is created for N-Element Potential

phi = sym('phi%d%d', [1 N])

%Define Normal Gradient N-element Matrix

%???

%Define A Matrix 
% Attempts with adding matrcies



%Define B Matrix



%Solution of linear Equation


