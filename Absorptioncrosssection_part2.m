    %Clear All Saved Parameters
clear all; close all
% Define Parameters 
%All Units are SI
h_bar = 1.054571817E-34;
omega_p = (8.95*1.602E-19)/ h_bar;
epsilon_water = 1.33;
epsilon_0 = 8.85418782E-12;
e_gamma = (65.8E-3*1.602E-19)/h_bar;
c = 2.9979E8;
mu_0 = 1.257E-6;
r = 15E-9;
lambda_0 = 450E-9 ;
omega_0 = (2*pi*c)/lambda_0;
omega_barp = (2.96*1.602E-19)/h_bar;
gamma = (0.59*1.602E-19)/h_bar;
gamma_lille = (0.59*1.602E-19)/h_bar;
n_1 = 1.333;

%For Loop Ensures a matrix is calculates for each step in the loop, in this
%case between 10 nm and 1 micrometer
%Nsteps can be increased to smooth the graph

i = sqrt(-1);
Nsteps=400;
lambdav=linspace(10e-9,1000e-9,Nsteps);
for j = 1:Nsteps;
  lambda(j) = lambdav(j);
  
  
  epsilon_drude(j) = 9 - (((omega_p)^2) / (((2*pi*c)/lambda(j))*(((2*pi*c)/lambda(j))+i*e_gamma))) - ((omega_barp)^2 / ((((2*pi*c)/lambda(j))^2)-((omega_0)^2)+(i*((2*pi*c)/lambda(j))*gamma)));
  
  absval = (3*epsilon_water)/((2*epsilon_water)+epsilon_drude(j));
  mysigma(j) = (((4/3*pi*r^3)*((2*pi*c)/lambda(j))*epsilon_0*imag(epsilon_drude(j))/(sqrt(epsilon_water)))*((sqrt(epsilon_0/mu_0))*((abs(absval))^2)))*10^18;
  
  %the E18 in the previous equation is to go from SI units of m^2 to nm^2
  %to allow for comparison 
  
end 

%Next The Data must be imported from Measured Value. Reference here
%Johnson.xlsx and OlmonPRB.xlsx, eventually also measured values!
%If you have files with problems importing look at whether a "," or "." is
%used, excel will in certain cases change this depending on language. Then
%1.000 becomes 1,0E3. Interchangably this is obviously an issue 

%Johnson And Christy Data % File is read, vectors are defined and plotted
% Eventual factors multiplied on are to account for unit differences in
% imported sheets

A = xlsread('Johnson.xlsx');
wl = (10^-6).*A(:,1);
njon = A(:,2);
kvecjon = A(:,3);
epsilonjon = (njon + (i.*kvecjon)).^2;

%Olmon et. Al Single Crystal Read
OlmonPRB2012SC = xlsread('CorrectedValuesOlmon.xlsx');
wlolmon = (1).*OlmonPRB2012SC(:,2);
epsilonolmon = OlmonPRB2012SC(:,4);

%Absorbance Seed
%The absorbance here must be scalled down such they can be compared with
%our nm^2 drude absorption function, hence the factor
%It would be ideal here to normalize the values for better comparison, food
%for thought
Seed = xlsread('10FoldDiluteSeed1.xlsx');
wlseed = (10^-9).*Seed(:,1);
absseed = (1.8*10^-2).*Seed(:,2);


%figure 1 is our Drude Comparioson with measured values for the imaginary
%part
%Units are arbitrary due to the drude function being dimensionless
figure(1)
plot((lambda.*10^9),imag(epsilon_drude),(wl.*10^9),imag(epsilonjon),(wlolmon.*10^9),epsilonolmon)
axis([200,1000,0,7])
grid
title('Permittivity of Au Comparison(Imaginary)');
xlabel('Wavelength [nm]')
ylabel('Arb. Units')

%Figure 2 is the Absorption cross section plot which is compared with
%out measured spectroscopy data

figure(2)
plot((lambda.*10^9), mysigma,(wlseed.*10^9),absseed)
grid
title('Absorption Cross Section - Spherical Au NPs');
xlabel('Wavelength [nm]')
ylabel('nm^2')


%Figure 3 compares the theoretical and measured real values for our
%permittivity 

figure(3)
plot((lambda.*10^9),real(epsilon_drude),(wl.*10^9),real(epsilonjon))
axis([200,1000,-50,10])
grid
title('Permittivity of Au Comparison(Real)');
xlabel('Wavelength [nm]')
ylabel('Arb. Units')


%Everything Below can be regarded as a dump for the Individual polarization
%parts

%x = i*((e_gamma*(omega_p)^2))/(((2*pi*c)/lambda(j))*((((2*pi*c)/lambda(j))^2)+e_gamma));
  %epsilon_drude(j) = 1-((omega_p)^2)/((((2*pi*c)/lambda(j))^2)+((e_gamma)^2)) + x;
  
  %andet forsøg:
  %x = i*(((gamma*(omega_barp^2)*(c/lambda(j))))/((((omega_0^2)-((c/lambda(j))^2))^2 + ((gamma^2)*((c/lambda(j))^2)))));
  %epsilon_drude(j) = 1 + ((((omega_barp^2)*((omega_0^2)-((c/lambda(j))^2))))/(((omega_0^2)-((c/lambda(j))^2))^2 + ((gamma^2)*((c/lambda(j))^2)))) + x;
  
  %tredje forsøg
  
  %fourth attempt
  
  %x1 = i*((e_gamma*(omega_p)^2))/(((2*pi*c)/lambda(j))*((((2*pi*c)/lambda(j))^2)+e_gamma));
  %epsilon_drudeinter(j) = -((omega_p)^2)/((((2*pi*c)/lambda(j))^2)+((e_gamma)^2)) + x1;
  
  %x2 = i*(((gamma*(omega_barp^2)*((2*pi*c)/lambda(j))))/((((omega_0^2)-(((c*2*pi)/lambda(j))^2))^2 + ((gamma^2)*(((c*2*pi)/lambda(j))^2)))));
  %epsilon_drudeintra(j) = -((((omega_barp^2)*((omega_0^2)-(((c*2*pi)/lambda(j))^2))))/(((omega_0^2)-(((c*2*pi)/lambda(j))^2))^2 + ((gamma^2)*(((c*2*pi)/lambda(j))^2)))) + x2;
 
 %epsilon_drude(j) = 1 + epsilon_drudeintra(j) + epsilon_drudeinter(j)
 
 %med refractive index i nævnerne af første led
  %absval = (3*epsilon_water)/((2*epsilon_water)+epsilon_drude(j));
  %mysigma(j) = ((4/3*pi*r^3)*((2*pi*c)/lambda(j))*epsilon_0*imag(epsilon_drude(j))/(n_1))*((sqrt(epsilon_0/mu_0))*((abs(absval))^2));
