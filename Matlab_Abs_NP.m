
%Load xlsx files
Wavelength = xlsread('AUNP_test.xlsx','A:A');
Abs_NP = xlsread('AUNP_test.xlsx','B:B');
Abs_Water = xlsread('WaterRefTest.xlsx','B:B');

%Calculate the Nanoparticles arbsorption
Abs_Total = Abs_NP - Abs_Water;

%Plotting the data into a graph
figure('name','xxx');
plot(Wavelength,Abs_Total);
xlabel('wavelength \lambda');
ylabel('arbsorption');
title('test')

