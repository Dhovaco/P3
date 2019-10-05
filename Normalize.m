
%Reads .xlsx file%
%Input file name%
A = xlsread('xxx.xlsx', 'A:A'); 
B = xlsread('xxx.xlsx', 'B:B');

%Normalizing data%
Np = normalize(B, 'norm')
Nr = normalize(B, 'range')

%Plots figure%
%Input name of graph%
figure('Name','XXX');
plot(A,Np)
xlabel('Wavelength \lambda')
ylabel('A.U.')
title('Normalized Graph')

%Plots figure%
%Input name of graph%
figure('Name','XXX');
plot(A,Nr)
xlabel('Wavelength \lambda')
ylabel('A.U.')
title('Normalized Height Graph')
