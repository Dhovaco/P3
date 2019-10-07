

%Parameterize area
%Define MatrixA MainA
%Define MatrixB MainB
%Define Other Variables




%Parameterize area
%Radius Definition
R =1;

N_theta = 8;
N_phi = 2*N_theta;
N = N_phi*N_theta;

%Position Vector Matrix

counter=1;
for jphi=1:N_phi,
   for jtheta=1:N_theta
       
      theta = ((pi)/N_theta)*(jtheta - 0.5);
        phi = ((2*pi)/N_phi)*(jphi - 0.5);
        
        
        x(counter)=sin(theta)*cos(phi);
        y(counter)=sin(theta)*sin(phi);
        z(counter)=cos(theta);
        counter=counter+1;
    end
end


vecR=[(x);(y);(z)];

normR = (R^-1).*vecR;

%Area for each element

counter=1;
    for j_theta=1:N_theta
        for j_phi = 1:N_phi
            
            theta2 = ((pi)/N_theta)*(j_theta - 0.5);
        phi2 = ((2*pi)/N_phi)*(j_phi - 0.5);
                
AreaR(counter) = (R^2)*((2*pi)/phi2)*(cos(theta2 - (pi/(2*N_theta)))-cos(theta2 + (pi/(2*N_theta))));

 counter=counter+1;
        end 
    end

    A_vec = AreaR;
   
 % Defining P vector 
 
 P = [vecR;A_vec;normR];
 
%Define MatrixA

counter = 1;
for i = 1:N;
    
    delta_i(counter)= P(4,i);
counter = counter + 1;
end

A =(sqrt((delta_i.')/(4*pi))).*eye(N);

 counter = 1;
for ig = 1:N
    i = ig

for jg = 1:N   
    
j=jg
    
g(i,j) = ((1/(4*pi)).*sqrt((P(1,j) - P(1,i))^2 + (P(2,j) - P(2,i))^2 + (P(3,j) - P(3,i))^2))*P(4,j);
    
end

counter = counter + 1;

end

MainA = g + A;


%Define MatrixB 

  counter = 1;
for ig = 1:N
    i = ig

for jg = 1:N   
    
j=jg
    
if not(i==j)

MainB(i,j) = -(1/(4*pi))*((P(1,j)-P(1,i))*(P(5,j))+(P(2,j)-P(2,i))*(P(6,j))+(P(3,j)-P(3,i))*(P(7,j)))*P(4,i)*(1/((sqrt((P(1,j) - P(1,i))^2 + (P(2,j) - P(2,i))^2 + (P(3,j) - P(3,i))^2))^3));

else 
    MainB(i,j)= 0;
end

end

counter = counter + 1;

end

%Define Other Variables
E_0 = 1;
vec0 = zeros(N,1);
phi_0 = E_0.*ones(N,1);
e_w = sym('e_w')
e_m = sym('e_m')

%phi = cell(N,1)
%psi = cell(N,1)



%System 
%Mainsys = [((1/2).*eye(N,N)-MainB),MainA;((1/2).*eye(N,N)+MainB),-(e_w/e_m).*MainA]

%ans = inv(Mainsys)*[psi_0;vec0]
