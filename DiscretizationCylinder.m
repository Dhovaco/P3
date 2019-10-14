clear all;
%Cylinder
%========
%Radius of cylinder
R=7;
%number of points along z in the cylinder
Nzc=14;

%number of points round the circumference
Nphi=20;

%remember to clear workspace after changing Nphi

%halv length of cylinder
L=10;

%z value along the centre axes of cylinder
%z_cyl=linspace(-L+(L/Nzc),L-(L/Nzc),Nzc);
%z=z_cyl;

dz_cyl=L/Nzc;

%pac is "Place Along the Centre"
%rac is "points Round About the Circumference" 

counterc=1;
for pac=1:Nzc
    for rac=1:Nphi
    
    phi=((2*pi)/Nphi)*(rac-0.5);
    
    %position z,x,y
    zc(counterc)=-L+(L/Nzc)+(2*L/Nzc)*(-1+pac);
    xc(counterc)=R*cos(phi);
    yc(counterc)=R*sin(phi);
    
    %normal vector
    z_n_cyl(counterc)=0;
    x_n_cyl(counterc)=cos(phi);
    y_n_cyl(counterc)=sin(phi);
    
    areaCyl(counterc)=R*(2*L/Nzc)*((2*pi)/Nphi);
    
    counterc=counterc+1;
    end
end

normVeccyl=[x_n_cyl;y_n_cyl;z_n_cyl];

%
N_theta = 10;
N_phi = 24;
N = N_phi*N_theta;

%Position Vector Matrix for the 4 quarter-spheres- 
%that makes the ends
countersph=1;
for jphi=1:N_phi/4
   for jtheta=1:N_theta
       
      theta = ((pi)/N_theta)*(jtheta - 0.5);
        phis = ((2*pi)/N_phi)*(jphi - 0.5);
        
        zs(countersph)=R*sin(theta)*cos(phis)+L;
        ys(countersph)=R*sin(theta)*sin(phis);
        xs(countersph)=R*cos(theta);
        
        normsph1z(countersph)=sin(theta)*cos(phis);
        normsph1y(countersph)=sin(theta)*sin(phis);
        normsph1x(countersph)=cos(theta);
        
        areasph(countersph) = (R^2)*((2*pi)/phi)*(cos(theta - (pi/(2*N_theta)))-cos(theta + (pi/(2*N_theta))));
        
        countersph=countersph+1;
    end
end

for jphi=1:N_phi/4
   for jtheta=1:N_theta
       
      theta = ((pi)/N_theta)*(jtheta - 0.5);
        phis = ((2*pi)/N_phi)*(jphi - 0.5);
        
        
        zs(countersph)=R*sin(theta)*cos(phis)+L;
        ys(countersph)=R*(-sin(theta))*sin(phis);
        xs(countersph)=R*cos(theta);
        
        normsph1z(countersph)=sin(theta)*cos(phis);
        normsph1y(countersph)=-sin(theta)*sin(phis);
        normsph1x(countersph)=cos(theta);
        
        areasph(countersph) = (R^2)*((2*pi)/phi)*(cos(theta - (pi/(2*N_theta)))-cos(theta + (pi/(2*N_theta))));
        
        countersph=countersph+1;
    end
end

for jphi=1:N_phi/4
   for jtheta=1:N_theta
       
      theta = ((pi)/N_theta)*(jtheta - 0.5);
        phis = ((2*pi)/N_phi)*(jphi - 0.5);
        
        
        zs(countersph)=R*(-sin(theta))*cos(phis)-L;
        ys(countersph)=R*sin(theta)*sin(phis);
        xs(countersph)=R*cos(theta);
        
        normsph1z(countersph)=-sin(theta)*cos(phis);
        normsph1y(countersph)=sin(theta)*sin(phis);
        normsph1x(countersph)=cos(theta);
        
        areasph(countersph) = (R^2)*((2*pi)/phi)*(cos(theta - (pi/(2*N_theta)))-cos(theta + (pi/(2*N_theta))));
        
        countersph=countersph+1;
    end
end

for jphi=1:N_phi/4
   for jtheta=1:N_theta
       
      theta = ((pi)/N_theta)*(jtheta - 0.5);
        phis = ((2*pi)/N_phi)*(jphi - 0.5);
        
        
        zs(countersph)=R*(-sin(theta))*cos(phis)-L;
        ys(countersph)=R*(-sin(theta))*sin(phis);
        xs(countersph)=R*cos(theta);
        
        normsph1z(countersph)=-sin(theta)*cos(phis);
        normsph1y(countersph)=-sin(theta)*sin(phis);
        normsph1x(countersph)=cos(theta);
        
        areasph(countersph) = (R^2)*((2*pi)/phi)*(cos(theta - (pi/(2*N_theta)))-cos(theta + (pi/(2*N_theta))));
        
        countersph=countersph+1;
    end
end

counter=counterc+countersph-2;
normVecSph=[;normsph1x;normsph1y;normsph1z];
normVec_cyl_sph=[normVeccyl normVecSph];
vecC=[(xc);(yc);(zc)];
vecsph=[(xs);(ys);(zs)];
x=[xc xs];
y=[yc ys];
z=[zc zs];
vecarea_cyl_sph=[areaCyl areasph];

vec_cyl_sph=[(x);(y);(z)];

plot3(y,z,x,'*')

AllInfo=[(vec_cyl_sph);(vecarea_cyl_sph);(normVec_cyl_sph)];


