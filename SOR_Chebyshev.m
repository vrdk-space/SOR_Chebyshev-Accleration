%% Generalized Program for Estimation of Electric Potential and Field in Thunderstorms
% In this code, we define a charge structure within a cloud.
% Using that information, we solve for Electric Potential and Electric
% field due to the charge structure
clc
clear all

% Note : In the current code, Boundary Potential function (BPI) is written for cases : 
% i) L and h is same for All 3 Axes and ,
% ii) Lx = Ly & hx = hy. Lz and hz can be more than or less than Lx and hx 
% respectively (i.e) need not be the same. 

type = 1; % Type = 1 if Simulation Domain is Cube (L and h same in all 3 axis)
          % Type = 2 if Simulation Domain is Cuboid (L and h same only in X
          % and Y axis)
          
% All other aspects of the code, including SOR solver, is implemented in
% general terms. 

%% Define the arbitary charge configuration and dimension of cube and boundaries

%% Simulation Input

dim = 3;     % dim = 2 for 2D and dim = 3 for 3D
PEC = 1;     % If Ground is assumed to be perfect conductor, PEC = 1 (else 0) 
             % (For Earth, PEC = 1)
CDI = 2;     % CDI = 1 for Gaussian Charge Distribution, 
             % CDI = 2 for Simple Constant Charge Distribution

% Simulation Cube Length Along Z-Axis,X-Axis and Y-Axis respectively
L = [15e3;15e3;15e3]; % L = [Lz ; Lx ; Ly]

% Step Size Along Z-Axis, X-axis and Y-Axis respectively
h = [71.9;71.9;71.9];    % h = [hz ; hx ; hy]

                                            
z = 0:h(1):L(1);     
z = z';             % Note : If PEC = 1, z = 0 plane is assumed to be PEC 
                    % and code is written accordingly.
nz = size(z,1);

x = (-L(2)/2):h(2):(L(2)/2);
x = x';
nx = size(x,1);

if (dim == 2)
    y = 0;
    n = [nz;nx];
else
    y = (-L(3)/2):h(3):(L(3)/2);
    y = y';
    ny = size(y,1);
    n = [nz;nx;ny];
end

zeros_2D = zeros(nz,nx);
zeros_3D = zeros(nz,nx,ny);

% Cloud Charge Parameters : This is to define charge within the simulation
% cube. I have defined for a charge configuration expected in couds during
% thunderstorms. This can be redefined for your appropriate application.

if(CDI == 1)
    % Iudin et al, 2015 Parameters
    z_ccp = [4.5e3;6.7e3;9.6e3;11e3];             % Altitude in [m]
    t_ccp = [1e3;0.6e3;0.6e3;0.4e3];              % Thickness in [m]
    r_ccp = [2.6e3;3e3;3.5e3;4e3];                % Radius in [m]
    rho0 = [0.33e-9;-1.96e-9;0.95e-9;-0.58e-9];   % Max Charge Density in [C/m]
    gamma = 4;                                    % Transition Factor
end

if(CDI == 2)
    
    % Benkhe et al, 2005 Parameters
    z_ccp = [3.5e3;6e3;9e3];                     % Altitude in [m]
    t_ccp = [2e3;2e3;2e3];                       % Thickness in [m]
    r_ccp = [2e3;2e3;2.5e3];                     % Radius in [m]
    q_tot = [4.5;-45;31.5];                      % Total Charge [C]
    rho0 = zeros(size(z_ccp,1),1);
    Vol = zeros(size(z_ccp,1),1);
    for i = 1:size(z_ccp,1)
        Vol = pi*(r_ccp(i)^2)*t_ccp(i);
        rho0(i) = q_tot(i)/Vol;
    end
end
%% Solving Poisson's Equation to Determine Electric Potential

% SOR with Chebyshev Accleration is used to solve Poisson's Equation
% For theory, refer : http://people.eecs.berkeley.edu/~demmel/cs267/lecture24/lecture24.html
% For further details refer Chapter 20 in the text " Numerical Recipes"
epsilon_0 = 8.854e-12;
c1 = (h(1)^2)/epsilon_0;
rhs = c1.*rho;

% Step 1: Initialising Red and Black Points

Black = Chess(dim,n); % Points = 1 are Black Points and has only Red Neignbours
                      % Points = 0 are Red Points and has only Black Neighbours
                      
% Step 2: Boundary Potential Initializer

if(dim == 2)
    V_old = BPI(rho,dim,n,h,PEC,x,y,z,z_ccp,t_ccp,r_ccp,L,type);
    V_init = V_old;
    V_new = V_old;
end

if(dim == 3)
    V_old = BPI(rho,dim,n,h,PEC,x,y,z,z_ccp,t_ccp,r_ccp,L,type);
    V_init = V_old;
    V_new = V_old;
end

%% Step 3: Estimating Potential In Simulation Domain
start = 1;
f = 1;
loop = 0;
w = 1;
r_jac = cos(pi/n(2));

% i) 2-D Grid
if(dim == 2)
    h0 = ((h(1))^2 + (h(2))^2);
    k0 = (1/(2*8.854e-12))*(((h(1)*h(2))^2)/h0);
    rhs = rho.*k0;
    cz = (((h(2))^2)/h0)*0.5;
    cx = (((h(1))^2)/h0)*0.5;
    while(start)
      tic  
      loop = loop + 1  
      % Updating Black Points
      for i = 2:n(1)-1
          for j = 2:n(2)-1
              if(Black(i,j) == 1)
                  term1 = cz*(V_old(i+1,j) + V_old(i-1,j));
                  term2 = cx*(V_old(i,j+1) + V_old(i,j-1));
                  term3 = -1*V_old(i,j);
                  Residue = term1 + term2 + term3 + rhs(i,j,k);
                  
                  V_new(i,j) = V_old(i,j) + (w*Residue);                 
              end
          end
      end
      
      % Updating Red Points
      for i = 2:n(1)-1
          for j = 2:n(2)-1
              if(Black(i,j) == 0)
                  term1 = cz*(V_new(i+1,j) + V_new(i-1,j));
                  term2 = cx*(V_new(i,j+1) + V_new(i,j-1));
                  term3 = -1*V_old(i,j);
                  Residue = term1 + term2 + term3 + rhs(i,j,k);
                  
                  V_new(i,j) = V_old(i,j) + (w*Residue);                 
              end
          end
      end
      
      % Updating SOR Weight Parameter using Chebyshev Accleration Scheme
      % (Chap 20 of Numerical Recipes Textbook)
      if(loop == 1)
          w = 1/(1-(0.5*((r_jac)^2)));         
      end
      if(loop > 1)
          w = 1/(1-(0.25*w*((r_jac)^2)));
      end
      
      % Checking if desired Tolerance level reached
      R_max(f) = max(max(abs(V_new-V_old)));
      f = f+1;

      if(max(max((abs(V_new-V_old))))<=0.01)
          start = 0;
      else
          V_old = V_new;
      end
      toc  
    end    
end

% ii) 3-D Grid
if(dim == 3)
    h0 = ((h(1)*h(2))^2 + (h(2)*h(3))^2 + (h(1)*h(3))^2);
    k0 = (1/(2*8.854e-12))*(((h(1)*h(2)*h(3))^2)/h0);
    rhs = rho.*k0;
    cz = (((h(2)*h(3))^2)/h0)*0.5;
    cx = (((h(1)*h(3))^2)/h0)*0.5;
    cy = (((h(1)*h(2))^2)/h0)*0.5;
    
    while(start)
        
      tic  
      loop = loop + 1  
      % Updating Black Points
      for i = 2:n(1)-1
          for j = 2:n(2)-1
              for k = 2:n(3)-1
                  if(Black(i,j,k) == 1)
                      term1 = cz*(V_old(i+1,j,k) + V_old(i-1,j,k));
                      term2 = cx*(V_old(i,j+1,k) + V_old(i,j-1,k));
                      term3 = cy*(V_old(i,j,k+1) + V_old(i,j,k-1));
                      term4 = -1*V_old(i,j,k);
                      Residue = term1 + term2 + term3 + rhs(i,j,k) + term4;
                  
                      V_new(i,j,k) = V_old(i,j,k) + (w*Residue);  
                  end
              end
          end
      end
      
      if(loop == 1)
          w = 1/(1-(0.5*((r_jac)^2)));         
      end
      if(loop > 1)
          w = 1/(1-(0.25*w*((r_jac)^2)));
      end
      
      % Updating Red Points
      for i = 2:n(1)-1
          for j = 2:n(2)-1
              for k = 2:n(3)-1
                  if(Black(i,j,k) == 0)
                      term1 = cz*(V_new(i+1,j,k) + V_new(i-1,j,k));
                      term2 = cx*(V_new(i,j+1,k) + V_new(i,j-1,k));
                      term3 = cy*(V_new(i,j,k+1) + V_new(i,j,k-1));
                      term4 = -1*V_old(i,j,k);
                      Residue = term1 + term2 + term3 + rhs(i,j,k) + term4;
                  
                      V_new(i,j,k) = V_old(i,j,k) + (w*Residue);  
                  end
              end
          end
      end
      
      % Updating SOR Weight Parameter using Chebyshev Accleration Scheme
      % (Chap 20 of Numerical Recipes Textbook)
      if(loop > 1)
          w = 1/(1-(0.25*w*((r_jac)^2)));
      end
      
      % Checking if desired Tolerance level reached
      R_max(f) = max(max(max(abs(V_new-V_old))));
      f = f+1;

      if(max(max(max(abs(V_new-V_old))))<=0.01)
          start = 0;
      else
          V_old = V_new;
      end
      toc  
        
    end
end

%% Calculating Electric Field
tic
if(dim == 2)
    Ex = zeros_2D;
    Ez = zeros_2D;
    E_tot = zeros_2D;  

else
    Ex = zeros_3D;
    Ey = zeros_3D;
    Ez = zeros_3D;
    E_tot = zeros_3D;  
end
count = 1;

% a) 2-D
if(dim == 2)
    for i = 2:n(1)-1
        for j = 2:n(2)-1
            Ez(i,j) = (V_new(i-1,j)-V_new(i+1,j))/(2*h(1));
            Ex(i,j) = (V_new(i,j-1)-V_new(i,j+1))/(2*h(2));
            E_tot(i,j) = sqrt((Ez(i,j)^2)+(Ex(i,j)^2));
            if(E_tot(i,j) >= E_be(i))
                z_start(count) = i;
                x_start(count) = j;
                count = count + 1;
            end
        end
    end
end
% b) 3-D
if(dim == 3)
    for i = 2:n(1)-1
        for j = 2:n(2)-1
            for k = 2:n(3)-1  
                
                Ez(i,j,k) = (V_new(i-1,j,k)-V_new(i+1,j,k))/(2*h(1));
                Ex(i,j,k) = (V_new(i,j-1,k)-V_new(i,j+1,k))/(2*h(2));
                Ey(i,j,k) = (V_new(i,j,k-1)-V_new(i,j,k+1))/(2*h(3));                
                E_tot(i,j,k) = sqrt((Ez(i,j,k)^2)+(Ex(i,j,k)^2)+(Ey(i,j,k)^2));
                
                if(E_tot(i,j,k) >= E_be(i))
                    z_start(count) = i;
                    x_start(count) = j;
                    y_start(count) = k;
                    count = count + 1;
                end
                               
            end
        end
    end
end 
toc

