function Vb = BPI(rho,dim,n,h,PEC,x,y,z,z_ccp,t_ccp,r_ccp,L,type) 
k_elec = (1/(4*pi*8.854e-12));
Vb = 0;
%% 2-D
if(dim == 2)
   Vb = zeros(n(1),n(2));
   Area = h(1)*h(2);
   q = rho.*Area;
   q_imag = -q;
   
   %% Top and Bottom Surface
   %% Side Surfaces
   
end

%% 3-D
if(dim == 3)
    if(type == 1)
        Vb = zeros(n(1),n(2),n(3)); 
        Vol = h(1)*h(2)*h(3);
        q = rho.*Vol;
        q_imag = -q;  
        loop_b = 0;
        for b1 = 1:n(1)
            for b2 = 1:ceil(n(1)/2)
                tic
                loop_b = loop_b + 1
                % z = Ground Surface
                z_z1 = z(1); x_z1 = x(b1); y_z1 = y(b2);
                Vz1 = 0;
            
                % z = +L surface
                z_zn = z(n(1)); x_zn = x(b1); y_zn = y(b2);
                Vzn = 0;
      
                % x = -L/2 Surface
                z_x1 = z(b1); x_x1 = x(1); y_x1 = y(b2);
                Vx1 = 0;
            
                for m = 1:size(z_ccp,1)                
                    temp1 = ceil(((z_ccp(m)/h(1))+1));
                    temp2 = ceil(((z_ccp(m)+t_ccp(m))/h(1))+1);
                    temp3 = ceil((r_ccp(m)/h(1))+1);
                    temp4 = n(1)-temp3+1;
                
                    for a = temp1:temp2
                        for b = temp3:temp4
                            for c =temp3:temp4
                                z2 = z(a);  x2 = x(b);  y2 = y(c);
                                % z = Ground Surface
                                if (z(1) ~= 0)
                                    z_temp1 = (z_z1-z2)^2; x_temp1 = (x_z1-x2)^2; y_temp1 = (y_z1-y2)^2;
                                    dist1 = sqrt(z_temp1 + x_temp1 + y_temp1);
                                    Vz1 = Vz1 + ((k_elec*q(a,b,c))/dist1);
                                    %Vz1 = Vz1 + ((k_elec*rho0(m)*h*h*h)/dist1);
                                end
                            
                                % z = +L Surface
                                z_temp2 = (z_zn-z2)^2; x_temp2 = (x_zn-x2)^2; y_temp2 = (y_zn-y2)^2;
                                dist2 = sqrt(z_temp2 + x_temp2 + y_temp2);
                                Vzn = Vzn + ((k_elec*q(a,b,c))/dist2);
                                %Vzn = Vzn + ((k_elec*rho0(m)*h*h*h)/dist2);
                            
                                % x = -L/2 Surface
                                z_temp3 = (z_x1-z2)^2; x_temp3 = (x_x1-x2)^2; y_temp3 = (y_x1-y2)^2;
                                dist3 = sqrt(z_temp3 + x_temp3 +y_temp3);
                                Vx1 = Vx1 + ((k_elec*q(a,b,c))/dist3);
                                %Vx1 = Vx1 + ((k_elec*rho0(m)*h*h*h)/dist3);
                            
                                if (PEC == 1)
                                    % z = Lower Limit
                                    if(z(1) ~=0)
                                        z_temp1_imag = (z_z1 + z2)^2;
                                        dist1_imag = sqrt(z_temp1_imag + x_temp1 + y_temp1);
                                        Vz1 = Vz1 + ((k_elec*q_imag(a,b,c))/dist1_imag);
                                        %Vz1 = Vz1 + ((k_elec*rho0(m)*h*h*h)/dist1_imag);
                                    end
                                    % z = +L Surface
                                    z_temp2_imag = (z_zn + z2)^2;
                                    dist2_imag = sqrt(z_temp2_imag + x_temp2 + y_temp2);
                                    Vzn = Vzn + ((k_elec*q_imag(a,b,c))/dist2_imag);
                                    %Vzn = Vzn + ((k_elec*rho0(m)*h*h*h)/dist2_imag);
                                
                                    % x = -L/2 Surface
                                    z_temp3_imag = (z_x1 + z2)^2;
                                    dist3_imag = sqrt(z_temp3_imag + x_temp3 +y_temp3);
                                    Vx1 = Vx1 + ((k_elec*q_imag(a,b,c))/dist3_imag);
                                    %Vx1 = Vx1 + ((k_elec*rho0(m)*h*h*h)/dist3_imag);
                                
                                end
                            
                            end
                        end
                    end                
                end
            
                % z = Ground
                Vb(1,b1,b2) = Vz1;
                Vb(1,b1,n(3)-b2+1) = Vz1;
            
                % z = +L surface
                Vb(n(1),b1,b2) = Vzn;
                Vb(n(1),b1,n(3)-b2+1) = Vzn;
            
                % x = -L/2 Surface
                Vb(b1,1,b2) = Vx1;
                Vb(b1,1,n(3)-b2+1) = Vx1;
                % x = L/2 Surface
                Vb(b1,n(2),b2) = Vx1;
                Vb(b1,n(2),n(3)-b2+1) = Vx1;
                % y = -L/2 surface
                Vb(b1,b2,1) = Vx1;
                Vb(b1,n(2)-b2+1,1) = Vx1;
                % y = L/2 Surface
                Vb(b1,b2,n(3)) = Vx1;
                Vb(b1,n(2)-b2+1,n(3)) = Vx1;
                toc
            end
        end
    end
    
    if(type == 2)
        Vb = zeros(n(1),n(2),n(3)); 
        Vol = h(1)*h(2)*h(3);
        q = rho.*Vol;
        q_imag = -q;  
        loop_b = 0;
        %% Top and Bottom Surface
        for b1 = 1:n(2)
            for b2 = 1:ceil(n(3)/2)
                tic
                loop_b = loop_b + 1
                % Bottom Surface
                z_z1 = z(1); x_z1 = x(b1); y_z1 = y(b2);
                Vz1 = 0;

                % Top Surface
                z_zn = z(n(1)); x_zn = x(b1); y_zn = y(b2);
                Vzn = 0;

                for m = 1:size(z_ccp,1)
                    temp1 = ceil(((z_ccp(m)/h(1))+1));
                    temp2 = ceil(((z_ccp(m)+t_ccp(m))/h(1))+1);
                    temp3 = ceil((((L(2)/2)-r_ccp(m))/h(2))+1);
                    temp4 = n(2)-temp3+1;

                    for a = temp1:temp2
                        for b = temp3:temp4
                            for c =temp3:temp4

                                z2 = z(a);  x2 = x(b);  y2 = y(c);

                                % Bottom Surface
                                if (z(1) ~= 0)
                                    z_temp1 = (z_z1-z2)^2; x_temp1 = (x_z1-x2)^2; y_temp1 = (y_z1-y2)^2;
                                    dist1 = sqrt(z_temp1 + x_temp1 + y_temp1);
                                    Vz1 = Vz1 + ((k_elec*q(a,b,c))/dist1);
                                end

                                % Top Surface
                                z_temp2 = (z_zn-z2)^2; x_temp2 = (x_zn-x2)^2; y_temp2 = (y_zn-y2)^2;
                                dist2 = sqrt(z_temp2 + x_temp2 + y_temp2);
                                Vzn = Vzn + ((k_elec*q(a,b,c))/dist2);

                                if (PEC == 1)
                                    % Bottom Surface
                                    if(z(1) ~=0)
                                        z_temp1_imag = (z_z1 + z2)^2;
                                        dist1_imag = sqrt(z_temp1_imag + x_temp1 + y_temp1);
                                        Vz1 = Vz1 + ((k_elec*q_imag(a,b,c))/dist1_imag);                                    
                                    end

                                    % Top Surface
                                    z_temp2_imag = (z_zn + z2)^2;
                                    dist2_imag = sqrt(z_temp2_imag + x_temp2 + y_temp2);
                                    Vzn = Vzn + ((k_elec*q_imag(a,b,c))/dist2_imag);
                                end
                           end
                       end
                   end                   
               end

               % Bottom Surface
               Vb(1,b1,b2) = Vz1;
               Vb(1,b1,n(3)-b2+1) = Vz1;

               % Top Surface
               Vb(n,b1,b2) = Vzn;
               Vb(n,b1,n(3)-b2+1) = Vzn;   
               toc
           end
       end

       %% Side Surface
       for b1 = 1:n(1)
           for b2 = 1:ceil(n(2)/2)
               tic
               loop_b = loop_b + 1
               % x = -L/2 Surface
               z_x1 = z(b1); x_x1 = x(1); y_x1 = y(b2);
               Vx1 = 0;
               for m = 1:size(z_ccp,1)
                   temp1 = ceil(((z_ccp(m)/h(1))+1));
                   temp2 = ceil(((z_ccp(m)+t_ccp(m))/h(1))+1);
                   temp3 = ceil((((L(2)/2)-r_ccp(m))/h(2))+1);
                   temp4 = n(2)-temp3+1;              
                   for a = temp1:temp2
                       for b = temp3:temp4
                           for c =temp3:temp4

                               z2 = z(a);  x2 = x(b);  y2 = y(c);
                               % x = -L/2 Surface
                               z_temp3 = (z_x1-z2)^2; x_temp3 = (x_x1-x2)^2; y_temp3 = (y_x1-y2)^2;
                               dist3 = sqrt(z_temp3 + x_temp3 +y_temp3);
                               Vx1 = Vx1 + ((k_elec*q(a,b,c))/dist3);
                               if(PEC == 1)
                                   % x = -L/2 Surface
                                   z_temp3_imag = (z_x1 + z2)^2;
                                   dist3_imag = sqrt(z_temp3_imag + x_temp3 +y_temp3);
                                   Vx1 = Vx1 + ((k_elec*q_imag(a,b,c))/dist3_imag);
                               end
                           end
                       end
                   end
               end
               % x = -L/2 Surface
               Vb(b1,1,b2) = Vx1;
               Vb(b1,1,n(3)-b2+1) = Vx1;
               % x = L/2 Surface
               Vb(b1,n(2),b2) = Vx1;
               Vb(b1,n(2),n(3)-b2+1) = Vx1;
               % y = -L/2 surface
               Vb(b1,b2,1) = Vx1;
               Vb(b1,n(2)-b2+1,1) = Vx1;
               % y = L/2 Surface
               Vb(b1,b2,n(3)) = Vx1;
               Vb(b1,n(2)-b2+1,n(3)) = Vx1;
               toc
           end
       end
    end
end

end