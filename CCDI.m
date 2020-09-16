function [rho_temp,q_temp] = CCDI(z_ccp,t_ccp,r_ccp,rho0,gamma,x,y,z,n,h,dim,CDI)

if(dim == 2)
    rho_temp = zeros(n(1),n(2));
    q_temp = zeros(size(z_ccp,1),1);
    var3 = 1;
else
    rho_temp = zeros(n(1),n(2),n(3));
    var3 = n(3);
    q_temp = zeros(size(z_ccp,1),1);
end

for u = 1:n(1)
    for v = 1:n(2)
        for w = 1:var3
            
            % Initialize a Dummy rho for each iteration
            rho_temp1 = 0;
            
            % Determine Axial Distance of Point Z - Axis
            if(dim == 2)
                r = abs(x(v));
            else
                r = (x(v)^2) + (y(w)^2);
                r = sqrt(r);
            end
            %% Gaussian Charge Distribution
            if (CDI == 1)
                for m = 1:size(z_ccp,1)
                    if(z(u) >= z_ccp(m) && z(u) <= (z_ccp(m) + t_ccp(m)))
                        if(r <= r_ccp(m))
                            term1 = (z(u) - z_ccp(m))/t_ccp(m);
                            term1 = term1^(2*gamma);
                            term2 = r/r_ccp(m);
                            term2 = term2^(2*gamma);
                            rho_temp1 = rho0(m)*exp(-term1-term2);
                            if(dim == 2)
                                q_temp(m) = q_temp(m)+ (rho_temp1*h(1)*h(2));
                            end
                            if(dim ==3)
                                q_temp(m) = q_temp(m)+ (rho_temp1*h(1)*h(2)*h(3));
                            end
                        end
                    end
                end
            end  
            
            %% Simple Constant Charge Distribution
            if (CDI == 2)
                for m = 1:size(z_ccp,1)
                    if(z(u) >= z_ccp(m) && z(u) <= (z_ccp(m) + t_ccp(m)))
                        if(r <= r_ccp(m))
                            rho_temp1 = rho0(m);
                            if(dim == 2)
                                q_temp(m) = q_temp(m)+ (rho_temp1*h(1)*h(2));
                            end
                            if(dim ==3)
                                q_temp(m) = q_temp(m)+ (rho_temp1*h(1)*h(2)*h(3));
                            end
                        end
                    end
                end
            end
            
            if(dim == 2)
                rho_temp(u,v) = rho_temp1;
            else
                rho_temp(u,v,w) = rho_temp1;
            end
            
        end
    end
end

end