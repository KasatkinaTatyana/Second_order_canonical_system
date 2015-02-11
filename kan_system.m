function dy = kan_system(t,y)
dy = zeros(2,1);

global control_arr
global count


flag = false;
for i=2:count
    y_0_i = control_arr(i,1);
    y_end_i = control_arr(i,2);
    
    if ((y(1) < y_end_i)&&(y(1) > y_i_0)) || (abs(y(1) - y_0_i) < 1e-6)
        coefs = control_arr(i,3:6);
        c0 = coefs(1);
        c1 = coefs(2);
        c2 = coefs(3);
        c3 = coefs(4);
        d = control_arr(i,7);
        
        if (abs(control_arr(i,8) - 1) < 1e-6)
            Psi = c0 + c1*(y(1) - y_0_i) + c2*(y(1) - y_0_i)^2 + c3*(y(1) - y_0_i)^3 + d*(y(1) - y_0_i)^2*(y(1) - y_end_i)^2;
            
            dPsi = c1 + 2*c2*(y(1) - y_0_i) + 3*c3*(y(1) - y_0_i)^2 + ...
                2*d*(y(1) - y_0_i)^2*(y(1) - y_end_i) + 2*d*(y(1) - y_0_i)*(y(1) - y_end_i)^2;
             
        elseif (abs(control_arr(i,8) - 2) < 1e-6)
            y_0_glob = control_arr(i,9);
            Psi = c0 + c1*(y(1) - y_0_glob) + c2*(y(1) - y_0_glob)^2 + c3*(y(1) - y_0_glob)^3 +...
                d*((y(1) - y_0_i)/(y_end_i - y_0_i))^2*(3 - 2*((y(1) - y_0_i)/(y_end_i - y_0_i) ));
            
            dPsi = c1 + 2*c2*(y(1) - y_0_glob) + 3*c3*(y(1) - y_0_glob)^2 + ...
                d*(6*((y(1) - y_0_i)/(y_end_i - y_0_i)) - 6*((y(1) - y_0_i)/(y_end_i - y_0_i))^2);
        end
        
        flag = true;
    end
    
    if (flag == true)
        break;
    end
end

if (flag == false)
    y_0 = control_arr(1,1);
    coefs = control_arr(1,3:6);
    c0 = coefs(1);
    c1 = coefs(2);
    c2 = coefs(3);
    c3 = coefs(4);
    Psi = c0 + c1*(y(1) - y_0) + c2*(y(1) - y_0)^2 + c3*(y(1) - y_0)^3;
    dPsi = c1 + 2*c2*(y(1) - y_0) + 3*c3*(y(1) - y_0)^2;
end

u = Psi*dPsi - sin(y(1));

dy(1) = y(2);
dy(2) = sin(y(1)) + u;