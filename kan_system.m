function dy = kan_system(t,y)
dy = zeros(2,1);

global control_arr
global count

global y_prog Psi_prog dPsi_prog
global dt


flag = false;
% Впоследствии надо сделать от count до 2
for i=2:count
    y_0_i = control_arr(i,1);
    y_end_i = control_arr(i,2);
    
    if ((y(1) < y_end_i)&&(y(1) > y_0_i)) || (abs(y(1) - y_0_i) < 1e-6)
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
            u = Psi*dPsi - sin(y(1));
             
        elseif (abs(control_arr(i,8) - 2) < 1e-6)
            y_0_glob = control_arr(i,9);
            t_0 = control_arr(i,10);
            Psi = c0 + c1*(y(1) - y_0_glob) + c2*(y(1) - y_0_glob)^2 + c3*(y(1) - y_0_glob)^3 +...
                d*((y(1) - y_0_i)/(y_end_i - y_0_i))^2*(3 - 2*((y(1) - y_0_i)/(y_end_i - y_0_i) ));
            
            dPsi = c1 + 2*c2*(y(1) - y_0_glob) + 3*c3*(y(1) - y_0_glob)^2 + ...
                d*(6*((y(1) - y_0_i)/(y_end_i - y_0_i)) - 6*((y(1) - y_0_i)/(y_end_i - y_0_i))^2);
            r1 = -0.01; r2 = -10;
            c1_st = -(r1 + r2);
            c0_st = r1*r2;

            ind = floor((t - t_0)/dt) + 1;
            
            if (ind > length(Psi_prog))
                ind = length(Psi_prog);
            end
            
            u = Psi_prog(ind)*dPsi_prog(ind) - sin(y(1)) - c1_st*(y(2) - Psi_prog(ind)) - c0_st*(y(1) - y_prog(ind));

        else
            Psi = c0 + c1*(y(1) - y_0_i) + c2*(y(1) - y_0_i)^2 + c3*(y(1) - y_0_i)^3;
            dPsi = c1 + 2*c2*(y(1) - y_0_i) + 3*c3*(y(1) - y_0_i)^2;
            u = Psi*dPsi - sin(y(1));
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
    u = Psi*dPsi - sin(y(1));
end


plot(y(1), u,'Marker','*');

dy(1) = y(2);
dy(2) = sin(y(1)) + u;