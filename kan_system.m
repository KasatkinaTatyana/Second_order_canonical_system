function dy = kan_system(t,y)
dy = zeros(2,1);

global control_arr
global count

global c0 c1 c2 c3 y_0_c y_0_i y_end_i d % объявлены глобальными, чтобы интегрировать систему уравнение dy /dt = Psi(y) по времени
global Mas_u

flag = false;

% Для каждого типа кривой (8 координата в строках массива control_arr)
% управление строится по своему

for i=1:count
    y_0_i = control_arr(i,1);
    y_end_i = control_arr(i,2);
    
    if ((y(1) < y_end_i)&&(y(1) > y_0_i)) || (abs(y(1) - y_0_i) < 1e-6)
        coefs = control_arr(i,3:6);
        c0 = coefs(1);
        c1 = coefs(2);
        c2 = coefs(3);
        c3 = coefs(4);
        d = control_arr(i,7);
        y_0_c = control_arr(i,9);
        t_0 = control_arr(i,10);
        
        if (abs(control_arr(i,8) - 1) < 1e-6)
            [T, Y] = ode45(@Trajectory_type_1,[t_0 t], y_0_i);
            y_t = Y(end);
            
            Psi_t = c0 + c1*(y_t - y_0_i) + c2*(y_t - y_0_i)^2 + c3*(y_t - y_0_i)^3 + d*(y_t - y_0_i)^2*(y_t - y_end_i)^2;
            
            dPsi_t = c1 + 2*c2*(y_t - y_0_i) + 3*c3*(y_t - y_0_i)^2 + ...
                2*d*(y_t - y_0_i)^2*(y_t - y_end_i) + 2*d*(y_t - y_0_i)*(y_t - y_end_i)^2;
            u = Psi_t*dPsi_t - sin(y(1));
             
        elseif (abs(control_arr(i,8) - 2) < 1e-6)           
            [T, Y] = ode45(@Trajectory_type_2,[t_0 t], y_0_i);
            y_t = Y(end);
            
            Psi_t = c0 + c1*(y_t - y_0_c) + c2*(y_t - y_0_c)^2 + c3*(y_t - y_0_c)^3 +...
                d*((y_t - y_0_i)/(y_end_i - y_0_i))^2*(3 - 2*((y_t - y_0_i)/(y_end_i - y_0_i) ));
            
            dPsi_t = c1 + 2*c2*(y_t - y_0_c) + 3*c3*(y_t - y_0_c)^2 + ...
                d*(6*(y_t - y_0_i)/(y_end_i - y_0_i)^2 - 6*(y_t - y_0_i)^2/(y_end_i - y_0_i)^3);
            r1 = -1; r2 = -0.1;
            c1_st = -(r1 + r2);
            c0_st = r1*r2;

            u = Psi_t*dPsi_t - sin(y(1)) - c1_st*(y(2) - Psi_t) - c0_st*(y(1) - y_t);
        else         
            if (t_0 == t)
                y_t = y_0_i;
            else
                [T, Y] = ode45(@Trajectory_type_0,[t_0 t], y_0_i);
                y_t = Y(end);
            end
            Psi_t = c0 + c1*(y_t - y_0_c) + c2*(y_t - y_0_c)^2 + c3*(y_t - y_0_c)^3;
            dPsi_t = c1 + 2*c2*(y_t - y_0_c) + 3*c3*(y_t - y_0_c)^2;
            u = Psi_t*dPsi_t - sin(y(1));
        end 
        flag = true;
        break;
    end
    
end

% последняя точка ?
if (flag == false)
    y_0 = control_arr(count,1);
    
    coefs = control_arr(count,3:6);
    c0 = coefs(1);
    c1 = coefs(2);
    c2 = coefs(3);
    c3 = coefs(4);

    y_0_c = control_arr(count,9);
    t_0 = control_arr(count,10);
    
    if (t_0 == t)
        y_t = y_0_i;
    else
        [T, Y] = ode45(@Trajectory_type_0,[t_0 t], y_0);
        y_t = Y(end);
    end
    Psi_t = c0 + c1*(y_t - y_0_c) + c2*(y_t - y_0_c)^2 + c3*(y_t - y_0_c)^3;
    dPsi_t = c1 + 2*c2*(y_t - y_0_c) + 3*c3*(y_t - y_0_c)^2;
    u = Psi_t*dPsi_t - sin(y(1));
end

plot(y(1), u,'Marker','*');

Mas_u = [Mas_u; y(1) u];

dy(1) = y(2);
dy(2) = sin(y(1)) + u;