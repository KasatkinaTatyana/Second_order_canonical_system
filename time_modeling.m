function F = time_modeling(Y_0)
global control_arr
%% Время движения по траектории
y_0 = control_arr(1,1);
y_end = control_arr(1,2);
coefs = control_arr(1,3:6);
c0_glob = coefs(1);
c1_glob = coefs(2);
c2_glob = coefs(3);
c3_glob = coefs(4);
F_quad = @(x)1./(c0_glob + c1_glob*(x - Y_0(1)) + c2_glob*(x - Y_0(1)).^2 + c3_glob*(x - Y_0(1)).^3);
tend = quad(F_quad,y_0,y_end);

global count
count = length(control_arr(:,1));

   
N = 1000000;    
global y_prog Psi_prog dPsi_prog
global dt
    

for i=2:count
    y_0_i = control_arr(i,1);
    y_end_i = control_arr(i,2);
    coefs = control_arr(i,3:6);
    c0 = coefs(1);
    c1 = coefs(2);
    c2 = coefs(3);
    c3 = coefs(4);
    d = control_arr(i,7);
    
    if (abs(control_arr(i,8) - 1) < 1e-6)
        F_quad_i = @(x)1./(c0 + c1*(x - y_0_i) + c2*(x - y_0_i).^2 + c3*(x - y_0_i).^3 + d*((x - y_0_i).^2).*(x - y_end_i).^2);
        t_i = quad(F_quad_i,y_0_i,y_end_i);
        t_old = quad(F_quad,y_0_i,y_end_i);
        tend = tend - t_old + t_i;
    elseif (abs(control_arr(i,8) - 2) < 1e-6)
        % Определяем интервал времени, который является начальным для этого
        % фрагмента кривой
        t_0 = quad(F_quad,y_0,y_0_i);
        control_arr(i,10) = t_0;
        for j=2:count
            if (i ~= j)&&(control_arr(j,2) < control_arr(i,1))
                if (abs(control_arr(j,8) - 1) < 1e-6)
                    y_0_j = control_arr(j,1);
                    y_end_j = control_arr(j,2);
                    coefs_j = control_arr(j,3:6);
                    c0_j = coefs_j(1);
                    c1_j = coefs_j(2);
                    c2_j = coefs_j(3);
                    c3_j = coefs_j(4);
                    d_j = control_arr(j,7);
                    
                    F_quad_j = @(x)1./(c0_j + c1_j*(x - y_0_j) + c2_j*(x - y_0_j).^2 + c3_j*(x - y_0_j).^3 + ...
                        d_j*((x - y_0_j).^2).*(x - y_end_j).^2);
                    t_j = quad(F_quad_j,y_0_j,y_end_j);
                    control_arr(i,10) = t_0 - quad(F_quad,y_0_j,y_end_j) + t_j;
                end
            end
        end
        
        y_0_glob = control_arr(i,9);
        F_quad_i = @(x)1./(c0 + c1*(x - y_0_glob) + c2*(x - y_0_glob).^2 + c3*(x - y_0_glob).^3 +...
            d*((x - y_0_i)/(y_end_i - y_0_i)).^2.*(3 - 2*((x - y_0_i)/(y_end_i - y_0_i) )) );
        t_i = quad(F_quad_i,y_0_i,y_end_i);
        t_old = quad(F_quad,y_0_i,y_end_i);
        tend = tend - t_old + t_i;
        
        dt = t_i/N;
        dy_prog = (y_end_i - y_0_i) / N;
        
        y_prog = y_0_i:dy_prog:y_end_i;
        Psi_prog = c0 + c1*(y_prog - y_0_glob) + c2*(y_prog - y_0_glob).^2 + c3*(y_prog - y_0_glob).^3 +...
            d*((y_prog - y_0_i)/(y_end_i - y_0_i)).^2.*(3 - 2*((y_prog - y_0_i)/(y_end_i - y_0_i) ));
        dPsi_prog = c1 + 2*c2*(y_prog - y_0_glob) + 3*c3*(y_prog - y_0_glob).^2 + ...
            d*(6*((y_prog - y_0_i)/(y_end_i - y_0_i)) - 6*((y_prog - y_0_i)/(y_end_i - y_0_i)).^2);
    else
        F_quad_i = @(x)1./(c0 + c1*(x - y_0_i) + c2*(x - y_0_i).^2 + c3*(x - y_0_i).^3);
        t_i = quad(F_quad_i,y_0_i,y_end_i);
        t_old = quad(F_quad,y_0_i,y_end_i);
        tend = tend - t_old + t_i;
    end
end
%%

options = odeset('RelTol',1e-6,'AbsTol',[1e-6 1e-6]);
[T, Y] = ode45(@kan_system, [0, tend], [Y_0(1) Y_0(2)], options);
F = Y;
