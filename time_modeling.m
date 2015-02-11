function F = time_modeling(Y_0)
global control_arr
%% Время движения по траектории
y_0 = control_arr(1,1);
y_end = control_arr(1,2);
coefs = control_arr(1,3:6);
c0 = coefs(1);
c1 = coefs(2);
c2 = coefs(3);
c3 = coefs(4);
F_quad = @(x)1./(c0 + c1*(x - Y_0(1)) + c2*(x - Y_0(1)).^2 + c3*(x - Y_0(1)).^3);
tend = quad(F_quad,y_0,y_end);

global count
count = length(control_arr(:,1));

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
        F_quad_i = @(x)1./(c0 + c1*(x - y_0_i) + c2*(x - y_0_i).^2 + c3*(x - y_0_i).^3 + d*(x - y_0_i).^2*(x - y_end_i).^2);
        t_i = quad(F_quad_i,y_0_i,y_end_i);
        t_old = quad(F_quad,y_0_i,y_end_i);
        tend = tend - t_old + t_i;
    elseif (abs(control_arr(i,8) - 2) < 1e-6)
        y_0_glob = control_arr(i,9);
        F_quad_i = @(x)1./(c0 + c1*(x - y_0_glob) + c2*(x - y_0_glob).^2 + c3*(x - y_0_glob).^3 +...
            d*((x - y_0_i)/(y_end_i - y_0_i)).^2.*(3 - 2*((x - y_0_i)/(y_end_i - y_0_i) )) );
        t_i = quad(F_quad_i,y_0_i,y_end_i);
        t_old = quad(F_quad,y_0_i,y_end_i);
        tend = tend - t_old + t_i;
    end
end
%%
[T, Y] = ode45(@kan_system, [0, tend], [Y_0(1) Y_0(2)]);
F = Y;
