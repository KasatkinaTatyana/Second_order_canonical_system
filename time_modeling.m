function [time, F] = time_modeling(Y_0)
% В данной функции производится интегрирование системы по времени
global control_arr

control_arr = sortrows(control_arr,1);
%% Время движения по траектории
y_0 = control_arr(1,1);
y_end = control_arr(1,2);
coefs = control_arr(1,3:6);
c0_glob = coefs(1);
c1_glob = coefs(2);
c2_glob = coefs(3);
c3_glob = coefs(4);
F_quad = @(x)1./(c0_glob + c1_glob*(x - Y_0(1)) + c2_glob*(x - Y_0(1)).^2 + c3_glob*(x - Y_0(1)).^3); % debug ?
tend = quad(F_quad,y_0,y_end);

global count
count = length(control_arr(:,1));
    
for i=2:count
    control_arr(i,10) = tend;
    
    y_0_i = control_arr(i,1);
    y_end_i = control_arr(i,2);
    coefs = control_arr(i,3:6);
    c0 = coefs(1);
    c1 = coefs(2);
    c2 = coefs(3);
    c3 = coefs(4);
    d = control_arr(i,7);
    y_0_c = control_arr(i,9);
    
    if (abs(control_arr(i,8) - 1) < 1e-6)
        F_quad_i = @(x)1./(c0 + c1*(x - y_0_i) + c2*(x - y_0_i).^2 + c3*(x - y_0_i).^3 + d*((x - y_0_i).^2).*(x - y_end_i).^2);
        t_i = quad(F_quad_i,y_0_i,y_end_i);
        tend = tend + t_i;
    elseif (abs(control_arr(i,8) - 2) < 1e-6)       
        F_quad_i = @(x)1./(c0 + c1*(x - y_0_c) + c2*(x - y_0_c).^2 + c3*(x - y_0_c).^3 +...
            d*((x - y_0_i)/(y_end_i - y_0_i)).^2.*(3 - 2*((x - y_0_i)/(y_end_i - y_0_i) )) );
        t_i = quad(F_quad_i,y_0_i,y_end_i);
        tend = tend + t_i;
    else
        F_quad_i = @(x)1./(c0 + c1*(x - y_0_c) + c2*(x - y_0_c).^2 + c3*(x - y_0_c).^3);
        t_i = quad(F_quad_i,y_0_i,y_end_i);
        tend = tend + t_i;
    end
end
%%

options = odeset('RelTol',1e-6,'AbsTol',[1e-6 1e-6]);
global Mas_u % массив управлений
[T, Y] = ode45(@kan_system, [0, tend], [Y_0(1) Y_0(2)], options);

figure();
hold on; grid on;
xlabel('y(t)'); ylabel('u(y(t))');
plot(Mas_u(:,1),Mas_u(:,2));
F = Y;
time = T;
