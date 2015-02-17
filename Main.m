clear all
close all
clc
%% Ограничения
%  constr1 < Psi(y) < constr2
%  dy_constr1 < dPsi(y) / dy < dy_constr2

y0 = pi / 4;
dy0 = 0.1;
yend = 3*pi/4;
dyend = 0.5;
u0 = 0;
uend = 0;
ddy0 = u0 + sin(y0);
ddyend = uend + sin(yend);


global constr1 constr2
constr1 = -0.85;
constr2 = 1.6;
dy_constr1 = -1.8;
dy_constr2 = 0.5;
%%
% global Y_0 Y_end c0 c1 c2 c3 % debug
Y_0 =   [y0    dy0    ddy0];
Y_end = [yend  dyend  ddyend];

N = 1000;

% Строится фазовая кривая вида
% Psi = c0 + c1*(y - y0) + c2*(y - y0)^2 + c3*(y - y0)^3 + d*(y - y0)^2*(y - yend)^2,
% соединяющая
% начальное положение Y_0 и конечное положение Y_end
delta = Y_end(1) - Y_0(1);

c0 = Y_0(2);
c1 = Y_0(3)/Y_0(2);

Ma = [delta^2     delta^3;
      2*delta     3*delta^2];
  
Mb = [Y_end(2) - c1*delta - c0;
      Y_end(3)/Y_end(2) - c1];
Mc = inv(Ma)*Mb;

c2 = Mc(1);
c3 = Mc(2);

Expr = c0 + c1*delta + c2*delta^2 + c3*delta^3;

Expr_1 = c1 + 2*c2*delta + 3*c3*delta^2;

d = 0;
%%
figure(1);
hold on; grid on;
title('Psi(y)');
xlabel('y');
ylabel('dy / d\tau');

dy = (Y_end(1) - Y_0(1)) / N;
y=Y_0(1):dy:Y_end(1);

Psi = c0 + c1*(y - Y_0(1)) + c2*(y - Y_0(1)).^2 + c3*(y - Y_0(1)).^3;
dPsi = c1 + 2*c2*(y - Y_0(1)) + 3*c3*(y - Y_0(1)).^2;

plot(y,Psi);

figure(2);
hold on; grid on;
xlabel('y');
ylabel('dy / d\tau');
title('dPsi(y) / dy');
plot(y,dPsi,'b');

figure(3);
hold on; grid on;
xlabel('y');
ylabel('dy / d\tau');
title('Управление u(y)');
u = Psi.*dPsi - sin(y);

u_t = u; % управление, которое используется для моделирования

plot(y,u);

%% debug
% r1 = -0.1; r2 = -0.1;
% c1_st = -(r1 + r2);
% c0_st = r1*r2;
% 
% [T, Y] = ode45(@kan_system, [0, tend], [y0 dy0]);
% 
% set(0,'CurrentFigure',1);
% plot(Y(:,1),Y(:,2));
%% Массив, в котором будут содержаться коэффициенты для управлений
global control_arr
control_arr = [Y_0(1) Y_end(1) c0 c1 c2 c3 0 0 Y_0(1) 0];
% четвертая цифра с конца - значение параметра d, третья цифра с конца - тип кривой
% вторая цифра  - y_0_c - это значение, которое вычитается из y это может быть Y_0(1) или Y_0(1) из
% метода lower/upper_constraint_dPsi; последняя цифра - если тип кривой
% равен 2, то вычисляется сколько времени прошло с начала, до момента
% прихода в начальную точку текущего фрагмента, это значение вычисляется в
% time_modeling.
% типы будут следующие:
% 0 - кривая вида Psi(y) = c0 + c1*(y - y_0) + c2*(y - y_0)^2 + c3*(y - y_0)^3;
% 1 - кривая вида Psi(y) = c0 + c1*(y - y_0) + c2*(y - y_0)^2 + c3*(y - y_0)^3 + d*(y - y0)^2*(y - y_end)^2
% 2 - кривая вида Psi(y) = c0 + c1*(y - y_0_glob) + c2*(y - y_0_glob)^2 + c3*(y -y_0_glob)^3 + d*y_norm.^2.*(3 - 2*y_norm), 
% где y_norm = (y - y_0) / (y_end - y_0);
% первый элемент этой структуры

%% Убираю выход за ограничение dPsi / dy < dy_constr2
% [data1, data2] = upper_constraint_dPsi(Y_0, Y_end, dy, dy_constr2);
% y_1 = data1(1,:);
% Psi_1 = data1(2,:);
% dPsi_1 = data1(3,:);
% y_2 = data2(1,:);
% Psi_2 = data2(2,:);
% dPsi_2 = data2(3,:);
% 
% set(0,'CurrentFigure',1);
% plot(y_1, Psi_1, 'g','LineWidth',2);
% plot(y_2, Psi_2, 'm','LineWidth',2);
% 
% set(0,'CurrentFigure',2);
% plot(y_1, dPsi_1, 'g','LineWidth',2);
% plot(y_2, dPsi_2, 'm','LineWidth',2);

% %% Убираю выход за ограничение dPsi / dy > dy_constr1
[data1, data2, par, coefs] = lower_constraint_dPsi(Y_0, Y_end, dy, dy_constr1);
y_1 = data1(1,:);
Psi_1 = data1(2,:);
dPsi_1 = data1(3,:);
u_1 = Psi_1.*dPsi_1 - sin(y_1);
y_2 = data2(1,:);
Psi_2 = data2(2,:);
dPsi_2 = data2(3,:);
u_2 = Psi_2.*dPsi_2 - sin(y_2);

set(0,'CurrentFigure',1);
plot(y_1, Psi_1, 'g');
plot(y_2, Psi_2, 'm');

set(0,'CurrentFigure',2);
plot(y_1, dPsi_1, 'g');
plot(y_2, dPsi_2, 'm');

set(0,'CurrentFigure',3);
plot(y_1, u_1, 'g');
plot(y_2, u_2, 'm');

%% параллельно формирую управление u_t, которое будет являться управлением
%  с учетом корректировок
arr_1 = abs(y - y_1(1));
[temp, ind] = min(arr_1);

for i = 1:length(y_1)
    u_t(ind + i - 1) = u_1(i);
end

arr_2 = abs(y - y_2(1));
[temp, ind] = min(arr_2);

for i = 1:length(y_2)
    u_t(ind + i - 1) = u_2(i);
end

row = [y_1(1) y_1(end) coefs(1,:) par(1) 2 par(2) 0];

add_control_branch(row);

row = [y_2(1) y_2(end) coefs(2,:) 0 0 par(3) 0];

add_control_branch(row);
%% Но производная все еще выходит за ограничение, поэтому провожу процедуру
%  коректировки второй раз

% Y_0_new = [data1(1,end) data1(2,end) data1(3,end)*data1(2,end)];
% [data1, data2] = lower_constraint_dPsi(Y_0_new, Y_end, dy, dy_constr1);
% y_1 = data1(1,:);
% Psi_1 = data1(2,:);
% dPsi_1 = data1(3,:);
% u_1 = Psi_1.*dPsi_1 - sin(y_1);
% y_2 = data2(1,:);
% Psi_2 = data2(2,:);
% dPsi_2 = data2(3,:);
% u_2 = Psi_2.*dPsi_2 - sin(y_2);
% 
% set(0,'CurrentFigure',1);
% plot(y_1, Psi_1, 'g');
% plot(y_2, Psi_2, 'm');
% 
% set(0,'CurrentFigure',2);
% plot(y_1, dPsi_1, 'g');
% plot(y_2, dPsi_2, 'm');
% 
% set(0,'CurrentFigure',3);
% plot(y_1, u_1, 'g');
% plot(y_2, u_2, 'm');
% 
% %% параллельно формирую управление u_t, которое будет являться управлением
% %  с учетом корректировок
% arr_1 = abs(y - y_1(1));
% [temp, ind] = min(arr_1);
% 
% for i = 1:length(y_1)
%     u_t(ind + i - 1) = u_1(i);
% end
% 
% arr_2 = abs(y - y_2(1));
% [temp, ind] = min(arr_2);
% 
% for i = 1:length(y_2)
%     u_t(ind + i - 1) = u_2(i);
% end
%% Убираю выход за ограничение на Psi(y) < constr2
% если кривая уже корректировалась 
% из-за выхода за ограничение по dPsi(y) / dy. Вместо точки Y_0 надо взять
% Y_0 = [data1(1,end) data1(2,end) data1(3,end)*data1(2,end)]; 
[data, par, coefs] = upper_constraint_Psi(Y_0, Y_end, dy); 
y_3 = data(1,:);
Psi_3 = data(2,:);
dPsi_3 = data(3,:);
u_3 = Psi_3.*dPsi_3 - sin(y_3);

set(0,'CurrentFigure',1);
plot(y_3, Psi_3, 'r');
set(0,'CurrentFigure',2);
plot(y_3, dPsi_3, 'r');
set(0,'CurrentFigure',3);
plot(y_3, u_3, 'r');

arr_3 = abs(y - y_3(1));
[temp, ind] = min(arr_3);

for i = 1:length(y_3)
    u_t(ind + i - 1) = u_3(i);
end

row = [y_3(1) y_3(end) coefs(1,:) par(1) 1 0 0];

add_control_branch(row);
% %% Убираю выход за ограничение на Psi(y) > constr1
% % % если кривая уже корректировалась 
% % % из-за выхода за ограничение по dPsi(y) / dy. Вместо точки Y_0 надо взять
% % Y_0 = [data1(1,end) data1(2,end) data1(3,end)*data1(2,end)]; 
% % data = lower_constraint_Psi(Y_0, Y_end, dy); 
% % y_3 = data(1,:);
% % Psi_3 = data(2,:);
% % dPsi_3 = data(3,:);
% % 
% % set(0,'CurrentFigure',1);
% % plot(y_3, Psi_3, 'r','LineWidth',2);
% % set(0,'CurrentFigure',2);
% % plot(y_3, dPsi_3, 'r');
% В значение params первого элемента массива структур записываю значимое количество элементов

%% ----debug
%     y_0_i = control_arr(2,1);
%     y_end_i = control_arr(2,2);
%     coefs = control_arr(2,3:6);
%     c0 = coefs(1);
%     c1 = coefs(2);
%     c2 = coefs(3);
%     c3 = coefs(4);
%     d = control_arr(2,7);
%     y_0_glob = control_arr(2,9);
%     F_quad_i = @(x)1./(c0 + c1*(x - y_0_glob) + c2*(x - y_0_glob).^2 + c3*(x - y_0_glob).^3 +...
%             d*((x - y_0_i)/(y_end_i - y_0_i)).^2.*(3 - 2*((x - y_0_i)/(y_end_i - y_0_i) )) );
%     t_i = quad(F_quad_i,y_0_i,y_end_i);
%     options = odeset('RelTol',1e-6,'AbsTol',[1e-6 1e-6]);
%     
%     N = N*10000;
%     dy_prog = (y_1(end) - y_1(1)) / N;
%    
%     
%     global y_prog Psi_prog dPsi_prog
%     global dt
%     dt = t_i/N;
%     
%     y_prog = y_1(1):dy_prog:y_1(end);
%     Psi_prog = c0 + c1*(y_prog - y_0_glob) + c2*(y_prog - y_0_glob).^2 + c3*(y_prog - y_0_glob).^3 +...
%         d*((y_prog - y_0_i)/(y_end_i - y_0_i)).^2.*(3 - 2*((y_prog - y_0_i)/(y_end_i - y_0_i) ));
%     dPsi_prog = c1 + 2*c2*(y_prog - y_0_glob) + 3*c3*(y_prog - y_0_glob).^2 + ...
%         d*(6*((y_prog - y_0_i)/(y_end_i - y_0_i)) - 6*((y_prog - y_0_i)/(y_end_i - y_0_i)).^2);
%     
%     [T, Y] = ode45(@kan_system_debug, [0, t_i], [y_1(1) Psi_1(1)], options);
    
%     Y = zeros(N+1,2);
%     Y(1,:) = [y_1(1) Psi_1(1)];
%    
%     i = 1;
%     
%     r1 = -6; r2 = -6;
%     c1_st = -(r1 + r2);
%     c0_st = r1*r2;
% 
%     
%     for t=0:dt:t_i
%         Psi = c0 + c1*(Y(i,1) - y_0_glob) + c2*(Y(i,1) - y_0_glob)^2 + c3*(Y(i,1) - y_0_glob)^3 +...
%         d*((Y(i,1) - y_0_i)/(y_end_i - y_0_i))^2*(3 - 2*((Y(i,1) - y_0_i)/(y_end_i - y_0_i) ));
% 
%         dPsi = c1 + 2*c2*(Y(i,1) - y_0_glob) + 3*c3*(Y(i,1) - y_0_glob)^2 + ...
%         d*(6*((Y(i,1) - y_0_i)/(y_end_i - y_0_i)) - 6*((Y(i,1) - y_0_i)/(y_end_i - y_0_i))^2);
% 
% 
%         % u = Psi*dPsi - sin(Y(i,1));
%         u = Psi_prog(i)*dPsi_prog(i) - sin(Y(i,1)) - c1_st*(Y(i,2) - Psi_prog(i)) - c0_st*(Y(i,1) - y_prog(i));
%         Y(i+1,1) = y_prog(i) + Psi_prog(i)*dt;
%         Y(i+1,2) = Psi_prog(i) + (Psi_prog(i)*dPsi_prog(i))*dt;
%         
%         Y(i+1,1) = Y(i,1) + Y(i,2)*dt;
%         Y(i+1,2) = Y(i,2) + (sin(Y(i,1)) + u)*dt;
%         
%         i = i+1;
%     end
%     set(0,'CurrentFigure',1);
%     plot(Y(:,1), Y(:,2),'y');
%%

% plot(y,u_t,'y');

traj = time_modeling(Y_0);

set(0,'CurrentFigure',1);
plot(traj(:,1), traj(:,2),'y');