clear all
close all
clc
%% �����������
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

% �������� ������� ������ ����
% Psi = c0 + c1*(y - y0) + c2*(y - y0)^2 + c3*(y - y0)^3 + d*(y - y0)^2*(y - yend)^2,
% �����������
% ��������� ��������� Y_0 � �������� ��������� Y_end
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
title('���������� u(y)');
u = Psi.*dPsi - sin(y);

u_t = u; % ����������, ������� ������������ ��� �������������

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
%% ������, � ������� ����� ����������� ������������ ��� ����������
global control_arr
control_arr = [Y_0(1) Y_end(1) c0 c1 c2 c3 0 0 Y_0(1) 0];
% 1 : y_0 - ��������� ����� ������;
% 2 : y_end - �������� ����� ������;
% 3-6 : ������������;
% 7 : d - ��������;
% 8 :  ��� ������
% 9 : y_0_c - ��� ��������, ������� ���������� �� y ��� ����� ���� Y_0(1) ��� Y_0(1) ��
% ������ lower/upper_constraint_dPsi; 
% 10 : ������� ������� ������ � ������, �� �������
% ������� � ��������� ����� �������� ���������, ��� �������� ����������� �
% time_modeling.
% ���� ����� ���������:
% 0 - ������ ���� Psi(y) = c0 + c1*(y - y_0_c) + c2*(y - y_0_c)^2 + c3*(y - y_0_c)^3;
% 1 - ������ ���� Psi(y) = c0 + c1*(y - y_0) + c2*(y - y_0)^2 + c3*(y - y_0)^3 + d*(y - y0)^2*(y - y_end)^2
% 2 - ������ ���� Psi(y) = c0 + c1*(y - y_0_c) + c2*(y - y_0_c)^2 + c3*(y -y_0_c)^3 + d*y_norm.^2.*(3 - 2*y_norm), 
% ��� y_norm = (y - y_0) / (y_end - y_0);

%% ������ ����� �� ����������� dPsi / dy < dy_constr2
% [data1, data2, par, coefs] = upper_constraint_dPsi(Y_0, Y_end, dy, dy_constr2);
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
%
% row = [y_1(1) y_1(end) coefs(1,:) par(1) 2 par(2) 0];
% 
% add_control_branch(row);
% 
% row = [y_2(1) y_2(end) coefs(2,:) 0 0 par(3) 0];
% 
% add_control_branch(row);

% %% ������ ����� �� ����������� dPsi / dy > dy_constr1
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

%% ����������� �������� ���������� u_t, ������� ����� ������ ��� ��������
% ����������
  
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
%% �� ����������� ��� ��� ������� �� �����������, ������� ������� ���������
%  ������������ ������ ���

% Y_0_new = [data1(1,end) data1(2,end) data1(3,end)*data1(2,end)];
% [data1, data2, par, coefs] = lower_constraint_dPsi(Y_0_new, Y_end, dy, dy_constr1);
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
% row = [y_1(1) y_1(end) coefs(1,:) par(1) 2 par(2) 0];
% 
% add_control_branch(row);
% 
% row = [y_2(1) y_2(end) coefs(2,:) 0 0 par(3) 0];
% 
% add_control_branch(row);
%% ������ ����� �� ����������� �� Psi(y) < constr2
% ���� ������ ��� ���������������� 
% ��-�� ������ �� ����������� �� dPsi(y) / dy. ������ ����� Y_0 ���� �����
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
%% ������ ����� �� ����������� �� Psi(y) > constr1
% % ���� ������ ��� ���������������� 
% % ��-�� ������ �� ����������� �� dPsi(y) / dy. ������ ����� Y_0 ���� �����
% Y_0 = [data1(1,end) data1(2,end) data1(3,end)*data1(2,end)]; 
% [data, par, coefs] = lower_constraint_Psi(Y_0, Y_end, dy); 
% y_3 = data(1,:);
% Psi_3 = data(2,:);
% dPsi_3 = data(3,:);
% row = [y_3(1) y_3(end) coefs(1,:) par(1) 1 0 0];
% add_control_branch(row);
% 
% set(0,'CurrentFigure',1);
% plot(y_3, Psi_3, 'r','LineWidth',2);
% set(0,'CurrentFigure',2);
% plot(y_3, dPsi_3, 'r');


% plot(y,u_t,'y');

% �������������� ������� �� �������
traj = time_modeling(Y_0);

set(0,'CurrentFigure',1);
plot(traj(:,1), traj(:,2),'y');