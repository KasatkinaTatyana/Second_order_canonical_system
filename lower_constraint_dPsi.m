function [F1, F2, params, coefs] = lower_constraint_dPsi(Y_0, Y_end, dy, dy_constr1)
% Функция генерирует кривую, отрезок, на котором она строится, и ее
% производную таким образом, чтобы выполнялось ограничение dPsi(y) / dy > dy_constr1

% F1 соответсвует первой кривой (левой),
% F2 - правой, они имеют следующую структуру
% массив y - массив функции Psi - массив dPsi / dy

% params и coefs нужны для того, чтобы затем однозначно восстановить кривую
% и построить нужное управление

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

y=Y_0(1):dy:Y_end(1);

Psi = c0 + c1*(y - Y_0(1)) + c2*(y - Y_0(1)).^2 + c3*(y - Y_0(1)).^3;
dPsi = c1 + 2*c2*(y - Y_0(1)) + 3*c3*(y - Y_0(1)).^2;

N_shift = 30;

N = (Y_end(1) - Y_0(1))/dy + 1;

flag = 0;
% Ищем границы интервала, на котором нарушается ограничение
for i=1:N
    % ------- -- dPsi / dy > constr -----------------------------------
    if ((flag == 0)&&(dPsi(i) <= dy_constr1))
    % -----------------------------------------------------------------
        if ((i - N_shift) < 1)
            y_left = Y_0;
        else
            y_left(1) = Y_0(1) + (i-1 - N_shift)*dy;  % отступили N_shift шагов назад от критической точки
            y_left(2) = Psi(i - N_shift);
            y_left(3) = (c1 + 2*c2*(y_left(1) - Y_0(1)) + 3*c3*(y_left(1) - Y_0(1))^2)*y_left(2);
        end
        flag = 1;
    end
    % --------- dPsi / dy > constr ------------------------------------
    if ((flag==1)&&(dPsi(i) > dy_constr1))
    % -----------------------------------------------------------------
    if ((i + N_shift) > N + 1)
        y_right = Y_end;
    else
        y_right(1) = Y_0(1) + (i-1 + N_shift)*dy;
        y_right(2) = Psi(i + N_shift);
        y_right(3) = (c1 + 2*c2*(y_right(1) - Y_0(1)) + 3*c3*(y_right(1) - Y_0(1))^2)*y_right(2);
    end
        break;
    end
end

y_1 = y_left(1):dy:y_right(1);
y_norm = (y_1 - y_left(1)) / (y_right(1) - y_left(1));     
h_d = 0.01;
for d = 0:h_d:10
    Psi_1 = c0 + c1*(y_1 - Y_0(1)) + c2*(y_1 - Y_0(1)).^2 + c3*(y_1 - Y_0(1)).^3 + d*y_norm.^2.*(3 - 2*y_norm);
    dPsi_1 = c1 + 2*c2*(y_1 - Y_0(1)) + 3*c3*(y_1 - Y_0(1)).^2 + d*(6* (y_1 - y_left(1)) / (y_right(1) - y_left(1))^2 ...
        - 6*(y_1 - y_left(1)).^2 / (y_right(1) - y_left(1))^3);
    if (min(dPsi_1) > dy_constr1)
        break;
    end
end

Y_1_end = [y_right(1) Psi_1(end) dPsi_1(end)*Psi_1(end)];

replace_part = curve_synthesis(Y_1_end, Y_end, dy, 0);

Psi_2 = replace_part(1,:);
dPsi_2 = replace_part(2,:);

y_2 = y_right(1):dy:Y_end(1);

F1 = [y_1; Psi_1; dPsi_1];
F2 = [y_2; Psi_2; dPsi_2];
params = [d; Y_0(1); y_right(1)];

coefs_2 = calc_coefs(Y_1_end, Y_end);

coefs = [c0 c1 c2 c3; coefs_2; 0 0 0 0];