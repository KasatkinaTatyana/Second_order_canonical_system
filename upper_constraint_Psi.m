function [F, params, coefs] = upper_constraint_Psi(Y_0, Y_end, dy)
% ������� ���������� ������, �������, �� ������� ��� ��������, � ��
% ����������� ����� �������, ����� ����������� ����������� Psi(y) < constr2

% F ����� ��������� ���������
% ������ y - ������ ������� Psi - ������ dPsi / dy

% params � coefs ����� ��� ����, ����� ����� ���������� ������������ ������
% � ��������� ������ ����������
global constr2

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

% dy = (Y_end(1) - Y_0(1)) / N;
N = (Y_end(1) - Y_0(1))/dy + 1;
y=Y_0(1):dy:Y_end(1);

Psi = c0 + c1*(y - Y_0(1)) + c2*(y - Y_0(1)).^2 + c3*(y - Y_0(1)).^3;

d_min = -100;
d_max = 100;
hd = 0.01;
N_shift = 30;

cond = false;


while not(cond)
    flag = 0;
    % ���� ������� ���������, �� ������� ���������� �����������
    for i=1:N
        % ------- 1) y' < constr ------------------------------------------
        if ((flag == 0)&&(Psi(i) >= constr2))
        % ------- 2) y' > constr ------------------------------------------
        % if ((flag == 0)&&(Psi(i) <= constr1))
        % -----------------------------------------------------------------
            y_left(1) = Y_0(1) + (i-1 - N_shift)*dy;  % ��������� N_shift ����� ����� �� ����������� �����
            y_left(2) = Psi(i - N_shift);
            y_left(3) = (c1 + 2*c2*(y_left(1) - Y_0(1)) + 3*c3*(y_left(1) - Y_0(1))^2)*y_left(2);
            flag = 1;
        end
        % -------1)  y' < constr ------------------------------------------
        if ((flag==1)&&(Psi(i) < constr2))
        % -------2) y' > constr -------------------------------------------
        % if ((flag==1)&&(Psi(i) > constr1))
        % -----------------------------------------------------------------
            y_right(1) = Y_0(1) + (i-1 + N_shift)*dy;
            y_right(2) = Psi(i + N_shift);
            y_right(3) = (c1 + 2*c2*(y_right(1) - Y_0(1)) + 3*c3*(y_right(1) - Y_0(1))^2)*y_right(2);
            break;
        end
    end
    
    for d=d_min:hd:d_max
        cond = IsCurveExist_up_constr(y_left, y_right, dy, d);
        if (cond == true)
            break;
        end
    end
    
    if not(cond)
        N_shift = N_shift + 1;
    end
end


if (cond)
    % ������� �������� �� d, �� ������� ����������� �� ����������
    d_lims = d_interval_up_constr(y_left, y_right, dy, d);
    % ������ �������������� ������� ��� ������������ d
    replace_part = curve_synthesis(y_left, y_right, dy, d_lims(2));
    Psi_2 = replace_part(1,:);
    dPsi_2 = replace_part(2,:);
    x_2 = y_left(1):dy:y_right(1);
    F = [x_2; Psi_2; dPsi_2];
else
    disp('�������� �������� ���� ������');
    F = false;
end

params = [d_lims(2); 0; 0];

coefs = [calc_coefs(y_left, y_right); 0 0 0 0; 0 0 0 0];