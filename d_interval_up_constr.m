function F = d_interval_up_constr(Y_0, Y_end, dy, d0)
% Определяются границы интервала по d, на котором не нарушается
% ограничение. При условии что при d0 ограничение не нарушается.

global constr1 constr2

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

N = (Y_end(1) - Y_0(1))/dy + 1;
% 
% N=N*100;
% 
% dy = (Y_end(1) - Y_0(1)) / N;



Psi = zeros(1,N);

if (d0 > 0)
    h_d = 0.01;
end
if (d0 < 0)
    h_d = -0.01;
end

flag = 0;
for d=d0:h_d:(10*d0)
    i = 1;
    for y = Y_0(1):dy:Y_end(1)
        Psi(i) = c0 + c1*(y - Y_0(1)) + c2*(y - Y_0(1))^2 + c3*(y - Y_0(1))^3 + ...
            d*(y - Y_0(1))^2*(y - Y_end(1))^2;
        %-----1) y' < constr-----------------------------------------------
        if ((Psi(i) >= constr2)||(Psi(i) <= 0))
        %------------------------------------------------------------------
            flag = 1;
            break;
        end
        i = i + 1;
    end
    if (flag == 1)
        break;
    end
end

d_max = d - h_d;

flag = 0;
for d=d0:(-h_d):(-10*d0)
    i = 1;
    for y = Y_0(1):dy:Y_end(1)
        Psi(i) = c0 + c1*(y - Y_0(1)) + c2*(y - Y_0(1))^2 + c3*(y - Y_0(1))^3 + ...
            d*(y - Y_0(1))^2*(y - Y_end(1))^2;
        %-----1) y' < constr-----------------------------------------------
        if ((Psi(i) >= constr2)||(Psi(i) <= 0))
        %------------------------------------------------------------------
            flag = 1;
            break;
        end
        i = i + 1;
    end
    if (flag == 1)
        break;
    end
end

d_min = d + h_d;

if (d0 >= 0)
    F = [d_min d_max];
else
    F = [d_max d_min];
end