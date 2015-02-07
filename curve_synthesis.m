function F = curve_synthesis(Y_0, Y_end, dy, d)
% Возвращается значение точек кривой, которая соединяет точки Y_0 и Y_end.
% dy - шаг разбиения. Данная кривая имеет вид
% y' = c0 + c1*(y - y0) + c2*(y - y0)^2 + c3*(y - y0)^3 + d*(y - y0)^2*(y -
% yend)^2

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
%----------------------Проверка--------------------------------------------
Expr = c0 + c1*delta + c2*delta^2 + c3*delta^3;

Expr_1 = c1 + 2*c2*delta + 3*c3*delta^2;
%--------------------------------------------------------------------------

N = (Y_end(1) - Y_0(1))/dy + 1;
% 
% N=N*100;
% 
% dy = (Y_end(1) - Y_0(1)) / N;

Psi = zeros(1,N);
dPsi = zeros(1,N);

i = 1;
for y = Y_0(1):dy:Y_end(1)
    Psi(i) = c0 + c1*(y - Y_0(1)) + c2*(y - Y_0(1))^2 + c3*(y - Y_0(1))^3 + ...
        d*(y - Y_0(1))^2*(y - Y_end(1))^2;
    dPsi(i) = c1 + 2*c2*(y - Y_0(1)) + 3*c3*(y - Y_0(1))^2 + ...
        2*d*(y - Y_0(1))^2*(y - Y_end(1)) + 2*d*(y - Y_0(1))*(y - Y_end(1))^2;
    i = i + 1;
end

F = [Psi; dPsi];