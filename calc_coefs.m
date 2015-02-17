function F = calc_coefs(Y_0, Y_end)
% Y_0, Y_end - состояния динамической системы вида [y, dy/dt, d^2y/dt^2],
% y, dy/dt - координаты фазовой плоскости
% Функция вычисляет коэффициенты c0 c1 c2 c3 функции вида
% Psi(y) = c0 + c1*(y - Y_0(1)) + c2*(y - Y_0(1))^2 + c3*(y - Y_0(1))^3;
% которая соединяет начальное положение Y_0 с конечным Y_end

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

F = [c0 c1 c2 c3];