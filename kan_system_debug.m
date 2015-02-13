function dy = kan_system_debug(t,y)
dy = zeros(2,1);

global control_arr
global count

global y_prog Psi_prog dPsi_prog
global t_i dt

    y_0_i = control_arr(2,1);
    y_end_i = control_arr(2,2);
    coefs = control_arr(2,3:6);
    c0 = coefs(1);
    c1 = coefs(2);
    c2 = coefs(3);
    c3 = coefs(4);
    d = control_arr(2,7);
    y_0_glob = control_arr(2,9);
    Psi = c0 + c1*(y(1) - y_0_glob) + c2*(y(1) - y_0_glob)^2 + c3*(y(1) - y_0_glob)^3 +...
        d*((y(1) - y_0_i)/(y_end_i - y_0_i))^2*(3 - 2*((y(1) - y_0_i)/(y_end_i - y_0_i) ));

    dPsi = c1 + 2*c2*(y(1) - y_0_glob) + 3*c3*(y(1) - y_0_glob)^2 + ...
        d*(6*((y(1) - y_0_i)/(y_end_i - y_0_i)) - 6*((y(1) - y_0_i)/(y_end_i - y_0_i))^2);


% u = Psi*dPsi - sin(y(1));
r1 = -0.01; r2 = -10;
c1_st = -(r1 + r2);
c0_st = r1*r2;

i = floor(t/dt) + 1;

u = Psi_prog(i)*dPsi_prog(i) - sin(y(1)) - c1_st*(y(2) - Psi_prog(i)) - c0_st*(y(1) - y_prog(i));


plot(y(1), u,'Marker','*');

dy(1) = y(2);
dy(2) = sin(y(1)) + u;