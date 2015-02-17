function dy=Trajectory_type_1(t,y)
global c0 c1 c2 c3 y_0_i y_end_i d
dy=zeros(1,1);

dy(1) = c0 + c1*(y - y_0_i) + c2*(y - y_0_i)^2 + c3*(y - y_0_i)^3 + d*(y - y_0_i)^2*(y - y_end_i)^2;