function dy=Trajectory_type_0(t,y)
global c0 c1 c2 c3 y_0_c 
dy=zeros(1,1);

dy(1) = c0 + c1*(y - y_0_c) + c2*(y - y_0_c)^2 + c3*(y - y_0_c)^3;