function F = add_control_branch(row)
global control_arr

N = length(control_arr(:,1));

y_0_add = row(1);
y_end_add = row(2);

for i=1:N
    y_0_i = control_arr(i,1);
    y_end_i = control_arr(i,2);
    
    if  (abs(y_0_add - y_0_i) < 1e-8)&&(abs(y_end_add - y_end_i) < 1e-8)
        control_arr(i,:) = row;
        break;
    end
    if  ( (abs(y_0_add - y_0_i) < 1e-8)&&(y_end_add < y_end_i) )
        control_arr(i,:) = row;
        control_arr = [control_arr; [y_end_add y_end_i control_arr(i,3:10)]];
        break;
    end
    if  ( (y_0_add > y_0_i)&&(abs(y_end_add - y_end_i) < 1e-8) )
        control_arr(i,:) = [y_0_i y_0_add control_arr(i,3:10)];
        control_arr = [control_arr; row];
        break;
    end
    if ( (y_0_add > y_0_i)&&(y_end_add < y_end_i) )
        control_arr(i,:) = [y_0_i y_0_add control_arr(i,3:10)];
        control_arr = [control_arr; row];
        control_arr = [control_arr; [y_end_add y_end_i control_arr(i,3:10)]];
        break;
    end
end