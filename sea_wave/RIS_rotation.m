function [R,quaternion,euler] = RIS_rotation(value)
%RIS_ROTATION Body-to-world ZYX rotation; input Euler [roll;pitch;yaw] or wxyz.
    arguments
        value (:,1) double {mustBeReal,mustBeFinite}
    end
    if numel(value)==3
        c = cos(value/2); s = sin(value/2);
        quaternion = [c(1)*c(2)*c(3)+s(1)*s(2)*s(3); ...
            s(1)*c(2)*c(3)-c(1)*s(2)*s(3); ...
            c(1)*s(2)*c(3)+s(1)*c(2)*s(3); ...
            c(1)*c(2)*s(3)-s(1)*s(2)*c(3)];
    elseif numel(value)==4 && norm(value)>0
        quaternion = value/norm(value);
    else
        error('RIS_rotation:Size','Supply three Euler angles or a nonzero wxyz quaternion.');
    end
    w = quaternion(1); x = quaternion(2); y = quaternion(3); z = quaternion(4);
    R = [1-2*(y*y+z*z),2*(x*y-w*z),2*(x*z+w*y); ...
        2*(x*y+w*z),1-2*(x*x+z*z),2*(y*z-w*x); ...
        2*(x*z-w*y),2*(y*z+w*x),1-2*(x*x+y*y)];
    euler = [atan2(R(3,2),R(3,3));asin(max(-1,min(1,-R(3,1))));atan2(R(2,1),R(1,1))];
end
