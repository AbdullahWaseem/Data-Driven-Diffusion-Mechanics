function [rot, k] = findrotationNR(gradC, gradC_str, tenJ, tenJ_str, Dref)
% This function provides the angle for which the internal state and the
% data minimizes the distance. It uses Newton-Raphson method to linearize
% the system. 
%
% Author:  Abdullah Waseem
% Created: 20-March-2020
% Contact: engineerabdullah@ymail.com
%
% Initialization
% If gradC_str=0 then the function becomes undefinied, so perturbing a bit.
if gradC_str(1)==0 
    gradC_str(1) = eps;
     tenJ_str(1) = -Dref(1,1)*gradC_str(1);
end
%
rot = 1;
tol = 1e-8;
ktol= 1;
k   = 0;
[f_p, df_p] = f_df(gradC, gradC_str,rot, tenJ, tenJ_str, Dref);
%
% Newton-Raphson Algorithm
while ktol>tol && k<1e3
    % Solution
    drot = -( f_p/df_p );
    rot = rot + drot;
    % Periodic angle
    if abs(rot)>2*pi
        rot = rem((rot),2*pi);
    end
    [f,df] = f_df(gradC, gradC_str, rot, tenJ, tenJ_str, Dref);
    ktol = norm(f-f_p);
    f_p = f; df_p = df;
    k = k+1;
end
%
% -------------------------------------------
% Function and its Derivative
%
function [f,df] = f_df(gradC, gradC_str, rot, tenJ, tenJ_str, Dref)
    % Rotation matrix and its derivative
    Q = [ cos(rot) -sin(rot); sin(rot)  cos(rot)];
    dQ= [-sin(rot) -cos(rot); cos(rot) -sin(rot)];
    % Function
    f= (dQ'*gradC_str)'*Dref*(gradC-Q'*gradC_str) + ...
        (dQ'*tenJ_str)'*(Dref\eye(2))*(tenJ-Q'*tenJ_str);
    % Derivative of the function
    df= -(Q'*gradC_str)'*Dref*(gradC-Q'*gradC_str) - ...
         (Q'*tenJ_str)'*(Dref\eye(2))*(tenJ-Q'*tenJ_str) - ...
         (dQ'*gradC_str)'*Dref*(dQ'*gradC_str) - ...
         (dQ'*tenJ_str)'*(Dref\eye(2))*(dQ'*tenJ_str);
    %
end
%
end