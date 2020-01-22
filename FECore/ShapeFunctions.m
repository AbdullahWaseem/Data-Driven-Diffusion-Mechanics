% In isogeometric FEM master element is constant. So, it's shape function
% and their gradients does not change. Hence, these quantities can be
% evalutated in the beginning of the FEM code once and for all.
%
% Author: Abdullah Waseem       
% Created: 19-August-2017       
% Last Modified: 11-March-2018   
% Contact: engineerabdullah@ymail.com
%
% INPUT: 
%       ngp   = number of gauss points
%       zta   = the zeta coordinates generated by gaussianquad function
%       elementtype = the type of the element as explained above
%       dim   = problem dimension
% 
% OUTPUT:
%       N     = shape functions
%       dphi    = gradient of shape functions
% 
%
% IT DOES NOT MATTER IN WHAT ORDER ENUMERATION OVER THE NODES IS WHILE PERFORMING THE ELEMENT INTEGRAL
% SO NODE NUMBERING WHILE DOING INTEGRATION IS NOT IMPORTANT i.e. 1 + 2 = 2 + 1
% 
% References:
%   For Quads: http://www.me.metu.edu.tr/courses/me582/files/Handouts/Shape_Functions.pdf
%   For Tri  : http://math2.uncc.edu/~shaodeng/TEACHING/math5172/Lectures/Lect_15.PDF
% 
% Quadrangle:            Quadrangle8:            Quadrangle9:
% 
%       e
%       ^
%       |
% 4-----------3          4-----7-----3           4-----7-----3 
% |     |     |          |           |           |           | 
% |     |     |          |           |           |           | 
% |     +---- | --> z    8           6           8     9     6 
% |           |          |           |           |           | 
% |           |          |           |           |           | 
% 1-----------2          1-----5-----2           1-----5-----2 
% elementnumber==3      elementnumber==16      elementnumber==10                  
% 
% 
% Initializing shape functions
PHI  = zeros(nne, ngp);
% Initializing gradient of shape functions.
DPHI = zeros(nne, dim, ngp);
% 
switch elementnumber
% For Quadratic 4 Node Element.
case 3
% The gauss point loop
for gs = 1 : ngp
    % Evaluation of zeta = z and eta = e at that gauss point
     z = glz(gs,1);
	 e = glz(gs,2);
     % Shape functions
     PHI(1, gs) = 0.25 * (1-z) * (1-e);
     PHI(2, gs) = 0.25 * (1+z) * (1-e);
     PHI(3, gs) = 0.25 * (1+z) * (1+e);
     PHI(4, gs) = 0.25 * (1-z) * (1+e);
     % Derivatives of shape functions.
     DPHI(1, 1, gs) = -0.25 * (1-e);
     DPHI(1, 2, gs) = -0.25 * (1-z);
     DPHI(2, 1, gs) =  0.25 * (1-e);
     DPHI(2, 2, gs) = -0.25 * (1+z);
     DPHI(3, 1, gs) =  0.25 * (1+e);
     DPHI(3, 2, gs) =  0.25 * (1+z);
     DPHI(4, 1, gs) = -0.25 * (1+e);
     DPHI(4, 2, gs) =  0.25 * (1-z);
end
% Q4 ends
% 
% For Quadratic 8 Node Element.
case 16
% The gauss point loop
for gs = 1 : ngp
    % Evaluation of zeta = z and eta = e at that gauss point
    z = glz(gs,1);       
	e = glz(gs,2);
    % Shape functions
    PHI(1, gs) = 0.25 * (1-z)   * (1-e)  * (-z-e-1);
    PHI(2, gs) = 0.25 * (1+z)   * (1-e)  * ( z-e-1);
    PHI(3, gs) = 0.25 * (1+z)   * (1+e)  * ( z+e-1);
    PHI(4, gs) = 0.25 * (1-z)   * (1+e)  * (-z+e-1);
    PHI(5, gs) = 0.50 * (1-z^2) * (1-e)  ;
    PHI(6, gs) = 0.50 * (1+z)   * (1-e^2);
    PHI(7, gs) = 0.50 * (1-z^2) * (1+e)  ;
    PHI(8, gs) = 0.50 * (1-z)   * (1-e^2);
    % Derivative of shape fucntions
    DPHI(1, 1, gs) = - 0.25 * (e-1)   * (2*z+e);
    DPHI(1, 2, gs) = - 0.25 * (z-1)   * (z+2*e);
    DPHI(2, 1, gs) = - 0.25 * (e-1)   * (2*z-e);
    DPHI(2, 2, gs) = - 0.25 * (z+1)   * (z-2*e);
    DPHI(3, 1, gs) =   0.25 * (e+1)   * (2*z+e);
    DPHI(3, 2, gs) =   0.25 * (z+1)   * (z+2*e);
    DPHI(4, 1, gs) =   0.25 * (e+1)   * (2*z-e);
    DPHI(4, 2, gs) =   0.25 * (z-1)   * (z-2*e);
    DPHI(5, 1, gs) =   0.50 * (-2*z)  *   (1-e);
    DPHI(5, 2, gs) =   0.50 * (1-z^2) *    (-1);
    DPHI(6, 1, gs) =   0.50 * ( 1)    * (1-e^2);
    DPHI(6, 2, gs) =   0.50 * (1+z)   *  (-2*e);
    DPHI(7, 1, gs) =   0.50 * (-2*z)  *   (1+e);
    DPHI(7, 2, gs) =   0.50 * (1-z^2) *     (1);
    DPHI(8, 1, gs) =   0.50 *  (-1)   * (1-e^2);
    DPHI(8, 2, gs) =   0.50 * (1-z)   *  (-2*e);
end
% Q8 Ends here.
%       
% For Quadratic 9 Node Element.
case 10
% The gauss point loop
for gs = 1 : ngp
    % Evaluation of zeta = z and eta = e at that gauss point
    z = glz(gs,1);       
	e = glz(gs,2);
    % Shape functions
    PHI(1, gs) = 0.25 * (z^2 - z) * (e^2 - e);
    PHI(2, gs) = 0.25 * (z^2 + z) * (e^2 - e);
    PHI(3, gs) = 0.25 * (z^2 + z) * (e^2 + e);
    PHI(4, gs) = 0.25 * (z^2 - z) * (e^2 + e);
    PHI(5, gs) = 0.50 * (1 - z^2) * (e^2 - e);
    PHI(6, gs) = 0.50 * (z^2 + z) * (1 - e^2);
    PHI(7, gs) = 0.50 * (1 - z^2) * (e^2 + e);
    PHI(8, gs) = 0.50 * (z^2 - z) * (1 - e^2);
    PHI(9, gs) =        (1 - z^2) * (1 - e^2);
    % Derivative of shape fucntions
    DPHI(1, 1, gs) = 0.25 * (2*z - 1) * (e^2 - e);
    DPHI(1, 2, gs) = 0.25 * (z^2 - z) * (2*e - 1);
    DPHI(2, 1, gs) = 0.25 * (2*z + 1) * (e^2 - e);
    DPHI(2, 2, gs) = 0.25 * (z^2 + z) * (2*e - 1);
    DPHI(3, 1, gs) = 0.25 * (2*z + 1) * (e^2 + e);
    DPHI(3, 2, gs) = 0.25 * (z^2 + z) * (2*e + 1);
    DPHI(4, 1, gs) = 0.25 * (2*z - 1) * (e^2 + e);
    DPHI(4, 2, gs) = 0.25 * (z^2 - z) * (2*e + 1);
    DPHI(5, 1, gs) = 0.50 * (-2 * z ) * (e^2 - e);
    DPHI(5, 2, gs) = 0.50 * (1 - z^2) * (2*e - 1);
    DPHI(6, 1, gs) = 0.50 * (2*z + 1) * (1 - e^2);
    DPHI(6, 2, gs) = 0.50 * (z^2 + z) * (-2 * e );
    DPHI(7, 1, gs) = 0.50 * (-2 * z ) * (e^2 + e);
    DPHI(7, 2, gs) = 0.50 * (1 - z^2) * (2*e + 1);
    DPHI(8, 1, gs) = 0.50 * (2*z - 1) * (1 - e^2);
    DPHI(8, 2, gs) = 0.50 * (z^2 - z) * (-2 * e );
    DPHI(9, 1, gs) =        (-2 * z ) * (1 - e^2);
    DPHI(9, 2, gs) =        (1 - z^2) * (-2 * e );     
end
% Q9 Ends
% 
% TRI3                   TRI6:                  TRI10:
% 
% e                                                              
% ^                                          
% |                                       
% 3                       3                   3            
% |`\                     |`\                 | \       
% |  `\                   |  `\               8   7
% |    `\                 6    `5             |     \
% |      `\               |      `\           9  10  6 
% |        `\             |        `\         |         \
% 1----------2 --> z      1-----4----2        1---4---5---2
% elementnumber=2        elementnumber=9     elementnumber=21
% 
% For Triangular 3 Node Element.
case 2
% The gauss point loop
for gs = 1 : ngp
    % Evaluation of zeta = z and eta = e at that gauss point
    z = glz(gs,1);
	e = glz(gs,2);
    % Shape functions
    PHI(1, gs) = 1-z-e;
    PHI(2, gs) = z;
    PHI(3, gs) = e;
    % Derivative of shape functions
    DPHI(1, 1, gs) = -1;
    DPHI(1, 2, gs) = -1;
    DPHI(2, 1, gs) =  1;
    DPHI(2, 2, gs) =  0;
    DPHI(3, 1, gs) =  0;
    DPHI(3, 2, gs) =  1;
end
% TRI3 Ends
% 
% For Triangular 6 Node Element.
case 9
% The gauss point loop
for gs = 1 : ngp
    % Evaluation of zeta = z and eta = e at that gauss point
    z = glz(gs,1);       
	e = glz(gs,2);
    % Shape functions
    PHI(1, gs) = 2 * (1-z-e) * (0.5-z-e);
    PHI(2, gs) = 2 *    z    *  (z-0.5);
    PHI(3, gs) = 2 *    e    *  (e-0.5);
    PHI(4, gs) = 4 * (1-z-e) *     z;
    PHI(5, gs) = 4 *    z    *     e;
    PHI(6, gs) = 4 * (1-z-e) *     e;
    % Derivative of shape functions.
    DPHI(1, 1, gs) =  4*z + 4*e - 3;
    DPHI(1, 2, gs) =  4*z + 4*e - 3;
    DPHI(2, 1, gs) =  4*z - 1;
    DPHI(2, 2, gs) =  0;
    DPHI(3, 1, gs) =  0;
    DPHI(3, 2, gs) =  4*e - 1;
    DPHI(4, 1, gs) =  4 - 8*z - 4*e;
    DPHI(4, 2, gs) = -4 * z;
    DPHI(5, 1, gs) =  4 * e;
    DPHI(5, 2, gs) =  4 * z;
    DPHI(6, 1, gs) = -4 * e;
    DPHI(6, 2, gs) =  4 - 4*z - 8*e; 
end
% TRI6 Ends
% 
% For Triangular 10 Node Element.
case 21
% The gauss point loop
for gs = 1 : ngp
    % Evaluation of zeta = z and eta = e at that gauss point
    z = glz(gs,1);       
	e = glz(gs,2);
    % Shape functions
    PHI(1, gs) = 1/2 *    z    *    (3*z-1)    *    (3*z-2);
    PHI(2, gs) = 1/2 *    e    *    (3*e-1)    *    (3*e-2);
    PHI(3, gs) = 1/2 * (1-z-e) * (3*(1-z-e)-1) * (3*(1-z-e)-2);
    PHI(4, gs) = 9/2 *    z    *       e       *    (3*z-1);
    PHI(5, gs) = 9/2 *    z    *       e       *    (3*e-1);
    PHI(6, gs) = 9/2 *    e    *    (1-z-e)    *    (3*e-1);
    PHI(7, gs) = 9/2 *    e    *    (1-z-e)    * (3*(1-z-e)-1);
    PHI(8, gs) = 9/2 * (1-z-e) *       z       * (3*(1-z-e)-1);
    PHI(9, gs) = 9/2 * (1-z-e) *       z       *    (3*z-1);
    PHI(10, gs)= 27  *    z    *       e       *    (1-z-e);
    % Derivative of shape functions.
    DPHI(1, 1, gs) = 27/2 *z^2 - 9*z + 1;
    DPHI(1, 2, gs) =     0;
    DPHI(2, 1, gs) =     0;
    DPHI(2, 2, gs) = 27/2 *    e^2   - 9*e + 1;
    DPHI(3, 1, gs) =  1/2 * (-27*z^2 + z*(36-54*e) - 27*e^2 + 36*e - 11);
    DPHI(3, 2, gs) =  1/2 * (-27*z^2 + z*(36-54*e) - 27*e^2 + 36*e - 11);
    DPHI(4, 1, gs) =  9/2 *  (6*z-1) * e;
    DPHI(4, 2, gs) =  9/2 *     z    * e * (3*z-1);
    DPHI(5, 1, gs) =  9/2 *     e    * (3*e-1);
    DPHI(5, 2, gs) =  9/2 *     z    * (6*e-1);
    DPHI(6, 1, gs) = -9/2 *     e    * (3*e-1);
    DPHI(6, 2, gs) = -9/2 * (z*(6*e-1) + z^2 - 8*e + 1);
    DPHI(7, 1, gs) =  9/2 *     e    * (6*z + 6*e - 5);
    DPHI(7, 2, gs) =  9/2 * (3*z^2   + z*(12*e-5) + 9*e^2 - 10*e + 2);
    DPHI(8, 1, gs) =  9/2 * (9*z^2   + 2*z*(6*e-5) + 3*e^2 - 5*e + 2);
    DPHI(8, 2, gs) =  9/2 *     z    * (6*z + 6*e -5);
    DPHI(9, 1, gs) = -9/2 * (9*z^2   + z*(6*e-8) - e + 1);
    DPHI(9, 2, gs) = -9/2 *     z    * (3*z-1);
    DPHI(10,1, gs) = -27  *     e    * (2*z + e + 1);
    DPHI(10,2, gs) = -27  *     z    * (z + 2*e + 1);
end
% TRI10 Ends

end
% 

% 