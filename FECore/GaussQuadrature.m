% This script provides gauss quadrature points and weights.
% 
% Author: Abdullah Waseem       
% Created: 19-August-2017       
% Last Modified: 11-March-2018
% Contact: engineerabdullah@ymail.com
% 
% References:
%   For Quads: http://www.me.metu.edu.tr/courses/me582/files/Handouts/Shape_Functions.pdf
%   For Tri  : http://math2.uncc.edu/~shaodeng/TEACHING/math5172/Lectures/Lect_15.PDF

%% Quadrangles 
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
if (ngp==1 && (elementnumber==3 || elementnumber==10 || elementnumber==16))
    glz(1,1) = 0;
    glz(1,2) = 0;
    glw(1) = 2;
    % Area denominator.
    Ad = 1;     
elseif (ngp==4 && (elementnumber==3 || elementnumber==10 || elementnumber==16))    
    % At Gauss point 1
    glz(1,1) = -0.5773502692; 
    glz(1,2) = -0.5773502692;
    glw(1)   =  1.0 * 1.0;    
    % At Gauss point 2
    glz(2,1) =  0.5773502692;
    glz(2,2) = -0.5773502692;
    glw(2)   =  1.0 * 1.0;    
    % At Gauss point 3
    glz(3,1) =  0.5773502692;
    glz(3,2) =  0.5773502692;
    glw(3)   =  1.0 * 1.0;    
    % At Gauss point 4
    glz(4,1) = -0.5773502692;
    glz(4,2) =  0.5773502692;
    glw(4)   =  1.0 * 1.0;    
    % Area denominator.
    Ad = 1;    
%     
elseif (ngp==9 && (elementnumber==3 || elementnumber==10 || elementnumber==16))    
    % At Gauss point 1
    glz(1,1) = -0.774596669241483;
    glz(1,2) = -0.774596669241483;
    glw(1)   =  0.555555555555555 * 0.555555555555555;    
    % At Gauss point 2
    glz(2,1) =  0.774596669241483;
    glz(2,2) = -0.774596669241483;
    glw(2)   =  0.555555555555555 * 0.555555555555555;    
    % At Gauss point 3
    glz(3,1) =  0.774596669241483;
    glz(3,2) =  0.774596669241483;
    glw(3)   =  0.555555555555555 * 0.555555555555555;     
    % At Gauss point 4
    glz(4,1) = -0.774596669241483;
    glz(4,2) =  0.774596669241483;
    glw(4)   =  0.555555555555555 * 0.555555555555555;         
    % At Gauss point 5
    glz(5,1) =  0.000000000000000;
    glz(5,2) = -0.774596669241483;
    glw(5)   =  0.888888888888888 * 0.555555555555555;    
    % At Gauss point 6
    glz(6,1) =  0.774596669241483;
    glz(6,2) =  0.000000000000000;
    glw(6)   =  0.555555555555555 * 0.888888888888888;        
    % At Gauss point 7
    glz(7,1) =  0.000000000000000;
    glz(7,2) =  0.774596669241483;
    glw(7)   =  0.888888888888888 * 0.555555555555555;    
    % At Gauss point 8
    glz(8,1) = -0.774596669241483;
    glz(8,2) =  0.000000000000000;
    glw(8)   =  0.555555555555555 * 0.888888888888888;   
    % At Gauss point 9
    glz(9,1) =  0.000000000000000;
    glz(9,2) =  0.000000000000000;
    glw(9)   =  0.888888888888888 * 0.888888888888888;    
    % Area denominator.
    Ad = 1;

%% Triangles 
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
% 
% For Triangular 3 Node Element.
elseif (ngp==3 && (elementnumber==2 || elementnumber==9 || elementnumber==21))
    % At Gauss point 1
    glz(1,1) = 1/6;
    glz(1,2) = 1/6;
    glw(1)   = 1/3;    
    % At Gauss point 2
    glz(2,1) = 2/3;
    glz(2,2) = 1/6;
    glw(2)   = 1/3;    
    % At Gauss point 3
    glz(3,1) = 1/6;
    glz(3,2) = 2/3;
    glw(3)   = 1/3;    
    % Area denominator.
    Ad = 2;   
elseif (ngp==6 && (elementnumber==2 || elementnumber==9 || elementnumber==21))
    % At Gauss point 1
    glz(1,1) = 0.091576213509771;
    glz(1,2) = 0.091576213509771;
    glw(1)   = 0.109951743655322;    
    % At Gauss point 2
    glz(2,1) = 0.816847572980459;
    glz(2,2) = 0.091576213509771;
    glw(2)   = 0.109951743655322;    
    % At Gauss point 3
    glz(3,1) = 0.091576213509771;
    glz(3,2) = 0.816847572980459;
    glw(3)   = 0.109951743655322;    
    % At Gauss point 4
    glz(4,1) = 0.445948490915965;
    glz(4,2) = 0.108103018168070;
    glw(4)   = 0.223381589678011;    
    % At Gauss point 5
    glz(5,1) = 0.445948490915965;
    glz(5,2) = 0.445948490915965;
    glw(5)   = 0.223381589678011;    
    % At Gauss point 6
    glz(6,1) = 0.108103018168070;
    glz(6,2) = 0.445948490915965;
    glw(6)   = 0.223381589678011;    
    % Area denominator.
    Ad = 2;    
elseif (ngp==12 && (elementnumber==2 || elementnumber==9 || elementnumber==21))
    % At Gauss point 1
    glz(1,1) = 0.06308901449150;
    glz(1,2) = 0.06308901449150;
    glw(1)    = 0.05084490637021;    
    % At Gauss point 2
    glz(2,1) = 0.87382197101700;
    glz(2,2) = 0.06308901449150;
    glw(2)   = 0.05084490637021;    
    % At Gauss point 3
    glz(3,1) = 0.06308901449150;
    glz(3,2) = 0.87382197101700;
    glw(3)   = 0.05084490637021;    
    % At Gauss point 4
    glz(4,1) = 0.31035245103378;
    glz(4,2) = 0.05314504984482;
    glw(4)   = 0.08285107561837;    
    % At Gauss point 5
    glz(5,1) = 0.63650249912140;
    glz(5,2) = 0.05314504984482;
    glw(5)   = 0.08285107561837;    
    % At Gauss point 6
    glz(6,1) = 0.63650249912140;
    glz(6,2) = 0.31035245103378;
    glw(6)   = 0.08285107561837; 
    % At Gauss point 7
    glz(7,1) = 0.31035245103378;
    glz(7,2) = 0.63650249912140;
    glw(7)   = 0.08285107561837;    
    % At Gauss point 8
    glz(8,1) = 0.05314504984482;
    glz(8,2) = 0.63650249912140;
    glw(8)   = 0.08285107561837;    
    % At Gauss point 9
    glz(9,1) = 0.05314504984482;
    glz(9,2) = 0.31035245103378;
    glw(9)   = 0.08285107561837;    
    % At Gauss point 10
    glz(10,1) = 0.24928674517091;
    glz(10,2) = 0.24928674517091;
    glw(10)   = 0.11678627572638;    
    % At Gauss point 11
    glz(11,1) = 0.50142650965818;
    glz(11,2) = 0.24928674517091;
    glw(11)   = 0.11678627572638;    
    % At Gauss point 12
    glz(12,1) = 0.24928674517091;
    glz(12,2) = 0.50142650965818;
    glw(12)   = 0.11678627572638;   
    
    % Area denominator.
    Ad = 2;
end

