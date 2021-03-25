clear all
clc
close all
%% Define the parameters
syms s real;% Scaling factor
syms dx dy real;% Normalized translation distance along x and y axis
syms alpha real;% Rotation angle
syms r0 theta real;% Normalized polar radius and polar angle of the transformed pupil
% Coordinates transformation
x=s*r0*cos(theta+alpha)+dx;
y=s*r0*sin(theta+alpha)+dy;
r=sqrt(x^2+y^2);
% Cartesian representation of 36 terms of standard Zernike polynomials
Z=[1
    x
    y
    2*r^2-1
    2*x^2-r^2
    2*x*y
    x*(3*r^2-2)
    y*(3*r^2-2)
    4*x^3-3*r^2*x
    -y*(r^2-4*x^2)
    6*r^4-6*r^2+1
    -(r^2-2*x^2)*(4*r^2-3)
    2*x*y*(4*r^2-3)
    r^4-8*r^2*x^2+8*x^4
    -4*x*y*(r^2-2*x^2)
    x*(10*r^4-12*r^2+3)
    y*(10*r^4-12*r^2+3)
    -x*(5*r^2-4)*(3*r^2-4*x^2)
    -y*(r^2-4*x^2)*(5*r^2-4)
    5*r^4*x-20*r^2*x^3+16*x^5
    y*(r^4-12*r^2*x^2+16*x^4)
    20*r^6-30*r^4+12*r^2-1
    -(r^2-2*x^2)*(15*r^4-20*r^2+6)
    2*x*y*(15*r^4-20*r^2+6)
    (6*r^2-5)*(r^4-8*r^2*x^2+8*x^4)
    -4*x*y*(r^2-2*x^2)*(6*r^2-5)
    -r^6+18*r^4*x^2-48*r^2*x^4+32*x^6
    2*x*y*(3*r^4-16*r^2*x^2+16*x^4)
    x*(35*r^6-60*r^4+30*r^2-4)
    y*(35*r^6-60*r^4+30*r^2-4)
    -x*(3*r^2-4*x^2)*(21*r^4-30*r^2+10)
    -y*(r^2-4*x^2)*(21*r^4-30*r^2+10)
    x*(7*r^2-6)*(5*r^4-20*r^2*x^2+16*x^4)
    y*(7*r^2-6)*(r^4-12*r^2*x^2+16*x^4)
    -7*r^6*x+56*r^4*x^3-112*r^2*x^5+64*x^7
    -y*(r^6-24*r^4*x^2+80*r^2*x^4-64*x^6)];
%Polar representation of 36 terms of standard Zernike polynomials
Z0=[1
    r0*cos(theta)
    r0*sin(theta)
    2*r0^2-1
    r0^2*cos(2*theta)
    r0^2*sin(2*theta)
    (3*r0^3-2*r0)*cos(theta)
    (3*r0^3-2*r0)*sin(theta)
    r0^3*cos(3*theta)
    r0^3*sin(3*theta)
    6*r0^4-6*r0^2+1
    (4*r0^4-3*r0^2)*cos(2*theta)
    (4*r0^4-3*r0^2)*sin(2*theta)
    r0^4*cos(4*theta)
    r0^4*sin(4*theta)
    (10*r0^5-12*r0^3+3*r0)*cos(theta)
    (10*r0^5-12*r0^3+3*r0)*sin(theta)
    (5*r0^5-4*r0^3)*cos(3*theta)
    (5*r0^5-4*r0^3)*sin(3*theta)
    r0^5*cos(5*theta)
    r0^5*sin(5*theta)
    20*r0^6-30*r0^4+12*r0^2-1
    (15*r0^6-20*r0^4+6*r0^2)*cos(2*theta)
    (15*r0^6-20*r0^4+6*r0^2)*sin(2*theta)
    (6*r0^6-5*r0^4)*cos(4*theta)
    (6*r0^6-5*r0^4)*sin(4*theta)
    r0^6*cos(6*theta)
    r0^6*sin(6*theta)
    (35*r0^7-60*r0^5+30*r0^3-4*r0)*cos(theta)
    (35*r0^7-60*r0^5+30*r0^3-4*r0)*sin(theta)
    (21*r0^7-30*r0^5+10*r0^3)*cos(3*theta)
    (21*r0^7-30*r0^5+10*r0^3)*sin(3*theta)
    (7*r0^7-6*r0^5)*cos(5*theta)
    (7*r0^7-6*r0^5)*sin(5*theta)
    r0^7*cos(7*theta)
    r0^7*sin(7*theta)];
%Norms of Zernike polynomials
N=[pi;pi/4;pi/4;pi/3;pi/6;pi/6;pi/8;pi/8;pi/8;pi/8;pi/5;pi/10;pi/10;pi/10;...
    pi/10;pi/12;pi/12;pi/12;pi/12;pi/12;pi/12;pi/7;pi/14;pi/14;pi/14;...
    pi/14;pi/14;pi/14;pi/16;pi/16;pi/16;pi/16;pi/16;pi/16;pi/16;pi/16];
%Transformation matrix
n=7;
totalN=n*(n+3)/2+1;
T1=sym(zeros(totalN,totalN));
tic
parfor i=1:totalN
    for j=1:totalN
        T1(i,j)=int(int(expand(Z0(i)*Z(j)*r0/N(i)),r0,0,1),theta,0,2*pi);
    end
end
toc
