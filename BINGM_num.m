function b=BINGM_num(a,n)
% This function is different from the function BINGM, which is for symbol
% cumputation;
% The input,a, and the output, b, are column vectors
% n is the power exponent of input a
MAG=sqrt(a(1)^2+a(2)^2);
theta=atan2(a(2),a(1));
b=[MAG^n*cos(n*theta);MAG^n*sin(n*theta)];
end