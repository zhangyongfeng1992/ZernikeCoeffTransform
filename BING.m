function c=BING(a,b)
% The input,a and b, and the output, c, are column vectors
c1=a(1)*b(1)-a(2)*b(2);
c2=a(1)*b(2)+a(2)*b(1);
c=[c1;c2];
end