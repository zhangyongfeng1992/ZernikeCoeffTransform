function b=BINGM(a,n)
% The input,a, and the output, b, are column vectors
% n is the power exponent of input a
if n>1
    b=BING(a,BINGM(a,n-1));
elseif n==1
    b=a;
else
    b=sym([1;0]);
end
% b=simplify(b);
end