function E=E_fun(np,m,k)
E=(-1)^(np-k)*factorial(np+abs(m)+k)/factorial(k)/factorial(np-k)/factorial(abs(m)+k);
end