function Ep=Ep_fun(np,m,k,l,p,q,t,r)
Ep=E_fun(np,m,k)/(1+(r==floor(q/2))*(mod(q,2)==0))*nchoosek(m,l)/...
    factorial(r)/factorial(q-r)*factorial(k)/factorial(p)/factorial(t);
end