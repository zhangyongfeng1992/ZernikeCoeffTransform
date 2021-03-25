function C=C_fun(Left_p,Right_U_q,Right_D_q)
t=(Right_D_q-Right_U_q)/2;
syms s real
C=symsum(E_fun(t,Right_U_q,s)*(Right_U_q+2*t+1)/(Left_p+Right_U_q+s+1),s,0,t);
end