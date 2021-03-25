function C=C_fun_num(Left_p,Right_U_q,Right_D_q)
% This function is different from the function C_fun, which is for symbol
% computation;
t=(Right_D_q-Right_U_q)/2;
C=0;
parfor s=0:t
    C=C+E_fun(t,Right_U_q,s)*(Right_U_q+2*t+1)/(Left_p+Right_U_q+s+1);
end
end