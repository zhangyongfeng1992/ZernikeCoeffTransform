clear all
clc
close all
%% Table for Zernike index transformation between j and n/m
n=2;
jnm=[];
counter=1;
for n_index=0:n
    for m_index=mod(n_index,2):2:n_index
        if m_index==0
            jnm=[jnm;counter,n_index,m_index];
            counter=counter+1;
        else
            jnm=[jnm;counter,n_index,m_index;counter+1,n_index,-m_index];
            counter=counter+2;
        end
    end
end
%% Computation of transformation matrix
syms dx dy alpha s real
d=[dx;dy];
% syms D beta alpha s real
% d=[D*cos(beta);D*sin(beta)];
e_a=[cos(alpha);sin(alpha)];
T=sym(zeros(jnm(end,1)));
syms A1 A2 real;
tic
for j_index=1:size(jnm,1)
    n_index=jnm(j_index,2);
    m_index=jnm(j_index,3);
    if m_index>=0 % for case of m_index>=0
        A=[A1;A2];
        np=(n_index-abs(m_index))/2;
        for k=0:1:np
            for l=0:1:abs(m_index) % expand item (A+B)^k
                for p=0:1:k % expand item [(A+B)¡¤(A+B)]^k
                    for q=0:1:k-p % expand item [(A+B)¡¤(A+B)]^k
                        t=k-p-q;
                        for r=0:floor(q/2) % expand item (A¡¤B)^q
                            Ep=Ep_fun(np,abs(m_index),k,l,p,q,t,r);
                            for f=0:1:p+r
                                mode_n=l+q-2*r+2*f;
                                mode_m=l+q-2*r;
                                j_mode=find((jnm(:,2)==mode_n).*(jnm(:,3)==mode_m));% index j for derived mode
                                coeff1=s^(2*p+q+l)*C_fun(p+r,mode_m,mode_n)*Ep*DOT(d,d)^(t+r)*...
                                    BING(A,BING(BINGM(CONJ(e_a),l+q-2*r),BING(BINGM(d,q-2*r),BINGM(CONJ(d),m_index-l))));
                                if m_index==0 % for symmetric original mode
                                    if mode_m==0 % for symmetric derived mode
                                        T(j_mode,j_index)=T(j_mode,j_index)+subs(coeff1(1),A2,0)/A1;
                                    else % for asymmetric derived mode
                                        T(j_mode:j_mode+1,j_index)=T(j_mode:j_mode+1,j_index)+subs(coeff1,A2,0)/A1;
                                    end
                                else % for asymmetric original mode
                                    if mode_m==0 % for symmetric derived mode
                                        T(j_mode,j_index)=T(j_mode,j_index)+subs(coeff1(1),A2,0)/A1;
                                        T(j_mode,j_index+1)=T(j_mode,j_index+1)+subs(coeff1(1),A1,0)/A2;
                                    else % for asymmetric derived mode
                                        T(j_mode:j_mode+1,j_index)=T(j_mode:j_mode+1,j_index)+subs(coeff1,A2,0)/A1;
                                        T(j_mode:j_mode+1,j_index+1)=T(j_mode:j_mode+1,j_index+1)+subs(coeff1,A1,0)/A2;
                                    end
                                end
                            end
                            if l>=q-2*r
                                for f=0:1:p+q-r
                                    mode_n=l-q+2*r+2*f;
                                    mode_m=l-q+2*r;
                                    j_mode=find((jnm(:,2)==mode_n).*(jnm(:,3)==mode_m));% index j for derived mode
                                    coeff2=s^(2*p+q+l)*C_fun(p+q-r,mode_m,mode_n)*Ep*DOT(d,d)^(t+r)*...
                                        BING(A,BING(BINGM(CONJ(e_a),l-q+2*r),BINGM(CONJ(d),m_index-l+q-2*r)));
                                    if m_index==0 % for symmetric original mode
                                        if mode_m==0 % for symmetric derived mode
                                            T(j_mode,j_index)=T(j_mode,j_index)+subs(coeff2(1),A2,0)/A1;
                                        else % for asymmetric derived mode
                                            T(j_mode:j_mode+1,j_index)=T(j_mode:j_mode+1,j_index)+subs(coeff2,A2,0)/A1;
                                        end
                                    else % for asymmetric original mode
                                        if mode_m==0 % for symmetric derived mode
                                            T(j_mode,j_index)=T(j_mode,j_index)+subs(coeff2(1),A2,0)/A1;
                                            T(j_mode,j_index+1)=T(j_mode,j_index+1)+subs(coeff2(1),A1,0)/A2;
                                        else % for asymmetric derived mode
                                            T(j_mode:j_mode+1,j_index)=T(j_mode:j_mode+1,j_index)+subs(coeff2,A2,0)/A1;
                                            T(j_mode:j_mode+1,j_index+1)=T(j_mode:j_mode+1,j_index+1)+subs(coeff2,A1,0)/A2;
                                        end
                                    end
                                end
                            else
                                for f=0:1:p+r+l
                                    mode_n=q-2*r-l+2*f;
                                    mode_m=q-2*r-l;
                                    j_mode=find((jnm(:,2)==mode_n).*(jnm(:,3)==mode_m));% index j for derived mode
                                    coeff2=s^(2*p+q+l)*C_fun(p+r+l,mode_m,mode_n)*Ep*DOT(d,d)^(t+r)*...
                                        BING(CONJ(A),BING(BINGM(CONJ(e_a),q-2*r-l),BINGM(d,m_index-l+q-2*r)));
                                    if m_index==0 % for symmetric original mode
                                        if mode_m==0 % for symmetric derived mode
                                            T(j_mode,j_index)=T(j_mode,j_index)+subs(coeff2(1),A2,0)/A1;
                                        else % for asymmetric derived mode
                                            T(j_mode:j_mode+1,j_index)=T(j_mode:j_mode+1,j_index)+subs(coeff2,A2,0)/A1;
                                        end
                                    else % for asymmetric original mode
                                        if mode_m==0 % for symmetric derived mode
                                            T(j_mode,j_index)=T(j_mode,j_index)+subs(coeff2(1),A2,0)/A1;
                                            T(j_mode,j_index+1)=T(j_mode,j_index+1)+subs(coeff2(1),A1,0)/A2;
                                        else % for asymmetric derived mode
                                            T(j_mode:j_mode+1,j_index)=T(j_mode:j_mode+1,j_index)+subs(coeff2,A2,0)/A1;
                                            T(j_mode:j_mode+1,j_index+1)=T(j_mode:j_mode+1,j_index+1)+subs(coeff2,A1,0)/A2;
                                        end
                                    end
                                end
                            end
                        end
                    end
                end
            end
        end
    else % for case of m_index<0. It could be skipped, because the Zernike terms for m_index>0 and m_index<0 are combined into a vector and dealted together!
        continue;
    end
end
toc
T
