clear all
clc
close all
%% Table for Zernike index transformation between j and n/m
n=300;
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
tic
n_index=2;% non-negative
m_index=0;% non-negative, m<=n, and mod(n-m,2)=0
np=(n_index-abs(m_index))/2;
Prediction=[];
for k=0:1:np
    for l=0:1:abs(m_index) % expand item (A+B)^k
        for p=0:1:k % expand item [(A+B)¡¤(A+B)]^k
            for q=0:1:k-p % expand item [(A+B)¡¤(A+B)]^k
                t=k-p-q;
                for r=0:floor(q/2) % expand item (A¡¤B)^q
                    for f=0:1:p+r
                        mode_n=l+q-2*r+2*f;
                        mode_m=l+q-2*r;
                        Prediction=[Prediction;mode_n mode_m];% Vector Zernike mode
                    end
                    if l>=q-2*r
                        for f=0:1:p+q-r
                            mode_n=l-q+2*r+2*f;
                            mode_m=l-q+2*r;
                            Prediction=[Prediction;mode_n mode_m];% Vector Zernike mode
                        end
                    else
                        for f=0:1:p+r+l
                            mode_n=q-2*r-l+2*f;
                            mode_m=q-2*r-l;
                            Prediction=[Prediction;mode_n mode_m];% Vector Zernike mode
                        end
                    end
                end
            end
        end
    end
end
[~,I]=sort(Prediction(:,2));
Prediction=Prediction(I,:);
[~,I]=sort(Prediction(:,1));
Prediction=Prediction(I,:);
[UniqueMode,~,~]=unique(Prediction,'rows','first')
