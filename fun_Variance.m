function variance= fun_Variance(I)
%FUN_VARIANCE �˴���ʾ�йش˺�����ժҪ
%���㷽��
    variance=0;
    [L,W]=size(I);
    L_m=0;
    W_m=0;

    N=0;%��ɫ�����ܺ�

    for i=1:L
        for j=1:W
            if I(i,j)==1
                N=N+1;
                L_m=L_m+i;
                W_m=W_m+j;
            end
        end
    end

    L_m=L_m/N;
    W_m=W_m/N;

    for i=1:L
        for j=1:W
            if I(i,j)==1
                variance=variance+(i-L_m)^2+(j-W_m)^2;
            end
        end
    end

    variance=variance/N;

end

