function N_wai = fun_bao1(I)
%FUN_Y �˴���ʾ�йش˺�����ժҪ
%   �˴���ʾ��ϸ˵��
%   �˴���ʾ��ϸ˵��
N_wai=0;
[N_shu,N_heng]=size(I);

for i=1:N_shu
    for j=1:N_heng
        if I(i,j)==0%����������
            N_wai=N_wai+1;
        else
            break
        end
    end
    
    for j=N_heng:-1:1
        if I(i,j)==0%����������
            N_wai=N_wai+1;
        else
            break
        end
    end
end


for i=1:N_heng
    for j=1:N_shu
        if I(j,i)==0%����������
            N_wai=N_wai+1;
        else
            break
        end
    end
    
    for j=N_shu:-1:1
        if I(j,i)==0%����������
            N_wai=N_wai+1;
        else
            break
        end
    end
end


N_wai=N_wai+l_x(I)/10000+l_y(I)/10000;

% N_wai=N_wai*1000+l_x(I)/25+l_y(I);

% N_wai=N_wai;



end

