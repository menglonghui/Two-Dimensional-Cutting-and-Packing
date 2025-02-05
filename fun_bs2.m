function N = fun_bs2(I,I2_t)
%FUN_Y 此处显示有关此函数的摘要
%   此处显示详细说明
%   此处显示详细说明
if ((l_x(I) > 1.85 * l_x(I2_t)) && (l_y(I) > 1.85 * l_y(I2_t))) || (l_x(I) > 2 * l_x(I2_t)) || (l_y(I) > 2 * l_y(I2_t))
    N = numel(I(:,1)) * numel(I(1,:));
else
    I_h=zeros(size(I));
    I_s=zeros(size(I));
    
    [N_shu,N_heng]=size(I);
    % S=N_shu * N_heng;
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %横向分析
    for i=1:N_shu
        for j=1:N_heng
            if I(i,j)==0%碰到黑像素
                I_h(i,j)=1;
            else
                break
            end
        end
        
        for j=N_heng:-1:1
            if I(i,j)==0%碰到黑像素
                I_h(i,j)=1;
            else
                break
            end
        end
    end
    
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %竖向分析
    for j=1:N_heng
        for i=1:N_shu
            if I(i,j)==0%碰到黑像素
                I_s(i,j)=1;
            else
                break
            end
        end
        
        for i=N_shu:-1:1
            if I(i,j)==0%碰到黑像素
                I_s(i,j)=1;
            else
                break
            end
        end
    end
    
    % I_ime=I_h | I_s;
    % imshow(I_ime)
    % N=-sum(I_s(:))-sum(I_h(:)) + l_x(I)/10000 + l_y(I)/10000;
    I1 = I_s | I;
    I2 = I_h | I;
    I_c = I1 & I2;
    I_c =~ I_c;
    N = sum(I_c(:));
end
end

