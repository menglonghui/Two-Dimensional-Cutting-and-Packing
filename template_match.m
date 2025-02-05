% function [M] = template_match(I_mix,I_single)
I_mix_re=~I_mix;

I=double(I_mix_re);
[Ir, Ic] = size(I);
figure(1),imshow(I)
T = trans2lefttop(I2_t_max_or);
T = imcrop(T,[0, 0, l_x(T), l_y(T)]);
I_base = zeros(size(I_mix));
T_tem = paste_img(I_base,T,1,1);
% T = imrotate(T,90,'bicubic');
[Tr, Tc] = size(T);
figure(2),imshow(T)
R = normxcorr2(T,I);
R = imcrop(R,[Tc Tr Ic Ir]);
R_sorted = sort(R(:), 'descend');
tic
a = 0;
% rv_notok = zeros(1,1);
rv_notok = cell(1, 1);
rv_notok{1}=[10000,10000];
rv_num = 0;
while a <= length(R_sorted)
    a = a + 100;
    [r, c, v] = find(R == R_sorted(a));
    for xulie1 = 1: length(c)
        for xulie2 =1:length(rv_notok)
            if (r(xulie1) == rv_notok{xulie2}(1)) && (c(xulie1) == rv_notok{xulie2}(2))
                ...
            else
            T_moved = imtranslate(T_tem,[c(xulie1),r(xulie1)]);
            I_com = I_mix | T_moved;
            
            %     for pp = 1:length(rv_notok)
            %
            %
            %     end
            
            if sum(I_com(:)) == sum(I_mix(:)) + sum(T_tem(:))
                
                a = 0;
                I_mix = I_mix | T_moved;
                I_mix_re = ~I_mix;
                I = double(I_mix_re);
                R = normxcorr2(T,I);
                R = imcrop(R,[Tc Tr Ic Ir]);
                R_sorted = sort(R(:), 'descend');
                imshow(I_com)
                
            else
                rv_num = rv_num + 1;
                rv_notok{rv_num} = [(c(xulie1)),r(xulie1)];
                %         rv_notok(rv_num,2) = r(1);
            end
            end
        end
    end
end