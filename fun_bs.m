function result = fun_bs(I)
%FUN_BS 此处显示有关此函数的摘要
%   此处显示详细说明
% I = trans2lefttop(I);
% I = imcrop(I,[0,0,l_x(I),l_y(I)]);
I1 = I;
I2 = I;
I3 = I;
I4 = I;

% S = numel(I(:,1)) * numel(I(1,:));
for i = 1:numel(I(:,1))
    for j = 1:numel(I(1,:))
        if I1(i,j) == 0
            I1(i,j) = 1;
        else
            break
        end
    end
end

for i = 1:numel(I(:,1))
    for j = numel(I(1,:)):(-1):1
        if I2(i,j) == 0
            I2(i,j) = 1;
        else
            break
        end
    end
end

for j = 1:numel(I(1,:))
    for i = 1:numel(I(:,1))
        if I3(i,j) == 0
            I3(i,j) = 1;
        else
            break
        end
    end
end

for j = 1:numel(I(1,:))
    for i = numel(I(:,1)):(-1):1
        if I4(i,j) == 0
            I4(i,j) = 1;
        else
            break
        end
    end
end

% N_bao = S-sum(I(:));
% result = N_bao + S;
% result = sum(I1(:)) + sum(I2(:)) + sum(I3(:)) + sum(I4(:));
I12 = (I1 | I2);
I34 = (I3 | I4);
I1234=I12 | I34;
% result = - sum(I1(:)) - sum(I2(:)) - sum(I3(:))- sum(I4(:)) + l_x(I)/10 + l_y(I)/10;
result = - sum(I1234(:))*100 + l_x(I)/100000 + l_y(I)/100000;
% result = - sum(I1234(:)) + l_x(I)/10 + l_y(I)/10;
end

