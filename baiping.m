function [dx1,dy1,dx2,dy2,dx3,dy3,sita1,ZP,YP,SP,XP,I1_hr,I2_hr,I3_hr] = baiping(I1_hr,I2_hr,I3_hr,I2_hr_3,I2_hr_3_10,I2_t_max,...
    I1_hr_10,I2_hr_10,I3_hr_10,I2_t_max_10,sita1,beta,gama,ZP,YP,SP,XP,dx1,dy1,dx2,dy2,dx3,dy3,imrf,imrf_special)
%BAIPING 此处显示有关此函数的摘要
%   此处显示详细说明
% gama为为了将图像进行水平排布为旋转的角度
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% A = logical(A);%
if dy3 ~= 0
    sita1 = sita1 + gama;
    I1_hr_10 = imrotate(I1_hr_10, beta + gama,'bicubic'); I1_hr_10(I1_hr_10 < 0.5) = 0; I1_hr_10 = logical(I1_hr_10);%(I1_hr_10 ~= 0) = 1;
    I2_hr_10 = imrotate(I2_hr_10, beta + gama,'bicubic'); I2_hr_10(I2_hr_10 < 0.5) = 0; I2_hr_10 = logical(I2_hr_10);%I2_hr_10(I2_hr_10 ~= 0) = 1;
    I3_hr_10 = imrotate(I3_hr_10, gama,'bicubic'); I3_hr_10(I3_hr_10 < 0.5) = 0; I3_hr_10 = logical(I3_hr_10);%I3_hr_10(I3_hr_10 ~= 0) = 1;
    I2_t_max_10 = imrotate(I2_t_max_10,sita1,'bicubic'); I2_t_max_10(I2_t_max_10 < 0.5) = 0; I2_t_max_10 = logical(I2_t_max_10);%I2_t_max_10(I2_t_max_10 ~= 0) = 1;
    
    I1_hr = imrotate(I1_hr, beta + gama,'bicubic'); I1_hr(I1_hr < 0.5) = 0; I1_hr = logical(I1_hr);%I1_hr(I1_hr ~= 0) = 1;
    I2_hr = imrotate(I2_hr, beta + gama,'bicubic'); I2_hr(I2_hr < 0.5) = 0; I2_hr = logical(I2_hr);%I2_hr(I2_hr ~= 0) = 1;
    I3_hr = imrotate(I3_hr, gama,'bicubic'); I3_hr(I3_hr < 0.5) = 0; I3_hr = logical(I3_hr);%I3_hr(I3_hr ~= 0) = 1;
    I2_t_max = imrotate(I2_t_max,sita1,'bicubic'); I2_t_max(I2_t_max < 0.5) = 0; I2_t_max = logical(I2_t_max);%I2_t_max(I2_t_max ~= 0) = 1;
    I2_hr_3 = imrotate(I2_hr_3, gama,'bicubic'); I2_hr_3(I2_hr_3 < 0.5) = 0; I2_hr_3 = logical(I2_hr_3);%I2_hr_3(I2_hr_3 ~= 0) = 1;
    I2_hr_3_10 = imrotate(I2_hr_3_10, gama,'bicubic'); I2_hr_3_10(I2_hr_3_10 < 0.5) = 0; I2_hr_3_10 = logical(I2_hr_3_10);%I2_hr_3_10(I2_hr_3_10 ~= 0) = 1;
    
%     I2_t_max=logical(I2_t_max);
%     I2_hr_3=logical(I2_hr_3);
%     I2_hr_3_10=logical(I2_hr_3_10);
    %下面进行判断图形旋转后属于什么性质
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %第一次判断
    I2_t_max_r = imrotate(I2_t_max,180,'bicubic'); I2_t_max_r(I2_t_max_r < 0.5) = 0; I2_t_max_r = logical(I2_t_max_r);%I2_t_max_r(I2_t_max_r ~= 0) = 1;
    dx1_abs_10 = (l_x(I1_hr_10) - l_x(I2_t_max_10)) / imrf_special;
    dy1_abs_10 = (l_y(I1_hr_10) - l_y(I2_t_max_10)) / imrf_special;
    dx1_abs = l_x(I1_hr) - l_x(I2_t_max);
    dy1_abs = l_y(I1_hr) - l_y(I2_t_max);
    
    I_base = zeros(round(1.5 * l_y(I1_hr)),round(1.5 * l_x(I1_hr)));
    %     I_base = zeros(round(2.8 * l_y(I2_t_max)),round(2.8 * l_x(I2_t_max)));
    I2_t_max = trans2lefttop(I2_t_max);
    I2_t_max = imcrop(I2_t_max,[0,0,l_x(I2_t_max),l_y(I2_t_max)]);
    I2_t_max_r = trans2lefttop(I2_t_max_r);
    I2_t_max_r = imcrop(I2_t_max_r,[0,0,l_x(I2_t_max_r),l_y(I2_t_max_r)]);
    
    I2_t_max = paste_img(I_base,I2_t_max,1,1);
    %     I2_t_max = trans2mid(I2_t_max);
    I2_t_max_r = paste_img(I_base,I2_t_max_r,1,1);
    %     I2_t_max_r = trans2mid(I2_t_max_r);
    I1_comparison = trans2lefttop(I1_hr);
    I1_comparison = imcrop(I1_comparison,[0,0,l_x(I1_comparison),l_y(I1_comparison)]);
    I1_comparison = paste_img(I_base,I1_comparison,1,1);
    I1_comparison =logical(I1_comparison);
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %1 左上
    I2_t_max = trans2lefttop(I2_t_max);
    I2_t_max_r = trans2lefttop(I2_t_max_r);
    I2_t_max_r = imtranslate(I2_t_max_r,[dx1_abs,dy1_abs]);
    img_combine1 = trans2lefttop(I2_t_max | I2_t_max_r);
    I1_C = img_combine1 | I1_comparison;
    diff1 = sum(I1_C(:)) - sum(I1_comparison(:));
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %2 左下
    I2_t_max = trans2leftbottom(I2_t_max);
    I2_t_max_r = trans2leftbottom(I2_t_max_r);
    I2_t_max_r = imtranslate(I2_t_max_r,[dx1_abs,-dy1_abs]);
    img_combine2 = trans2lefttop(I2_t_max | I2_t_max_r);
    I2_C = img_combine2 | I1_comparison;
    diff2 = sum(I2_C(:)) - sum(I1_comparison(:));
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %3 右下
    I2_t_max = trans2rightbottom(I2_t_max);
    I2_t_max_r = trans2rightbottom(I2_t_max_r);
    I2_t_max_r = imtranslate(I2_t_max_r,[-dx1_abs,-dy1_abs]);
    img_combine3 = trans2lefttop(I2_t_max | I2_t_max_r);
    I3_C = img_combine3 | I1_comparison;
    diff3 = sum(I3_C(:)) - sum(I1_comparison(:));
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %4 右上
    I2_t_max = trans2righttop(I2_t_max);
    I2_t_max_r = trans2righttop(I2_t_max_r);
    I2_t_max_r = imtranslate(I2_t_max_r,[-dx1_abs,dy1_abs]);
    img_combine4 = trans2lefttop(I2_t_max | I2_t_max_r);
    I4_C = img_combine4 | I1_comparison;
    diff4 = sum(I4_C(:)) - sum(I1_comparison(:));
    
    diff = [diff1,diff2,diff3,diff4];
    if diff1 == min(diff)
        %imshow(img_combine1);pause(1);
        dx1 =  dx1_abs_10 / imrf;
        dy1 =  dy1_abs_10 / imrf;
    elseif diff2 == min(diff)
        %imshow(img_combine2);pause(1);
        dx1 =  dx1_abs_10 / imrf;
        dy1 = -dy1_abs_10 / imrf;
    elseif diff3 == min(diff)
        %imshow(img_combine3);pause(1);
        dx1 = -dx1_abs_10 / imrf;
        dy1 = -dy1_abs_10 / imrf;
    else%diff4 == min(diff)
        %imshow(img_combine4);pause(1);
        dx1 = -dx1_abs_10 / imrf;
        dy1 =  dy1_abs_10 / imrf;
    end
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %第二次判断
    dx2_abs_10 = (l_x(I2_hr_10) - l_x(I1_hr_10)) / imrf_special;
    dy2_abs_10 = (l_y(I2_hr_10) - l_y(I1_hr_10)) / imrf_special;
    dx2_abs = l_x(I2_hr) - l_x(I1_hr);
    dy2_abs = l_y(I2_hr) - l_y(I1_hr);
    
    I_base = zeros(round(1.5 * l_y(I2_hr)),round(1.5 * l_x(I2_hr)));
    %     I_base = zeros(round(4*l_y(I1_hr)),round(4*l_x(I1_hr)));
    I1_hr = trans2lefttop(I1_hr);
    I1_hr = imcrop(I1_hr,[0,0,l_x(I1_hr),l_y(I1_hr)]);
    I1_hr = paste_img(I_base,I1_hr,1,1);
    I1_hr = trans2mid(I1_hr,0);
    I1_hr_r = I1_hr;
    I2_comparison = trans2lefttop(I2_hr);
    I2_comparison = imcrop(I2_comparison,[0,0,l_x(I2_comparison),l_y(I2_comparison)]);
    I2_comparison = paste_img(I_base,I2_comparison,1,1);
    I2_comparison =logical(I2_comparison);
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %1 左上
    I1_hr = trans2lefttop(I1_hr);
    I1_hr_r = trans2lefttop(I1_hr_r);
    I1_hr_r = imtranslate(I1_hr_r,[dx2_abs,dy2_abs]);
    img_combine1 = trans2lefttop(I1_hr | I1_hr_r);
    I1_C = img_combine1 | I2_comparison;
    diff1 = sum(I1_C(:)) - sum(I2_comparison(:));
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %2 左下
    I1_hr = trans2leftbottom(I1_hr);
    I1_hr_r = trans2leftbottom(I1_hr_r);
    I1_hr_r = imtranslate(I1_hr_r,[dx2_abs,-dy2_abs]);
    img_combine2 = trans2lefttop(I1_hr | I1_hr_r);
    I2_C = img_combine2 | I2_comparison;
    diff2 = sum(I2_C(:)) - sum(I2_comparison(:));
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %3 右下
    I1_hr = trans2rightbottom(I1_hr);
    I1_hr_r = trans2rightbottom(I1_hr_r);
    I1_hr_r = imtranslate(I1_hr_r,[-dx2_abs,-dy2_abs]);
    img_combine3 = trans2lefttop(I1_hr | I1_hr_r);
    I3_C = img_combine3 | I2_comparison;
    diff3 = sum(I3_C(:)) - sum(I2_comparison(:));
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %4 右上
    I1_hr = trans2righttop(I1_hr);
    I1_hr_r = trans2righttop(I1_hr_r);
    I1_hr_r = imtranslate(I1_hr_r,[-dx2_abs,dy2_abs]);
    img_combine4 = trans2lefttop(I1_hr | I1_hr_r);
    I4_C = img_combine4 | I2_comparison;
    diff4 = sum(I4_C(:)) - sum(I2_comparison(:));
    
    diff = [diff1,diff2,diff3,diff4];
    if diff1 == min(diff)
        %imshow(img_combine1);pause(1);
        dx2 =  dx2_abs_10 / imrf;
        dy2 =  dy2_abs_10 / imrf;
    elseif diff2 == min(diff)
        %imshow(img_combine2);pause(1);
        dx2 =  dx2_abs_10 / imrf;
        dy2 = -dy2_abs_10 / imrf;
    elseif diff3 == min(diff)
        %imshow(img_combine3);pause(1);
        dx2 = -dx2_abs_10 / imrf;
        dy2 = -dy2_abs_10 / imrf;
    else%diff4 == min(diff)
        %imshow(img_combine4);pause(1);
        dx2 = -dx2_abs_10 / imrf;
        dy2 =  dy2_abs_10 / imrf;
    end
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %第三次判断
    dx3_abs_10 = (l_x(I3_hr_10) - l_x(I2_hr_3_10)) / imrf_special;
    dy3_abs_10 = (l_y(I3_hr_10) - l_y(I2_hr_3_10)) / imrf_special;
    dx3_abs = l_x(I3_hr) - l_x(I2_hr_3);
    dy3_abs = l_y(I3_hr) - l_y(I2_hr_3);
    
    lymax = max(l_y(I2_hr_3),l_y(I3_hr));
    lxmax = max(l_x(I2_hr_3),l_x(I3_hr));
    I_base = zeros(round(1.5 * lymax),round(1.5 * lxmax));
    I2_hr_3 = trans2lefttop(I2_hr_3);
    I2_hr_3 = imcrop(I2_hr_3,[0,0,l_x(I2_hr_3),l_y(I2_hr_3)]);
    I2_hr_3 = paste_img(I_base,I2_hr_3,1,1);
    I2_hr_3 = trans2mid(I2_hr_3,0);
    I2_hr_r = I2_hr_3;
    I2_hr_r(I2_hr_r<0.5)=0;
    I2_hr_r=logical(I2_hr_r);
    I3_comparison = trans2lefttop(I3_hr);
    I3_comparison = imcrop(I3_comparison,[0,0,l_x(I3_comparison),l_y(I3_comparison)]);
    I3_comparison = paste_img(I_base,I3_comparison,1,1);
    I3_comparison =logical(I3_comparison);
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %1 左上
    I2_hr = trans2lefttop(I2_hr_r);
    I2_hr_r = trans2lefttop(I2_hr_r);
    I2_hr_r = imtranslate(I2_hr_r,[dx3_abs,dy3_abs]);
    img_combine1 = trans2lefttop(I2_hr | I2_hr_r);
    I1_C = img_combine1 | I3_comparison;
    diff1 = sum(I1_C(:)) - sum(I3_comparison(:));
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %2 左下
    I2_hr = trans2leftbottom(I2_hr);
    I2_hr_r = trans2leftbottom(I2_hr_r);
    I2_hr_r = imtranslate(I2_hr_r,[dx3_abs,-dy3_abs]);
    img_combine2 = trans2lefttop(I2_hr | I2_hr_r);
    I2_C = img_combine2 | I3_comparison;
    diff2 = sum(I2_C(:)) - sum(I3_comparison(:));
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %3 右下
    I2_hr = trans2rightbottom(I2_hr);
    I2_hr_r = trans2rightbottom(I2_hr_r);
    I2_hr_r = imtranslate(I2_hr_r,[-dx3_abs,-dy3_abs]);
    img_combine3 = trans2lefttop(I2_hr | I2_hr_r);
    I3_C = img_combine3 | I3_comparison;
    diff3 = sum(I3_C(:)) - sum(I3_comparison(:));
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %4 右上
    I2_hr = trans2righttop(I2_hr);
    I2_hr_r = trans2righttop(I2_hr_r);
    I2_hr_r = imtranslate(I2_hr_r,[-dx3_abs,dy3_abs]);
    img_combine4 = trans2lefttop(I2_hr | I2_hr_r);
    I4_C = img_combine4 | I3_comparison;
    diff4 = sum(I4_C(:)) - sum(I3_comparison(:));
    
    diff = [diff1,diff2,diff3,diff4];
    if diff1 == min(diff)
        %imshow(img_combine1);pause(1);
        dx3 =  dx3_abs_10 / imrf;
        dy3 =  dy3_abs_10 / imrf;
    elseif diff2 == min(diff)
        %imshow(img_combine2);pause(1);
        dx3 =  dx3_abs_10 / imrf;
        dy3 = -dy3_abs_10 / imrf;
    elseif diff3 == min(diff)
        %imshow(img_combine3);pause(1);
        dx3 = -dx3_abs_10 / imrf;
        dy3 = -dy3_abs_10 / imrf;
    else%diff4 == min(diff)
        %imshow(img_combine4);pause(1);
        dx3 = -dx3_abs_10 / imrf;
        dy3 =  dy3_abs_10 / imrf;
    end
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %下面判断方向
    ZP = 0;%左偏
    YP = 0;%右偏
    
    if dy2 >= 0 && dx2 >= 0
        ZP = 1;
    elseif dy2 >= 0 && dx2 < 0
        YP = 1;
    end
    
    if dy2 < 0 && dx2 >= 0
        YP = 1;
    elseif dy2<0 && dx2 < 0
        ZP = 1;
    end
    
    if ZP == 1
        dx2 = -abs(dx2);
        dy2 = -abs(dy2);
    else
        dx2 =  abs(dx2);
        dy2 = -abs(dy2);
    end
    
    SP = 0;%上偏
    XP = 0;%下偏
    
    %     if dx3 >= 0 && dy3 >= 0
    %         XP = 1;
    %     elseif dx3 >= 0 && dy3 < 0
    %         SP = 1;
    %     end
    %
    %     if dx3<0 && dy3 >= 0
    %         SP = 1;
    %     elseif dx3<0 && dy3 < 0
    %         XP = 1;
    %     end
    
    %     if SP == 1
    %         dx3 = abs(dx3);
    %         dy3 = -abs(dy3);
    %     else
    %         dx3 = abs(dx3);
    %         dy3 = abs(dy3);
    %     end
    beishu = round(dy3/dy2);%dy2始终是负数，不可能是0
    dx3 = abs(abs(dx3) - beishu * dx2);
    %     dy3 = round(dy3/dy2) * dy2;
    dy3 = 0;
else
    
    I1_hr = imrotate(I1_hr, beta + gama,'bicubic'); I1_hr(I1_hr < 0.5) = 0; I1_hr(I1_hr ~= 0) = 1;
    I2_hr = imrotate(I2_hr, beta + gama,'bicubic'); I2_hr(I2_hr < 0.5) = 0; I2_hr(I2_hr ~= 0) = 1;
    I3_hr = imrotate(I3_hr, gama,'bicubic'); I3_hr(I3_hr < 0.5) = 0; I3_hr(I3_hr ~= 0) = 1;
end
end

