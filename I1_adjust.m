function [dx1,dy1] = I1_adjust(I2_t_max_10,I1_hr_10,I2_t_max,I1_hr,imrf_special,imrf,angle)
%I1_ADJUST 此处显示有关此函数的摘要
%   此处显示详细说明
I2_t_max_10 = double(logical(I2_t_max_10));
I2_t_max = double(logical(I2_t_max));
I2_t_max_r = I2_t_max;
I1_hr_10 = double(logical(I1_hr_10));
I1_hr = double(logical(I1_hr));

I2_t_max_10 = trans2mid(I2_t_max_10,1);
I2_t_max = trans2mid(I2_t_max,1);
I2_t_max_r = trans2mid(I2_t_max_r,1);
I1_hr_10 = trans2mid(I1_hr_10,1);
I1_hr = trans2mid(I1_hr,1);

I2_t_max_r = imrotate(I2_t_max_r,180,'bicubic');

I2_t_max_10_ro = imrotate(I2_t_max_10,angle,'bicubic');
I2_t_max_ro = imrotate(I2_t_max,angle,'bicubic');
I2_t_max_r_ro = imrotate(I2_t_max_r,angle,'bicubic');
I1_hr_10_ro = imrotate(I1_hr_10,angle,'bicubic');
I1_hr_ro = imrotate(I1_hr,angle,'bicubic');

I2_t_max_10_ro(I2_t_max_10_ro < 0.5) = 0;
I2_t_max_ro(I2_t_max_ro < 0.5) = 0;
I2_t_max_r_ro(I2_t_max_r_ro < 0.5) = 0;
I1_hr_10_ro(I1_hr_10_ro < 0.5) = 0;
I1_hr_ro(I1_hr_ro < 0.5) = 0;

I2_t_max_10_ro = logical(I2_t_max_10_ro);
I2_t_max_ro = logical(I2_t_max_ro);
I2_t_max_r_ro = logical(I2_t_max_r_ro);
I1_hr_10_ro = logical(I1_hr_10_ro);
I1_hr_ro = logical(I1_hr_ro);

dx1_abs_10 = (l_x(I1_hr_10_ro) - l_x(I2_t_max_10_ro)) / imrf_special;
dy1_abs_10 = (l_y(I1_hr_10_ro) - l_y(I2_t_max_10_ro)) / imrf_special;
dx1_abs = l_x(I1_hr_ro) - l_x(I2_t_max_ro);
dy1_abs = l_y(I1_hr_ro) - l_y(I2_t_max_ro);

I_base = zeros(round(2.8 * l_y(I2_t_max_ro)),round(2.8 * l_x(I2_t_max_ro)));
I2_t_max_ro = trans2lefttop(I2_t_max_ro);
I2_t_max_ro = imcrop(I2_t_max_ro,[0,0,l_x(I2_t_max_ro),l_y(I2_t_max_ro)]);
I2_t_max_r_ro = trans2lefttop(I2_t_max_r_ro);
I2_t_max_r_ro = imcrop(I2_t_max_r_ro,[0,0,l_x(I2_t_max_r_ro),l_y(I2_t_max_r_ro)]);

I2_t_max_ro = paste_img(I_base,I2_t_max_ro,1,1);
I2_t_max_r_ro = paste_img(I_base,I2_t_max_r_ro,1,1);

I1_comparison = trans2lefttop(I1_hr_ro);
I1_comparison = imcrop(I1_comparison,[0,0,l_x(I1_comparison),l_y(I1_comparison)]);
I1_comparison = paste_img(I_base,I1_comparison,1,1);
I1_comparison =logical(I1_comparison);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%1 左上
I2_t_max_ro = trans2lefttop(I2_t_max_ro);
I2_t_max_r_ro = trans2lefttop(I2_t_max_r_ro);
I2_t_max_r_ro = imtranslate(I2_t_max_r_ro,[dx1_abs,dy1_abs]);
img_combine1 = trans2lefttop(I2_t_max_ro | I2_t_max_r_ro);
I1_C = img_combine1 | I1_comparison;
diff1 = sum(I1_C(:)) - sum(I1_comparison(:));
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%2 左下
I2_t_max_ro = trans2leftbottom(I2_t_max_ro);
I2_t_max_r_ro = trans2leftbottom(I2_t_max_r_ro);
I2_t_max_r_ro = imtranslate(I2_t_max_r_ro,[dx1_abs,-dy1_abs]);
img_combine2 = trans2lefttop(I2_t_max_ro | I2_t_max_r_ro);
I2_C = img_combine2 | I1_comparison;
diff2 = sum(I2_C(:)) - sum(I1_comparison(:));
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%3 右下
I2_t_max_ro = trans2rightbottom(I2_t_max_ro);
I2_t_max_r_ro = trans2rightbottom(I2_t_max_r_ro);
I2_t_max_r_ro = imtranslate(I2_t_max_r_ro,[-dx1_abs,-dy1_abs]);
img_combine3 = trans2lefttop(I2_t_max_ro | I2_t_max_r_ro);
I3_C = img_combine3 | I1_comparison;
diff3 = sum(I3_C(:)) - sum(I1_comparison(:));
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%4 右上
I2_t_max_ro = trans2righttop(I2_t_max_ro);
I2_t_max_r_ro = trans2righttop(I2_t_max_r_ro);
I2_t_max_r_ro = imtranslate(I2_t_max_r_ro,[-dx1_abs,dy1_abs]);
img_combine4 = trans2lefttop(I2_t_max_ro | I2_t_max_r_ro);
I4_C = img_combine4 | I1_comparison;
diff4 = sum(I4_C(:)) - sum(I1_comparison(:));

diff = [diff1,diff2,diff3,diff4];
if diff1 == min(diff)
    dx1 =  dx1_abs_10 / imrf;
    dy1 =  dy1_abs_10 / imrf;
elseif diff2 == min(diff)
    dx1 =  dx1_abs_10 / imrf;
    dy1 = -dy1_abs_10 / imrf;
elseif diff3 == min(diff)
    dx1 = -dx1_abs_10 / imrf;
    dy1 = -dy1_abs_10 / imrf;
else%diff4 == min(diff)
    dx1 = -dx1_abs_10 / imrf;
    dy1 =  dy1_abs_10 / imrf;
end

end

