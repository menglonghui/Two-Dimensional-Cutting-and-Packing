function [I_mix,coordinate_addi_top_single_by_single,N_left] = addi_top_single_by_single(I_mix,I1,I2_t_max_or,heiban,dx1,dy1,...
    sita1,N_left)
%ADDI_TOP 此处显示有关此函数的摘要
%   此处显示详细说明
%先进行预处理，截取上面的部分进行分析
if sita1 == 90 || sita1 == 270
    xz = 90;
else
    xz = 0;
end

I1 = trans2lefttop(I1);
I1 = imcrop(I1,[0,0,l_x(I1),l_y(I1)]);
I1 = paste_img(heiban,I1,1,1);
I1_temp = imtranslate(I1,[l_x(I1),0]);
I1_temp0 = I1_temp;
I_com = I1 | I1_temp;
ii = 0;
while sum(I_com(:)) == 2 * sum(I1(:))
    ii = ii - 1;
    I1_temp = imtranslate(I1_temp0,[ii,0]);
    I_com = I1 | I1_temp;
end

dx2_t = l_x(I1) + ii + 1;

if sita1 == 180 || sita1 == 270
    dx1 = - dx1;
    dy1 = - dy1;
end

%dx1和dy1是翻转后的图像相对于未翻转图像的距离

%coordinate for I2_t_max_or
if dx1 >= 0 && dy1 >= 0
    
    dx0 = l_x(I2_t_max_or) / 2;
    dy0 = numel(heiban(:,1)) - l_y(I2_t_max_or) / 2;
    
    dx = dx0;
    dy = dy0;
    ci = 0;
    while dx <= (numel(heiban(1,:)) - l_x(I2_t_max_or) / 2)
        ci = ci + 1;
        coordinate_zheng(ci,1) = dx;
        coordinate_zheng(ci,2) = dy;
        coordinate_zheng(ci,3) = 0;
        dx = dx + dx2_t;
    end
    
    dx = dx0 + dx1;
    dy = dy0 - dy1;
    ci = 0;
    while dx <= (numel(heiban(1,:)) - l_x(I2_t_max_or) / 2)
        ci = ci + 1;
        coordinate_fan(ci,1) = dx;
        coordinate_fan(ci,2) = dy;
        coordinate_fan(ci,3) = 180;
        dx = dx + dx2_t;
    end
    
elseif dx1 >= 0 && dy1 < 0
    
    dx0 = l_x(I2_t_max_or) / 2;
    dy0 = numel(heiban(:,1)) - l_y(I2_t_max_or) / 2 + dy1;
    
    dx = dx0;
    dy = dy0;
    ci = 0;
    while dx <= (numel(heiban(1,:)) - l_x(I2_t_max_or) / 2)
        ci = ci + 1;
        coordinate_zheng(ci,1) = dx;
        coordinate_zheng(ci,2) = dy;
        coordinate_zheng(ci,3) = 0;
        dx = dx + dx2_t;
    end
    
    dx = dx0 + dx1;
    dy = dy0 - dy1;
    ci = 0;
    while dx <= (numel(heiban(1,:)) - l_x(I2_t_max_or) / 2)
        ci = ci + 1;
        coordinate_fan(ci,1) = dx;
        coordinate_fan(ci,2) = dy;
        coordinate_fan(ci,3) = 180;
        dx = dx + dx2_t;
    end
    
elseif dx1 < 0 && dy1 >= 0
    
    dx0 = l_x(I2_t_max_or) / 2 - dx1;
    dy0 = numel(heiban(:,1)) - l_y(I2_t_max_or) / 2;
    
    dx = dx0;
    dy = dy0;
    ci = 0;
    while dx <= (numel(heiban(1,:)) - l_x(I2_t_max_or) / 2)
        ci = ci + 1;
        coordinate_zheng(ci,1) = dx;
        coordinate_zheng(ci,2) = dy;
        coordinate_zheng(ci,3) = 0;
        dx = dx + dx2_t;
    end
    
    dx = dx0 + dx1;
    dy = dy0 - dy1;
    ci = 0;
    while dx <= (numel(heiban(1,:)) - l_x(I2_t_max_or) / 2)
        ci = ci + 1;
        coordinate_fan(ci,1) = dx;
        coordinate_fan(ci,2) = dy;
        coordinate_fan(ci,3) = 180;
        dx = dx + dx2_t;
    end
    
elseif dx1 < 0 && dy1 < 0
    
    dx0 = l_x(I2_t_max_or) / 2 - dx1;
    dy0 = numel(heiban(:,1)) - l_y(I2_t_max_or) / 2 + dy1;
    
    dx = dx0;
    dy = dy0;
    ci = 0;
    while dx <= (numel(heiban(1,:)) - l_x(I2_t_max_or) / 2)
        ci = ci + 1;
        coordinate_zheng(ci,1) = dx;
        coordinate_zheng(ci,2) = dy;
        coordinate_zheng(ci,3) = 0;
        dx = dx + dx2_t;
    end
    
    dx = dx0 + dx1;
    dy = dy0 - dy1;
    ci = 0;
    while dx <= (numel(heiban(1,:)) - l_x(I2_t_max_or) / 2)
        ci = ci + 1;
        coordinate_fan(ci,1) = dx;
        coordinate_fan(ci,2) = dy;
        coordinate_fan(ci,3) = 180;
        dx = dx + dx2_t;
    end
end

coordinate_addi_top_single_by_single = [coordinate_zheng; coordinate_fan];
[coordinate_addi_top_single_by_single,~] = sortrows(coordinate_addi_top_single_by_single, [1 2]);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if exist('coordinate_addi_top_single_by_single','var')==0%如果该值不存在
    coordinate_addi_top_single_by_single = [];
end

if isempty(coordinate_addi_top_single_by_single) ==  0%说明其不为空
    if length(coordinate_addi_top_single_by_single) >= N_left%可以通过这个方式把零件排完
        coordinate_addi_top_single_by_single = coordinate_addi_top_single_by_single(1:N_left,:);
        N_left = 0;
    else %不能通过这个方式把零件排完
        N_left = N_left - length(coordinate_addi_top_single_by_single);
        %         coordinate_addi_top_single_by_single = coordinate_addi_top_single_by_single;
    end
    
    I2_t_max_or = trans2lefttop(I2_t_max_or);
    T_0 = paste_img(heiban,I2_t_max_or,1,1);
    T_0 = trans2mid(T_0,0);
    
    for kk = 1 : length(coordinate_addi_top_single_by_single(:,1))
        T_moved = imrotate(T_0,coordinate_addi_top_single_by_single(kk,3));
        T_moved = trans2leftbottom(T_moved);
        
        dx_move = coordinate_addi_top_single_by_single(kk,1) - l_x(I2_t_max_or) / 2;
        dy_move = -(coordinate_addi_top_single_by_single(kk,2) - l_y(I2_t_max_or) / 2);
        
        T_moved = imtranslate(T_moved,[dx_move, dy_move]);
        I_mix = I_mix | T_moved;
    end
    %     imshow(I_mix)
    I_mix_show = imresize(I_mix,0.5);
    imshow(I_mix_show,'Border','tight')
    set(gcf,'ToolBar','none','MenuBar','none','NumberTitle','off');
    pause(0.0001)
end

coordinate_addi_top_single_by_single(:,3) = coordinate_addi_top_single_by_single(:,3) + xz;

end