function [coordinate_addi_top] = addi_top(dx3,dy2,coordinate1_pre,coordinate2_pre,middle_coordinate_x,middle_coordinate_y,...
    heiban,I_s_m0,I_s_m_r0,I2_t_temp0,I2_t_max_t_10_noexpand,I1,lx_single,ly_single,imrf,imrf_special,I2_t_max_10,I1_hr_10,I2_t_max,I1_hr)
%ADDI_TOP 此处显示有关此函数的摘要
%   此处显示详细说明
%先进行预处理，截取上面的部分进行分析
I_s_m0 = double(I_s_m0);
I_s_m_r0 = double(I_s_m_r0);
I_mix_pre = heiban;

for i =1 : length(coordinate1_pre(:,1))
    if  coordinate1_pre(i,2) <= -numel(heiban(:,1)) + ly_single + 2 * abs(dy2)
        I_s_m = imtranslate(I_s_m0,[coordinate1_pre(i,1),coordinate1_pre(i,2)]);
        I_s_m(I_s_m < 0.5) = 0;
        I_s_m = logical(I_s_m);
        I_mix_pre = I_mix_pre | I_s_m;
    end
end

for i =1 : length(coordinate2_pre(:,1))
    if  coordinate2_pre(i,2) <= -numel(heiban(:,1)) + ly_single + 2 * abs(dy2)
        I_s_m_r = imtranslate(I_s_m_r0,[coordinate2_pre(i,1),coordinate2_pre(i,2)]);
        I_s_m_r(I_s_m_r < 0.5) = 0;
        I_s_m_r = logical(I_s_m_r);
        I_mix_pre = I_mix_pre | I_s_m_r;
    end
end
[~,~,nt_pre,nb_pre] = margin(I_mix_pre);
I_mix_pre_crop = imcrop(I_mix_pre,[0,0,numel(I_mix_pre(1,:)),numel(I_mix_pre(:,1)) - nb_pre]);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%首先，对I1进行旋转判断
LXY = zeros(180,3);
I1 = double(I1);

for angle = 1:180
    I1_ro = imrotate(I1,angle,'bicubic');
    I1_ro(I1_ro < 0.5) = 0;
    I1_ro(I1_ro ~= 0) = 1;
    lx_t = l_x(I1_ro);
    ly_t = l_y(I1_ro);
    LXY(angle,1) = lx_t;
    LXY(angle,2) = ly_t;
    LXY(angle,3) = angle;
end
LXY_y_ascend = sortrows(LXY,2);
LXY_y_min = LXY_y_ascend(1,2);

if LXY_y_min <= nt_pre
    LXY_y_suitable = LXY;
    index = (LXY_y_suitable(:,2) > nt_pre);
    LXY_y_suitable(index,:) = [];
    LXY_y_suitable = sortrows(LXY_y_suitable,2);
    LXY_y_suitable = LXY_y_suitable(end,:);
end

if exist('LXY_y_suitable','var')%如果该值存在,说明上面的空隙可以容下I1
    %先求dx1_top,dy1_top
    [dx1_top,dy1_top] = I1_adjust(I2_t_max_10,I1_hr_10,I2_t_max,I1_hr,imrf_special,imrf,LXY_y_suitable(3));
    
    %下面开始求dx2_top
    I1_hr = double(logical(I1_hr));
    I1_hr_top_ro = imrotate(I1_hr,LXY_y_suitable(3),'bicubic');
    I1_hr_top_ro(I1_hr_top_ro < 0.5) = 0;
    I1_hr_top_ro = logical(I1_hr_top_ro);
    I1_hr_top_ro = trans2lefttop(I1_hr_top_ro);
    I1_hr_top_ro = imcrop(I1_hr_top_ro,[0,0,l_x(I1_hr_top_ro),l_y(I1_hr_top_ro)]);
    I_base = zeros(round(1.2 * l_y(I1_hr_top_ro)), round(2.2 * l_x(I1_hr_top_ro)));
    I1_hr_top_ro = paste_img(I_base,I1_hr_top_ro,1,1);
    I1_hr_top_ro1 = imtranslate(I1_hr_top_ro,[l_x(I1_hr_top_ro),0]);
    I1_hr_top_ro_com = I1_hr_top_ro | I1_hr_top_ro1;
    N_single = sum(I1_hr_top_ro(:));
    j = 0;
    while (sum(I1_hr_top_ro_com(:)) / N_single) == 2
        j = j - 1;
        I1_hr_top_ro_com = I1_hr_top_ro | imtranslate(I1_hr_top_ro1,[j,0]);
    end
    j = j + 1;
    x_top_jiange = (l_x(I1_hr_top_ro) + j) / imrf;
    
    %     I2_t_max_t_10_noexpand = trans2mid(I2_t_max_t_10_noexpand,1);
    I2_t_max_t_10_noexpand_ro_t = imrotate(I2_t_max_t_10_noexpand,LXY_y_suitable(3),'bicubic','crop');
    I2_t_max_t_10_noexpand_ro_t(I2_t_max_t_10_noexpand_ro_t < 0.5) = 0;
    I2_t_max_t_10_noexpand_ro_t(I2_t_max_t_10_noexpand_ro_t ~= 0) = 1;
    
    [nl,~,~,nb] = margin(I2_t_max_t_10_noexpand_ro_t);
    middle_coordinate_x_ro_t = (nl + l_x(I2_t_max_t_10_noexpand_ro_t) / 2) / (imrf * imrf_special);
    middle_coordinate_y_ro_t = (nb + l_y(I2_t_max_t_10_noexpand_ro_t) / 2) / (imrf * imrf_special);
    
    middle_dx_t = middle_coordinate_x - middle_coordinate_x_ro_t;
    middle_dy_t = middle_coordinate_y - middle_coordinate_y_ro_t;
    
    I_base = zeros(numel(I_mix_pre_crop(:,1)),numel(I_mix_pre_crop(1,:)));
    
    I_try = double(I2_t_temp0);
    I_try = imrotate(I_try,LXY_y_suitable(3),'bicubic');
    I_try(I_try < 0.5) = 0;
    I_try(I_try ~= 0) = 1;
    
    lx_I_try = l_x(I_try);
    ly_I_try = l_y(I_try);
    
    x_first = lx_I_try / 2;
    y_first = numel(I_mix_pre(:,1)) - ly_I_try / 2;
    
    if dx1_top < 0
        x_first = x_first - dx1_top;
    end
    
    if dy1_top < 0
        y_first = y_first + dy1_top;
    end
    
    x_first_r = x_first + dx1_top;
    y_first_r = y_first - dy1_top;
    
    xmm = x_first;
    k = 0;
    coordinate_addi1 = zeros(1,3);%完全为了消除警告
    coordinate_addi2 = zeros(1,3);%完全为了消除警告
    while xmm <= numel(I_base(1,:)) - lx_I_try / 2
        k = k + 1;
        coordinate_addi1(k,1) = xmm + middle_dx_t;
        coordinate_addi1(k,2) = y_first + middle_dy_t;
        coordinate_addi1(k,3) = LXY_y_suitable(3);
        
        coordinate_addi1(k,4) = round(coordinate_addi1(k,1) / (lx_single / 2));
        coordinate_addi1(k,5) = coordinate_addi1(k,2);
        xmm = xmm + x_top_jiange;
    end
    
    xmm_r = x_first_r;
    k = 0;
    while xmm_r <= numel(I_base(1,:)) - lx_I_try / 2
        k = k + 1;
        coordinate_addi2(k,1) = xmm_r - middle_dx_t;
        coordinate_addi2(k,2) = y_first_r - middle_dy_t;
        coordinate_addi2(k,3) = LXY_y_suitable(3) + 180;
        
        coordinate_addi2(k,4) = round(coordinate_addi2(k,1) / (lx_single / 2));
        coordinate_addi2(k,5) = coordinate_addi2(k,2);
        xmm_r = xmm_r + x_top_jiange;
    end
    coordinate_addi_top = [coordinate_addi1; coordinate_addi2];
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%如果I1容不下，则对I2_t_temp0进行旋转判断
if exist('coordinate_addi_top','var')==0%如果该值不存在，说明上面的空间未能容下I1
    LXY= zeros(180,3);
    I2_t_temp0_double = double(I2_t_temp0);
    
    %     I2_t_temp0_double=double(I2_t_temp0);
    for angle = 1:180
        I2_t_temp0_ro=imrotate(I2_t_temp0_double,angle,'bicubic');
        I2_t_temp0_ro(I2_t_temp0_ro < 0.5)=0;
        I2_t_temp0_ro(I2_t_temp0_ro ~= 0)=1;
        lx_t=l_x(I2_t_temp0_ro);
        ly_t=l_y(I2_t_temp0_ro);
        LXY(angle,1) = lx_t;
        LXY(angle,2) = ly_t;
        LXY(angle,3) = angle;
    end
    LXY_y_ascend = sortrows(LXY,2);
    LXY_y_min = LXY_y_ascend(1,2);
    
    if LXY_y_min <= nt_pre
        LXY_y_suitable = LXY;
        index = (LXY_y_suitable(:,2) > nt_pre);
        LXY_y_suitable(index,:) = [];
        LXY_y_suitable = sortrows(LXY_y_suitable,2);
        LXY_y_suitable = LXY_y_suitable(end,:);
    else
        LXY_y_suitable1 = LXY_y_ascend(1,:);
        LXY_y_suitable2 = LXY_y_ascend(1,:);
        LXY_y_suitable2(3) = LXY_y_suitable2(3) + 180;
    end
    
    if LXY_y_min <= nt_pre %说明上面的空隙完全可以将零件(旋转后)放进去
        %         I2_t_max_t_10_noexpand = trans2mid(I2_t_max_t_10_noexpand,1);
        I2_t_max_t_10_noexpand_ro_t = imrotate(I2_t_max_t_10_noexpand,LXY_y_suitable(3),'bicubic','crop');
        I2_t_max_t_10_noexpand_ro_t(I2_t_max_t_10_noexpand_ro_t < 0.5) = 0;
        I2_t_max_t_10_noexpand_ro_t(I2_t_max_t_10_noexpand_ro_t ~= 0) = 1;
        
        [nl,~,~,nb] = margin(I2_t_max_t_10_noexpand_ro_t);
        middle_coordinate_x_ro_t = (nl + l_x(I2_t_max_t_10_noexpand_ro_t) / 2) / (imrf * imrf_special);
        middle_coordinate_y_ro_t = (nb + l_y(I2_t_max_t_10_noexpand_ro_t) / 2) / (imrf * imrf_special);
        
        middle_dx_t = middle_coordinate_x - middle_coordinate_x_ro_t;
        middle_dy_t = middle_coordinate_y - middle_coordinate_y_ro_t;
        
        %首先，将I2_t_temp0旋转合适的角度，然后移动到左上角
        I_try = double(I2_t_temp0);
        I_try = imrotate(I_try,LXY_y_suitable(3),'bicubic');
        I_try(I_try < 0.5) = 0;
        I_try(I_try ~= 0) = 1;
        I_try = trans2lefttop(I_try);
        I_try = imcrop(I_try,[0,0,l_x(I_try),l_y(I_try)]);
        I_base = zeros(numel(I_mix_pre_crop(:,1)),numel(I_mix_pre_crop(1,:)));
        
        lx_I_try = l_x(I_try);
        ly_I_try = l_y(I_try);
        I_try = paste_img(I_base,I_try,1,1);
        N_I_try = sum(I_try(:));
        
        for i=l_x(I_try) : (-1) : 1
            I_try_t = imtranslate(I_try, [i,0]);
            I_try_double = I_try | I_try_t;
            if sum(I_try_double(:)) ~= 2 * N_I_try
                x_jiange = i+1;
                break
            end
        end
        n_addi = floor(numel(I_mix_pre_crop(1,:)) / x_jiange) - 1;
        I_addi_com = I_try;
        coordinate_addi_top = zeros(n_addi,3);
        for i = 0 : n_addi
            I_try_t = imtranslate(I_try,[i * x_jiange,0]);
            I_addi_com = I_addi_com | I_try_t;
            %             imshow(I_addi_com)
            coordinate_addi_top(i+1,1) = i * x_jiange + lx_I_try/2;
            coordinate_addi_top(i+1,2) = numel(I_mix_pre(:,1)) - ly_I_try/2;
            coordinate_addi_top(i+1,3) = LXY_y_suitable(3);
            
            coordinate_addi_top(i+1,1) = coordinate_addi_top(i+1,1) + middle_dx_t;
            coordinate_addi_top(i+1,2) = coordinate_addi_top(i+1,2) + middle_dy_t;
            coordinate_addi_top(i+1,4) = round(coordinate_addi_top(i+1,1) / (lx_single / 2));
            coordinate_addi_top(i+1,5) = coordinate_addi_top(i+1,2);
        end
    end
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if exist('coordinate_addi_top','var')==0%如果该值不存在，说明上面的空间未能容下I2_t_temp0
    %如果I2_t_temp0也容不下，则强制性的对I2_t_temp0进行旋转，并取l_y最小值的角度进行向右移动，看看能不能容下，如果容不下，则旋转180°再进行右移判断
    %  LXY_y_suitable1和LXY_y_suitable2 在前面已经求得
    if nt_pre >= LXY_y_min / 4
        
        %         I2_t_max_t_10_noexpand = trans2mid(I2_t_max_t_10_noexpand,1);
        I2_t_max_t_10_noexpand_ro_t = imrotate(I2_t_max_t_10_noexpand,LXY_y_suitable1(3),'bicubic','crop');
        I2_t_max_t_10_noexpand_ro_t(I2_t_max_t_10_noexpand_ro_t < 0.5) = 0;
        I2_t_max_t_10_noexpand_ro_t(I2_t_max_t_10_noexpand_ro_t ~= 0) = 1;
        
        [nl,~,~,nb] =margin(I2_t_max_t_10_noexpand_ro_t);
        middle_coordinate_x_ro_t = (nl + l_x(I2_t_max_t_10_noexpand_ro_t)/2)/(imrf * imrf_special);
        middle_coordinate_y_ro_t = (nb + l_y(I2_t_max_t_10_noexpand_ro_t)/2)/(imrf * imrf_special);
        
        middle_dx_t = middle_coordinate_x - middle_coordinate_x_ro_t;
        middle_dy_t = middle_coordinate_y - middle_coordinate_y_ro_t;
        
        coordinate_addi1 = zeros(1,3);%完全为了消除警告
        coordinate_addi2 = zeros(1,3);%完全为了消除警告
        p1 = 0;
        p2 = 0;
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %先尝试LXY_y_suitable1
        %首先，将I2_t_temp0旋转合适的角度，然后移动到左上角
        I_try = double(I2_t_temp0);
        I_try = imrotate(I_try,LXY_y_suitable1(3),'bicubic');
        I_try(I_try < 0.5) = 0;
        I_try(I_try ~= 0) = 1;
        I_try = trans2lefttop(I_try);
        I_try = imcrop(I_try,[0,0,l_x(I_try),l_y(I_try)]);
        I_base = zeros(numel(I_mix_pre_crop(:,1)),numel(I_mix_pre_crop(1,:)));
        I_try = paste_img(I_base,I_try,1,1); I_try = logical(I_try);
        %     N_I_mix_pre_crop = sum(I_mix_pre_crop(:));
        
        lx_I_try = l_x(I_try);
        ly_I_try = l_y(I_try);
        
        N_I_try = sum(I_try(:));
        
        I_try_t = imtranslate(I_try, [lx_I_try,0]);
        for k = 0 : -1 : -lx_I_try
            I_try_t = imtranslate(I_try_t, [k,0]);
            I_com = I_try | I_try_t;
            if sum(I_com(:)) ~= 2 * N_I_try
                x_top_jiange_min = lx_I_try + k + 1;
                break
            end
        end
        
        x_first = lx_I_try / 2;
        y_first = numel(I_mix_pre(:,1)) - ly_I_try / 2;
        %     I_addi_com = I_try;
        x_trans = 0;
        tiaoyue = x_top_jiange_min;
        y_t_try = round(min(max(3 * dx3, 3 * (lx_I_try)),(numel(I_try(1,:)) - lx_I_try)));
        N_I_mix_pre_crop = sum(I_mix_pre_crop(:));
        
        while x_trans <= y_t_try
            I_try_t = imtranslate(I_try,[x_trans, 0]);
            I_com = I_mix_pre_crop | I_try_t;
            %             imshow(I_com)
            if sum(I_com(:)) - N_I_mix_pre_crop == N_I_try%说明图形并未发生重叠
                p1 = p1 + 1;%多了一个合适的位置
                
                coordinate_addi1(p1,1) = x_first + x_trans + middle_dx_t;
                coordinate_addi1(p1,2) = y_first + middle_dy_t;
                coordinate_addi1(p1,3) = LXY_y_suitable1(3);
                
                coordinate_addi1(p1,4) = round(coordinate_addi1(p1,1) / (lx_single / 2));
                coordinate_addi1(p1,5) = coordinate_addi1(p1,2);
                
                x_trans = x_trans + tiaoyue;
                %                 imshow(I_com)
            else
                x_trans = x_trans + 1;
            end
        end
        
        if p1 >= 2
            while x_trans <= (numel(I_try(1,:)) - lx_I_try)
                I_try_t = imtranslate(I_try,[x_trans, 0]);
                I_com = I_mix_pre_crop | I_try_t;
                %                 imshow(I_com)
                if sum(I_com(:)) - N_I_mix_pre_crop == N_I_try%说明图形并未发生重叠
                    p1 = p1 + 1;%多了一个合适的位置
                    
                    coordinate_addi1(p1,1) = x_first + x_trans + middle_dx_t;
                    coordinate_addi1(p1,2) = y_first + middle_dy_t;
                    coordinate_addi1(p1,3) = LXY_y_suitable1(3);
                    
                    coordinate_addi1(p1,4) = round(coordinate_addi1(p1,1) / (lx_single / 2));
                    coordinate_addi1(p1,5) = coordinate_addi1(p1,2);
                    
                    x_trans = x_trans + tiaoyue;
                    %                     imshow(I_com)
                else
                    x_trans = x_trans + 1;
                end
            end
        end
        
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %再尝试LXY_y_suitable2
        I_try = double(I2_t_temp0);
        I_try = imrotate(I_try,LXY_y_suitable2(3),'bicubic');
        I_try(I_try < 0.5) = 0;
        I_try(I_try ~= 0) = 1;
        I_try = trans2lefttop(I_try);
        I_try = imcrop(I_try,[0,0,l_x(I_try),l_y(I_try)]);
        I_base = zeros(numel(I_mix_pre_crop(:,1)),numel(I_mix_pre_crop(1,:)));
        I_try = paste_img(I_base,I_try,1,1);I_try = logical(I_try);
        
        x_trans = 0;
        tiaoyue = x_top_jiange_min;
        
        while x_trans <= y_t_try
            I_try_t=imtranslate(I_try,[x_trans, 0]);
            I_com = I_mix_pre_crop | I_try_t;
            %             imshow(I_com)
            if sum(I_com(:)) - N_I_mix_pre_crop == N_I_try%说明图形并未发生重叠
                p2 = p2 + 1;%多了一个合适的位置
                
                coordinate_addi2(p2,1) = x_first + x_trans - middle_dx_t;
                coordinate_addi2(p2,2) = y_first - middle_dy_t;
                coordinate_addi2(p2,3) = LXY_y_suitable2(3);
                
                coordinate_addi2(p2,4) = round(coordinate_addi2(p2,1) / (lx_single / 2));
                coordinate_addi2(p2,5) = coordinate_addi2(p2,2);
                
                x_trans = x_trans + tiaoyue;
                %                 imshow(I_com)
            else
                x_trans = x_trans + 1;
            end
        end
        
        if p2 >= 2
            while x_trans <= (numel(I_try(1,:)) - lx_I_try)
                I_try_t=imtranslate(I_try,[x_trans, 0]);
                I_com = I_mix_pre_crop | I_try_t;
                %                 imshow(I_com)
                if sum(I_com(:)) - N_I_mix_pre_crop == N_I_try%说明图形并未发生重叠
                    p2 = p2 + 1;%多了一个合适的位置
                    
                    coordinate_addi2(p2,1) = x_first + x_trans - middle_dx_t;
                    coordinate_addi2(p2,2) = y_first - middle_dy_t;
                    coordinate_addi2(p2,3) = LXY_y_suitable2(3);
                    
                    coordinate_addi2(p2,4) = round(coordinate_addi2(p2,1) / (lx_single / 2));
                    coordinate_addi2(p2,5) = coordinate_addi2(p2,2);
                    
                    x_trans = x_trans + tiaoyue;
                    %                     imshow(I_com)
                else
                    x_trans = x_trans + 1;
                end
            end
        end
        
        if p1 >= p2 && p1~=0
            coordinate_addi_top = coordinate_addi1;
        elseif p1 < p2
            coordinate_addi_top = coordinate_addi2;
        elseif p1 == p2 && p1 == 0
            coordinate_addi_top = [];
        end
        
    end
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if exist('coordinate_addi_top','var')==0%如果该值不存在
    coordinate_addi_top = [];
end
end

