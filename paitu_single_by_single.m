function [I_mix,coordinate,board_is_full,N_left] = paitu_single_by_single(I_mix,I2_t_max_or,N_left,sita_xz)
%   此处显示详细说明
tic
isrectangle = IsRectangle(I2_t_max_or);
Rectangle_used = 0;
I2_t_max_or_filled = imfill(I2_t_max_or, 'holes');
N_threshold = sum(I2_t_max_or_filled(:));

X_max = numel(I_mix(1,:));
lxymax = max(l_x(I2_t_max_or), l_y(I2_t_max_or));
% I2_t_max_or_filled = imfill(I2_t_max_or, 'holes');
% if sum(I_mix(:)) == 0%空板子，还没有排任何的零件
%     % X_initial = min(round(1.5*lxymax),X_max);
%     X_initial = min(round(0.5 * lxymax),X_max);
% else
%     [~,nr,~,~] = margin(I_mix);
%     X_initial = max(numel(I_mix(1,:)) - nr - lxymax, lxymax);
% end

[~,nr,~,~] = margin(I_mix);
% X_initial = max(floor(l_x(I_mix)/(1.5 * lxymax)) * lxymax, 2 * lxymax);
if nr == 0
    X_initial =  X_max;
else
    X_initial = max(numel(I_mix(1,:)) - nr - lxymax, round(1.5 * lxymax));
end
X_initial = min(X_initial,X_max);

if l_x(I_mix) > numel(I_mix(1,:)) * 1 / 2
    X_initial = numel(I_mix(1,:));
end
X_extend_to = X_initial ;

board_is_full = 0;
coordinate = [];
coordinate_N = 0;
f_times = 1;%寻找次数
I_mix_re =~ I_mix;
I=double(I_mix_re);

cc = bwconncomp(I);%其找的是白洞，不是黑洞
% 获取空心区域的像素索引
holePixelIdxList = cc.PixelIdxList;
% 获取空心区域的数量
numHoles = numel(holePixelIdxList);
for i = 1:numHoles
    
    % 获取当前空心区域的像素索引
    pixelIdxList = holePixelIdxList{i};
    % 根据像素索引创建空心区域的掩模
    holeMask = true(size(I));
    
    if length(pixelIdxList) < N_threshold
        holeMask(pixelIdxList) = false;
    end
    
    I = I & holeMask;
end

% [Ir, Ic] = size(I);
% figure(1),imshow(I)
T = trans2lefttop(I2_t_max_or);
T = imcrop(T,[0, 0, l_x(T), l_y(T)]);
I_base = zeros(l_y(T)+2,l_x(T)+2);
% I_base = zeros(l_y(T)+0,l_x(T)+1);
T = paste_img(I_base,T,1,1);
T = trans2mid(T,0);
% T = trans2lefttop(T);

I_base = zeros(size(I_mix));

T_0 = imrotate(T,0,'bicubic');
T_90 = imrotate(T,90,'bicubic');
T_180 = imrotate(T,180,'bicubic');
T_270 = imrotate(T,270,'bicubic');


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% T_0 = trans2lefttop(T_0);
% T_90 = trans2lefttop(T_90);
% T_180 = trans2lefttop(T_180);
% T_270 = trans2lefttop(T_270);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


T_tem0 = paste_img(I_base,T_0,1,1);
T_tem90 = paste_img(I_base,T_90,1,1);
T_tem180 = paste_img(I_base,T_180,1,1);
T_tem270 = paste_img(I_base,T_270,1,1);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% T_tem0 = trans2lefttop(T_tem0);
% T_tem90 = trans2lefttop(T_tem90);
% T_tem180 = trans2lefttop(T_tem180);
% T_tem270 = trans2lefttop(T_tem270);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%






T_0 = imfill(T_0, 'holes');
T_90 = imfill(T_90, 'holes');
T_180 = imfill(T_180, 'holes');
T_270 = imfill(T_270, 'holes');

[margin_left_0,margin_right_0,margin_top_0,margin_bottom_0] = margin(T_tem0);
[margin_left_90,margin_right_90,margin_top_90,margin_bottom_90] = margin(T_tem90);
[margin_left_180,margin_right_180,margin_top_180,margin_bottom_180] = margin(T_tem180);
[margin_left_270,margin_right_270,margin_top_270,margin_bottom_270] = margin(T_tem270);

% [Tr0, Tc0] = size(T_0);
% [Tr90, Tc90] = size(T_90);
% [Tr180, Tc180] = size(T_180);
% [Tr270, Tc270] = size(T_270);
% figure(2),imshow(T)

R0 = normxcorr2(T_0,I);
columnsToDelete = [1 : (numel(T_0(1,:)) - margin_left_0 - 1), ((numel(T_0(1,:)) + margin_right_0+1)) : size(R0,2)];
rowsToDelete = [1:(numel(T_0(:,1)) - margin_top_0 - 1), ((numel(T_0(:,1)) + margin_bottom_0+1)):size(R0,1)];
R0(:,columnsToDelete) = [];
R0(rowsToDelete,:) = [];
%%%%%%%%%%%%%%%%
R0=R0(:,1:min(X_initial,numel(R0(1,:))));
%%%%%%%%%%%%%%%%
R0_sorted = R0(:);
R0_sorted(:,2) = 0;

R90 = normxcorr2(T_90,I);
columnsToDelete = [1:(numel(T_90(1,:)) - margin_left_90 - 1), ((numel(T_90(1,:)) + margin_right_90+1)):size(R90,2)];
rowsToDelete = [1:(numel(T_90(:,1)) - margin_top_90 - 1), ((numel(T_90(:,1)) + margin_bottom_90+1)):size(R90,1)];
R90(:,columnsToDelete) = [];
R90(rowsToDelete,:) = [];
%%%%%%%%%%%%%%%%
R90=R90(:,1:min(X_initial,numel(R90(1,:))));
%%%%%%%%%%%%%%%%
R90_sorted = R90(:);
R90_sorted(:,2)=90;

R180 = normxcorr2(T_180,I);
columnsToDelete = [1:(numel(T_180(1,:)) - margin_left_180 - 1), ((numel(T_180(1,:)) + margin_right_180+1)):size(R180,2)];
rowsToDelete = [1:(numel(T_180(:,1)) - margin_top_180 - 1), ((numel(T_180(:,1)) + margin_bottom_180+1)):size(R180,1)];
R180(:,columnsToDelete) = [];
R180(rowsToDelete,:) = [];
%%%%%%%%%%%%%%%%
R180=R180(:,1:min(X_initial,numel(R180(1,:))));
%%%%%%%%%%%%%%%%
R180_sorted = R180(:);
R180_sorted(:,2)=180;

R270 = normxcorr2(T_270,I);
columnsToDelete = [1:(numel(T_270(1,:)) - margin_left_270 - 1), ((numel(T_270(1,:)) + margin_right_270+1)):size(R270,2)];
rowsToDelete = [1:(numel(T_270(:,1)) - margin_top_270 - 1), ((numel(T_270(:,1)) + margin_bottom_270+1)):size(R270,1)];
R270(:,columnsToDelete) = [];
R270(rowsToDelete,:) = [];
%%%%%%%%%%%%%%%%
R270=R270(:,1:min(X_initial,numel(R270(1,:))));
%%%%%%%%%%%%%%%%
R270_sorted = R270(:);
R270_sorted(:,2)=270;

Rtotal = [R0_sorted;R90_sorted;R180_sorted;R270_sorted];
Rtotal_sorted = sortrows(Rtotal,1,'descend');
Rtotal_sorted= unique(Rtotal_sorted,'rows','stable');
l_cut = ceil(length(Rtotal_sorted)/1);
Rtotal_sorted = Rtotal_sorted(1:l_cut,:);


cumlate = 0;
a = randi([1, 4]);
a=4;
rc_notok0(1,1)=24400;
rc_notok0(1,2)=24400;
rc_notok90(1,1)=24400;
rc_notok90(1,2)=24400;
rc_notok180(1,1)=24400;
rc_notok180(1,2)=24400;
rc_notok270(1,1)=24400;
rc_notok270(1,2)=24400;
rc_num0 = 0;
rc_num90 = 0;
rc_num180 = 0;
rc_num270 = 0;
f_compare = 1;
% chushu = 3000 / lxymax;
chushu = randi([2500, 3500]) / lxymax;
%  X_extend_to =1220;

% if sum(I_mix(:)) == 0%板子上什么都没有
%     a_compare = floor(length(Rtotal) / 2);
% else
%     a_compare = floor(length(Rtotal_sorted) / 2);
% end

a_compare = floor(length(Rtotal_sorted) / 2);

while (f_times <= 10000) && (a < a_compare) && (N_left ~= 0) && ~((a > 4) && (sum(I_mix(:)) == 0))   %如果不满足条件，则认为图已经排满了，或者某个零件排完了
    
    %     if  f_times <= 50
    %         bujin = 1;
    %     elseif f_times > 50 && f_times <= 100
    %         bujin = 100;
    %     elseif f_times > 100 && f_times <= 1000
    %         bujin = 100;
    %     elseif f_times > 1000 && f_times <= 2000
    %         bujin = 100;
    %     elseif f_times > 2000 && f_times <= 3000
    %         bujin = 100;
    %     elseif f_times > 3000 && f_times <= 4000
    %         bujin = 100;
    %     else
    %         bujin = 500;
    %     end
    
        if  f_times <= 110
            bujin = 1;
        elseif f_times > 100 && f_times <= 1000
%             bujin=10;
            bujin = randi([5, 20]);
%             bujin = max(round(f_times / chushu),1);
%             bujin = min(bujin, 100);
        elseif f_times > 1000 && f_times <= 2000
            bujin = 25;
        elseif f_times > 2000 && f_times <= 3000
            bujin = 30;
        elseif f_times > 3000 && f_times <= 4000
            bujin = 40;
        else
            bujin = 50;
        end
    
%     if  f_times <= 100
%         bujin = 1;
%     elseif f_times > 100 && f_times <= 2000
%         bujin = 50;
% %         bujin = max(round(f_times / chushu),1);
% %         bujin = min(bujin, 100);
%     elseif f_times > 2000 && f_times <= 3000
%         bujin = 20;
%     elseif f_times > 3000 && f_times <= 4000
%         bujin = 30;
%     else
%         bujin = 50;
%     end
%     
% %     bujin = 1;
    
    if floor(f_times / 1000) == f_compare
        
        f_compare = f_compare + 1;
        pp = ceil(f_times / 2000) - 1;
        
        if pp > 0
            se = strel('rectangle',[3*pp 3*pp]);
            I_mix_re = ~(imdilate(I_mix,se));%是白色的区域膨胀
            I = double(I_mix_re);
            %         I = imdilate(I,se); %是白色的区域膨胀
            cc = bwconncomp(I);%其找的是白洞，不是黑洞
            % 获取空心区域的像素索引
            holePixelIdxList = cc.PixelIdxList;
            % 获取空心区域的数量
            numHoles = numel(holePixelIdxList);
            for i = 1:numHoles
                
                % 获取当前空心区域的像素索引
                pixelIdxList = holePixelIdxList{i};
                % 根据像素索引创建空心区域的掩模
                holeMask = true(size(I));
                
                if length(pixelIdxList) < (1 + (f_compare-1) / 10) * N_threshold
                    holeMask(pixelIdxList) = false;
                end
                
                I = I & holeMask;
            end
        end
        
        %         cumlate = 1;
        
        if isrectangle == 0%如果不是矩形
            
            if Rectangle_used == 0
                if f_times < 1000
                    T_0_t = T_0;
                    T_90_t = T_90;
                    T_180_t = T_180;
                    T_270_t = T_270;
                else
                    T_0_t = shape_to_rec(T_0);
                    T_90_t = shape_to_rec(T_90);
                    T_180_t = shape_to_rec(T_180);
                    T_270_t = shape_to_rec(T_270);
                    a = 1;
                    f_times =1;
                    bujin = 1;
                    Rectangle_used = 1;
                    f_compare = 1;
                end
            end
            
        else
            
            T_0_t = T_0;
            T_90_t = T_90;
            T_180_t = T_180;
            T_270_t = T_270;
            
        end
        
        R0 = normxcorr2(T_0_t,I);
        columnsToDelete = [1 : (numel(T_0_t(1,:)) - margin_left_0-1), ((numel(T_0_t(1,:)) + margin_right_0+1)) : size(R0,2)];
        rowsToDelete = [1 : (numel(T_0_t(:,1)) - margin_top_0-1), ((numel(T_0_t(:,1)) + margin_bottom_0+1)):size(R0,1)];
        R0(:,columnsToDelete) = [];
        R0(rowsToDelete,:) = [];
        %%%%%%%%%%%%%%%%min(X_initial,numel(R0(1,:))));
        R0=R0(:,1:min(X_extend_to,numel(R0(1,:))));
        %%%%%%%%%%%%%%%%
        R0_sorted = R0(:);
        R0_sorted(:,2)=0;
        
        R90 = normxcorr2(T_90_t,I);
        columnsToDelete = [1:(numel(T_90_t(1,:)) - margin_left_90-1), ((numel(T_90_t(1,:)) + margin_right_90+1)):size(R90,2)];
        rowsToDelete = [1:(numel(T_90_t(:,1)) - margin_top_90-1), ((numel(T_90_t(:,1)) + margin_bottom_90+1)):size(R90,1)];
        R90(:,columnsToDelete) = [];
        R90(rowsToDelete,:) = [];
        %%%%%%%%%%%%%%%%
        R90=R90(:,1:min(X_extend_to,numel(R90(1,:))));
        %%%%%%%%%%%%%%%%
        R90_sorted = R90(:);
        R90_sorted(:,2)=90;
        
        R180 = normxcorr2(T_180_t,I);
        columnsToDelete = [1:(numel(T_180_t(1,:)) - margin_left_180-1), ((numel(T_180_t(1,:)) + margin_right_180+1)):size(R180,2)];
        rowsToDelete = [1:(numel(T_180_t(:,1)) - margin_top_180-1), ((numel(T_180_t(:,1)) + margin_bottom_180+1)):size(R180,1)];
        R180(:,columnsToDelete) = [];
        R180(rowsToDelete,:) = [];
        %%%%%%%%%%%%%%%%
        R180=R180(:,1:min(X_extend_to,numel(R180(1,:))));
        %%%%%%%%%%%%%%%%
        R180_sorted = R180(:);
        R180_sorted(:,2)=180;
        
        R270 = normxcorr2(T_270_t,I);
        columnsToDelete = [1:(numel(T_270_t(1,:)) - margin_left_270-1), ((numel(T_270_t(1,:)) + margin_right_270+1)):size(R270,2)];
        rowsToDelete = [1:(numel(T_270_t(:,1)) - margin_top_270-1), ((numel(T_270_t(:,1)) + margin_bottom_270+1)):size(R270,1)];
        R270(:,columnsToDelete) = [];
        R270(rowsToDelete,:) = [];
        %%%%%%%%%%%%%%%%
        R270=R270(:,1:min(X_extend_to,numel(R270(1,:))));
        %%%%%%%%%%%%%%%%
        R270_sorted = R270(:);
        R270_sorted(:,2)=270;
        
        Rtotal = [R0_sorted;R90_sorted;R180_sorted;R270_sorted];
        Rtotal_sorted = sortrows(Rtotal,1,'descend');
        Rtotal_sorted= unique(Rtotal_sorted,'rows','stable');
        l_cut = ceil(length(Rtotal_sorted)/1);
        Rtotal_sorted = Rtotal_sorted(1:l_cut,:);
        a_compare = floor(length(Rtotal_sorted) / 2);
    end
    
    if X_extend_to < X_max
        if (f_times > 5000) || (a > (length(Rtotal) / 30)) %如果满足条件，则认为图需要向右延伸
            X_extend_to = X_extend_to + round(2 * lxymax);
            X_extend_to = min(X_extend_to, X_max);
            
            if min(X_extend_to, X_max) == X_max %说明板子已经最大了
                cumlate = 0;
            else %板子没有达到最大
                cumlate = 1;
            end
        end
    end
    
    if cumlate >= 1
        Rectangle_used = 0;
        f_times = 0;%寻找次数
        cumlate = 0;
        a = 1;
        
        I_mix_re = ~I_mix;
        I = double(I_mix_re);
        
        cc = bwconncomp(I);%其找的是白洞，不是黑洞
        % 获取空心区域的像素索引
        holePixelIdxList = cc.PixelIdxList;
        % 获取空心区域的数量
        numHoles = numel(holePixelIdxList);
        for i = 1:numHoles
            
            % 获取当前空心区域的像素索引
            pixelIdxList = holePixelIdxList{i};
            % 根据像素索引创建空心区域的掩模
            holeMask = true(size(I));
            
            if length(pixelIdxList) < N_threshold
                holeMask(pixelIdxList) = false;
            end
            
            I = I & holeMask;
        end
        
        R0 = normxcorr2(T_0,I);
        columnsToDelete = [1 : (numel(T_0(1,:)) - margin_left_0-1), ((numel(T_0(1,:)) + margin_right_0+1)) : size(R0,2)];
        rowsToDelete = [1 : (numel(T_0(:,1)) - margin_top_0-1), ((numel(T_0(:,1)) + margin_bottom_0+1)):size(R0,1)];
        R0(:,columnsToDelete) = [];
        R0(rowsToDelete,:) = [];
        %%%%%%%%%%%%%%%%min(X_initial,numel(R0(1,:))));
        R0=R0(:,1:min(X_extend_to,numel(R0(1,:))));
        %%%%%%%%%%%%%%%%
        R0_sorted = R0(:);
        R0_sorted(:,2) = 0;
        
        R90 = normxcorr2(T_90,I);
        columnsToDelete = [1:(numel(T_90(1,:)) - margin_left_90-1), ((numel(T_90(1,:)) + margin_right_90+1)):size(R90,2)];
        rowsToDelete = [1:(numel(T_90(:,1)) - margin_top_90-1), ((numel(T_90(:,1)) + margin_bottom_90+1)):size(R90,1)];
        R90(:,columnsToDelete) = [];
        R90(rowsToDelete,:) = [];
        %%%%%%%%%%%%%%%%
        R90=R90(:,1:min(X_extend_to,numel(R90(1,:))));
        %%%%%%%%%%%%%%%%
        R90_sorted = R90(:);
        R90_sorted(:,2)=90;
        
        R180 = normxcorr2(T_180,I);
        columnsToDelete = [1:(numel(T_180(1,:)) - margin_left_180-1), ((numel(T_180(1,:)) + margin_right_180+1)):size(R180,2)];
        rowsToDelete = [1:(numel(T_180(:,1)) - margin_top_180-1), ((numel(T_180(:,1)) + margin_bottom_180+1)):size(R180,1)];
        R180(:,columnsToDelete) = [];
        R180(rowsToDelete,:) = [];
        %%%%%%%%%%%%%%%%
        R180=R180(:,1:min(X_extend_to,numel(R180(1,:))));
        %%%%%%%%%%%%%%%%
        R180_sorted = R180(:);
        R180_sorted(:,2) = 180;
        
        R270 = normxcorr2(T_270,I);
        columnsToDelete = [1:(numel(T_270(1,:)) - margin_left_270-1), ((numel(T_270(1,:)) + margin_right_270+1)):size(R270,2)];
        rowsToDelete = [1:(numel(T_270(:,1)) - margin_top_270-1), ((numel(T_270(:,1)) + margin_bottom_270+1)):size(R270,1)];
        R270(:,columnsToDelete) = [];
        R270(rowsToDelete,:) = [];
        %%%%%%%%%%%%%%%%
        R270=R270(:,1:min(X_extend_to,numel(R270(1,:))));
        %%%%%%%%%%%%%%%%
        R270_sorted = R270(:);
        R270_sorted(:,2)=270;
        
        Rtotal = [R0_sorted;R90_sorted;R180_sorted;R270_sorted];
        Rtotal_sorted = sortrows(Rtotal,1,'descend');
        Rtotal_sorted= unique(Rtotal_sorted,'rows','stable');
        l_cut = ceil(length(Rtotal_sorted)/1);
        Rtotal_sorted = Rtotal_sorted(1:l_cut,:);
        a_compare = floor(length(Rtotal_sorted) / 2);
    end
    
    if a > length(Rtotal_sorted)
        imshow(I2_t_max_or)
        pause(1000)
    end
    
    if Rtotal_sorted(a,2) == 0
        [r, c, ~] = find(R0 == Rtotal_sorted(a,1)); 
        isMatch = any(rc_notok0(:, 1) == r(1) & rc_notok0(:, 2) == c(1));
    elseif Rtotal_sorted(a,2) == 90
        [r, c, ~] = find(R90 == Rtotal_sorted(a,1)); 
        isMatch = any(rc_notok90(:, 1) == r(1) & rc_notok90(:, 2) == c(1));
    elseif Rtotal_sorted(a,2) == 180
        [r, c, ~] = find(R180 == Rtotal_sorted(a,1)); 
        isMatch = any(rc_notok180(:, 1) == r(1) & rc_notok180(:, 2) == c(1));
    elseif Rtotal_sorted(a,2) == 270
        [r, c, ~] = find(R270 == Rtotal_sorted(a,1)); 
        isMatch = any(rc_notok270(:, 1) == r(1) & rc_notok270(:, 2) == c(1));
    end
    
    while isMatch%如果isMatch值为1，说明其遇到了已经判断过的不合适的值
        a = a + 1;
        if Rtotal_sorted(a,2) == 0
            [r, c, ~] = find(R0 == Rtotal_sorted(a,1)); 
            isMatch = any(rc_notok0(:, 1) == r(1) & rc_notok0(:, 2) == c(1));
        elseif Rtotal_sorted(a,2) == 90
            [r, c, ~] = find(R90 == Rtotal_sorted(a,1)); 
            isMatch = any(rc_notok90(:, 1) == r(1) & rc_notok90(:, 2) == c(1));
        elseif Rtotal_sorted(a,2) == 180
            [r, c, ~] = find(R180 == Rtotal_sorted(a,1)); 
            isMatch = any(rc_notok180(:, 1) == r(1) & rc_notok180(:, 2) == c(1));
        elseif Rtotal_sorted(a,2) == 270
            [r, c, ~] = find(R270 == Rtotal_sorted(a,1)); 
            isMatch = any(rc_notok270(:, 1) == r(1) & rc_notok270(:, 2) == c(1));
        end
    end
    
    if Rtotal_sorted(a,2) == 0
        T_tem = T_tem0;
        T_moved = imtranslate(T_tem,[c(1)-margin_left_0-1,r(1)-margin_top_0-1]);
    elseif Rtotal_sorted(a,2) == 90
        T_tem = T_tem90;
        T_moved = imtranslate(T_tem,[c(1)-margin_left_90-1,r(1)-margin_top_90-1]);
    elseif Rtotal_sorted(a,2) == 180
        T_tem = T_tem180;
        T_moved = imtranslate(T_tem,[c(1)-margin_left_180-1,r(1)-margin_top_180-1]);
    elseif Rtotal_sorted(a,2) == 270
        T_tem = T_tem270;
        T_moved = imtranslate(T_tem,[c(1)-margin_left_270-1,r(1)-margin_top_270-1]);
    end
    
    I_com = I_mix | T_moved;
    %          imshow(I_com)
    
    if abs(sum(I_com(:)) - sum(I_mix(:)) - sum(T_tem(:))) == 0 %说明零件的位置可行，并未发生重叠
        %这里需要优化，因为bujin可能较大，需要退回
        if a - bujin + 1 > 0
            b_initial = a - bujin + 1;
        else
            b_initial = 1;
        end
        
        for b = b_initial : a
            
            if Rtotal_sorted(b,2) == 0
                [r, c, ~] = find(R0 == Rtotal_sorted(b,1)); 
                isMatch = any(rc_notok0(:, 1) == r(1) & rc_notok0(:, 2) == c(1));
            elseif Rtotal_sorted(b,2) == 90
                [r, c, ~] = find(R90 == Rtotal_sorted(b,1)); 
                isMatch = any(rc_notok90(:, 1) == r(1) & rc_notok90(:, 2) == c(1));
            elseif Rtotal_sorted(b,2) == 180
                [r, c, ~] = find(R180 == Rtotal_sorted(b,1)); 
                isMatch = any(rc_notok180(:, 1) == r(1) & rc_notok180(:, 2) == c(1));
            elseif Rtotal_sorted(b,2) == 270
                [r, c, ~] = find(R270 == Rtotal_sorted(b,1)); 
                isMatch = any(rc_notok270(:, 1) == r(1) & rc_notok270(:, 2) == c(1));
            end
            
            if isMatch
                disp(['a返回值：',num2str(b)]);
                continue;
            end
            
            if Rtotal_sorted(b,2) == 0
                T_tem = T_tem0;
                T_moved = imtranslate(T_tem,[c(1)-margin_left_0-1,r(1)-margin_top_0-1]);
            elseif Rtotal_sorted(b,2) == 90
                T_tem = T_tem90;
                T_moved = imtranslate(T_tem,[c(1)-margin_left_90-1,r(1)-margin_top_90-1]);
            elseif Rtotal_sorted(b,2) == 180
                T_tem = T_tem180;
                T_moved = imtranslate(T_tem,[c(1)-margin_left_180-1,r(1)-margin_top_180-1]);
            elseif Rtotal_sorted(b,2) == 270
                T_tem = T_tem270;
                T_moved = imtranslate(T_tem,[c(1)-margin_left_270-1,r(1)-margin_top_270-1]);
            end
            
            I_com = I_mix | T_moved;
            
            if abs(sum(I_com(:)) - sum(I_mix(:)) - sum(T_tem(:))) == 0 %满足条件
                disp(['bujin：',num2str(bujin)]);
                disp(['f_times：',num2str(f_times)]);                
                
                
                cumlate = cumlate + 1;
                N_left = N_left - 1;
                coordinate_N = coordinate_N + 1;
                [x_left,~,~,y_bottom] = margin(T_moved);
                coordinate(coordinate_N,1) = x_left + l_x(T_moved) / 2; %#ok<AGROW>
                coordinate(coordinate_N,2) = y_bottom + l_y(T_moved) / 2; %#ok<AGROW>
                coordinate(coordinate_N,3) = Rtotal_sorted(b,2) + sita_xz; %#ok<AGROW>
                
                I_mix = I_mix | T_moved;
                I_mix_show = imresize(I_mix,0.5);
                imshow(I_mix_show,'Border','tight')
                set(gcf,'ToolBar','none','MenuBar','none','NumberTitle','off');
                %set(gcf, 'CloseRequestFcn', @(src,event) myCloseFcn(src,event,h));
                pause(0.0001)
                toc
                tic
                break
                
            else
                
                if Rtotal_sorted(b,2) == 0
                    rc_num0 = rc_num0 + 1;
                    rc_notok0(rc_num0,1) = c(1);
                    rc_notok0(rc_num0,2) = r(1);
                elseif Rtotal_sorted(b,2) == 90
                    rc_num90 = rc_num90 + 1;
                    rc_notok90(rc_num90,1) = c(1);
                    rc_notok90(rc_num90,2) = r(1);
                elseif Rtotal_sorted(b,2) == 180
                    rc_num180 = rc_num180 + 1;
                    rc_notok180(rc_num180,1) = c(1);
                    rc_notok180(rc_num180,2) = r(1);
                elseif Rtotal_sorted(b,2) == 270
                    rc_num270 = rc_num270 + 1;
                    rc_notok270(rc_num270,1) = c(1);
                    rc_notok270(rc_num270,2) = r(1);
                end
            end
        end
        
    else
        
        if Rtotal_sorted(a,2) == 0
            rc_num0 = rc_num0 + 1;
            rc_notok0(rc_num0,1) = c(1);
            rc_notok0(rc_num0,2) = r(1);
        elseif Rtotal_sorted(a,2) == 90
            rc_num90 = rc_num90 + 1;
            rc_notok90(rc_num90,1) = c(1);
            rc_notok90(rc_num90,2) = r(1);
        elseif Rtotal_sorted(a,2) == 180
            rc_num180 = rc_num180 + 1;
            rc_notok180(rc_num180,1) = c(1);
            rc_notok180(rc_num180,2) = r(1);
        elseif Rtotal_sorted(a,2) == 270
            rc_num270 = rc_num270 + 1;
            rc_notok270(rc_num270,1) = c(1);
            rc_notok270(rc_num270,2) = r(1);
        end
    end
    
    a = a + bujin;
    f_times = f_times + 1;%寻找次数
    
    %     if sum(I_mix(:)) == 0%板子上什么都没有
    %         a_compare = floor(length(Rtotal) / 2);
    %     else
    %         a_compare = floor(length(Rtotal_sorted) / 2);
    %     end
end

if ((f_times > 10000) || (a >= floor(length(Rtotal_sorted) / 2))) && (N_left ~= 0) && ~((a > 4) && (sum(I_mix(:)) == 0))%这种情况，说明板子满了
    board_is_full = 1;
end

N_left = N_left - length(coordinate);

end
% delete(gcf);delete(gcf),delete(gcf);
% I_mix = ~I_mix_re;