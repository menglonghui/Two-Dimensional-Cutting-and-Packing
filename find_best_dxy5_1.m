function [coordinate,N] = find_best_dxy5_1(lx_heiban,ly_heiban,X_interval1,Y_interval1,X_interval2,Y_interval2,lx1,ly1,lx2,ly2)
%FIND_BEST_DXY5 此处显示有关此函数的摘要
%find_best_dxy5函数对图像完全是横平竖直的排列
%   此处显示详细说明
x_along = 0;
y_along = 0;
% if isempty(heiban)
if lx_heiban == 0 || ly_heiban == 0
    coordinate = [];
    N = 0;
else
    %     lx_heiban = numel(heiban(1,:));
    %     ly_heiban = numel(heiban(:,1));
    
    %首先，图形不旋转的情况下，求出X_interval，Y_interval
    %     I_base = zeros(2 * l_y(I2_t_max_or), 2 * l_x(I2_t_max_or));
    %     I_t = paste_img(I_base,trans2lefttop(I2_t_max_or),1,1);
    %     I_t1 = trans2lefttop(I_t);
    %     I_t2 = imtranslate(I_t1,[l_x(I_t1),0]);
    %
    %     I = I_t1 | I_t2;
    %     j = 0;
    %     while sum(I(:)) == 2 * sum(I_t(:))
    %         j = j - 1;
    %         I_t2 = imtranslate(I_t2,[-1,0]);
    %         I = I_t1 | I_t2;
    %     end
    %     X_interval = l_x(I_t1) + j + 1;
    %
    %
    %     I_t1 = I_t1 | imtranslate(I_t1, [X_interval, 0]);
    %     I_t2 =imtranslate(I_t1,[0,l_y(I_t1)]);
    %     I = I_t1 | I_t2;
    %     j = 0;
    %     while sum(I(:)) == 4 * sum(I_t(:))
    %         j = j - 1;
    %         I_t2 = imtranslate(I_t2,[0,-1]);
    %         I = I_t1 | I_t2;
    %     end
    %     Y_interval = l_y(I_t1) + j + 1;
    %
    %     %如果图形不旋转
    %     lx1 = l_x(I_t);
    %     ly1 = l_y(I_t);
    %     X_interval1 = X_interval;
    %     Y_interval1 = Y_interval;
    %
    %     %如果图形发生旋转
    %     lx2 = l_y(I_t);
    %     ly2 = l_x(I_t);
    %     X_interval2 = Y_interval;
    %     Y_interval2 = X_interval;
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %首先，依靠x轴进行铺图，即所有的图保持水平,那么heiban可以分为heiban_top,heiban_bo
    %如果图形不旋转，则X,Y方向最多能排的图形的数目N1x1,N1y1分别为：
    if lx_heiban < lx1
        N1x1 = 0;
    else
        N1x1 = floor((lx_heiban -  lx1) / X_interval1) + 1;
    end
    
    if ly_heiban < ly1
        N1y1 = 0;
    else
        N1y1 = floor((ly_heiban -  ly1) / Y_interval1) + 1;
    end
    %如果图形发生旋转，则X,Y方向最多能排的图形的数目N1x2,N1y2分别为：
    if lx_heiban < lx2
        N1x2 = 0;
    else
        N1x2 = floor((lx_heiban -  lx2) / X_interval2) + 1;
    end
    
    if ly_heiban < ly2
        N1y2 = 0;
    else
        N1y2 = floor((ly_heiban -  ly2) / Y_interval2) + 1;
    end
    
    N1 = zeros(1,3);
    %下面揣测怎么铺图录用率更高
    %m,n分别为不旋转图形和旋转后的图形所排的行数
    for m = 0 : N1y1
        if m == 0
            n = N1y2;
        else
            if (ly_heiban - ly1 - (m - 1) * Y_interval1 - ly2) < 0
                n = 0;
            else
                n = floor(((ly_heiban - ly1 - (m - 1) * Y_interval1) - ly2) / Y_interval2) + 1;
            end
        end
        N1(m+1,1) = m;
        N1(m+1,2) = n;
        N1(m+1,3) = m * N1x1 + n * N1x2;
    end
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %其次，依靠y轴进行铺图，即所有的图保持竖直,那么heiban可以分为heiban_lef,heiban_rig
    %如果图形不旋转，则X,Y方向最多能排的图形的数目N2x1,N2y1分别为：
    if lx_heiban <  lx1
        N2x1 = 0;
    else
        N2x1 = floor((lx_heiban -  lx1) / X_interval1) + 1;
    end
    
    if ly_heiban <  ly1
        N2y1 = 0;
    else
        N2y1 = floor((ly_heiban -  ly1) / Y_interval1) + 1;
    end
    
    %如果图形发生旋转，则X,Y方向最多能排的图形的数目N2x2,N2y2分别为：
    if lx_heiban <  lx2
        N2x2 = 0;
    else
        N2x2 = floor((lx_heiban -  lx2) / X_interval2) + 1;
    end
    
    if ly_heiban <  ly2
        N2y2 = 0;
    else
        N2y2 = floor((ly_heiban -  ly2) / Y_interval2) + 1;
    end
    
    N2 = zeros(1,3);
    %下面揣测怎么铺图录用率更高
    %p,q分别为不旋转图形和旋转后的图形所排的列数
    for p = 0 : N2x1
        if p == 0
            q = N2x2;
        else
            if (lx_heiban - lx1 - (p - 1) * X_interval1 - lx2) < 0
                q = 0;
            else
                q = floor(((lx_heiban - lx1 - (p - 1) * X_interval1) - lx2) / X_interval2) + 1;
            end
        end
        N2(p+1,1) = p;
        N2(p+1,2) = q;
        N2(p+1,3) = p * N2y1 + q * N2y2;
    end
    
    N1_ascend = sortrows(N1,3);
    N2_ascend = sortrows(N2,3);
    N1max = N1_ascend(end,:);
    N2max = N2_ascend(end,:);
    
    if N1max(1,3) >= N2max(1,3)
        Nmax = N1max;
    else
        Nmax = N2max;
    end
    
    if Nmax(1,3) > 0
        if N1max(1,3) >= N2max(1,3)
            x_along = 1;
        else
            y_along = 1;
        end
        
        coordinate = zeros(1,3);
        
        %暂时先放N1,再放N2
        mm = 0;
        if x_along == 1
            for j = 1 : Nmax(1,1)%未旋转的图形
                for i = 1 : N1x1
                    mm = mm + 1;
                    coordinate(mm,1) = (i - 1) * X_interval1 + lx1 / 2;
                    coordinate(mm,2) = (j - 1) * Y_interval1 + ly1 / 2;
                    coordinate(mm,3) = 0;
                end
            end
            
            if Nmax(1,1) == 0
                Y_occupied = 0;
            else
                Y_occupied = (Nmax(1,1) - 1) * Y_interval1 + ly1;
            end
            
            for j = 1 : Nmax(1,2)%旋转后的图形
                for i = 1 : N1x2
                    mm = mm + 1;
                    coordinate(mm,1) = (i - 1) * X_interval2 + lx2 / 2;
                    coordinate(mm,2) = (j - 1) * Y_interval2 + ly2 / 2 + Y_occupied;
                    coordinate(mm,3) = 90;
                end
            end
            
        elseif y_along == 1
            
            for j = 1 : Nmax(1,1) %Nmax(1,1)表示不旋转图形排的列数
                for i = 1 : N1y1%未旋转的图形
                    mm = mm + 1;
                    coordinate(mm,1) = (j - 1) * X_interval1 + lx1 / 2;
                    coordinate(mm,2) = (i - 1) * Y_interval1 + ly1 / 2;
                    coordinate(mm,3) = 0;
                end
            end
            
            if Nmax(1,1) == 0
                X_occupied = 0;
            else
                X_occupied = (Nmax(1,1) - 1) * X_interval1 + lx1;
            end
            
            for j = 1 : Nmax(1,2)%这里排旋转后的图形排的列数
                for i = 1 : N1y2
                    mm = mm + 1;
                    coordinate(mm,1) = (j - 1) * X_interval2 + lx2 / 2 + X_occupied;
                    coordinate(mm,2) = (i - 1) * Y_interval2 + ly2 / 2;
                    coordinate(mm,3) = 90;
                end
            end
        end
        
        N = length(coordinate(:,1));
    else
        coordinate = [];
        N=0;
        
        % if N1max(1,3) >= N2max(1,3)
        %     sita1 = 0;
        %     N = N1max(1,3);
        % else
        %     sita1 = 90;
        %     N = N2max(1,3);
        % end
        
        % struct5.alg1 = [];
        % struct5.alg2 = [];
        % struct5.dx1 = [];
        % struct5.dx2 = [];
        % struct5.dx3 = [];
        % struct5.dy1 = [];
        % struct5.dy2 = [];
        % struct5.dy3 = [];
        % struct5.sita10 = [];
        % struct5.sita1 = sita1;
        % struct5.beta = [];
        % struct5.I1 = [];
        % struct5.I2 = [];
        % struct5.I3 = [];
        % struct5.SP = [];
        % struct5.XP = [];
        % struct5.ZP = [];
        % struct5.YP = [];
        % struct5.N = N;
        % struct5.ii = [];
    end
end

