function [struct5] = find_best_dxy5(lx_heiban,ly_heiban,I2_t_max,imrf,sita_xz)
%FIND_BEST_DXY5 此处显示有关此函数的摘要
%find_best_dxy5函数对图像完全是横平竖直的排列
%   此处显示详细说明
% imrf = imrf * imrf_special;
S.number = [];
S.coordinate = [];
% x_along = 0;
% y_along = 0;
% lx_heiban = numel(heiban(1,:));
% ly_heiban = numel(heiban(:,1));

%首先，图形不旋转的情况下，求出X_interval，Y_interval
I2_t_max = logical(I2_t_max);
I_base = zeros(2 * l_y(I2_t_max), 2 * l_x(I2_t_max));
% I_base = logical(I_base);
I_t = logical(paste_img(I_base,trans2lefttop(I2_t_max),1,1));
I_t1 = trans2lefttop(I_t);
I_t2 = imtranslate(I_t1,[l_x(I_t1),0]);

I = I_t1 | I_t2;
j = 0;
while sum(I(:)) == 2 * sum(I_t(:))
    j = j - 1;
    I_t2 = imtranslate(I_t2,[-1,0]);
    I = I_t1 | I_t2;
end
X_interval = l_x(I_t1) + j + 1;


I_t1 = I_t1 | imtranslate(I_t1, [X_interval, 0]);
I_t2 =imtranslate(I_t1,[0,l_y(I_t1)]);
I = I_t1 | I_t2;
j = 0;
while sum(I(:)) == 4 * sum(I_t(:))
    j = j - 1;
    I_t2 = imtranslate(I_t2,[0,-1]);
    I = I_t1 | I_t2;
end
Y_interval = l_y(I_t1) + j + 1;

%如果图形不旋转
lx1 = l_x(I_t) / imrf;
ly1 = l_y(I_t) / imrf;
X_interval1 = X_interval / imrf;
Y_interval1 = Y_interval / imrf;

%如果图形发生旋转
lx2 = l_y(I_t) / imrf;
ly2 = l_x(I_t) / imrf;
X_interval2 = Y_interval / imrf;
Y_interval2 = X_interval / imrf;

%这里求长度系列
Nxmax = ceil(max(lx_heiban,ly_heiban) / lx1);
Nymax = ceil(max(lx_heiban,ly_heiban) / ly1);

L_ser = zeros(1);
kk = 0;
for i = 0 : Nxmax
    for j = 0 : Nymax
        kk = kk + 1;
        L_ser(kk,1) = i * lx1 + j * ly1;
        
    end
end
L_ser = [L_ser;lx_heiban;ly_heiban];
[L_ser,~,~] = unique(L_ser);
L_ser = sort(L_ser);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%首先，依靠x轴进行铺图，即所有的图保持水平,那么heiban可以分为heiban_top,heiban_bot
%如果图形不旋转，则X,Y方向最多能排的图形的数目N1x1,N1y1分别为：
% if lx_heiban < lx1
%     N1x1 = 0;
% else
%     N1x1 = floor((lx_heiban -  lx1) / X_interval1) + 1;
% end

% if ly_heiban < ly1
%     N1y1 = 0;
% else
%     N1y1 = floor((ly_heiban -  ly1) / Y_interval1) + 1;
% end
%如果图形发生旋转，则X,Y方向最多能排的图形的数目N1x2,N1y2分别为：
% if lx_heiban < lx2
%     N1x2 = 0;
% else
%     N1x2 = floor((lx_heiban -  lx2) / X_interval2) + 1;
% end
%
% if ly_heiban < ly2
%     N1y2 = 0;
% else
%     N1y2 = floor((ly_heiban -  ly2) / Y_interval2) + 1;
% end

% N1 = zeros(1,3);
%下面揣测怎么铺图录用率更高
%m,n分别为不旋转图形和旋转后的图形所排的行数
N1 = repmat(S, 1, 1);
for a = 1 : length(L_ser)
    L_try = L_ser(a);
    
    if L_try > ly_heiban
        break
    end
    
    if L_try == 0
        %         heiban_bot = [];
        %         heiban_top =heiban;
        lx_heiban_bot = 0;
        ly_heiban_bot = 0;
        
        lx_heiban_top = lx_heiban;
        ly_heiban_top = ly_heiban;
    else
        %         heiban_bot = zeros(ceil(L_try),lx_heiban);
        %         heiban_top = zeros(ly_heiban - ceil(L_try),lx_heiban);
        
        ly_heiban_bot = L_try;
        lx_heiban_bot = lx_heiban;
        
        ly_heiban_top = ly_heiban - L_try;
        lx_heiban_top = lx_heiban;
        
    end
    
    [coordinate_top,N_top] = find_best_dxy5_1(lx_heiban_top,ly_heiban_top,X_interval1,Y_interval1,X_interval2,Y_interval2,lx1,ly1,lx2,ly2);
    [coordinate_bot,N_bot] = find_best_dxy5_1(lx_heiban_bot,ly_heiban_bot,X_interval1,Y_interval1,X_interval2,Y_interval2,lx1,ly1,lx2,ly2);
    
    if isempty(coordinate_top) == 0
        coordinate_top(:,2) = coordinate_top(:,2) + L_try;
    end
    
    N1(a).number = N_top + N_bot;
    N1(a).coordinate = [coordinate_bot;coordinate_top];
end
% for m = 0 : N1y1
%     if m == 0
%         heiban_bot = [];
%         heiban_top =heiban;
%     else
%         heiban_bot = zeros((m - 1) * Y_interval1 + ly1,lx_heiban);
%         heiban_top = zeros(ly_heiban - ((m - 1) * Y_interval1 + ly1),lx_heiban);
%     end
%
%     [coordinate_top,N_top] = find_best_dxy5_1(heiban_top,X_interval1,Y_interval1,X_interval2,Y_interval2,lx1,ly1,lx2,ly2);
%     [coordinate_bot,N_bot] = find_best_dxy5_1(heiban_bot,X_interval1,Y_interval1,X_interval2,Y_interval2,lx1,ly1,lx2,ly2);
%     if isempty(coordinate_top) == 0 && m > 0
%         coordinate_top(:,2) = coordinate_top(:,2) + (m - 1) * Y_interval1 + ly1;
%     end
%     N1(m+1).number = N_top + N_bot;
%     N1(m+1).coordinate = [coordinate_bot;coordinate_top];
% end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%其次，依靠y轴进行铺图，即所有的图保持竖直,那么heiban可以分为heiban_lef,heiban_rig
%如果图形不旋转，则X,Y方向最多能排的图形的数目N2x1,N2y1分别为：
% if lx_heiban <  lx1
%     N2x1 = 0;
% else
%     N2x1 = floor((lx_heiban -  lx1) / X_interval1) + 1;
% end

% if ly_heiban <  ly1
%     N2y1 = 0;
% else
%     N2y1 = floor((ly_heiban -  ly1) / Y_interval1) + 1;
% end

%如果图形发生旋转，则X,Y方向最多能排的图形的数目N2x2,N2y2分别为：
% if lx_heiban <  lx2
%     N2x2 = 0;
% else
%     N2x2 = floor((lx_heiban -  lx2) / X_interval2) + 1;
% end
%
% if ly_heiban <  ly2
%     N2y2 = 0;
% else
%     N2y2 = floor((ly_heiban -  ly2) / Y_interval2) + 1;
% end

% N2 = zeros(1,3);
%下面揣测怎么铺图录用率更高
%p,q分别为不旋转图形和旋转后的图形所排的列数
N2 = repmat(S, 1, 1);
for a = 1 : length(L_ser)
    L_try = L_ser(a);
    
    if L_try > lx_heiban
        break
    end
    
    if L_try == 0
        %         heiban_lef = [];
        %         heiban_rig =heiban;
        
        lx_heiban_lef = 0;
        ly_heiban_lef = 0;
        
        lx_heiban_rig =lx_heiban;
        ly_heiban_rig =ly_heiban;
    else
        %         heiban_lef = zeros(ly_heiban,ceil(L_try));
        %         heiban_rig = zeros(ly_heiban,lx_heiban - ceil(L_try));
        
        lx_heiban_lef = L_try;
        ly_heiban_lef = ly_heiban;
        
        lx_heiban_rig = lx_heiban - L_try;
        ly_heiban_rig = ly_heiban;
        
    end
    
    [coordinate_lef,N_lef] = find_best_dxy5_1(lx_heiban_lef,ly_heiban_lef,X_interval1,Y_interval1,X_interval2,Y_interval2,lx1,ly1,lx2,ly2);
    [coordinate_rig,N_rig] = find_best_dxy5_1(lx_heiban_rig,ly_heiban_rig,X_interval1,Y_interval1,X_interval2,Y_interval2,lx1,ly1,lx2,ly2);
    
    if isempty(coordinate_rig) == 0
        coordinate_rig(:,1) = coordinate_rig(:,1) + L_try;
    end
    
    N2(a).number = N_lef + N_rig;
    N2(a).coordinate = [coordinate_lef;coordinate_rig];
end
% for p = 0 : N2x1
%     if p == 0
%         heiban_lef = [];
%         heiban_rig =heiban;
%     else
%         heiban_lef = zeros(ly_heiban,(p - 1) * X_interval1 + lx1);
%         heiban_rig = zeros(ly_heiban,lx_heiban - ((p - 1) * X_interval1 + lx1));
%     end
%
%     [coordinate_lef,N_lef] = find_best_dxy5_1(heiban_lef,X_interval1,Y_interval1,X_interval2,Y_interval2,lx1,ly1,lx2,ly2);
%     [coordinate_rig,N_rig] = find_best_dxy5_1(heiban_rig,X_interval1,Y_interval1,X_interval2,Y_interval2,lx1,ly1,lx2,ly2);
%     if isempty(coordinate_rig) == 0 && p > 0
%         coordinate_rig(:,1) = coordinate_rig(:,1) + (p - 1) * X_interval1 + lx1;
%     end
%     N2(p+1).number = N_lef + N_rig;
%     N2(p+1).coordinate = [coordinate_lef;coordinate_rig];
% end

N = [N1,N2];
[~, so]=sort([N.number], 'descend'); %so为降序后的序号排序
N_paixu=zeros(length(so),1);
for k = 1 : length(so)
    N_paixu(k) =  N(so(k)).number;
end
coordinate_best = N((so(1))).coordinate;
coordinate_best = sortrows(coordinate_best,[1 2]);
% N1_ascend = sortrows(N1,3);
% N2_ascend = sortrows(N2,3);
% N1max = N1_ascend(end,:);
% N2max = N2_ascend(end,:);
%
% if N1max(1,3) >= N2max(1,3)
%     Nmax = N1max;
% else
%     Nmax = N2max;
% end
%
% if N1max(1,3) >= N2max(1,3)
%     x_along = 1;
% else
%     y_along = 1;
% end
%
% coordinate = zeros(1,3);
%
% %暂时先放N1,再放N2
% mm = 0;
% if x_along == 1
%     for j = 1 : Nmax(1,1)%未旋转的图形
%         for i = 1 : N1x1
%             mm = mm + 1;
%             coordinate(mm,1) = (i - 1) * X_interval1 + lx1 / 2;
%             coordinate(mm,2) = (j - 1) * Y_interval1 + ly1 / 2;
%             coordinate(mm,3) = 0;
%         end
%     end
%
%     if Nmax(1,1) == 0
%         Y_occupied = 0;
%     else
%         Y_occupied = (Nmax(1,1) - 1) * Y_interval1 + ly1;
%     end
%
%     for j = 1 : Nmax(1,2)%旋转后的图形
%         for i = 1 : N1x2
%             mm = mm + 1;
%             coordinate(mm,1) = (i - 1) * X_interval2 + lx2 / 2;
%             coordinate(mm,2) = (j - 1) * Y_interval2 + ly2 / 2 + Y_occupied;
%             coordinate(mm,3) = 90;
%         end
%     end
%
% elseif y_along == 1
%
%     for j = 1 : Nmax(1,1) %Nmax(1,1)表示不旋转图形排的列数
%         for i = 1 : N1y1%未旋转的图形
%             mm = mm + 1;
%             coordinate(mm,1) = (j - 1) * X_interval1 + lx1 / 2;
%             coordinate(mm,2) = (i - 1) * Y_interval1 + ly1 / 2;
%             coordinate(mm,3) = 0;
%         end
%     end
%
%     if Nmax(1,1) == 0
%         X_occupied = 0;
%     else
%         X_occupied = (Nmax(1,1) - 1) * X_interval1 + lx1;
%     end
%
%     for j = 1 : Nmax(1,2)%这里排旋转后的图形排的列数
%         for i = 1 : N1y2
%             mm = mm + 1;
%             coordinate(mm,1) = (j - 1) * X_interval2 + lx2 / 2 + X_occupied;
%             coordinate(mm,2) = (i - 1) * Y_interval2 + ly2 / 2;
%             coordinate(mm,3) = 90;
%         end
%     end
% end

% if N1max(1,3) >= N2max(1,3)
%     sita1 = 0;
%     N = N1max(1,3);
% else
%     sita1 = 90;
%     N = N2max(1,3);
% end
coordinate_best(:,3) = coordinate_best(:,3) + sita_xz;

struct5.alg1 = [];
struct5.alg2 = [];
struct5.dx1 = [];
struct5.dx2 = [];
struct5.dx3 = [];
struct5.dy1 = [];
struct5.dy2 = [];
struct5.dy3 = [];
struct5.sita10 = [];
struct5.sita1 = [];
struct5.beta = [];
struct5.I1 = [];
struct5.I2 = [];
struct5.I3 = [];
struct5.SP = [];
struct5.XP = [];
struct5.ZP = [];
struct5.YP = [];
struct5.N = length(coordinate_best);
struct5.ii = [];
struct5.coordinate = coordinate_best;
disp('find_best_dxy5完成');
end

