function [struct1] = find_best_dxy1(heiban,I2_t0,I2_t,I2_t_max,Parent0,template,...
    key,maxIt_va,nPop,nVar,n,lx1,ly1,nC,maxIt,W0,L0,kb,partgap,N_times,imrf_special,sita_xz)
%FIND_BEST_DXY 此处显示有关此函数的摘要
%   此处显示详细说明
if max(lx1,ly1) > 500
    imrf = 0.5;
else
    imrf = 1;
end
heiban = imresize(heiban,imrf);
I2_t0 = imresize(I2_t0,imrf);
I2_t = imresize(I2_t,imrf);
I2_t_max = imresize(I2_t_max,imrf);
lx1 = lx1 * imrf;
ly1 = ly1 * imrf;
partgap = partgap * imrf;

tx = 0;
ijmax1 = 6;
ijmax2 = 8;
ijmax3 = 10;
alg1=3;
alg2=3;

struct.alg1=3;
struct.alg2=3;

struct.dx1=[];
struct.dy1=[];

struct.dx2=[];
struct.dy2=[];

struct.dx3=[];
struct.dy3=[];

struct.sita10=[];
struct.sita1=[];
struct.beta=[];

struct.I1=[];
struct.I2=[];
struct.I3=[];

struct.SP=[];
struct.XP=[];
struct.ZP=[];
struct.YP=[];
struct.N=[];

struct1=repmat(struct, 1, 1);

[dx1,dx2,dx3,dy1,dy2,dy3,sita10,sita1,beta,~,I1,I2,~,SP,XP,ZP,YP,ii,RF] = ...%beta的值只有两种可能，要么是0，要么是90
    find_position(I2_t,Parent0,template,key,alg1,alg2,maxIt_va,nPop,nVar,n,lx1,ly1,nC,maxIt,W0,L0,kb,N_times,tx,ijmax1,ijmax2,ijmax3);
%寻找第1次排图的高精度位置
I2_t_max = trans2lefttop(I2_t_max);
I2_t_max = imcrop(I2_t_max,[0,0,l_x(I2_t_max),l_y(I2_t_max)]);
I_base = zeros(ceil(l_y(I2_t_max)*(l_y(I1)/l_y(I2_t)))+round(8*5)+2,ceil(l_x(I2_t_max)*(l_x(I1)/l_x(I2_t))+round(8*5))+2);
I2_t_max = paste_img(I_base,I2_t_max,1,1);
I2_t_max = trans2mid(I2_t_max,0);
I2_t_max_10 = imresize(I2_t_max,imrf_special);
I2_t_max_10(I2_t_max_10<0.5)=0;
I2_t_max_10 = logical(I2_t_max_10);
if dx1 ~= 0
    dx1 = abs(dx1) / dx1 * (l_x(I1) - l_x(I2_t_max));
end
if dy1 ~= 0
    dy1 = abs(dy1) / dy1 * (l_y(I1) - l_y(I2_t_max));
end
[dx111,dy111,I1_hr] = high_accuracy1(I2_t_max,alg1,key,dx1,dy1,sita10,1,tx,4,RF);
%寻找第2次排图的高精度位置
I1_hr_t = trans2lefttop(I1_hr);
I1_hr_t = imcrop(I1_hr_t,[0,0,l_x(I1_hr_t),l_y(I1_hr_t)]);
I1_hr_10 = imresize(I1_hr_t,imrf_special);
I_base = zeros(ceil(l_y(I1_hr_t)*(l_y(I2)/l_y(I1)))+round(9*5)+2,ceil(l_x(I1_hr_t)*(l_x(I2)/l_x(I1))+round(9*5))+2);
I1_hr_t = paste_img(I_base,I1_hr_t,1,1);
I1_hr_t = trans2mid(I1_hr_t,0);
if dx2 ~= 0
    dx2 = abs(dx2) / dx2 * (l_x(I2) - l_x(I1));
end
if dy2 ~= 0
    dy2 = abs(dy2) / dy2 * (l_y(I2) - l_y(I1));
end
[dx222,dy222,I2_hr] = high_accuracy2(I1_hr_t,alg2,key,dx2,dy2,1,5,RF);%这里的dx2要换

if abs(dx222) >= abs(dy222)
    beta1 = 90;
else
    beta1 = 0;
end

if beta == beta1
    ...
else %beta ~= beta1
if beta == 0 && beta1 == 90
    beta = 90;
    dx3_t = -dy3;
    dy3_t = dx3;
    dx3 = dx3_t;
    dy3 = dy3_t;
else %beta == 90 && beta1 == 0
    beta = 0;
    dx3_t = dy3;
    dy3_t = -dx3;
    dx3 = dx3_t;
    dy3 = dy3_t;
end
end

sita1 = sita1 + beta;
if beta == 90
    dx111_ro =  dy111;
    dy111_ro = -dx111;
else
    dx111_ro = dx111;
    dy111_ro = dy111;
end

if beta == 90
    dx222_ro =  dy222;
    dy222_ro = -dx222;
else
    dx222_ro = dx222;
    dy222_ro = dy222;
end
%寻找第3次排图的高精度位置
kk_end = 1 + floor(abs(40/dy222_ro));
sita0 = atan(-dx222_ro/dy222_ro)*180/pi;%最终结果180度数制
% if dx222_ro ~= 0 && dy222_ro ~= 0
gaodu = zeros(1,2);

I3_t = imrotate(I1_hr,beta,'bicubic');
I3_t = trans2lefttop(I3_t);
I3_t = imcrop(I3_t,[0,0,l_x(I3_t),l_y(I3_t)]);
I_base = zeros(l_y(I3_t) + (2 * kk_end + 1) * abs(dy222_ro),3 * l_x(I3_t) + (2 * kk_end + 1) * abs(dx222_ro));
I3_t = paste_img(I_base,I3_t,1,1); I3_t = logical(I3_t);
I3_t = trans2mid(I3_t,0);
I3_tt = I3_t;
for kk = 1 : kk_end
    I3_tt = imtranslate(I3_tt,[dx222_ro,dy222_ro]);
    I3_t = I3_t | I3_tt;
end
I3_t = trans2lefttop(I3_t);
I3_t = imcrop(I3_t,[0,0,l_x(I3_t),l_y(I3_t)]);

Ymax = round(1.1*(1 * l_y(I3_t) + 2 * abs(dy3) + (2 * ii * abs(dy222_ro))));
Xmax = round(1.1*(1 * l_x(I3_t) + 2 * abs(dx3) + (2 * ii * abs(dx222_ro))));
N_I3_t_pix = sum(I3_t(:));

I_base = zeros(Ymax,Xmax);
I3_t = paste_img(I_base,I3_t,1,1); I3_t = logical(I3_t);
I3_t = trans2mid(I3_t,0);
I3_t_r = imtranslate(I3_t,[round(dx3 - ii * dx222_ro), round(dy3 - ii * dy222_ro)]);
I3_t_r1 = imtranslate(I3_t,[round(dx3 + ii * dx222_ro), round(dy3 + ii * dy222_ro)]);

while sum(I3_t_r(:)) ~= N_I3_t_pix || sum(I3_t_r1(:)) ~= N_I3_t_pix
    Ymax = round(1.2 * Ymax);
    Xmax = round(1.2 * Xmax);
    I_base = zeros(Ymax,Xmax);
    I3_t = paste_img(I_base,I3_t,1,1); I3_t = logical(I3_t);
    I3_t = trans2mid(I3_t,0);
    I3_t_r = imtranslate(I3_t,[round(dx3 - ii * dx222_ro), round(dy3 - ii * dy222_ro)]);
    I3_t_r1 = imtranslate(I3_t,[round(dx3 + ii * dx222_ro), round(dy3 + ii * dy222_ro)]);
end

I3_t_combine = I3_t | I3_t_r;
I3_t_combine_ro = imrotate(I3_t_combine,sita0,'bicubic');
x_del = min(round(9.8 * l_x(I3_t_combine_ro) / 20),round(l_x(I3_t_combine_ro)/2)-4);
m = 0;

for i=-ii:ii
    m = m + 1;
    I3_t_r = imtranslate(I3_t,[dx3 + i * dx222_ro, dy3 + i * dy222_ro]);
    %imshow(I3_t_r);
    I3_t_combine = I3_t | I3_t_r;
    I3_t_combine_ro = imrotate(I3_t_combine,sita0,'bicubic');
    %x_del = round(9 * l_x(I3_t_combine_ro) / 20);
    I3_t_combine_ro_t = trans2lefttop(I3_t_combine_ro);
    I3_t_combine_ro_t = imtranslate(I3_t_combine_ro_t,[-x_del,0]);
    
    while (sum(I3_t_combine_ro_t(:))) <= (sum(I3_t_combine_ro(:)) * 0.505)
        x_del = min(round(x_del * 0.8), x_del-2);
        I3_t_combine_ro_t = trans2lefttop(I3_t_combine_ro);
        I3_t_combine_ro_t = imtranslate(I3_t_combine_ro_t,[-x_del,0]);
    end
    
    I3_t_combine_ro_t = trans2righttop(I3_t_combine_ro_t);
    I3_t_combine_ro_t = imtranslate(I3_t_combine_ro_t,[ x_del,0]);
    I3_t_combine_ro_t = trans2mid(I3_t_combine_ro_t,0);
    %imshow(I3_t_combine_ro_t),pause(0.5)
    gaodu(m,1) = i;
    l_y_now = l_y(I3_t_combine_ro_t);
    gaodu(m,2) = l_y_now;
end


[~, so] = sort(gaodu(:,2), 'ascend');
i_best = gaodu(so(1),1);

while i_best == abs(ii)
    ii = ii + 1;
    Ymax = round(1.1*(1 * l_y(I3_t) + 2 * abs(dy3) + (2 * ii * abs(dy222_ro))));
    Xmax = round(1.1*(1 * l_x(I3_t) + 2 * abs(dx3) + (2 * ii * abs(dx222_ro))));
    I_base = zeros(Ymax,Xmax);
    I3_t = paste_img(I_base,I3_t,1,1); I3_t = logical(I3_t);
    I3_t = trans2mid(I3_t,0);
    I3_t_r = imtranslate(I3_t,[round(dx3 - ii * dx222_ro), round(dy3 - ii * dy222_ro)]);
    I3_t_r1 = imtranslate(I3_t,[round(dx3 + ii * dx222_ro), round(dy3 + ii * dy222_ro)]);
    
    while sum(I3_t_r(:)) ~= N_I3_t_pix || sum(I3_t_r1(:)) ~= N_I3_t_pix
        Ymax = round(1.2 * Ymax);
        Xmax = round(1.2 * Xmax);
        I_base = zeros(Ymax,Xmax);
        I3_t = paste_img(I_base,I3_t,1,1); I3_t = logical(I3_t);
        I3_t = trans2mid(I3_t,0);
        I3_t_r = imtranslate(I3_t,[round(dx3 - ii * dx222_ro), round(dy3 - ii * dy222_ro)]);
        I3_t_r1 = imtranslate(I3_t,[round(dx3 + ii * dx222_ro), round(dy3 + ii * dy222_ro)]);
    end
    
    m = 0;
    
    for i=-ii:ii
        m = m + 1;
        I3_t_r = imtranslate(I3_t,[dx3 + i * dx222_ro, dy3 + i * dy222_ro]);
        %imshow(I3_t_r);
        I3_t_combine = I3_t | I3_t_r;
        I3_t_combine_ro = imrotate(I3_t_combine,sita0,'bicubic');
        %x_del = round(9 * l_x(I3_t_combine_ro) / 20);
        I3_t_combine_ro_t = trans2lefttop(I3_t_combine_ro);
        I3_t_combine_ro_t = imtranslate(I3_t_combine_ro_t,[-x_del,0]);
        
        while (sum(I3_t_combine_ro_t(:))) <= (sum(I3_t_combine_ro(:)) * 0.505)
            x_del = min(round(x_del * 0.8), x_del-2);
            I3_t_combine_ro_t = trans2lefttop(I3_t_combine_ro);
            I3_t_combine_ro_t = imtranslate(I3_t_combine_ro_t,[-x_del,0]);
        end
        
        I3_t_combine_ro_t = trans2righttop(I3_t_combine_ro_t);
        I3_t_combine_ro_t = imtranslate(I3_t_combine_ro_t,[ x_del,0]);
        I3_t_combine_ro_t = trans2mid(I3_t_combine_ro_t,0);
        %imshow(I3_t_combine_ro_t),pause(0.5)
        gaodu(m,1) = i;
        l_y_now = l_y(I3_t_combine_ro_t);
        gaodu(m,2) = l_y_now;
    end
    
    [~, so] = sort(gaodu(:,2), 'ascend');
    i_best = gaodu(so(1),1);
end

dx3 = dx3 + i_best * dx222_ro;
dy3 = dy3 + i_best * dy222_ro;

I2_hr_3 = I3_t;

I2_hr_3 = trans2lefttop(I2_hr_3);
I2_hr_3 = imcrop(I2_hr_3,[0,0,l_x(I2_hr_3),l_y(I2_hr_3)]);

I2_hr_t = I3_t;
I2_hr_3_10 = imresize(I2_hr_3,imrf_special);

I2_hr = trans2lefttop(I2_hr);
I2_hr = imcrop(I2_hr,[0,0,l_x(I2_hr),l_y(I2_hr)]);
I2_hr_10 = imresize(I2_hr,imrf_special);
I_base = zeros(1*l_y(I2_hr_t)+round(abs(dy3))+60,1*l_x(I2_hr_t)+round(abs(dx3))+60);
I2_hr_t=trans2lefttop(I2_hr_t);
I2_hr_t = imcrop(I2_hr_t,[0,0,l_x(I2_hr_t),l_y(I2_hr_t)]);

I2_hr_t = paste_img(I_base,I2_hr_t,1,1);
I2_hr_t = trans2mid(I2_hr_t,0);
[dx333_ro,dy333_ro,I3_hr] = high_accuracy3(I2_hr_t,dx3,dy3,sita0,1,6,RF);
I3_hr = trans2lefttop(I3_hr);
I3_hr = imcrop(I3_hr,[0,0,l_x(I3_hr),l_y(I3_hr)]);
I3_hr_10 = imresize(I3_hr,imrf_special);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%A%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%以上已经找到了高精度的位置，下面进行排图
if beta == 90
    ZP = 1 - ZP;
    YP = 1 - YP;
end

dx1=dx111_ro;
dy1=dy111_ro;

dx2=dx222_ro;
dy2=dy222_ro;

if dy2 > 0
    dx2=-dx2;
    dy2=-dy2;
end

dx3=dx333_ro;
dy3=dy333_ro;
%[dx1,dy1,dx2,dy2,dx3,dy3,sita1] = find_suit_dxy(ZP,dx1,dy1,dx2,dy2,dx3,dy3,sita1);%最终不能确定dx3的符号，但是可以确定dy3的符号为非负
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%接下来将图排成水平，即将dy3变为0
if exist('i_best','var')
    gama = atan((dy333_ro - i_best * dy222_ro) / (dx333_ro - i_best * dx222_ro)) * 180 / pi;
    %gama为为了将图像进行水平排布为旋转的角度,有问题，因为有可能dx3和dy3不是相邻的两列的差值
else
    gama = atan(dy333_ro / dx333_ro) * 180 / pi;
end%经过验证，结果正确

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%下面将零件摆平，并寻找合适的间距值
[dx11,dy11,dx22,dy22,dx33,dy33,sita1,ZP,YP,~,~,~,~,~] = baiping(I1_hr,I2_hr,I3_hr,I2_hr_3,I2_hr_3_10,I2_t_max,I1_hr_10,I2_hr_10,I3_hr_10,...
    I2_t_max_10,sita1,beta,gama,ZP,YP,SP,XP,dx1,dy1,dx2,dy2,dx3,dy3,1,imrf_special);
%         imshow(I3_hr);
%以上所所得到的dx1,dx2,dx3,dy1,dy2,dy3均为摆平后的值

[dx11,dy11,dx22,dy22,dx33,dy33,sita1] = find_suit_dxy(ZP,dx11,dy11,dx22,dy22,dx33,dy33,sita1);%最终不能确定dx3的符号，但是可以确定dy3的符号为非负
%上面的表达式主要是调整dx3和dy3的值，并未对dx2和dy2做调整
sita2 = atan(-dx22/dy22);
sita0 = sita2*180/pi;%最终结果180度数制

I2_shu = imrotate(I1,sita1,'bicubic');
lmax = round(2*(max(l_y(I1)*2+abs(dy2)+abs(dy3),l_x(I1)*2+abs(dx2)+abs(dx3))));
I_base = zeros(lmax,lmax);
I2_shu = trans2lefttop(I2_shu);
I2_shu = imcrop(I2_shu,[0,0,l_x(I2_shu),l_y(I2_shu)]);
I2_shu = paste_img(I_base,I2_shu,1,1);
I2_shu = trans2mid(I2_shu,0);
%     I2_shu = imrotate(I2_shu,beta);
%     I2_shu = trans2rightbottom(I2_shu);
I2_shu_r = I2_shu;
I2_shu_r = imtranslate(I2_shu_r,[dx2,dy2]);
I2_shu = I2_shu | I2_shu_r;%正确

I3_shu = I2_shu;
I3_shu = trans2mid(I3_shu,0);
I3_shu_r = I3_shu;
I3_shu_r = imtranslate(I3_shu_r,[dx3,dy3]);
I3_shu = I3_shu | I3_shu_r;

I2_shu = trans2mid(I2_shu,0);
I3_shu = trans2mid(I3_shu,0);
I2_shu = imrotate(I2_shu,sita0,'bicubic');
I3_shu = imrotate(I3_shu,sita0,'bicubic');
%         imshow(I3_shu)

d3_shu = l_x(I3_shu)-l_x(I2_shu); %这里有问题
d3_shu = floor(d3_shu / cos(sita2));%/imrf;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
N_shu = ceil(numel(heiban(:,1)) / abs(dy2));  %结果正确 其表示纵向最多能排多少个图案
delta_xx = abs(numel(heiban(:,1)) * tan(sita2));%斜铺后横向多出来得到部分像素（实际上不存在）
x_total = delta_xx + numel(heiban(1,:));
N_dx3 = ceil(numel(heiban(1,:)) / d3_shu);%仅仅在板子上，横向走能排多少个
N_dx3_max = ceil(x_total / d3_shu);%+3000;%加上扩展的区域，横向走能排多少个

I2_t0_ex = logical(expand(I2_t0,partgap));
I_s_m = imrotate(I2_t0_ex,sita1,'bicubic');
lx_single = l_x(I_s_m);
ly_single = l_y(I_s_m);

N_to = zeros(abs(round(dx33) * round(dy22)),1);
q=0;
for x_py = 0:-2:-abs(round(dx33))
    for y_py = 0:2:abs(round(dy22))
        q = q + 1;
        [~,~,N] = paitu_pre(ZP,dx11,dx22,dx33,dy11,dy22,dy33,sita1,N_dx3,N_dx3_max,N_shu,heiban,lx_single,ly_single,x_py,y_py,sita_xz);
        N_to(q,1) = N;
    end
end
N = max(N_to);

dx1 = dx111/imrf;
dy1 = dy111/imrf;
dx2 = dx222/imrf;
dy2 = dy222/imrf;
dx3 = dx3/imrf;
dy3 = dy3/imrf;

struct1.alg1 = alg1;
struct1.alg2 = alg2;
struct1.dx1 = dx1;
struct1.dx2 = dx2;
struct1.dx3 = dx3;
struct1.dy1 = dy1;
struct1.dy2 = dy2;
struct1.dy3 = dy3;
struct1.sita10 = sita10;
struct1.sita1 = sita10;
struct1.beta = beta;
struct1.I1 = imresize(I1_hr,1/imrf);
struct1.I2 = imresize(I2_hr,1/imrf);
struct1.I3 = imresize(I3_hr,1/imrf);
struct1.SP = SP;
struct1.XP = XP;
struct1.ZP = ZP;
struct1.YP = YP;
struct1.N = N;
struct1.ii = ii;
struct1.coordinate = [];
% dxy_struct=[struct1,struct2,struct3,struct4];
disp('find_best_dxy1完成');
end
