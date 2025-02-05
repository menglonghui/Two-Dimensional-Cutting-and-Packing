function [dx1,dx2,dx3,dy1,dy2,dy3,sita10,sita1,beta,sita2,I1,I2,I3,SP,XP,ZP,YP,ii,RF] = ...
    find_position(I2_t,Parent0,template,key,alg1,alg2,maxIt_va,nPop,nVar,n,lx1,ly1,nC,maxIt,W0,L0,kb,N_times,tx,ijmax1,ijmax2,ijmax3) %#ok<*INUSD,*INUSL>
%FIND_POSITION 此处显示有关此函数的摘要
%   此处显示详细说明
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% I0=I2_t;
Lmax = max(l_x(I2_t),l_y(I2_t));

%RF图像缩小倍数

if Lmax <= 100
    RF = 1;
elseif Lmax > 100 && Lmax <= 400
    RF = 1/2;
else %Lmax > 200 && Lmax <= 500
    RF = 1/5;
end

I2_t = imresize(I2_t,RF);
I2_t = tuomao(I2_t);
lx1 = lx1 * RF;
ly1 = ly1 * RF;
W0 = ceil(W0 * RF);
L0 = ceil(L0 * RF);

%第一次迭代寻优
% newRc = zeros(max(round(2.1*l_y(I2_t)),numel(I2_t(:,1))),max(round(2.1*l_x(I2_t)),numel(I2_t(1,:))));
% I2_t_temp = paste_img(newRc,I2_t,1,1);
newRc = zeros(round(2.5 * l_y(I2_t)),round(2.5 * l_x(I2_t)));


I_base = logical(newRc);
I2_t = trans2lefttop(I2_t);
I2_t = imcrop(I2_t,[0,0,l_x(I2_t),l_y(I2_t)]);
I2_t = paste_img(I_base,I2_t,1,1);
I2_t = trans2mid(I2_t,0);
I2_t = logical(I2_t);

% I2_t = logical(I2_t_temp);%裁剪成newRc尺寸，newRc = zeros(min(2*lw,numel(I2_t(:,1))),min(2*lw,numel(I2_t(1,:))));
I2_t_r = imrotate(I2_t,180,'crop','bicubic');

% imshow(I2_t,'Border','tight');
% set(gcf,'ToolBar','none','MenuBar','none','NumberTitle','off');
% set(gcf, 'CloseRequestFcn', @(src,event) myCloseFcn(src,event,h));

Parent = Parent0;
flag = 0;

if isstruct(alg1)
    Alg1 = eval(['alg1.' key]);
else
    Alg1 = alg1;
end

while flag == 0
    Grandma = template;
    for It = 1 : maxIt_va %迭代次数
        disp(['外层迭代次数：',num2str(It)]);
        Mother = repmat(template, 1, 1);
        %         Parent = Initial_pop(Parent,nPop,nVar,tx);%
        Offspring = repmat(Parent, N_times * 3, 1);
        Parent = Initial_pop(Parent,nPop,nVar,tx);%
        Offspring = Initial_pop(Offspring,length(Offspring),nVar,tx);%
        
        Offspring = Variation(Offspring,nC,n,newRc,tx);%变异操作
        [Parent,Mother] = Up_pop_Mo1(Mother,Parent,Offspring,nPop,I2_t,I2_t_r,newRc,lx1,ly1,n,template,Alg1);%0.5s 更新种群,或者说更新父辈
        
        Mother_last.x = 0;
        Mother_last.y = 0;
        Mother_last.z = 0;
        
        It_Mo = 0;%Mother迭代次数
        %flag_y1 = 0;
        while It_Mo<maxIt%如果Mother超过maxIt次没有找到更好的解，则认为达到最优
            if isempty(Mother(1).x) == 0%Mother被赋值
                Parent(:) = Mother(end);%因为目前看的话Mother中的元素最优，所以直接Mother中的元素给Parent和Offspring进行赋值
                Offspring = repmat(Parent, N_times, 1);
                Offspring(:) = Mother(end);
                Offspring = Variation(Offspring,nC,n,newRc,tx);%变异操作
                [Parent,Mother] = Up_pop_Mo1(Mother,Parent(1),Offspring,nPop,I2_t,I2_t_r,newRc,lx1,ly1,n,template,Alg1);%2s多 更新种群,或者说更新父辈
            else
                while isempty(Mother(1).x)%循环，直至Mother被赋值，即出现完全分离的情况
                    Offspring = repmat(Parent, round(N_times * 5), 1);
                    Offspring = Variation(Offspring,nC,n,newRc,tx);%变异操作
                    [Parent,Mother] = Up_pop_Mo1(Mother,Parent,Offspring,nPop,I2_t,I2_t_r,newRc,lx1,ly1,n,template,Alg1);%2s多 更新种群,或者说更新父辈
                end
            end
            
            if Mother(end).z ~= Mother_last.z%Mother有更新
                Mother_last = Mother(end);
                %                 I = Mother_last.tu;
                %                 imshow(I,'Border','tight');
                disp(['Mother中的元素数目：',num2str(numel(Mother))]);
                It_Mo = 0;
            else%如果Mother没有更新
                It_Mo = It_Mo+1;
                disp(['Mother未更新次数：',num2str(It_Mo)]);
            end
        end
        
        if isempty(Grandma.x) == 1 %如果Grandma是空的
            Grandma = Mother(end);
        else
            if Grandma.z<Mother(end).z
                Grandma = Mother(end);
            end
        end
    end
    
    if It == maxIt_va
        flag = 1;
    end
    
    dx1 = Grandma.ddx;
    dy1 = Grandma.ddy;
    sita1 = Grandma.sita;
    I1 = Grandma.tu;
    %     I = trans2mid(I,0);
    I2_t = trans2lefttop(I2_t);
    I_base = zeros(ceil(l_y(I2_t)*(l_y(I1)/l_y(I2_t)))+20,ceil(l_x(I2_t)*(l_x(I1)/l_x(I2_t))+20));
    I2_t = paste_img(I_base,I2_t,1,1);
    [dx1,dy1,I1] = high_accuracy1(I2_t,alg1,key,dx1,dy1,sita1,1,0,ijmax1,1);
    %         imshow(I);
    disp('奶奶出现,第一次寻优完毕！');%得到目前最好的排布
    %     I1 = I;%I1为第一次组合成功后的结果
end

I2_t_temp = zeros(W0,L0);%W0和L0为空白板的尺寸

for i = 1:numel(I1(:,1))
    for j = 1:numel(I1(1,:))
        I2_t_temp(i,j) = I1(i,j);
    end
end

I2_t_temp = logical(I2_t_temp);

I2_t = I2_t_temp;%这里的I2_t是已经排好的两个图元的组合
I2_t = trans2leftbottom(I2_t);%将其移动到左下角
I2_t_double = I2_t;%I2_t_double为两个69式排好的图元，为后期的X方向的寻优做准备。

I2_t2 = I2_t;
I2_t2 = trans2mid(I2_t2,0);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%下面开始高度方向上的自循环，这是第二次迭代寻优
% eval(['newRc = zeros(round(2.2*l_y(I2_t2)),round(kb.' key '*l_x(I2_t2)));'])
newRc = zeros(round(2.3 * l_y(I2_t2)),round(2.3 * l_x(I2_t2)));
I_base = newRc;
I2_t2 = trans2lefttop(I2_t2);
I2_t2 = imcrop(I2_t2,[0,0,l_x(I2_t2),l_y(I2_t2)]);
I2_t2 = paste_img(I_base,I2_t2,1,1);
I2_t2 = logical(I2_t2);

lx2 = l_x(I2_t2);
ly2 = l_y(I2_t2);

I2_t2 = logical(I2_t2);
I2_t2 = trans2mid(I2_t2,0);
I2_t2_r = logical(trans2mid(imrotate(I2_t2,180,'crop','bicubic'),0));
% imshow(I2_t2,'Border','tight')

Parent = Parent0;%Parent0为空的
flag = 0;
% Alg2 = eval(['alg2.' key]);

if isstruct(alg1)
    Alg2 = eval(['alg2.' key]);
else
    Alg2 = alg2;
end
% N_times = N_times + 1;
while flag == 0
    Grandma = template;
    for It = 1 : maxIt_va %迭代次数
        disp(['外层迭代次数：',num2str(It)]);
        Mother = repmat(template, 1, 1);
        %         Parent = Initial_pop(Parent,nPop,nVar,0);%2s多   生成子代种群（包含交叉操作）
        Offspring = repmat(Parent, N_times, 1);
        Parent = Initial_pop(Parent,nPop,nVar,0);%2s多   生成子代种群（包含交叉操作）
        Offspring = Initial_pop(Offspring,length(Offspring),nVar,tx);%
        %         Offspring = Variation(Offspring,nC,n,newRc,tx);%变异操作
        [Parent,Mother] = Up_pop_Mo2(Mother,Parent,Offspring,nPop,I2_t2,I2_t2_r,newRc,lx2,ly2,n,template,Alg2);%0.5s 更新种群,或者说更新父辈
        
        Mother_last.x = 0;
        Mother_last.y = 0;
        Mother_last.z = 0;
        
        It_Mo = 0;%Mother迭代次数
        %flag_y1 = 0;
        while It_Mo<maxIt%如果Mother超过maxIt次没有找到更好的解，则认为达到最优
            if isempty(Mother(1).x) == 0%Mother被赋值
                Offspring = repmat(Parent, N_times, 1);
                Offspring(:) = Mother(end);
                Parent(:) = Mother(end);%因为目前看的话Mother中的元素最优，所以直接Mother中的元素给Parent和Offspring进行赋值
                Offspring = Variation(Offspring,nC,n,newRc,tx);%变异操作
                [Parent,Mother] = Up_pop_Mo2(Mother,Parent(1),Offspring,nPop,I2_t2,I2_t2_r,newRc,lx2,ly2,n,template,Alg2);%2s多 更新种群,或者说更新父辈
            else
                while isempty(Mother(end).x)%循环，直至Mother被赋值，即出现完全分离的情况
                    Offspring = repmat(Parent, N_times * 5, 1);
                    Offspring = Variation(Offspring,nC,n,newRc,tx);%变异操作
                    [Parent,Mother] = Up_pop_Mo2(Mother,Parent,Offspring,nPop,I2_t2,I2_t2_r,newRc,lx2,ly2,n,template,Alg2);%2s多 更新种群,或者说更新父辈
                end
            end
            
            if Mother(end).z ~= Mother_last.z%Mother有更新
                Mother_last = Mother(end);
                %                 I = Mother_last.tu;
                %                 imshow(I,'Border','tight');
                disp(['Mother中的元素数目：',num2str(numel(Mother))]);
                It_Mo = 0;
            else%如果Mother没有更新
                It_Mo = It_Mo+1;
                disp(['Mother未更新次数：',num2str(It_Mo)]);
            end
        end
        
        if isempty(Grandma(1).x) == 1 %如果Grandma是空的
            Grandma = Mother(end);
        else
            if Grandma(1).z<Mother(end).z
                Grandma = Mother(end);
            end
        end
    end
    
    if It == maxIt_va
        flag = 1;
    end
    
    I2 = Grandma.tu;
    %     I2 = trans2mid(I2,0);
    dx2 = Grandma.ddx;
    dy2 = Grandma.ddy;
    I1 = trans2lefttop(I1);
    I_base = zeros(ceil(l_y(I1) * (l_y(I2) / l_y(I1))) + 50, ceil(l_x(I1) * (l_x(I2) / l_x(I1)) + 50));
    I1 = paste_img(I_base,I1,1,1);
    [dx2,dy2,I2] = high_accuracy2(I1,alg2,key,dx2,dy2,1,ijmax2,1);
    %     imshow(I2);
    disp('奶奶出现,第二次寻优完毕！');%I2为第二次组合成功后的结果
end

ZP = 0;%左偏
YP = 0;%右偏

if dy2 >= 0 && dx2 >= 0
    ZP = 1;
elseif dy2 >= 0 && dx2<0
    YP = 1;
end

if dy2<0 && dx2 >= 0
    YP = 1;
elseif dy2<0 && dx2<0
    ZP = 1;
end

if ZP == 1
    dx2 = -abs(dx2);
    dy2 = -abs(dy2);
else
    dx2 = abs(dx2);
    dy2 = -abs(dy2);
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%以上完成了Y方向的排布，下面需要进行X方向的排布,这是第三次迭代
sita10 =sita1;

if abs(dx2) >= abs(dy2)
    beta = 90;
else
    beta = 0;
end

I2_t_double = trans2mid(I2_t_double,0);
I2_t_double = imrotate(I2_t_double,beta,'bicubic');
I2_t_double = trans2leftbottom(I2_t_double);

I_base = zeros(5*l_y(I2_t_double),4*l_x(I2_t_double));
% I_base = logical(I_base);
I_y5 = I2_t_double;
I_y5 = trans2lefttop(I_y5);
I_y5 = imcrop(I_y5,[0,0,l_x(I_y5),l_y(I_y5)]);
I_y5 = paste_img(I_base,I_y5,1,1);
I_y5 = trans2mid(I_y5,0);

if beta == 0
    
    dx2_t = dx2;
    dy2_t = dy2;
    
    if ZP == 1
        sita2 = -atan(abs(dx2_t/dy2_t));
        I_y5 = trans2rightbottom(I_y5);
    else
        sita2 = atan(abs(dx2_t/dy2_t));
        I_y5 = trans2leftbottom(I_y5);
    end
    
else
    
    dx2_t = -dy2;
    dy2_t =  dx2;
    if dy2_t > 0
        dx2_t = - dx2_t;
        dy2_t = - dy2_t;
    end
    
    if ZP == 0
        sita2 = -atan(abs(dx2_t/dy2_t));
        I_y5 = trans2rightbottom(I_y5);
    else
        sita2 = atan(abs(dx2_t/dy2_t));
        I_y5 = trans2leftbottom(I_y5);
    end
    
end

sita0 = atan(-dx2_t/dy2_t)*180/pi;%最终结果180度数制

I_y0 = I_y5;
for ii = 1:1
    I_y5 = I_y5 | imtranslate(I_y0,[ii*dx2_t, ii*dy2_t]);
end

% I_y5 = trans2mid(I_y5,0);
LY_I_y5 = l_y(I_y5);
LX_I_y5 = l_x(I_y5);
I_y5_s = imrotate(I_y5,sita0,'bicubic');%I_y5_temp为将I_y5搬竖直了
lx_I_y5_s = l_x(I_y5_s);%搬竖直后横向尺寸
lx_I_y5 = lx_I_y5_s/cos(sita2);%倾斜后的横向尺寸（注意是同一个y值下的横向尺寸）

ls = round(LY_I_y5+ceil(abs(dy2_t)/1.5));
lh = round(LX_I_y5 + lx_I_y5 * 1.2);

l_ver_add = 2 * abs(dy2_t);
while l_ver_add <= (ls / 2 + 25)
    ii = ii + 1;
    l_ver_add = l_ver_add + abs(dy2_t);
    I_y5 = I_y5 | imtranslate(I_y0,[ii*dx2_t, ii*dy2_t]);
    ls = ls + abs(dy2_t);
    lh = lh + abs(dx2_t);
end
ii = ii + 2;

newRc = zeros(ls,lh);%横竖正确
I_y5_guang = I_y5;
I_y5_guang = trans2lefttop(I_y5_guang);
I_y5_guang = imcrop(I_y5_guang,[0,0,l_x(I_y5),l_y(I_y5)]);
I_y5 = paste_img(newRc,I_y5_guang,1,1);
I_y5 = trans2mid(I_y5,0);

lx3 = l_x(I_y5);
ly3 = l_y(I_y5);
I_y5 = logical(I_y5);
I_y5_r = logical(trans2mid(imrotate(I_y5,180,'crop','bicubic'),0));
% imshow(I_y5,'Border','tight')

Parent = Parent0;
flag = 0;
% N_times = N_times + 1;
while flag == 0
    Grandma = repmat(template, 1, 1);
    for It = 1 : maxIt_va %迭代次数
        disp(['外层迭代次数：',num2str(It)]);
        Mother = template;
        %         Parent = Initial_pop(Parent,nPop,nVar,0);%2s多   生成子代种群（包含交叉操作）
        Offspring = repmat(Parent, N_times, 1);
        Parent = Initial_pop(Parent,nPop,nVar,0);%2s多   生成子代种群（包含交叉操作）
        Offspring = Initial_pop(Offspring,length(Offspring),nVar,tx);%
        %         Offspring = Variation(Offspring,nC,n,newRc,tx);%变异操作
        [Parent,Mother] = Up_pop_Mo3(Mother,Parent,Offspring,nPop,I_y5,I_y5_r,newRc,lx3,ly3,n,template,sita0);%0.5s 更新种群,或者说更新父辈
        
        Mother_last.x = 0;
        Mother_last.y = 0;
        Mother_last.z = 0;
        
        It_Mo = 0;%Mother迭代次数
        %flag_y1 = 0;
        while It_Mo<maxIt%如果Mother超过maxIt次没有找到更好的解，则认为达到最优
            if isempty(Mother(1).x) == 0%Mother被赋值
                Offspring = repmat(Parent, N_times, 1);
                Offspring(:) = Mother(end);
                Parent(:) = Mother(end);%因为目前看的话Mother中的元素最优，所以直接Mother中的元素给Parent和Offspring进行赋值
                Offspring = Variation(Offspring, nC, n, newRc,tx);%变异操作
                [Parent,Mother] = Up_pop_Mo3(Mother,Parent(1),Offspring,nPop,I_y5,I_y5_r,newRc,lx3,ly3,n,template,sita0);%2s多 更新种群,或者说更新父辈
            else
                while isempty(Mother(1).x)%循环，直至Mother被赋值，即出现完全分离的情况
                    Offspring = repmat(Parent, N_times * 5, 1);
                    Offspring = Variation(Offspring,nC,n,newRc,tx);%变异操作
                    [Parent,Mother] = Up_pop_Mo3(Mother,Parent,Offspring,nPop,I_y5,I_y5_r,newRc,lx3,ly3,n,template,sita0);%2s多 更新种群,或者说更新父辈
                end
            end
            
            if Mother(end).z ~= Mother_last.z%Mother有更新
                Mother_last = Mother(end);
                %                 I = Mother_last.tu;
                %                 imshow(I,'Border','tight');
                disp(['Mother中的元素数目：',num2str(numel(Mother))]);
                It_Mo = 0;
            else%如果Mother没有更新
                It_Mo = It_Mo+1;
                disp(['Mother未更新次数：',num2str(It_Mo)]);
            end
        end
        
        if isempty(Grandma(1).x) == 1 %如果Grandma是空的
            Grandma = Mother(end);
        else
            if Grandma(1).z<Mother(end).z
                Grandma = Mother(end);
            end
        end
    end
    
    if It == maxIt_va
        flag = 1;
    end
    
    I3 = Grandma.tu;
    dx3 = Grandma.ddx;
    dy3 = Grandma.ddy;
    
    I_y5 = trans2lefttop(I_y5);
    I_base = zeros(ceil(l_y(I_y5) * (l_y(I3) / l_y(I_y5))) + 50, ceil(l_x(I_y5) * (l_x(I3) / l_x(I_y5)) + 50));
    I_y5 = paste_img(I_base,I_y5,1,1);
    [dx3,dy3,I3] = high_accuracy3(I_y5,dx3,dy3,sita0,1,ijmax3,1);
    
    disp('奶奶出现,第三次寻优完毕！');%第三次寻优结果
end

SP = 0;%上偏
XP = 0;%下偏

if dx3 >= 0 && dy3 >= 0
    XP = 1;
elseif dx3 >= 0 && dy3<0
    SP = 1;
end

if dx3<0 && dy3 >= 0
    SP = 1;
elseif dx3<0 && dy3<0
    XP = 1;
end

if SP == 1
    dx3 = abs(dx3);
    dy3 = -abs(dy3);
else
    dx3 = abs(dx3);
    dy3 = abs(dy3);
end

%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%以上完成了第三次迭代
dx1 = round(dx1 / RF);
dy1 = round(dy1 / RF);
dx2 = round(dx2 / RF);
dy2 = round(dy2 / RF);
dx3 = round(dx3 / RF);
dy3 = round(dy3 / RF);

I1 = imresize(I1,1/RF);
I1 = tuomao(I1);
I2 = imresize(I2,1/RF);
I2 = tuomao(I2);
I3 = imresize(I3,1/RF);
I3 = tuomao(I3);
end
