function [I1,dx1,dy1,sita1] = find_I1_t(I2_t,I2_t_max_or,Parent0,alg1,template,key,maxIt_va,N_times,nPop,nVar,nC,n,maxIt)
%FIND 此处显示有关此函数的摘要
%   此处显示详细说明
% I0=I2_t;
lx1 = l_x(I2_t);
ly1 = l_y(I2_t);
% W0 = min(numel(I2_t(:,1)),numel(I2_t(1,:)));%空白板的宽度
% L0 = max(numel(I2_t(:,1)),numel(I2_t(1,:)));%空白板的长度
Lmax = max(l_x(I2_t),l_y(I2_t));
tx = 0;
ijmax1 = 6;
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
% W0 = ceil(W0 * RF);
% L0 = ceil(L0 * RF);

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
    [dx1,dy1,~] = high_accuracy1(I2_t,alg1,key,dx1,dy1,sita1,1,0,ijmax1,1);
    %         imshow(I);
    disp('奶奶出现,第一次寻优完毕！');%得到目前最好的排布
    %     I1 = I;%I1为第一次组合成功后的结果
    
end

% I2_t_temp = zeros(W0,L0);%W0和L0为空白板的尺寸
%
% for i = 1:numel(I1(:,1))
%     for j = 1:numel(I1(1,:))
%         I2_t_temp(i,j) = I1(i,j);
%     end
% end
%
% I2_t_temp = logical(I2_t_temp);
%
% I2_t = I2_t_temp;%这里的I2_t是已经排好的两个图元的组合
% I2_t = trans2leftbottom(I2_t);%将其移动到左下角
% I2_t_double = I2_t;%I2_t_double为两个69式排好的图元，为后期的X方向的寻优做准备。
%
% I2_t2 = I2_t;
% I2_t2 = trans2mid(I2_t2,0);

dx1 = round(dx1 / RF);
dy1 = round(dy1 / RF);

% I1 = imresize(I1,1/RF);
% I1 = tuomao(I1);

[dx1,dy1,I1] = high_accuracy1(I2_t_max_or,alg1,key,dx1,dy1,sita1,1,0,ijmax1,1);
end

