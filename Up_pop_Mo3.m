function [Parent,Mother] = Up_pop_Mo3(Mother,Parent,Offspring,nPop,I2_t,I2_t_r,newRc,lx,ly,n,template,sita0)
%UP_POP 此处显示有关此函数的摘要
%   后期的迭代组合
Nw=sum(I2_t(:))*n;
%I=newRc;
%I=logical(I);
Nx=numel(newRc(1,:));
Ny=numel(newRc(:,1));

dx=zeros(2,1);
dy=zeros(2,1);

Nm=numel(Mother);%Nm为Mother中元素的个数
newPop = [Parent; Offspring];

N_y1=0;%newPop中y值为1的元素的个数

% I2_t =logical(trans2mid(I2_t,1));
% I2_t_r =logical(imrotate(I2_t,180));
% I2_t_r =logical(trans2mid(I2_t_r,1));
for i=1:numel(newPop)%如果图元完全分开，则计算方差，否则，将方差置0
    
    j=1;
    [dx(j), dy(j)] = pop_decode_xy(newPop(i).x((1+25*(j-1)):(25*j)),Nx,Ny,lx,ly);
    I11=imtranslate(I2_t,[dx(j), dy(j)]);
    
    j=2;
    [dx(j), dy(j)] = pop_decode_xy(newPop(i).x((1+25*(j-1)):(25*j)),Nx,Ny,lx,ly);
    I22=imtranslate(I2_t_r,[dx(j), dy(j)]);
    
    I=I11 | I22;
    Nw_imerged=sum(I(:));
    
    if Nw_imerged==Nw
        %         imshow(I)
        N_y1=N_y1+1;
        newPop(i).y=1;
        newPop(i).z=(-1)*fun_XY_ro(I,sita0);
        newPop(i).tu=I;
        newPop(i).ddx=dx(2)-dx(1);
        newPop(i).ddy=dy(2)-dy(1);
    else
        newPop(i).y=0;
        newPop(i).z=0;
    end
end

[~, so] = sort([newPop.y], 'descend'); %so为降序后的序号排序
newPop = newPop(so);

newPop_t=repmat(template, 1, 1);
if N_y1==0
    ...
else
newPop_t=newPop(1:N_y1);

[~, so]=sort([newPop_t.z], 'descend'); %将newPop_t中的元素按照方差进行排序
newPop_t=newPop_t(so);
end

newPop(1:N_y1)=newPop_t;
Parent = newPop(1 : nPop);%更新种群，或者说更新父代

if Parent(1).y==1 && isempty(Mother(1).x)%只将y值为1的元素赋给Mother
    Mother(1)=Parent(1);
elseif (Parent(1).y==1) && (isempty(Mother(1).x)==0) && (Parent(1).z > Mother(end).z)%如果Mother被赋值过，则直接赋给Mother的下一个元素
    Nm=Nm+1;
    Mother(Nm)=Parent(1);
end
end