function [Parent,Mother] = Up_pop_Mo3(Mother,Parent,Offspring,nPop,I2_t,I2_t_r,newRc,lx,ly,n,template,sita0)
%UP_POP �˴���ʾ�йش˺�����ժҪ
%   ���ڵĵ������
Nw=sum(I2_t(:))*n;
%I=newRc;
%I=logical(I);
Nx=numel(newRc(1,:));
Ny=numel(newRc(:,1));

dx=zeros(2,1);
dy=zeros(2,1);

Nm=numel(Mother);%NmΪMother��Ԫ�صĸ���
newPop = [Parent; Offspring];

N_y1=0;%newPop��yֵΪ1��Ԫ�صĸ���

% I2_t =logical(trans2mid(I2_t,1));
% I2_t_r =logical(imrotate(I2_t,180));
% I2_t_r =logical(trans2mid(I2_t_r,1));
for i=1:numel(newPop)%���ͼԪ��ȫ�ֿ�������㷽����򣬽�������0
    
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

[~, so] = sort([newPop.y], 'descend'); %soΪ�������������
newPop = newPop(so);

newPop_t=repmat(template, 1, 1);
if N_y1==0
    ...
else
newPop_t=newPop(1:N_y1);

[~, so]=sort([newPop_t.z], 'descend'); %��newPop_t�е�Ԫ�ذ��շ����������
newPop_t=newPop_t(so);
end

newPop(1:N_y1)=newPop_t;
Parent = newPop(1 : nPop);%������Ⱥ������˵���¸���

if Parent(1).y==1 && isempty(Mother(1).x)%ֻ��yֵΪ1��Ԫ�ظ���Mother
    Mother(1)=Parent(1);
elseif (Parent(1).y==1) && (isempty(Mother(1).x)==0) && (Parent(1).z > Mother(end).z)%���Mother����ֵ������ֱ�Ӹ���Mother����һ��Ԫ��
    Nm=Nm+1;
    Mother(Nm)=Parent(1);
end
end