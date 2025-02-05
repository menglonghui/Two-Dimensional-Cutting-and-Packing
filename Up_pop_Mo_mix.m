function [Parent,Mother,z] = Up_pop_Mo_mix(Mother,Parent,Offspring,nPop,I_t,I_o,template)
%UP_POP_MO_MIX �˴���ʾ�йش˺�����ժҪ
%   �˴���ʾ��ϸ˵��
I_t = imresize(I_t, 0.25);
I_o = imresize(I_o, 0.25);
% Nw = sum(I_t(:)) + sum(I_o(:));
Nx=numel(I_o(1,:));
Ny=numel(I_o(:,1));
Nm=numel(Mother);%NmΪMother��Ԫ�صĸ���
% I_C=I_o;
I_t = logical(I_t);
I_t = trans2mid(I_t,0);
I_t_0 = I_t;
I_t_90 = imrotate(I_t,90,'bicubic','crop');
I_t_180 = imrotate(I_t,180,'bicubic','crop');
I_t_270 = imrotate(I_t,270,'bicubic','crop');

Nw = sum(I_t(:)) + sum(I_o(:));


lx = min(l_x(I_t),l_y(I_t));
ly = lx;

newPop = [Parent; Offspring];
N_y1=0;%newPop��yֵΪ1��Ԫ�صĸ���

for i=1:numel(newPop)%���ͼԪ��ȫ�ֿ�������㷽����򣬽�������0
    
    [dx, dy, afa] = pop_decode_xya_mix(newPop(i).x,Nx,Ny,lx,ly);
    
    if afa == 0
        I_t1 = I_t_0;
    elseif afa == 90
        I_t1 = I_t_90;
    elseif afa == 180
        I_t1 = I_t_180;
    elseif afa == 270
        I_t1 = I_t_270;
    end
    
    I_t1 = imtranslate(I_t1,[dx, dy]);
    
    I = I_o | I_t1;
    Nw_imerged = sum(I(:));
    
    if Nw_imerged == Nw %˵��ͼ����ȫ����
        %         I = expand_mix(I,100);
        %                         imshow(I)
        %         I_C = I_C | I;
        %         pause(1)
        N_y1 = N_y1+1;
        newPop(i).y = 1;
        
        I_ex = expand_mix(I,100);
        %         imshow(I_ex)
        newPop(i).z = - sum(I_ex(:));
        
        newPop(i).tu=I;
        newPop(i).sita=afa;
        
        [nl,~,~,nb] = margin(I);
        
        dx = nl + l_x(I)/2;
        dy = nb + l_y(I)/2;
        
        newPop(i).ddx=dx;
        newPop(i).ddy=dy;
        
    else
        newPop(i).y=0;
        newPop(i).z=0;
    end
end

[~, so] = sort([newPop.y], 'descend'); %soΪ�������������
newPop = newPop(so);

newPop_t = repmat(template, 1, 1);
if N_y1 == 0
    ...
else
newPop_t = newPop(1:N_y1);

[~, so] = sort([newPop_t.z], 'descend'); %��newPop_t�е�Ԫ�ذ��շ����������
newPop_t = newPop_t(so);
end

newPop(1:N_y1) = newPop_t;
Parent = newPop(1 : nPop);%������Ⱥ������˵���¸���

if Parent(1).y == 1 && isempty(Mother(1).x)%ֻ��yֵΪ1��Ԫ�ظ���Mother
    Mother(1) = Parent(1);
elseif (Parent(1).y == 1) && (isempty(Mother(1).x) == 0) && (Parent(1).z > Mother(end).z)%���Mother����ֵ������ֱ�Ӹ���Mother����һ��Ԫ��
    Nm = Nm + 1;
    Mother(Nm) = Parent(1);
end
I1=Mother(end).tu;
I2=I_o;
I1=imresize(I1,4);
I2=imresize(I2,4);
% imshow(I1 | I2)
z = - sum(I_ex(:));
end

