function [Parent,Mother] = Up_pop_Mo1(Mother,Parent,Offspring,nPop,I2_t,I2_t_r,newRc,lx,ly,n,template,Alg1)
%UP_POP �˴���ʾ�йش˺�����ժҪ
%   ǰ�ڵĵ������
Nw=sum(I2_t(:))*n;
I=newRc;
I=logical(I);
Nx=numel(I(1,:));
Ny=numel(I(:,1));

dx=zeros(2,1);
dy=zeros(2,1);
afa=zeros(2,1);

Nm=numel(Mother);%NmΪMother��Ԫ�صĸ���
newPop = [Parent; Offspring];

N_y1=0;%newPop��yֵΪ1��Ԫ�صĸ���
for i=1:numel(newPop)%���ͼԪ��ȫ�ֿ�������㷽����򣬽�������0
    
    [dx(1), dy(1), afa(1)] = pop_decode_xya(newPop(i).x(1:25),Nx,Ny,lx,ly);
    if afa(1) == 0
        I11=I2_t;%���﷭ת���Ǽ��к��ͼ��ߴ�
    else
        I11=I2_t_r;
    end
    I11=imtranslate(I11,[dx(1), dy(1)]);
    
    [dx(2), dy(2), afa(2)] = pop_decode_xya(newPop(i).x(26:50),Nx,Ny,lx,ly);
    if afa(2) == 0
        I22=I2_t;%���﷭ת���Ǽ��к��ͼ��ߴ�
    else
        I22=I2_t_r;
    end
    I22=imtranslate(I22,[dx(2), dy(2)]);
    
    I=I11 | I22;
    I_lt=trans2lefttop(I);
    Nw_imerged=sum(I(:));
    
    if Nw_imerged==Nw %˵��ͼ����ȫ����
        %         imshow(I)
        N_y1=N_y1+1;
        newPop(i).y=1;
        
        if Alg1==1
            newPop(i).z=(-1)*fun_XY(I);
        elseif Alg1==2
            newPop(i).z=(-1)*fun_YX(I);
        elseif Alg1==3
            newPop(i).z=(-1)*fun_bs2(I,I11)+(-1)*fun_Variance(I)/100000000;
        elseif Alg1==4
            newPop(i).z=(-1)*fun_Variance(I)+(-1)*fun_bs2(I,I11)/100000000;
        end
        
        %         if Alg1==0
        %             newPop(i).z=(-1)*fun_bs(I);
        %         elseif Alg1==1
        %             newPop(i).z=(-1)*fun_bao2(I,I2_t)+(-1)*fun_Variance(I)/1000;
        %         elseif Alg1==2
        %             newPop(i).z=(-1)*fun_YX(I);
        %         elseif Alg1==3
        %             newPop(i).z=(-1)*fun_XY(I);
        %         elseif Alg1==4
        %             newPop(i).z=(-1)*fun_Variance(I);
        %         elseif Alg1==5
        %             newPop(i).z=(-1)*fun_S(I);
        %         else
        %             newPop(i).z=(-1)*fun_Lxy(I);
        %         end
        
        newPop(i).tu=I;
        newPop(i).sita=afa(1);
        
        dx=abs(l_x(I)-l_x(I2_t));
        dy=abs(l_y(I)-l_y(I2_t));
        
        flag=0;
        
        I11_t=trans2leftbottom(I11);
        I22_t=trans2leftbottom(I22);
        I22_t=imtranslate(I22_t,[dx,-dy]);
        I_t=I11_t | I22_t;
        I_t=trans2lefttop(I_t);
        I_compare=(I_lt ~= I_t);
        N_compare=sum(I_compare(:));
        
        if N_compare==0
            flag=1;
            dx=abs(dx);
            dy=-abs(dy);
        end
        
        if flag==0
            I11_t=trans2rightbottom(I11);
            I22_t=trans2rightbottom(I22);
            I22_t=imtranslate(I22_t,[-dx,-dy]);
            I_t=I11_t | I22_t;
            I_t=trans2lefttop(I_t);
            I_compare=(I_lt ~= I_t);
            N_compare=sum(I_compare(:));
            if N_compare==0
                flag=1;
                dx=-abs(dx);
                dy=-abs(dy);
            end
        end
        
        if flag==0
            I11_t=trans2righttop(I11);
            I22_t=trans2righttop(I22);
            I22_t=imtranslate(I22_t,[-dx,dy]);
            I_t=I11_t | I22_t;
            I_t=trans2lefttop(I_t);
            I_compare=(I_lt ~= I_t);
            N_compare=sum(I_compare(:));
            if N_compare==0
                flag=1;
                dx=-abs(dx);
                dy=abs(dy);
            end
        end
        
        if flag==0
            I11_t=trans2lefttop(I11);
            I22_t=trans2lefttop(I22);
            I22_t=imtranslate(I22_t,[dx,dy]);
            I_t=I11_t | I22_t;
            I_t=trans2lefttop(I_t);
            I_compare=(I_lt ~= I_t);
            N_compare=sum(I_compare(:));
            if N_compare==0
                dx=abs(dx);
                dy=abs(dy);
            end
        end
        
        newPop(i).ddx=dx;
        newPop(i).ddy=dy;
        
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