function [I,ddx,ddy] = Nw_cal_xy(newRc,Mother,I2_t,~,lx,ly)
%NW_CAL �˴���ʾ�йش˺�����ժҪ
%   �˴���ʾ��ϸ˵��
    nl=zeros(2,1);
    nt=zeros(2,1);

%     I=newRc;
%     I=logical(I);
    Nx=numel(newRc(1,:));
    Ny=numel(newRc(:,1));
    dx=zeros(2,1);
    dy=zeros(2,1);

    j=1;
    [dx(j), dy(j)] = pop_decode_xy(Mother(end).x((1+25*(j-1)):(25*j)),Nx,Ny,lx,ly);
    I11=imtranslate(I2_t,[dx(j), dy(j)]);

    j=2;
    [dx(j), dy(j)] = pop_decode_xy(Mother(end).x((1+25*(j-1)):(25*j)),Nx,Ny,lx,ly);
    I22=imtranslate(I2_t,[dx(j), dy(j)]);

    I=I11 | I22;

    for i=1:numel(I11(1,:)) %newRec�ĳ�ʼֵΪI1
        if sum(I11(:,i))==0
            nl(1)=nl(1)+1;
        else
            break
        end
    end

    for i=1:numel(I22(1,:)) %newRec�ĳ�ʼֵΪI1
        if sum(I22(:,i))==0
            nl(2)=nl(2)+1;
        else
            break
        end
    end

    for i=1:numel(I11(:,1)) %newRec�ĳ�ʼֵΪI1
        if sum(I11(i,:))==0
            nt(1)=nt(1)+1;
        else
            break
        end
    end

    for i=1:numel(I22(:,1)) %newRec�ĳ�ʼֵΪI1
        if sum(I22(i,:))==0
            nt(2)=nt(2)+1;
        else
            break
        end
    end

    ddx=nl(2)-nl(1);
    ddy=-(nt(2)-nt(1));

end

