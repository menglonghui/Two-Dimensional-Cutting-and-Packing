function [I,ddx,ddy,sita] = Nw_cal_xya(newRc,Mother,I2_t,lx,ly)
%NW_CAL 此处显示有关此函数的摘要
%   此处显示详细说明
    nl=zeros(2,1);
    nt=zeros(2,1);
    %     I=newRc;
    Nx=numel(newRc(1,:));
    Ny=numel(newRc(:,1));
    dx=zeros(2,1);
    dy=zeros(2,1);
    sita=zeros(2,1);

    j=1;
    [dx(j), dy(j), sita(j)] = pop_decode_xya(Mother(end).x((1+25*(j-1)):(25*j)),Nx,Ny,lx,ly);
    I11=imrotate(I2_t,sita(j),'crop','bicubic');%其翻转的是原始图像大小
    % I11=trans2mid(I11);
    I11=imtranslate(I11,[dx(j), dy(j)]);

    j=2;
    [dx(j), dy(j), sita(j)] = pop_decode_xya(Mother(end).x((1+25*(j-1)):(25*j)),Nx,Ny,lx,ly);
    I22=imrotate(I2_t,sita(j),'crop','bicubic');%其翻转的是原始图像大小
    % I22=trans2mid(I22);
    I22=imtranslate(I22,[dx(j), dy(j)]);

    I=I11 | I22;

    for i=1:numel(I11(1,:)) %newRec的初始值为I1
        if sum(I11(:,i))==0
            nl(1)=nl(1)+1;
        else
            break
        end
    end

    for i=1:numel(I22(1,:)) %newRec的初始值为I1
        if sum(I22(:,i))==0
            nl(2)=nl(2)+1;
        else
            break
        end
    end

    for i=1:numel(I11(:,1)) %newRec的初始值为I1
        if sum(I11(i,:))==0
            nt(1)=nt(1)+1;
        else
            break
        end
    end

    for i=1:numel(I22(:,1)) %newRec的初始值为I1
        if sum(I22(i,:))==0
            nt(2)=nt(2)+1;
        else
            break
        end
    end

    ddx=nl(2)-nl(1);
    ddy=-(nt(2)-nt(1));
    sita=sita(1);

end

