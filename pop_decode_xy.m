function [dx, dy] = pop_decode_xy(bi,Nx,Ny,lx,ly)
%该函数将二进制编码转换为十进制编码
    dx=(Nx-lx)*(bin2dec(num2str(bi(:,1:12)))/4095);%范围为0~(Nx-sw)
    dy=(Ny-ly)*(bin2dec(num2str(bi(:,13:24)))/4095);%范围为0~(Ny-sw)
    %     sita=bin2dec(num2str(bi(:,25)));

    dx=round((Nx-lx)/2-dx);%dx的取值范围
    dy=round((Ny-ly)/2-dy);%dy的取值范围
    %     sita=sita*180;
end