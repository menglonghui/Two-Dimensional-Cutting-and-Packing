function [dx, dy, sita] = pop_decode_xya_mix(bi,Nx,Ny,lx,ly)
%该函数将二进制编码转换为十进制编码
    dx = (Nx - lx) * (bin2dec(num2str(bi(01 : 12))) / 4095);%范围为0~(Nx-sw)
    dy = (Ny - ly) * (bin2dec(num2str(bi(15 : 26))) / 4095);%范围为0~(Ny-sw)
    sita = bin2dec(num2str(bi(13:14)));

    dx = round(1.01 * ((Nx - lx) / 2 - dx));%dx的取值范围
    dy = round(1.01 * ((Ny - ly) / 2 - dy));%dy的取值范围
    sita = sita * 90;
end
