function [x_ave,y_ave] = cor_average(I)
I = logical(I);
x_sum = 0;
y_sum = 0;

N = 0;
for i=1:numel(I(:,1))
    for j=1:numel(I(1,:))
        if I(i,j) == 1
            x_sum = x_sum + j;
            y_sum = y_sum + i;
            N = N + 1;
        end
    end
end

x_ave = x_sum / N;
y_ave = numel(I(:,1)) - y_sum / N;
