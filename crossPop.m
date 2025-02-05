function [y1, y2] = crossPop(x1, x2) %交叉操作

    n = numel(x1);%n = numel(A) 返回数组A中元素个数。
    s = randi([2, n - 1]);% r = randi([imin,imax],...)  返回一个在[imin,imax]范围内的伪随机整数

    y1 = [x1(1 : s) x2(s + 1: end)];
    y2 = [x2(1 : s) x1(s + 1: end)];

end