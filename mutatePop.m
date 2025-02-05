function p = mutatePop(x,n,L_nR)  %变异
max_end=max(2,18-floor(L_nR/220));
min=1;
max_value=n*randi([1,max_end]);
byws=randi([min, max_value])-randi([0,1]);%变异位数
m = numel(x);
for j=1:byws
    s = randi([1, m]);
    x(s)=1-x(s);
end
p = x;
end