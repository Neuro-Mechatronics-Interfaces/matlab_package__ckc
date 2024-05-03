function xp = approxPrctile(x,p)
x = sort(x);
L = length(x);
xp=[];
for k=1:length(p)
    xp(k)=x(round(p(k)/100*L));
end