function X=panduan(A,p)
B=sum(sum(A))/16;
s=0;
for i1=1:4
    for j1=1:4
        s=s+(B-A(i1,j1))^2;
    end
end
if(s<p)
    X=1;
else
    X=0;
end
end