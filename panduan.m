function X=panduan(A,p)
B=sum(sum(A))/9;
s=0;
for i1=1:3
    for j1=1:3
        s=s+(B-A(i1,j1))^2;
    end
end
if(s<p)
    X=1;
else
    X=0;
end
end