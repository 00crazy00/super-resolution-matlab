 function A=tuxiangchazhi1(Y,apha)%进行图像的一维插值运算，基于thiele和netwon两种算法。X为灰色图像矩阵图像，apha为缩放系数，暂时假设aoha为2;
[x,y]=size(Y);
X=im2double(Y);
x1=x*apha;y1=y*apha;
A=zeros(x1,y1);
for i=apha:apha:x1
    for j=apha:apha:y1
        A(i,j)=X(i/apha,j/apha);
    end
end%将原图像的值对应赋予新图像的位置
a=floor(y1/6);a1=mod(y1,6);
b=round(x1/6);b1=mod(x1,6);
for i=apha:apha:x1%先对行进行变化，使用newton方法,采取每八个点为一段，若段内变化率大于30%，直接赋值相邻两点中值
    k=0;
    for j=1:6:y1
       k=k+1;
       if (k==a+1)%检验图像的每段中是否为8个点，每推进一段，重合两点
           if(a1~=0)   
          for l=j:j+a1-1
              A(i,l)=A(i,j);
          end 
           end
       elseif(k==1)
               if(bianhualv(A(i,apha),A(i,2*apha))&&bianhualv(A(i,2*apha),A(i,3*apha))&&bianhualv(A(i,3*apha),A(i,4*apha)))%如果变化率小于30%，进行插值
                   s=thiele([apha,2*apha,3*apha,4*apha],[A(i,apha),A(i,2*apha),A(i,3*apha),A(i,4*apha)],[1,3,5,7]);
                   A(i,1)=s(1,1);A(i,3)=s(2,1);A(i,5)=s(3,1);A(i,7)=s(4,1);
               else
                   A(i,1)=A(i,2);A(i,3)=(A(i,2)+A(i,4))/2;A(i,5)=(A(i,4)+A(i,6))/2;A(i,7)=(A(i,6)+A(i,8))/2;
               end
       elseif(k==a)
               break;
       else
               if(bianhualv(A(i,j),A(i,j+1))&&bianhualv(A(i,j+1),A(i,j+3))&&bianhualv(A(i,j+3),A(i,j+5))&&bianhualv(A(i,j+5),A(i,j+7)))%如果变化率小于30%，进行插值
                   s=thiele([j,j+1,j+3,j+5,j+7],[A(i,j),A(i,j+1),A(i,j+3),A(i,j+5),A(i,j+7)],[j+2,j+4,j+6]);
                   A(i,j+2)=s(1,1);A(i,j+4)=s(2,1);A(i,j+6)=s(3,1);
               else
                   A(i,j+2)=(A(i,j+1)+A(i,j+3))/2;A(i,j+4)=(A(i,j+3)+A(i,j+5))/2;A(i,j+6)=(A(i,j+5)+A(i,j+7))/2;
               end
        end
     end
end
for j=apha:apha:y1%先对行进行变化，使用newton方法,采取每八个点为一段，若段内变化率大于30%，直接赋值相邻两点中值
    k=0;
    for i=1:6:x1
       k=k+1;
       if (k==b+1)%检验图像的每段中是否为8个点，每推进一段，重合两点
           if(b1~=0)   
          for l=j:j+b1-1
              A(i,j)=A(i+l,j);
          end 
           end
       elseif(k==1)
               if(bianhualv(A(apha,i),A(2*apha,i))&&bianhualv(A(2*apha,i),A(3*apha,i))&&bianhualv(A(3*apha,i),A(4*apha,i)))%如果变化率小于30%，进行插值
                   s=thiele([apha,2*apha,3*apha,4*apha],[A(apha,i),A(2*apha,i),A(3*apha,i),A(4*apha,i)],[1,3,5,7]);
                   A(1,j)=s(1,1);A(3,j)=s(2,1);A(5,j)=s(3,1);A(7,j)=s(4,1);
               else
                   A(1,j)=A(2,j);A(3,j)=(A(2,j)+A(4,j))/2;A(5,j)=(A(4,j)+A(6,j))/2;A(7,j)=(A(6,j)+A(8,j))/2;
               end
       elseif(k==a)
               break;
       else
               if(bianhualv(A(i,j),A(i+1,j))&&bianhualv(A(i+1,j),A(i+3,j))&&bianhualv(A(i+3,j),A(i+5,j))&&bianhualv(A(i+5,j),A(i+7,j)))%如果变化率小于30%，进行插值
                   s=thiele([i,i+1,i+3,i+5,i+7],[A(i,j),A(i+1,j),A(i+3,j),A(i+5,j),A(i+7,j)],[i+2,i+4,i+6]);
                   A(i+2,j)=s(1,1);A(i+4,j)=s(2,1);A(i+6,j)=s(3,1);
               else
                   A(i+2,j)=(A(i+1,j)+A(i+3,j))/2;A(j+4,i)=(A(i+3,j)+A(i+5,j))/2;A(i+6,j)=(A(i+5,j)+A(i+7,j))/2;
               end
       end
   end
end
A=im2uint8(A);
end

function d=bianhualv(a,b)%计算相邻像素变化率
c=abs(a-b);
if(c/a<=0.3)
    d=1;
else
    d=0;
end
end
000000000000000000