 function [B1,X]=tuxiangcj(X,apha,p)%X为待处理图片，apha为放大倍数，先假设为2,区分平坦块系数
%预处理，先将X规范化，能进行4*4的块运算
X=rgb2gray(X);
X=im2double(X);
[x,y]=size(X);
A=X;
if(mod(x,2)~=0)
    A=[A;A(x,:)];
end
if(mod(y,2)~=0)
    A=[A,A(:,y)];
end
%对图像进行放大
[x1,y1]=size(A);
x2=x1*apha;y2=y1*apha;
B=zeros(x2,y2);%生成插值后矩阵模板
for i=apha:apha:x2
    for j=apha:apha:y2
        B(i,j)=A(i/apha,j/apha);
    end
end
%将原图像的值对应赋予新图像的位置
B1=B;
%对图像进行纹理块或平坦块判断;
[xc,yc]=meshgrid(linspace(0.1,0.4,4),linspace(0.1,0.4,4));
[xxc,yyc]=meshgrid(linspace(0.05,0.4,8),linspace(0.05,0.4,8));
for i=1:2:x1-2
    for j=1:2:y1-2
        L=[A(i,j),A(i,j+1),A(i,j+2),A(i,j+3);
            A(i+1,j),A(i+1,j+1),A(i+1,j+2),A(i+1,j+3);
            A(i+2,j),A(i+2,j+1),A(i+2,j+2),A(i+2,j+3);
            A(i+3,j),A(i+3,j+1),A(i+3,j+2),A(i+3,j+3)];
        x=panduan(L,p);
        if(x==1)%对平坦块使用牛顿牛顿法
            
           L1=netnetwon(xc,yc,L,xxc,yyc);%X,Y,Z为插值点坐标，T为待插值点坐标横纵坐标，P为T的插值结果
        else
           L1=NT(xc,yc,L,xxc,yyc);%X,Y,Z为插值点坐标，T为待插值点坐标横纵坐标，P为T的插值结果  
        end
        if(i==1&&j==1)
        B1(2*i-1:2*i+6,2*j-1:2*j+6)=L1;
        end
        if(i==1&&j~=1)
        L2=(L1(1:8,1:4)+B1(2*i-1:2*i+6,2*j-1:2*j+2))./2;
        L1=[L2,L1(1:8,5:8)];
        B1(2*i-1:2*i+6,2*j-1:2*j+6)=L1;
        end 
        if(i~=1&&j==1)
        L2=(L1(1:4,1:8)+B1(2*i-1:2*i+2,2*j-1:2*j+6))./2;
        L1=[L2;L1(5:8,1:8)];
        B1(2*i-1:2*i+6,2*j-1:2*j+6)=L1;
        end
        if(i~=1&&j~=1)
        L2=(L1(1:4,1:8)+B1(2*i-1:2*i+2,2*j-1:2*j+6))./2;
        L3=(L1(1:8,1:4)+B1(2*i-1:2*i+6,2*j-1:2*j+2))./2;
        L4=(L2(1:4,1:4)+L3(1:4,1:4))./2;
        L1=[L4,L2(1:4,5:8);L3(5:8,1:4),L1(5:8,5:8)];
        B1(2*i-1:2*i+6,2*j-1:2*j+6)=L1;
        end
    end
end
    B1=im2uint8(B1);
  imwrite(B1,'E:\test\new002.bmp');
  imwrite(X,'E:\test\001.bmp');
end