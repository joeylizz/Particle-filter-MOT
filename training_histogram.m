%==========================================================================
clc;clear all;
close all;
%==========================================================================
image=imread('115.jpg');
imshow(image)
a=imcrop(image,round([265-40/2 90-40/2 2*0.5*40 2.4*0.5*40]));
b=double(a);
s=size(a);
b1(:,:,1)=b(:,:,1)./(b(:,:,1)+b(:,:,2)+b(:,:,3));
b1(:,:,2)=b(:,:,2)./(b(:,:,1)+b(:,:,2)+b(:,:,3));
b1(:,:,3)=b(:,:,3)./(b(:,:,1)+b(:,:,2)+b(:,:,3));
ahistogram=zeros(8,8,8);
a=(s(1)^2+s(2)^2)^(0.5);
for i=1:s(1)
     for j=1:s(2)
        f_index=min(fix(b1(i,j,1)*8)+1,8);
         s_index=min(fix(b1(i,j,2)*8)+1,8);
          t_index=min(fix(b1(i,j,3)*8)+1,8);
          r=(((i-(s(1)/2))^2+(j-(s(2)/2))^2)^(1/2))/a;
       ahistogram(f_index,s_index,t_index)= ahistogram(f_index,s_index,t_index)+1-r^2;
    end
end
 ahistogram3=ahistogram/sum(sum(sum(ahistogram)));
save ahistogram3.mat ahistogram3 
% % ==========================================================================
