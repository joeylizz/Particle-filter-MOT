function [ahistogram_head_image1] = Update1( ahistogram_head_image1,I1)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here
        I1_size=size(I1); 
        I1=double(I1);
        beta=0.99;
        
        Ia(:,:,1)=I1(:,:,1)./(I1(:,:,1)+I1(:,:,2)+I1(:,:,3));
        Ia(:,:,2)=I1(:,:,2)./(I1(:,:,1)+I1(:,:,2)+I1(:,:,3));
        Ia(:,:,3)=I1(:,:,3)./(I1(:,:,1)+I1(:,:,2)+I1(:,:,3));
        
        a=sqrt(I1_size(1)^2+I1_size(2)^2);
        histogram=zeros(8,8,8);
            for k2=1:I1_size(1)
                for k3=1:I1_size(2)
                     r_index=min(fix(Ia(k2,k3,1)*8)+1,8);
                     g_index=min(fix(Ia(k2,k3,2)*8)+1,8);
                     b_index=min(fix(Ia(k2,k3,3)*8)+1,8);
                     r=(((k2-I1_size(1)/2)^2+(k3-I1_size(2)/2)^2)^(1/2))/a;
                     histogram(r_index,g_index,b_index)= histogram(r_index,g_index,b_index)+1-r^2;
                end
            end
         histogram=histogram/sum(sum(sum(histogram)));  
        ahistogram_head_image1=beta*ahistogram_head_image1+(1-beta)*histogram;
end

