function [p2] = histogram2(ahistogram_head_image2,I2)
%HISTOGRAM Summary of this function goes here

%   Detailed explanation goes here
        I1_size=size(I2); 
        I2=double(I2);
        
        Ia(:,:,1)=I2(:,:,1)./(I2(:,:,1)+I2(:,:,2)+I2(:,:,3));
        Ia(:,:,2)=I2(:,:,2)./(I2(:,:,1)+I2(:,:,2)+I2(:,:,3));
        Ia(:,:,3)=I2(:,:,3)./(I2(:,:,1)+I2(:,:,2)+I2(:,:,3));
        
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
        p2=0;
for i=1:8
    for j=1:8
       for k=1:8
            p2=p2+(histogram(i,j,k)*ahistogram_head_image2(i,j,k))^(1/2);
       end
    end
end


