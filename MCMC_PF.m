clc; clear all
close all
%==========================================================================
load multi_position;
load ahistogram_head_image1;
load ahistogram_head_image2;
%==========================================================================
%% 1.Initialisation--------------------------------------------------------
V = VideoWriter('newfile.avi');
open(V)
N_particle = 100;   % number of particles
N = 100;
N_pixel = 30; % number of pixel
theta = linspace(0,2*pi,30); % 
A = eye(4);    % state transition matrix
v_g=0.00005;   % variance of intensity gradient
v_c=0.00075;   % variance of colour histogram
s1=[268 96 17 22]';  %true position of first target
s2=[179 79 13 18]';  %true position of second target
alpha = 0.5; % Scaling parameter
B=30; % Burn period
p1=0 ; p2=0; p3=0; p4=0;
% prior as importance function
for n=1:N
particle1(:,n)= s1;
particle2(:,n)= s2;
end
%==========================================================================
%% 2.MCMC Particle filter and collect the gradient and histogram-----------
for i=115:225
    HEAD=imread([num2str(i),'.jpg']); % images read in
    I=double(rgb2gray(HEAD)); % obtain the grayscale intensity 
    s_particle1(:,1)=s1; % target 1 initial state
    s_particle2(:,1)=s2; % target 2 initial state
% 2.1 prediction
    for k=2:N_particle+B 
        s_particle1(:,k)= s_particle1(:,k-1) +[2.5*randn;1*randn;0.1*randn;0.1*randn]; % target 1
        s_particle2(:,k)= s_particle2(:,k-1) +[2.5*randn;1*randn;0.1*randn;0.1*randn]; % target 2
        % 2.2 obtain the intensity gradient
        for k1=1:N_pixel
            % target 1
            pixels1 = [s_particle1(1,k)+s_particle1(3,k)*cos(theta(k1)),s_particle1(2,k)+s_particle1(4,k)*sin(theta(k1))];
            a_p1=round(pixels1);
            gx1 = I(a_p1(2)-2,a_p1(1))-I(a_p1(2)+2,a_p1(1))+2*I(a_p1(2)-1,a_p1(1))-2*I(a_p1(2)+1,a_p1(1));% sobel opearator
            gy1 = I(a_p1(2),a_p1(1)-2)-I(a_p1(2),a_p1(1)+2)+2*I(a_p1(2),a_p1(1)-1)-2*I(a_p1(2),a_p1(1)+1);
            g1(k1)= sqrt(gx1^2+gy1^2); % magnitude
            % target 2
            pixels2 = [s_particle2(1,k)+s_particle2(3,k)*cos(theta(k1)),s_particle2(2,k)+s_particle2(4,k)*sin(theta(k1))];
            a_p2=round(pixels2);
            gx2 = I(a_p2(2)-2,a_p2(1))-I(a_p2(2)+2,a_p2(1))+2*I(a_p2(2)-1,a_p2(1))-2*I(a_p2(2)+1,a_p2(1));% sobel opearator
            gy2 = I(a_p2(2),a_p2(1)-2)-I(a_p2(2),a_p2(1)+2)+2*I(a_p2(2),a_p2(1)-1)-2*I(a_p2(2),a_p2(1)+1);
            g2(k1)= sqrt(gx2^2+gy2^2); % magnitude 
        end
        % target 1
        sum_g1(k)=sum(g1)/30; % normalization
        weight_g1(k)=(1/sqrt(2*pi*v_g))*exp(-(1/sum_g1(k))^2/(2*v_g)); % gradient weight
        I1=imcrop(HEAD,round([s_particle1(1,k)-s_particle1(3,k) s_particle1(2,k)-s_particle1(4,k) 2*s_particle1(3,k) 2*s_particle1(4,k)])); % target head   
        % target 2
        sum_g2(k)=sum(g2)/30; % normalization
        weight_g2(k)=(1/sqrt(2*pi*v_g))*exp(-(1/sum_g2(k))^2/(2*v_g)); % gradient weight
        I2=imcrop(HEAD,round([s_particle2(1,k)-s_particle2(3,k) s_particle2(2,k)-s_particle2(4,k) 2*s_particle2(3,k) 2*s_particle2(4,k)])); % target head
        % 2.3 obtain the histogram
        % target 1
        [p1] = histogram1(ahistogram_head_image1,I1); 
         d1  = 1-p1; % Bhattacharyya distance
        weight_c1(k)=(1/sqrt(2*pi*v_c))*exp(-d1/(2*v_c)); % weight on color histogram
        % target 2
        [p2] = histogram2(ahistogram_head_image2,I2);
         d2  = 1-p2; % Bhattacharyya distance
        weight_c2(k)=(1/sqrt(2*pi*v_c))*exp(-d2/(2*v_c)); % weight on color histogram
        % 2.4 Average weights
        weight_1(k) = weight_g1(k)* weight_c1(k); % Target 1
        weight_2(k) = weight_g2(k)* weight_c2(k); % Target 2
        % 2.5 Compute the acceptance ratio
        p1=0 ; p2=0; p3=0; p4=0;
        for n=1:N
            % target 1
            p1 = p1 + exp(-0.5*(s_particle1(:,k)-particle1(:,n))'*(s_particle1(:,k)-particle1(:,n)));
            p2 = p2 + exp(-0.5*(s_particle1(:,k-1)-particle1(:,n))'*(s_particle1(:,k-1)-particle1(:,n)));
            % target 2
            p3 = p3 + exp(-0.5*(s_particle2(:,k)-particle2(:,n))'*(s_particle2(:,k)-particle2(:,n)));
            p4 = p4 + exp(-0.5*(s_particle2(:,k-1)-particle2(:,n))'*(s_particle2(:,k-1)-particle2(:,n)));
        end
        a1=min(1,(weight_1(k)*p1/(weight_1(k-1)*p2)));
        if a1< rand
            s_particle1(:,k)= s_particle1(:,k-1); % acceptance ratio for target 1
        end
        a2=min(1,(weight_2(k)*p3/(weight_2(k-1)*p4)));
        if a2 < rand                              % acceptance ratio for target 2
            s_particle2(:,k)= s_particle2(:,k-1);
        end
    end   
    % 2.6 Discard the first B samples to form new particle sets
    % Target 1
    s_particle1 = s_particle1(:,B+1:N+B);
    particle1= s_particle1;
    s_estimated1= mean(s_particle1,2); % mean estimated position
    I1=imcrop(HEAD,round([s_estimated1(1)-s_estimated1(3) s_estimated1(2)-s_estimated1(4) 2*s_estimated1(3) 2*s_estimated1(4)]));
    [ahistogram_head_image1] = Update1(ahistogram_head_image1,I1); % update the target model (histogram)
    % Target 2
    s_particle2 = s_particle2(:,B+1:N+B);
    particle2= s_particle2;
    s_estimated2= mean(s_particle2,2); % mean estimated position
    I2=imcrop(HEAD,round([s_estimated2(1)-s_estimated2(3) s_estimated2(2)-s_estimated2(4) 2*s_estimated2(3) 2*s_estimated2(4)]));
    [ahistogram_head_image2] = Update2(ahistogram_head_image2,I2); % update the target model (histogram) 
    s1=s_estimated1; % reset the initial state
    s2=s_estimated2;
    % 2.7 Draw points to match the ellipse
    % target 1
    for k4=1:N_pixel
    pixels1 = [s_estimated1(1)+s_estimated1(3)*cos(theta(k4)),s_estimated1(2)+s_estimated1(4)*sin(theta(k4))]; 
    point_x1(1,k4)= pixels1(1);
    point_y1(1,k4)= pixels1(2);
    end
    % target 2
    for k4=1:N_pixel
    pixels2 = [s_estimated2(1)+s_estimated2(3)*cos(theta(k4)),s_estimated2(2)+s_estimated2(4)*sin(theta(k4))]; 
    point_x2(1,k4)= pixels2(1);
    point_y2(1,k4)= pixels2(2);
    end
    figure(1)   
    imshow(HEAD)
    hold on
    plot(point_x1,point_y1,'r')
    hold on
    plot(point_x2,point_y2,'r')
%   pause (0.1)
    % Record
    % target 1
    x1(1,i-114)=s_estimated1(1);
    y1(1,i-114)=s_estimated1(2); 
    % target 2
    x2(1,i-114)=s_estimated2(1);
    y2(1,i-114)=s_estimated2(2); 
    M(i-114) = getframe(gcf);
end
writeVideo(V,M);
%==========================================================================
%% 3.Euclidean Error-------------------------------------------------------
for t=1:111
    % Target 1
EUerror1(t)= sqrt((x1(t)-multi_position(1,t))^2+(y1(t)-multi_position(2,t))^2);
    % Target 2
EUerror2(t)= sqrt((x2(t)-multi_position(3,t))^2+(y1(t)-multi_position(4,t))^2);
end
% Target 1
MSE_error_x1 = 1/111*(x1-multi_position(1,:))*(x1-multi_position(1,:))'
MSE_error_y1 = 1/111*(y1-multi_position(2,:))*(y1-multi_position(2,:))'
% Target 2
MSE_error_x2 = 1/111*(x2-multi_position(3,:))*(x1-multi_position(3,:))'
MSE_error_y2 = 1/111*(y2-multi_position(4,:))*(y1-multi_position(4,:))'
%==========================================================================
%% 4.Figures---------------------------------------------------------------
figure(2) % compare the true and estimated x position 
plot(multi_position(1,:),'r');
hold on 
plot(x1,'b')
title('true position of x (target one) and MCMC estimate');
xlabel('Frames');ylabel('position x');
legend('True position x (target one)','MCMC');
%--------------------------------------------------------------------------
figure(3) % compare the true and estimated y position 
plot(multi_position(2,:),'r')
hold on
plot(y1,'b')
title('true position of y (target one) and MCMC estimate');
xlabel('Frames');ylabel('position y');
legend('True position y (target one)','MCMC');
%--------------------------------------------------------------------------
figure(4) % compare the true track and estimates 
plot(multi_position(1,:),multi_position(2,:),'r')
hold on 
plot(x1,y1,'b')
title('true track (target one) and MCMC estimate');
xlabel('x position of head');ylabel('y position of head');
legend('True track (target one)','MCMC');
%--------------------------------------------------------------------------
figure(5) % plot the error on x position
plot(EUerror1)
title('Euclidean error of target one');
xlabel('Frames');ylabel('magnitude(pixels)');
%--------------------------------------------------------------------------
figure(6) % compare the true and estimated x position 
plot(multi_position(3,:),'r');
hold on 
plot(x2,'b')
title('true position of x (target two) and MCMC estimate');
xlabel('Frames');ylabel('position x');
legend('True position x (target two)','MCMC');
%--------------------------------------------------------------------------
figure(7) % compare the true and estimated y position 
plot(multi_position(4,:),'r')
hold on
plot(y2,'b')
title('true position of y (target two) and MCMC estimate');
xlabel('Frames');ylabel('position y');
legend('True position y (target two)','MCMC');
%--------------------------------------------------------------------------
figure(8) % compare the true track and estimates 
plot(multi_position(3,:),multi_position(4,:),'r')
hold on 
plot(x2,y2,'b')
title('true track (target two) and MCMC estimate');
xlabel('x position of head');ylabel('y position of head');
legend('True track (target two)','MCMC');
%--------------------------------------------------------------------------
figure(9) % plot the error on x position
plot(EUerror2)
title('Euclidean error of target two)');
xlabel('Frames');ylabel('magnitude(pixels)');
%--------------------------------------------------------------------------