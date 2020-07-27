%% Head Tracking based on Particle filter (Two-target)
%These codes are provided for academic/research purposes only. 
% They are designed for readability and demonstration thus they are not optimized for speed.
% Author: Zeyu Fu
% PGR, Newcastle University, UK
%==========================================================================
clc; clear all
close all
%==========================================================================
load multi_position;
load ahistogram_head_image1;
load ahistogram_head_image2;
%==========================================================================
%% 1.Initialisation--------------------------------------------------------
N_particle = 100;   % number of particles
N = 100;
N_pixel = 30; % number of pixel
theta = linspace(0,2*pi,30); % 
A = eye(4);    % state transition matrix
v_g=0.00001;   % variance of intensity gradient
v_c=0.00005;   % variance of colour histogram
s1=[268 96 17 22]';  %true position of first target
s2=[179 79 13 18]';  %true position of second target
alpha = 0.5; %scaling parameter
% prior as importance function
for i=1:N_particle;
    s_particle01(:,i)=s1; % target 1
    s_particle02(:,i)=s2; % target 2
end
%==========================================================================
%% 2.Particle filter and collect the gradient and histogram----------------
for i=115:225
    HEAD=imread([ num2str(i),'.jpg']); % images read in
    I=double(rgb2gray(HEAD)); % obtain the grayscale intensity 
    % 2.1 prediction
    for k=1:N_particle 
        s_particle1(:,k)= A*s_particle01(:,k) +[2.4*randn;1*randn;0.1*randn;0.1*randn]; % target 1
        s_particle2(:,k)= A*s_particle02(:,k) +[2.4*randn;1*randn;0.1*randn;0.1*randn]; % target 2
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
    end
    % target 1
    weight_g1 = weight_g1/(sum(weight_g1));   % normalization of weight on gradient
    weight_c1 = weight_c1/(sum(weight_c1));      % normalization of weight on color histogram
    w1 = alpha*weight_g1 +(1-alpha)*weight_c1;  % combined weight
    w1 = w1/(sum(w1));                           % normalization
    % target 2
    weight_g2 = weight_g2/(sum(weight_g2));   % normalization of weight on gradient
    weight_c2 = weight_c2/(sum(weight_c2));      % normalization of weight on color histogram
    w2 = alpha*weight_g2 +(1-alpha)*weight_c2;  % combined weight
    w2 = w2/(sum(w2));                           % normalization
    
    % Resampling
    % target 1
    c(1)=w1(1); % intialize the CSW
    for k=2:N
        c(k)=c(k-1)+w1(k);  % construct CSW
    end
    k=1; % Start at the bottom of the CSW
    u(1)=rand*(N^-1); % starting point
    for j=1:N
        u(j)=u(1)+(N^-1)*(j-1); % Move along the CSW
        while u(j)>c(k)
            k=k+1;
        end
        s_particle1(:,j)=s_particle1(:,k);% assign new samples
    end
    % target 2
    c(1)=w2(1); % intialize the CSW
    for k=2:N
        c(k)=c(k-1)+w2(k);  % construct CSW
    end
    k=1; % Start at the bottom of the CSW
    u(1)=rand*(N^-1); % starting point
    for j=1:N
        u(j)=u(1)+(N^-1)*(j-1); % Move along the CSW
        while u(j)>c(k)
            k=k+1;
        end
        s_particle2(:,j)=s_particle2(:,k);% assign new samples
    end
    
    % Update
    % target 1
    s_estimated1 = s_particle1*w1'; % mean estimated position
    I1=imcrop(HEAD,round([s_estimated1(1)-s_estimated1(3) s_estimated1(2)-s_estimated1(4) 2*s_estimated1(3) 2*s_estimated1(4)]));
    [ahistogram_head_image1] = Update1(ahistogram_head_image1,I1); % update the target model (histogram)
    % target 2
    s_estimated2 = s_particle2*w2'; % mean estimated position
    I2=imcrop(HEAD,round([s_estimated2(1)-s_estimated2(3) s_estimated2(2)-s_estimated2(4) 2*s_estimated2(3) 2*s_estimated2(4)]));
    [ahistogram_head_image2] = Update2(ahistogram_head_image2,I2); % update the target model (histogram)
    for m=1:N_particle;
    s_particle01(:,m)=s_estimated1;
    s_particle02(:,m)=s_estimated2;% prior as importance function
    end
    % draw points to match the ellipse
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
    % plot points 
    figure(1)   
    imshow(HEAD)
    hold on
    plot(point_x1,point_y1,'r')
    hold on
    plot(point_x2,point_y2,'r')
    pause (0.1)
    % Record
    % target 1
    x1(1,i-114)=s_estimated1(1);
    y1(1,i-114)=s_estimated1(2); 
    % target 2
    x2(1,i-114)=s_estimated2(1);
    y2(1,i-114)=s_estimated2(2); 
end
%% 4.Calculate the Mean square error-----------------------------------------------------------------
for t=1:111
    % target 1
error_x1(t) = (x1(t)-multi_position(1,t))*(x1(t)-multi_position(1,t))';
error_y1(t) = (y1(t)-multi_position(2,t))*(y1(t)-multi_position(2,t))';
    % target 2
error_x2(t) = (x2(t)-multi_position(3,t))*(x2(t)-multi_position(3,t))';
error_y2(t) = (y2(t)-multi_position(4,t))*(y2(t)-multi_position(4,t))';
end
    % target 1
MSE_error_x1 = 1/111*(sum(error_x1)) %MSE
MSE_error_y1 = 1/111*(sum(error_y1))
    % target 2
MSE_error_x2 = 1/111*(sum(error_x2))  %MSE
MSE_error_y2 = 1/111*(sum(error_y2))
%==========================================================================
%% 5.Figures---------------------------------------------------------------
figure(2) % compare the true and estimated x position 
plot(multi_position(1,:),'r');
hold on 
plot(x1,'b')
title('true position of x (target one) and PF estimate(SIR)');
xlabel('Frames');ylabel('position x');
legend('True position x (target one)','Particle filter(SIR)');
%--------------------------------------------------------------------------
figure(3) % compare the true and estimated y position 
plot(multi_position(2,:),'r')
hold on
plot(y1,'b')
title('true position of y (target one) and PF estimate(SIR)');
xlabel('Frames');ylabel('position y');
legend('True position y (target one)','Particle filter(SIR)');
%--------------------------------------------------------------------------
figure(4) % compare the true track and estimates 
plot(multi_position(1,:),multi_position(2,:),'r')
hold on 
plot(x1,y1,'b')
title('true track (target one) and PF estimate(SIR)');
xlabel('x position of head');ylabel('y position of head');
legend('True track (target one)','Particle filter(SIR)');
%--------------------------------------------------------------------------
figure(5) % plot the error on x position
plot(MSE_error_x1)
title('Estimated error of x position (target one)');
xlabel('Frames');ylabel('magnitude');
%--------------------------------------------------------------------------
figure(6) % plot the error on y position
plot(MSE_error_y1)
title('Estimated error of y position (target one)');
xlabel('Frames');ylabel('magnitude');
%==========================================================================
figure(7) % compare the true and estimated x position 
plot(multi_position(3,:),'r');
hold on 
plot(x2,'b')
title('true position of x (target two) and PF estimate(SIR)');
xlabel('Frames');ylabel('position x');
legend('True position x (target two)','Particle filter(SIR)');
%--------------------------------------------------------------------------
figure(8) % compare the true and estimated y position 
plot(multi_position(4,:),'r')
hold on
plot(y2,'b')
title('true position of y (target two) and PF estimate(SIR)');
xlabel('Frames');ylabel('position y');
legend('True position y (target two)','Particle filter(SIR)');
%--------------------------------------------------------------------------
figure(9) % compare the true track and estimates 
plot(multi_position(3,:),multi_position(4,:),'r')
hold on 
plot(x2,y2,'b')
title('true track (target two) and PF estimate(SIR)');
xlabel('x position of head');ylabel('y position of head');
legend('True track (target two)','Particle filter(SIR)');
%--------------------------------------------------------------------------
figure(10) % plot the error on x position
plot(MSE_error_x2)
title('Estimated error of x position (target two)');
xlabel('Frames');ylabel('magnitude');
%--------------------------------------------------------------------------
figure(11) % plot the error on y position
plot(MSE_error_y2)
title('Estimated error of y position (target two)');
xlabel('Frames');ylabel('magnitude');
%==========================================================================