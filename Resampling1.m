function [s_particle1]=Resampling1(w1)
 c(1)=w1(1); % intialize the CSW
 N=100;
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