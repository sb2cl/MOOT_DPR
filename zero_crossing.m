v_mean=mean(sigma);
zero_mean=sigma-v_mean;
zero_cross=zeros(size(sigma));
counter=0;
for j=2:length(zero_mean)
    if zero_mean(j)>=0
        zero_cross(j)=1;
    end
    counter=counter+(zero_cross(j)-zero_cross(j-1))^2;
end
oscillation=counter/2;
