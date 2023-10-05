% ssta_now is 360*120*10593*20 and calculated by temperature minus seasonal
% varying climatology
depth = 20;
time = datenum(2020,12,31)-datenum(1992,1,1)+1;
features_train = zeros(11*11*depth,time);

% The next procedure calculates the corrcoef in a 10¡ãx10¡ãx200m (11x11x20) region
% for each point in the interior

for m = 6:360-5
    for n = 6:120-5
        ssta1 = ssta_now(m-5:m+5,n-5:n+5,:,:);
        
        for k = 1:depth
           for i = 1:11
              for j = 1:11
                  features_train(j + (i + (k-1)*11 -1)*11,:) = ssta1(i,j,:,k);
              end
           end
        end
        r = corrcoef(features_train');
        save(['./r_' num2str(m) '_' num2str(n) '.mat'],'r')
        
    end
end

%% for the boundary

for m = 1:5
    for n = 1:5
        ssta1 = ssta_now(1:11,1:11,:,:);
        
        for k = 1:depth
           for i = 1:11
              for j = 1:11
                  features_train(j + (i + (k-1)*11 -1)*11,:) = ssta1(i,j,:,k);
              end
           end
        end
        r = corrcoef(features_train');
        save(['./r_' num2str(m) '_' num2str(n) '.mat'],'r')
        
    end
end

for m = 1:5
    for n = 6:120-5
        ssta1 = ssta_now(1:11,n-5:n+5,:,:);
        
        for k = 1:depth
           for i = 1:11
              for j = 1:11
                  features_train(j + (i + (k-1)*11 -1)*11,:) = ssta1(i,j,:,k);
              end
           end
        end
        r = corrcoef(features_train');
        save(['./r_' num2str(m) '_' num2str(n) '.mat'],'r')
        
    end
end

for m = 1:5
    for n = 116:120
        ssta1 = ssta_now(1:11,end-10:end,:,:);
        
        for k = 1:depth
           for i = 1:11
              for j = 1:11
                  features_train(j + (i + (k-1)*11 -1)*11,:) = ssta1(i,j,:,k);
              end
           end
        end
        r = corrcoef(features_train');
        save(['./r_' num2str(m) '_' num2str(n) '.mat'],'r')
        
    end
end

for m = 6:360-5
    for n = 1:5
        ssta1 = ssta_now(m-5:m+5,1:11,:,:);
        
        for k = 1:depth
           for i = 1:11
              for j = 1:11
                  features_train(j + (i + (k-1)*11 -1)*11,:) = ssta1(i,j,:,k);
              end
           end
        end
        r = corrcoef(features_train');
        save(['./r_' num2str(m) '_' num2str(n) '.mat'],'r')
        
    end
end

for m = 356:360
    for n = 1:5
        ssta1 = ssta_now(end-10:end,1:11,:,:);
        
        for k = 1:depth
           for i = 1:11
              for j = 1:11
                  features_train(j + (i + (k-1)*11 -1)*11,:) = ssta1(i,j,:,k);
              end
           end
        end
        r = corrcoef(features_train');
        save(['./r_' num2str(m) '_' num2str(n) '.mat'],'r')
        
    end
end

for m = 356:360
    for n = 6:120-5
        ssta1 = ssta_now(end-10:end,n-5:n+5,:,:);
        
        for k = 1:depth
           for i = 1:11
              for j = 1:11
                  features_train(j + (i + (k-1)*11 -1)*11,:) = ssta1(i,j,:,k);
              end
           end
        end
        r = corrcoef(features_train');
        save(['./r_' num2str(m) '_' num2str(n) '.mat'],'r')
        
    end
end

for m = 356:360
    for n = 116:120
        ssta1 = ssta_now(end-10:end,end-10:end,:,:);
        
        for k = 1:depth
           for i = 1:11
              for j = 1:11
                  features_train(j + (i + (k-1)*11 -1)*11,:) = ssta1(i,j,:,k);
              end
           end
        end
        r = corrcoef(features_train');
        save(['./r_' num2str(m) '_' num2str(n) '.mat'],'r')
        
    end
end

for m = 6:360-5
    for n = 116:120
        ssta1 = ssta_now(m-5:m+5,end-10:end,:,:);
        
        for k = 1:depth
           for i = 1:11
              for j = 1:11
                  features_train(j + (i + (k-1)*11 -1)*11,:) = ssta1(i,j,:,k);
              end
           end
        end
        r = corrcoef(features_train');
        save(['./r_' num2str(m) '_' num2str(n) '.mat'],'r')
        
    end
end

