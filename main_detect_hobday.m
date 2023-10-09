%% Detect point-wise MHWs following Hobday et al.,2016 at each depth 
% Load 4-D temperature data
% sst_full 360*120*days*depth_layers e.g. 360*120*10593*20 
% 1992.1.1-2020.12.31 total days:10593d
% 20 is for 20 layers (e.g. 10m per layer for total 200m)

% Here we detect marine heatwaves based on the traditional definition of MHWs (Hobday et al. 2016). 
% We detected MHWs for each depth during 1992 to 2020 for climatologies and thresholds in 1992 to 2020.
for i = 1:20
[MHW,mclim,m90,mhw_ts]=detect(sst_full(:,:,:,i),datenum(1992,1,1):datenum(2020,12,31),datenum(1992,1,1),datenum(2020,12,31),datenum(1992,1,1),datenum(2020,12,31));
disp('saving')
save(['./depth' num2str(i) '/mclim.mat'],'mclim');
save(['./depth' num2str(i) '/MHW.mat'],'MHW');
save(['./depth' num2str(i) '/mhw_ts.mat'],'mhw_ts','-v7.3');
save(['./depth' num2str(i) '/m90.mat'],'m90','-v7.3');

% mhw metrics based on Hobday et al., 2016
% if i==1
%     metric_used={'Frequency','MeanInt','MaxInt','Duration'};
% 
%     for j=1:4
%         eval(['[mean_' metric_used{j} ',annual_' metric_used{j} ',trend_' metric_used{j} ',p_' metric_used{j} ']=mean_and_trend_test(MHW,mhw_ts,1993,' '''' 'Metric' '''' ',' 'metric_used{j}' ');'])
%     end
% 
%     save('./depth1/mean_Frequency.mat','mean_Frequency');
%     save('./depth1/mean_MeanInt.mat','mean_MeanInt');
%     save('./depth1/mean_MaxInt.mat','mean_MaxInt');
%     save('./depth1/mean_Duration.mat','mean_Duration');
% 
% end
end
%% Make masks following MHWs as 1 and non-MHW as 0
for i = 1:20
    load(['./depth' num2str(i) '/mhw_ts.mat'])
    mask = zeros(360,120,datenum(2020,12,31)-datenum(1992,1,1)+1);
    mask(~isnan(mhw_ts) & mhw_ts ~= 0) = 1;
    mask_now(:,:,:,i) = mask;
    disp(i)
end
%% Split into single file for different years
for i = 1992:2020
    mask_day = mask_now(:,:,datenum(i,1,1)-datenum(1992,1,1)+1:datenum(i,12,31)-datenum(1992,1,1)+1,:);
    save(['./mask_day_' num2str(i) '.mat'],'mask_day','-v7.3')
    disp(i)
end