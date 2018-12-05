function [LS_USA_45m, LS_USA_filt_45m] = gen_fine_res_LS(LandScan_data, LandScan_data_filt)
LS_USA_45m = cell2mat(arrayfun(@(LandScan_data) repmat(LandScan_data,2,2), LandScan_data, 'uniform',0));
LS_USA_filt_45m = cell2mat(arrayfun(@(LandScan_data_filt) repmat(LandScan_data_filt,2,2), LandScan_data_filt, 'uniform',0));

end