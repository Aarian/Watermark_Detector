function indx = var_max(cell_arr, cell_arr_length )
	% cell array contains subband images i.e cell{1}{1} means the first image 
	% cell_arr_length :  means number of subband images in the cell array
	% indx : index of the maximum variance image in the cell array 
	
	max_var_arr =[];
	for ii  = [1 : cell_arr_length] 
		max_var_arr =  [max_var_arr var(cell_arr{1}{ii}(:)) ];
	end
	[val ,idx] = max(max_var_arr);
	indx = idx;
end