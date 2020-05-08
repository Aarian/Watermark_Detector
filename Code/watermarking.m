function watermarked_img = watermarking( img ,  key, max_var_idx_subband , Gamma)
% watermarked_img : contains the concealed key value in the max-variance subband of the 2th level of the pyramid filtered image  
... (watermarked_img ==whole image)

% img :  initial img

% key : [-1,1] array with the size of the max-variance subband  ~ Unif()

%max_var_idx_subband : maximum variance subband index 

% Gamma :  0.2

addpath('D:\MSC\Term2\Random process\HWs\HW3\contourlet_toolbox') % toolbox

coeffs = pdfbdec(double(img), '9-7', 'pkva', [ 2, 3]) ;
pyramid_2th = coeffs (1,3) ; % 2th level of pyramid filter
copy = coeffs (1,3) ;
%size(key)

%%% Watermarking
watermarked_subband = pyramid_2th{1}{max_var_idx_subband}+ Gamma * key ; 
pyramid_2th{1}{max_var_idx_subband} = watermarked_subband;
coeffs (1,3) = pyramid_2th ;


%% reconstructing Water marked img

watermarked_img = pdfbrec( coeffs, '9-7', 'pkva') ;

%%----------------------- Just check for embbeding the data  
%xx = coeffs (1,3) ;
%yy = copy ;
%  yy{1}{max_var_idx_subband}== xx{1}{max_var_idx_subband}
end