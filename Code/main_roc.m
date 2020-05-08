imgFolder = 'D:\MSC\Term2\Random process\HWs\HW3\DataSet\g512_001';
%imgFolder = 'D:\MSC\Term2\Random process\HWs\HW3\DataSet\All';
imgPattern = fullfile(imgFolder, '*.pgm');
pgmFiles = dir(imgPattern);
display('	FileName    mu    alpha    beta') ; 
addpath('D:\MSC\Term2\Random process\HWs\HW3\Code\latexTable-master')
MLE_G=[] ; 
sz = [100 4];
ggamma = 0.2 ; 
varTypes = {'string' ,'double' ,'double' ,'double' , 'string' , 'double' , 'double'};
%Paramet_table = table('Size',sz,'VariableTypes',varTypes);
Paramet_table  = table({ ' ' }, [0] , [0] , [0],{' '} ,'VariableNames',{'FileName' 'mu' 'alpha' 'beta' 'Comment'  })
LKH_table = table({ ' ' }, [0] , [0] , [0],'VariableNames',{'FileName' 'WaterMarked_LKH' 'WithoutWatermaking_LKH' 'Hypo_Result' })

	for k = 1:length(pgmFiles)
	  baseFileName = pgmFiles(k).name;
	  fullFileName = fullfile(imgFolder, baseFileName);
	  fprintf(1, ' %s\n', baseFileName);
	  imageArray = imread(fullFileName);
	  %display (fullFileName)
	  %imshow(imageArray);  % Display image.
	  %drawnow; % Force display to update immediately.
	

	addpath('D:\MSC\Term2\Random process\HWs\HW3\contourlet_toolbox') % toolbox
	%addpath('D:\MSC\Term2\Random process\HWs\HW3\DataSet\g512_001') % Data set
	addpath(imgFolder)
	
	img = imread(baseFileName);
	coeffs = pdfbdec(double(img), '9-7', 'pkva', [ 2, 3]);

	%pyramid_2th = coeffs (1,3) ; % 2th level of pyramid filter
	
	%cll{1} 
	%var(cll{1}{1}(:))

	% ---------------------------------------------------------------key value generation 

	pyramid_2th = coeffs (1,3) ; % 2th level of pyramid filter
	idx =var_max(pyramid_2th , 8)  ;
	max_var_subband = pyramid_2th{1}{idx};
	max_var_subband_db = mat2gray(max_var_subband) ;
	
	size_pyramid = size(pyramid_2th{1}{idx}) ; 
	key = randi([0 1],size_pyramid(1),size_pyramid(2))*2 - 1;



	%var_max(cll,8)
	watermarked_img = watermarking(img , key,idx , ggamma);
	watermarked_img_db = mat2gray(watermarked_img);
	
	coeffs_water = pdfbdec(double(watermarked_img), '9-7', 'pkva', [ 2, 3]);
	pyramid_2th_water = coeffs_water (1,3);
	idx_water =var_max(pyramid_2th , 8)  ;
	max_var_subband_water = pyramid_2th_water{1}{idx_water};
	max_var_subband_water_db = mat2gray(max_var_subband_water) ;	
	
	% Show the reconstructed image and the original image

	%subplot(1,2,1), imagesc( img, [0, 255] ); 
	%title('Original image' ) ;
	%axis image off;
	%subplot(1,2,2), imagesc( watermarked_img, [0, 255] );
	%title('Watermerked image' ) ;
	%axis image off;
	GG = @(data,mu,alpha,beta)( ( exp ( - ( (abs ( data - mu ) / alpha) .^ beta ) )  * beta/ ( 2 * alpha * gamma ( 1 / beta ) )  )  ) ;
	GG_h0 = @(data,mu,alpha,beta)( log( exp ( - ( (abs ( data - mu ) / alpha) .^ beta ) )  * beta/ ( 2 * alpha * gamma ( 1 / beta ) )  )  ) ;
	GG_h1 = @(data,mu,alpha,beta,w,Gamma)( log( exp ( - ( (abs ( data -(mu+Gamma*w) ) / alpha) .^ beta ) )  * beta/ ( 2 * alpha * gamma ( 1 / beta ) )   ) ) ;
	
	GG_h00 = @(data,mu,alpha,beta)( ( exp ( - ( (abs ( data - mu ) / alpha) .^ beta ) )  * beta/ ( 2 * alpha * gamma ( 1 / beta ) )  )  ) ;
	GG_h11 = @(data,mu,alpha,beta,w,Gamma)( ( exp ( - ( (abs ( data -(mu+Gamma*w) ) / alpha) .^ beta ) )  * beta/ ( 2 * alpha * gamma ( 1 / beta ) )   ) ) ;	
	
	
	%Gu =  @(data,mu,sigma)( 1/(sigma*sqrt(2*pi))*exp(-(abs(data-mu)/(2*sigma^2).^2)) )
	%normImage = mat2gray(img);
	normImage_db = mat2gray(img);
	%imshow(normImage )


	%phat = mle(MPG ,'pdf',GG,'start',[.1,.8,1]) 
	%imshow(normImage_db )
	%phat = mle(normImage_db(:),'distribution','Normal')

	%phat = mle(normImage_db(:),'distribution','Normal');
	%display('in order:  mu, alpha, beta')
	try
		phat = mle(normImage_db(:) ,'pdf',GG,'start',[1,1,1]);
	catch
		Paramet_table = [Paramet_table ; {baseFileName,0 ,0 ,0},' '];
		LKH_table = [LKH_table ; {baseFileName, 0 , 0 ,0}];
		continue
	end
	%MLE_G = [MLE_G ; phat(1) ,phat(2) ,phat(3)];
	muu = phat(1) ; 
	alphaa = phat(2) ; 
	betta = phat(3) ;
	
	tereshold =  0 ; 
	P_detection = integral( @(data)GG_h11(data,muu,alphaa , betta, 1, 0.2) , treshold , inf ) + ...
						   %integral( @(data)GG_h11(data,muu,alphaa , betta, -1, 0.2) , treshold , inf )	
	P_fa = integral( @(data)GG_h00(data,muu,alphaa , betta) , treshold , inf )
	
	
	LK_h1_water = (GG_h1(max_var_subband_water_db(:) , muu , alphaa ,betta , key(:),ggamma));
	LK_h1 = (GG_h1 (max_var_subband_db(:) , muu , alphaa ,betta , key(:),ggamma)) ;  
	%LK_h1  = LK_h1  (LK_h1>=0.01)  ;
	%size(LK_h1)
	LK_h0_water = (GG_h0(max_var_subband_water_db(:) , muu , alphaa ,betta)) ;
	LK_h0 = (GG_h0(max_var_subband_db(:) , muu , alphaa ,betta)) ; 
	%LK_h0  = LK_h0 (LK_h0>=0.01)  ;
	%size(LK_h0)
	
	%likelihood = sum(LK_h1)/ sum(LK_h0)  ;
	%likelihood_water =  sum(LK_h1_water)/ sum(LK_h0_water)  ; 
	likelihood = sum(LK_h1)  - sum(LK_h0)  ;	
	likelihood_water =  sum(LK_h1_water) -  sum(LK_h0_water)  ;	
	
	Paramet_table = [Paramet_table ; { baseFileName,phat(1) ,phat(2) ,phat(3)},' ' ]; % append to end of table
	LKH_table = [LKH_table ; {baseFileName, likelihood_water , likelihood , likelihood_water > 0}] ; 
	
	
	
	
	%likelihood = prod(LK_h1) > 0.2*prod(LK_h0 );
	%h=lratiotest(LK_h0,LK_h1,1) ; 
	%size(h)
	%prod(LK_h0 )
	%prod(LK_h1 )
	
	
	%MLE_G = [MLE_G ; phat(1) ,phat(2) ];
	
	%imshow(cll{1}{1})
	%imshow(img )
	%imcoeff = showpdfb( coeffs) ;

	%img = imread('zoneplate.png')
	%decdemo( img , 'auto')
end
%display(Paramet_table)
display(LKH_table)
 %writetable(Paramet_table, 'LKH.csv','QuoteStrings',true);
%writetable(LKH_table, 'LKH.csv','QuoteStrings',true);
