% Load required packages
%ltfat version 2.4.0+ required
pkg load ltfat;
pkg load image;
pkg load geometry;

% Specify filter variant
% Consider only negative deviation from expected intensity
typemax = 0;      
% Consider only negative deviations from expected pair corr
typemaxcorr = 0;  
% Factor for correlation mask before thresholding
fact_corr = 1;  
% Factor for intensity mask before thresholding  
fact_intens = 1;  
%subtract this from all mask values in pair correlation filtering
% corresponds to  0.99 quantile in our setup
%threshold_corr = 3.58; #3.64, 3.54
%subtract this from all mask values in intensity filtering
%threshold_intens = 3.15; #3.12, 3.16
%subtract this from all mask values in combined filtering
%threshold_comb = 3.79; #3.79, 3.74
%0.999 quantile values over three noise realizations:
threshold_corr = 5.42; #5.94, 5.32
threshold_intens = 4.32; #4.26, 4.49
threshold_comb = 5.43; #5.53, 5.17
% Flag: Artificially add imaginary noise
imnoise = 0;   
% Flag: no signal desired  
onlynoise = 0;     
% negative lower limit in log plots
disp_lowlim = 100;
% negative upper limit in log plots
disp_uplim = 20;  

% Set parameters
% Number of scales
filtNos = 600;		
% Cauchy wavelet order
alphas  = 300;
% Vector of Scales		
y_scaling = 2.^(linspace(6, -3.3, filtNos));  
% Decimation factor in time
as = 1;		
% Load sample based gcorr standard deviations	
load 'stdgcorr.txt';

% Desired signal length
L = 2*44100;		
% Starting from sample 'start'		
start = 0.4*44100;	
% Signal to noise ratio for additive noise
SNR = 5;    


% Prepare to load test signal
% Place your test signal wav files here
path = './input/'; 
listing = dir([path,'*.flac']);
allwavsCell = arrayfun(@(wav) fullfile(path,wav.name), listing, ...
              'UniformOutput',0);
% The experiment output will be saved here, make sure the folder exists
writepath = './output/'; 

% Loop over files
% for ii=1:numel(allwavsCell) 
% or consider only one file
ii=1; 

% Select file and load
wavfile = allwavsCell{ii};
[~,filename,ext] = fileparts(wavfile);
[f,fs] = audioread(wavfile);

% Truncate signal and normalize conservatively
f = f(start+1:start+L,1);   
f = 0.3*normalize(f,'inf');

% Calculate Filterbank corresponding to continuous wavelet transform
[gs,info] = ltfatnote053_waveletfilters(L, y_scaling,{'cauchy',alphas});
gd = filterbankrealdual(gs, as, L);
fprintf('Filterbank calculated \n')

% Keep noiseless signal and its continuous wavelet transform
f_noiseless = f;
c_noiseless = ufilterbank(f_noiseless,gs,as);

% Add noise to the signal such that SNR is as specified
noise=randn(L,1);
l2f = sqrt(mean(abs(f_noiseless).^2));
noise=noise/sqrt(mean(abs(noise).^2))*sqrt(mean(abs(f_noiseless).^2))/SNR;
f = f_noiseless + noise;

% If 'onlynoise' is set, ignore signal and use pure noise
if onlynoise
  f=noise;
  filename='noise';
end;


% Write audio file of noisy input signal
filenameW = [writepath,filename,'_noisy.flac'];
audiowrite(filenameW,real(f),fs);   
    
% If 'imnoise' is set, add also imaginary noise
if imnoise
  inoise=1i*randn(L,1);
  inoise=inoise/sqrt(mean(abs(inoise).^2))*sqrt(mean(abs(f_noiseless).^2))/SNR;
  f = f + inoise;
end;


% Compute continuous wavelet coefficients of noisy signal
c_orig = ufilterbank(f,gs,as);
abss = abs(c_orig);
[N,M] = size(abss);
fprintf('Wavelet Transform Complete \n');

figure(1);
plotfilterbank(max(10^(-disp_lowlim/20), abs(c_orig)), as, ...
               'fc', info.fc*fs/2, 'fs', fs, ...
               'clim', [-disp_lowlim, -disp_uplim], 'audtick');

% Scale correction: due to dimensionless construction of the filterbank, 
% this is the factor that transforms the y_scaling values to correct y values
scale_corr = 5*(alphas-1)/(pi*fs);
% Grid specification:
% vectors specifying correct (x,y) values for (xind,yind) index pairs of c_orig
xgrid = (1:L)'*as/fs;
ygrid = y_scaling'*scale_corr;


% Find local minima in the matrix abss ignoring first row (corresponding to
%  low pass filter)
% These minima correspond to zeros of the CWT
locminmat = locminima(abss(:,2:end));
[locminxind,locminyind]=find(locminmat);

% Calculate correct coordinate positions of the zeros from the indices
locminx = xgrid(locminxind);
locminy = ygrid(locminyind);

% Sort local minima by their x coordinate for more efficient processing
[locmin,mysort]=sortrows([locminx,locminy],1);
locminx = locmin(:,1);
locminy = locmin(:,2);
[locminind,mysortind]=sortrows([locminxind,locminyind],1);
locminxind = locminind(:,1);
locminyind = locminind(:,2);


% Scatter plot of local minima / zeros
figure(2);
scatter(locminx, locminy,9);

% Calculate number of neighbors at a certain distance
%   On average we want to average over 5 zeros
%   Since the average number of zeros in a disk of radius r_base is 
%   (alphas*r_base^2)/(1-r_base^2)
%   we choose 
r_base = sqrt(5/(5+alphas));
% Number of different radii
nrr0 = 4;
nrr1 = 5;
hstep = r_base/3;
hwin = r_base/3;
r_1 = (3:(nrr1+2))*r_base/5;


% Initialize variables
% Number of neighbors at nrr0 different hyperbolic radii
nrofneig=zeros(nrr0,length(locminx));
% Local estimate of (corrected) pair correlation function
%   based on zeros in nrr1 different radii and their neighbors in nrr0 different
%   radii
gcorr=zeros(N,M-1,nrr0,nrr1);
% Local estimate of (unnormalized)intensity function based on 
%   zeros in nrr1 different radii
rhoest=zeros(N,M-1,nrr1);

% Separate zeros into regions that can be processed individually
maxrx = 2*r_1(nrr1)/(1-r_1(nrr1)^2) + 2*nrr0*hstep/(1-(nrr0*hstep)^2);
maxryfac = (1-r_1(nrr1))/(1+r_1(nrr1))*(1-nrr0*hstep)/(1+nrr0*hstep);
max_y = max(locminy);
min_y = min(ygrid);
nrofyshifts = ceil(log(min_y/max_y)/log(maxryfac));
total_processed = 0;
for ky = 1:nrofyshifts
  max_y_loc = max_y*maxryfac^(ky-1);
  min_y_loc = max_y*maxryfac^(ky);
  max_y_close = max_y*maxryfac^(max(0,ky-2));
  min_y_close = max_y*maxryfac^(ky+1);
  x_shift = maxrx*max_y_close;
  nrofshifts = ceil(L/fs/x_shift);
%consider iteratively x-ranges of width 3*x_shift shifted by x_shift

  for kx = 1:nrofshifts
    max_x_loc = (kx)*x_shift;
    min_x_loc = (kx-1)*x_shift;
    max_x_close = (kx+1)*x_shift;
    min_x_close = (kx-2)*x_shift;
    indinvestigated = find((locminx>=min_x_loc)&(locminx<max_x_loc)&...
                      (locminy>min_y_loc)&(locminy<=max_y_loc));
    indcloseext = find((locminx>=min_x_close)&(locminx<max_x_close)&...
                      (locminy>min_y_close)&(locminy<=max_y_close));
    total_processed = total_processed + length(indinvestigated);
    %calculate pseudo-hyperbolic distance between zeros in vicinity
    tempdists = zeros(length(indinvestigated),length(indcloseext));
    tempdists = ((locminx(indcloseext)'-locminx(indinvestigated)).^2 + ...
        (locminy(indcloseext)'-locminy(indinvestigated)).^2)...
        ./(locminy(indinvestigated)*locminy(indcloseext)');
    %pseudo-hyperbolic circle of radius r_base around (x_1, y_1) is the same as 
    % Euclidean circle of radius r_euclid around (x_euclid1,y_euclid1)
    x_euclid1 = zeros(length(indinvestigated),1);
    x_euclid1 = locminx(indinvestigated);
    y_euclid1 = zeros(length(indinvestigated),nrr1);
    y_euclid1 = locminy(indinvestigated)*((1+r_1.^2)./(1-r_1.^2));
    r_euclid = zeros(length(indinvestigated),nrr1);
    r_euclid = 2*locminy(indinvestigated)*(r_1./(1-r_1.^2));
    %Calculate which indices are within the circle 
    % of radius r_euclid around (x_euclid1,y_euclid1)
    yfrom = zeros(length(indinvestigated),nrr1);
    yto = zeros(length(indinvestigated),nrr1);
    yfrom = lookup(ygrid',y_euclid1+r_euclid)+1;
    yto = lookup(ygrid',y_euclid1-r_euclid);
    xfrom = zeros(length(indinvestigated),filtNos,nrr1);
    xto = zeros(length(indinvestigated),filtNos,nrr1);
    for rind=1:nrr1
      xfrom(:,:,rind) = max(ceil(fs/as*(x_euclid1-...
            sqrt(max(r_euclid(:,rind).^2-(ygrid'-y_euclid1(:,rind)).^2,0)))),1);
      xto(:,:,rind) = min(floor(fs/as*(x_euclid1+...
            sqrt(max(r_euclid(:,rind).^2-(ygrid'-y_euclid1(:,rind)).^2,0)))),N);
    end;
  
    %Calculate the number of neighbors at various distances r_0
    for rind=1:nrr0
      r_0 = rind*hstep;
      nrofneig(rind,indinvestigated) = sum(tempdists< 4*r_0^2/(1-r_0^2),2)-1;
    end;
    
    %Loop through zeros to calculate estimates of pair correlation gcorr 
    % and intensity rhoest for various radii
    for k = 1:length(indinvestigated)
      for rind=1:nrr1
        for yind = yfrom(k,rind):yto(k,rind)
          gcorr(xfrom(k,yind,rind):xto(k,yind,rind),yind,1:nrr0,rind) ...
            = gcorr(xfrom(k,yind,rind):xto(k,yind,rind),yind,1:nrr0,rind) + ...
            permute(ones(xto(k,yind,rind)-xfrom(k,yind,rind)+1,1) ...
                  *nrofneig(1:nrr0,indinvestigated(k))',[1,3,2]);
            
          rhoest(xfrom(k,yind,rind):xto(k,yind,rind),yind,rind) ...
            = rhoest(xfrom(k,yind,rind):xto(k,yind,rind),yind,rind) + 1;
        end;
        
      end;
      
    end;
    
  end;
  fprintf( 'Processed %d%% of data \n', 100*total_processed/length(locminx));
end;
fprintf( 'Processed 100%% of data, begin filtering \n');


%Normalize pair-correlation function correctly
const_fact = 1/(4*hwin*alphas);
for rind = nrr0:-1:(2*hwin/hstep+1)
  gcorr(:,:,rind,1:nrr1)=gcorr(:,:,rind,1:nrr1)...
                          -gcorr(:,:,rind-(2*hwin/hstep),1:nrr1);
  r_t = hstep*rind-hwin;
  gcorr(:,:, rind,1:nrr1)=gcorr(:,:, rind,1:nrr1)./...
        permute(max(rhoest,0.01),[1,2,4,3])...
        *const_fact*(1-r_t^2)^2/r_t;
end;
r_t = hwin;
gcorr(:,:,(2*hwin/hstep),1:nrr1) = gcorr(:,:,(2*hwin/hstep),1:nrr1)./...
      permute(max(rhoest,0.01),[1,2,4,3])...
      *const_fact*(1-r_t^2)^2/r_t;

%For some indices and radii rind1 there are no zeros in the vicinity and thus
% the estimator of gcorr results in a div by 0
% To avoid this we set those indices to the correct value gintendedalt but 
% keep track of the affected indices and account for them when averaging the
% mask (posgintsum)

for rind = (2*hwin/hstep):nrr0
  r_t = hstep*rind-hwin;
  r_low = hstep*rind-2*hwin;
  r_up = hstep*rind;
  s_low = 1-r_low^2;
  s_t = 1-r_t^2;
  s_up = 1-r_up^2;
  gintended = (s_t^alphas*(alphas*(1-s_t)-s_t*(1-s_t^alphas))^2 ...
          +(alphas*s_t^alphas*(1-s_t)-1+s_t^alphas)^2)/(1-s_t^alphas)^3;
  if r_low > 0
     gintendedalt = 1/(4*hwin)*(1-r_t^2)^2/r_t * ...
        (((alphas+1)*s_low^alphas*(1-s_low)^2 - (1-s_low^(alphas+1))^2 )/...
        (s_low *(1-s_low^alphas)^2)-...
        ((alphas+1)*s_up^alphas*(1-s_up)^2 - (1-s_up^(alphas+1))^2 )/...
        (s_up *(1-s_up^alphas)^2));
  else
     gintendedalt = 1/(4*hwin)*(1-r_t^2)^2/r_t * ...
        (-(alphas+1)/alphas-...
        ((alphas+1)*s_up^alphas*(1-s_up)^2 - (1-s_up^(alphas+1))^2 )/...
        (s_up *(1-s_up^alphas)^2));
  end
  addgint = zeros(N,M-1,1,nrr1);
  %addgint(gcorr(:,:,rind,:)==0) = gintended;
  addgint(gcorr(:,:,rind,:)==0) = gintendedalt;
  gcorr(:,:,rind,:) = gcorr(:,:,rind,:) + addgint; 
end;
posgint = true(N,M-1,1,nrr1);
%posgint(gcorr(:,:,nrr0,:)==gintended) = 0;
posgint(gcorr(:,:,nrr0,:)==gintendedalt) = 0;
posgintsum = sum(posgint,4);

%We begin with the 5 filtering procedures of the signal:
% 1. we filter the noiseless signal to the region of interest
% 2. we filter based on the deviation of gcorr from the expected (corrected)
%     pair correlation function
% 3. we filter based on the deviation of rhoest from the expected 
%     intensity function
% 4. we filter based on the deviation of gcorr from the expected (corrected)
%     pair correlation function and the deviation of rhoest from the expected 
%     intensity function
% 5. we filter the noisy signal to the region of interest

%We first define the region of interest
interestregion = false(N,M-1);

yupper = lookup(y_scaling*scale_corr,...
          y_scaling(filtNos)*scale_corr/maxryfac);
ylower = lookup(y_scaling*scale_corr,...
          y_scaling(1)*scale_corr*maxryfac)+1;
xlower = zeros(yupper - ylower + 1,1);
xupper = zeros(yupper - ylower + 1,1);


for yind = ylower:yupper
  xlower(yind-ylower+1) = ceil(y_scaling(yind)*scale_corr*maxrx*fs);
  xupper(yind-ylower+1) = N-xlower(yind-ylower+1)+1;
  interestregion(xlower(yind-ylower+1):xupper(yind-ylower+1),yind) = true;
end;


%ad 1.
% Filter noiseless signal to region of interest
maxryfac = (1-r_1(nrr1))/(1+r_1(nrr1))*(1-nrr0*hstep)/(1+nrr0*hstep);

masksm = interestregion;

c_manip = zeros(N,M);
c_manip(:,2:601) = c_noiseless(:,2:601).*masksm;

f_noiseless_filt = ifilterbank(c_manip,gd,as,'real');
% Write output file
filenameW = [writepath,filename,'_OrigFilt.flac'];
audiowrite(filenameW,real(f_noiseless_filt),fs);

c_recons = ufilterbank(f_noiseless_filt,gs,as);
figure(3);
plotfilterbank(max(10^(-disp_lowlim/20),abs(c_recons)),as,...
      'fc',info.fc*fs/2,'fs',fs,'clim',[-disp_lowlim,-disp_uplim],'audtick');




%ad 2.
% Filter noisy signal based on modified pair-correlation function
myfinalmask=zeros(N,M-1);
for rind = (2*hwin/hstep):nrr0
  r_t = hstep*rind-hwin;
  r_low = hstep*rind-2*hwin;
  r_up = hstep*rind;
  s_low = 1-r_low^2;
  s_t = 1-r_t^2;
  s_up = 1-r_up^2;
  %gintended = (s_t^alphas*(alphas*(1-s_t)-s_t*(1-s_t^alphas))^2 ...
  %        +(alphas*s_t^alphas*(1-s_t)-1+s_t^alphas)^2)/(1-s_t^alphas)^3;
  if r_low > 0
     gintendedalt = 1/(4*hwin)*(1-r_t^2)^2/r_t * ...
        (((alphas+1)*s_low^alphas*(1-s_low)^2 - (1-s_low^(alphas+1))^2 )/...
        (s_low *(1-s_low^alphas)^2)-...
        ((alphas+1)*s_up^alphas*(1-s_up)^2 - (1-s_up^(alphas+1))^2 )/...
        (s_up *(1-s_up^alphas)^2));
  else
     gintendedalt = 1/(4*hwin)*(1-r_t^2)^2/r_t * ...
        (-(alphas+1)/alphas-...
        ((alphas+1)*s_up^alphas*(1-s_up)^2 - (1-s_up^(alphas+1))^2 )/...
        (s_up *(1-s_up^alphas)^2));
  end
  for rind1 = 1:nrr1
    gcorrtemp = gcorr(:,:,rind,rind1);
      if typemaxcorr
        myfinalmask(interestregion) = myfinalmask(interestregion) + ...
            (max((gintendedalt) ...
            - (gcorrtemp(interestregion)),0)...
            /stdgcorr(rind,rind1)).^2;
      else
        myfinalmask(interestregion) = myfinalmask(interestregion) + ...
            (abs((gintendedalt) ...
            - (gcorrtemp(interestregion)))...
            /stdgcorr(rind,rind1)).^2;
      end;
    
  end;

end;

myfinalmask = myfinalmask/(nrr0 - 2*hwin/hstep+1)./posgintsum;

masksm = min(max(fact_corr*myfinalmask - threshold_corr,0),1);

figure(4);
plotfilterbank(min(max(masksm,0),1),as,...
      'fc',info.fc(2:601)*fs/2,'fs',fs,'clim',[0,1],'audtick','lin');

c_manip = zeros(N,M);
c_manip(:,2:601) = c_orig(:,2:601).*masksm;
figure(5);
plotfilterbank(max(10^(-disp_lowlim/20),abs(c_manip)),as,...
      'fc',info.fc*fs/2,'fs',fs,'clim',[-disp_lowlim,-disp_uplim],'audtick');

    
frecons = ifilterbank(c_manip,gd,as,'real');
snr_pair = sqrt(mean(f_noiseless_filt.^2))/...
            sqrt(mean((frecons-f_noiseless_filt).^2));
% Write output file
filenameW = [writepath,filename,'_reconsPaircor.flac'];
audiowrite(filenameW,real(frecons),fs);

c_recons = ufilterbank(frecons,gs,as);
figure(6);
plotfilterbank(max(10^(-disp_lowlim/20),abs(c_recons)),as,...
      'fc',info.fc*fs/2,'fs',fs,'clim',[-disp_lowlim,-disp_uplim],'audtick');



%ad 3. 
% Filter noisy signal based on intensity function
myfinalmask=zeros(N,M-1);
for rind = 1:nrr1
  rhointended = (alphas*r_1(rind)^2)/(1-r_1(rind)^2);
  myvarrho = alphas^2*r_1(rind)^4/(2*pi*(1-r_1(rind)^2)^2)*...
        integral(@(t) intinvarrho(t, r_1(rind), alphas), -pi, pi);
  rhoesttemp = rhoest(:,:,rind);
    if typemax
      myfinalmask(interestregion) = myfinalmask(interestregion) +  ...
          max(rhointended - rhoesttemp(interestregion),0).^2/myvarrho;
    else
      myfinalmask(interestregion) = myfinalmask(interestregion) +  ...
          abs(rhointended - rhoesttemp(interestregion)).^2/myvarrho;
    end;
  
end;

myfinalmask = myfinalmask/nrr1;
masksm = min(max(fact_intens*myfinalmask - threshold_intens,0),1);

figure(7);
plotfilterbank(min(max(masksm,0),1),as,...
      'fc',info.fc(2:601)*fs/2,'fs',fs,'clim',[0,1],'audtick','lin');

c_manip = zeros(N,M);
c_manip(:,2:601) = c_orig(:,2:601).*masksm;
figure(8);
plotfilterbank(max(10^(-disp_lowlim/20),abs(c_manip)),as,...
      'fc',info.fc*fs/2,'fs',fs,'clim',[-disp_lowlim,-disp_uplim],'audtick');


frecons = ifilterbank(c_manip,gd,as,'real');
snr_intens = sqrt(mean(f_noiseless_filt.^2))/...
              sqrt(mean((frecons-f_noiseless_filt).^2));
% Write output file
filenameW = [writepath,filename,'_reconsIntens.flac'];
audiowrite(filenameW,real(frecons),fs);

c_recons = ufilterbank(frecons,gs,as);
figure(9);
plotfilterbank(max(10^(-disp_lowlim/20),abs(c_recons)),as,...
      'fc',info.fc*fs/2,'fs',fs,'clim',[-disp_lowlim,-disp_uplim],'audtick');




%ad 4. 
% Filter noisy signal based on combination of intensity 
% and modified pair correlation function
myfinalmask=zeros(N,M-1);
for rind = (2*hwin/hstep):nrr0
  r_t = hstep*rind-hwin;
  r_low = hstep*rind-2*hwin;
  r_up = hstep*rind;
  s_low = 1-r_low^2;
  s_t = 1-r_t^2;
  s_up = 1-r_up^2;
  gintended = (s_t^alphas*(alphas*(1-s_t)-s_t*(1-s_t^alphas))^2 ...
          +(alphas*s_t^alphas*(1-s_t)-1+s_t^alphas)^2)/(1-s_t^alphas)^3;
  if r_low > 0
     gintendedalt = 1/(4*hwin)*(1-r_t^2)^2/r_t * ...
        (((alphas+1)*s_low^alphas*(1-s_low)^2 - (1-s_low^(alphas+1))^2 )/...
        (s_low *(1-s_low^alphas)^2)-...
        ((alphas+1)*s_up^alphas*(1-s_up)^2 - (1-s_up^(alphas+1))^2 )/...
        (s_up *(1-s_up^alphas)^2));
  else
     gintendedalt = 1/(4*hwin)*(1-r_t^2)^2/r_t * ...
        (-(alphas+1)/alphas-...
        ((alphas+1)*s_up^alphas*(1-s_up)^2 - (1-s_up^(alphas+1))^2 )/...
        (s_up *(1-s_up^alphas)^2));
  end
  for rind1 = 1:nrr1
    gcorrtemp = gcorr(:,:,rind,rind1);
      if typemaxcorr
        myfinalmask(interestregion) = myfinalmask(interestregion) + ...
            (max((gintendedalt) - (gcorrtemp(interestregion)),0)/...
                stdgcorr(rind,rind1)).^2;
      else
        myfinalmask(interestregion) = myfinalmask(interestregion) + ...
            (abs((gintendedalt) - (gcorrtemp(interestregion)))/...
                stdgcorr(rind,rind1)).^2;
      end;
    
  end;

end;
myfinalmask = fact_corr*myfinalmask/(nrr0 - 2*hwin/hstep+1)./posgintsum;

for rind = 1:nrr1
  rhointended = (alphas*r_1(rind)^2)/(1-r_1(rind)^2);
  myvarrho = alphas^2*r_1(rind)^4/(2*pi*(1-r_1(rind)^2)^2)*...
        integral(@(t) intinvarrho(t, r_1(rind), alphas), -pi, pi);
  rhoesttemp = rhoest(:,:,rind);
    if typemax
      myfinalmask(interestregion) ...
        = myfinalmask(interestregion) + fact_intens*...
          max(rhointended - rhoesttemp(interestregion),0).^2/myvarrho/nrr1;
    else
      myfinalmask(interestregion) ...
        = myfinalmask(interestregion) + fact_intens*...
          abs(rhointended - rhoesttemp(interestregion)).^2/myvarrho/nrr1;
    end;
  
end;

masksm = min(max(myfinalmask/sqrt(2) - threshold_comb,0),1);

figure(10);
plotfilterbank(min(max(masksm,0),1),as,...
      'fc',info.fc(2:601)*fs/2,'fs',fs,'clim',[0,1],'audtick','lin');

c_manip = zeros(N,M);
c_manip(:,2:601) = c_orig(:,2:601).*masksm;
figure(11);
plotfilterbank(max(10^(-disp_lowlim/20),abs(c_manip)),as,...
      'fc',info.fc*fs/2,'fs',fs,'clim',[-disp_lowlim,-disp_uplim],'audtick');



frecons = ifilterbank(c_manip,gd,as,'real');
snr_comb = sqrt(mean(f_noiseless.^2))/sqrt(mean((frecons-f_noiseless).^2));
% Write output file
filenameW = [writepath,filename,'_reconsComb.flac'];
audiowrite(filenameW,real(frecons),fs);

c_recons = ufilterbank(frecons,gs,as);
figure(12);
plotfilterbank(max(10^(-disp_lowlim/20),abs(c_recons)),as,...
      'fc',info.fc*fs/2,'fs',fs,'clim',[-disp_lowlim,-disp_uplim],'audtick');


%ad 5.
% Filter noisy signal to region of interest

masksm = interestregion;

%figure(13);
%plotfilterbank(min(max(masksm,0),1),as,...
%      'fc',info.fc(2:601)*fs/2,'fs',fs,'clim',[0,1],'audtick','lin');

c_manip = zeros(N,M);
c_manip(:,2:601) = c_orig(:,2:601).*masksm;
figure(14);
plotfilterbank(max(10^(-disp_lowlim/20),abs(c_manip)),as,...
      'fc',info.fc*fs/2,'fs',fs,'clim',[-disp_lowlim,-disp_uplim],'audtick');


frecons = ifilterbank(c_manip,gd,as,'real');
snr_filt = sqrt(mean(f_noiseless_filt.^2))/...
      sqrt(mean((frecons-f_noiseless_filt).^2));
% Write output file
filenameW = [writepath,filename,'_reconsFilt.flac'];
audiowrite(filenameW,real(frecons),fs);

c_recons = ufilterbank(frecons,gs,as);
figure(15);
plotfilterbank(max(10^(-disp_lowlim/20),abs(c_recons)),as,...
      'fc',info.fc*fs/2,'fs',fs,'clim',[-disp_lowlim,-disp_uplim],'audtick');



%Optional plotting of single masks:
%
%
%maskN = size(rhoest,1)/20;
%maskM = size(rhoest,2);
%%
%mymask=zeros(maskN,maskM,nrr1);
%for rind = 1:nrr1
%  rhointended = (alphas*r_1(rind)^2)/(1-r_1(rind)^2);
%  myvarrho = alphas^2*r_1(rind)^4/(2*pi*(1-r_1(rind)^2)^2)*...
%        integral(@(t) intinvarrho(t, r_1(rind), alphas), -pi, pi);
%  yupper = lookup(y_scaling*scale_corr,...
%            y_scaling(filtNos)*scale_corr/maxryfac);
%  ylower = lookup(y_scaling*scale_corr,...
%            y_scaling(1)*scale_corr*maxryfac)+1;
%  for yind = ylower:yupper
%    xlower = ceil(y_scaling(yind)*scale_corr*maxrx*fs);
%    xupper = N-xlower+1;
%    if typemax
%      mymask(ceil(xlower/20):ceil(xupper/20),yind,rind) = ...
%            2/9*max(rhointended - rhoest(ceil(xlower/20)*20:20:ceil(xupper/20)*20,yind,rind),0).^2/(myvarrho);
%    else
%      mymask(ceil(xlower/20):ceil(xupper/20),yind,rind) = ...
%            1/9*abs(rhointended - rhoest(ceil(xlower/20)*20:20:ceil(xupper/20)*20,yind,rind)).^2/(myvarrho);
%    end;
%      
%  end;
%  figure(20+rind);
%  
%  imagesc([0,2], [1,0],mymask(:,:,rind)', climits = [0, 1])
%end;
%%
%mymask=zeros(maskN,maskM,nrr0,nrr1);
%for rind = (2*hwin/hstep):nrr0
%  r_t = hstep*rind-hwin;
%  r_low = hstep*rind-2*hwin;
%  r_up = hstep*rind;
%  s_low = 1-r_low^2;
%  s_t = 1-r_t^2;
%  s_up = 1-r_up^2;
%  gintended = (s_t^alphas*(alphas*(1-s_t)-s_t*(1-s_t^alphas))^2 ...
%          +(alphas*s_t^alphas*(1-s_t)-1+s_t^alphas)^2)/(1-s_t^alphas)^3;
%  if r_low > 0
%     gintendedalt = 1/(4*hwin)*(1-r_t^2)^2/r_t * ...
%        (((alphas+1)*s_low^alphas*(1-s_low)^2 - (1-s_low^(alphas+1))^2 )/...
%        (s_low *(1-s_low^alphas)^2)-...
%        ((alphas+1)*s_up^alphas*(1-s_up)^2 - (1-s_up^(alphas+1))^2 )/...
%        (s_up *(1-s_up^alphas)^2));
%  else
%     gintendedalt = 1/(4*hwin)*(1-r_t^2)^2/r_t * ...
%        (-(alphas+1)/alphas-...
%        ((alphas+1)*s_up^alphas*(1-s_up)^2 - (1-s_up^(alphas+1))^2 )/...
%        (s_up *(1-s_up^alphas)^2));
%  end;
%  
%  yupper = lookup(y_scaling*scale_corr,...
%            y_scaling(filtNos)*scale_corr/maxryfac);
%  ylower = lookup(y_scaling*scale_corr,...
%            y_scaling(1)*scale_corr*maxryfac)+1;
%  for yind = ylower:yupper
%    xlower = ceil(y_scaling(yind)*scale_corr*maxrx*fs);
%    xupper = N-xlower+1;
%    if typemaxcorr
%      mymask(ceil(xlower/20):ceil(xupper/20),yind,rind,:)= ...
%        2/9*(max(((gintendedalt) - ...
%        (gcorr(ceil(xlower/20)*20:20:ceil(xupper/20)*20,yind,rind,:))),0)/...
%        stdgcorr(rind,rind1)).^2;
%    else
%      mymask(ceil(xlower/20):ceil(xupper/20),yind,rind,:)= ...
%        1/9*(abs((gintendedalt) - ...
%        (gcorr(ceil(xlower/20)*20:20:ceil(xupper/20)*20,yind,rind,:)))/...
%        stdgcorr(rind,rind1)).^2;
%    end;
%      
%  end;
%  
%  
%  for rind1=1:nrr1
%    figure(50+(rind-(2*hwin/hstep))*nrr1+rind1)
%    imagesc([0,2], [1,0],mymask(:,:,rind,rind1)', climits = [0, 1])
%  end;
%  
%end;

%Optional saving images

%hf = figure(2);
%filenameW = [writepath,'noiseppp.pdf'];
%print (hf, filenameW)

hf = figure(3);
filenameW = [writepath,filename,'scalfiltered.pdf'];
print (hf, filenameW)

hf = figure(15);
filenameW = [writepath,filename,'scalfilterednoisy.pdf'];
print (hf, filenameW)

hf = figure(4);
filenameW = [writepath,filename,'maskcorr.pdf'];
print (hf, filenameW)

hf = figure(5);
filenameW = [writepath,filename,'maskedscalcorr.pdf'];
print (hf, filenameW)

hf = figure(6);
filenameW = [writepath,filename,'scalcorr.pdf'];
print (hf, filenameW)

hf = figure(7);
filenameW = [writepath,filename,'maskintens.pdf'];
print (hf, filenameW)

hf = figure(8);
filenameW = [writepath,filename,'maskedscalintens.pdf'];
print (hf, filenameW)

hf = figure(9);
filenameW = [writepath,filename,'scalintens.pdf'];
print (hf, filenameW)

hf = figure(10);
filenameW = [writepath,filename,'maskcomb.pdf'];
print (hf, filenameW)

hf = figure(11);
filenameW = [writepath,filename,'maskedscalcomb.pdf'];
print (hf, filenameW)

hf = figure(12);
filenameW = [writepath,filename,'scalcomb.pdf'];
print (hf, filenameW)

