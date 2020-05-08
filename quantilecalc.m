%Script to calculate statistical properties of gcorr and rhoest 
% for one realization
x_size=size(gcorr,1);
y_size=size(gcorr,2);

meangcorr=zeros(nrr0,nrr1);
stdgcorrnew=zeros(nrr0,nrr1);
skewgcorr = zeros(nrr0,nrr1);
sqrtgintended=zeros(nrr0,1);
sqrtgintendedalt=zeros(nrr0,1);

interestregion = false(x_size,y_size);

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


myfinalmask=zeros(x_size,y_size);
for rind=2:4
  r_t = hstep*rind-hwin;
  r_low = hstep*rind-2*hwin;
  r_up = hstep*rind;
  s_low = 1-r_low^2;
  s_t = 1-r_t^2;
  s_up = 1-r_up^2;
  gintended(rind) = ((s_t^alphas*(alphas*(1-s_t)-s_t*(1-s_t^alphas))^2 ...
          +(alphas*s_t^alphas*(1-s_t)-1+s_t^alphas)^2)/(1-s_t^alphas)^3);
  if r_low > 0
     gintendedalt(rind) = (1/(4*hwin)*(1-r_t^2)^2/r_t * ...
        (((alphas+1)*s_low^alphas*(1-s_low)^2 - (1-s_low^(alphas+1))^2 )/...
        (s_low *(1-s_low^alphas)^2)-...
        ((alphas+1)*s_up^alphas*(1-s_up)^2 - (1-s_up^(alphas+1))^2 )/...
        (s_up *(1-s_up^alphas)^2)));
  else
     gintendedalt(rind) = (1/(4*hwin)*(1-r_t^2)^2/r_t * ...
        (-(alphas+1)/alphas-...
        ((alphas+1)*s_up^alphas*(1-s_up)^2 - (1-s_up^(alphas+1))^2 )/...
        (s_up *(1-s_up^alphas)^2)));
  end;
  for rind1 = 1:nrr1
    gcorrtemp=gcorr(:,:,rind,rind1);
    meansqgcorr(rind,rind1) = mean(sqrt(gcorrtemp(interestregion &posgint(:,:,1,rind1))));
    stdsqgcorr(rind,rind1) = std(sqrt(gcorrtemp(interestregion &posgint(:,:,1,rind1))));
    skewsqgcorr(rind,rind1) = skewness(sqrt(gcorrtemp(interestregion &posgint(:,:,1,rind1))));
    meangcorr(rind,rind1) = mean((gcorrtemp(interestregion &posgint(:,:,1,rind1))));
    stdgcorrnew(rind,rind1) = std((gcorrtemp(interestregion &posgint(:,:,1,rind1))));
    skewgcorr(rind,rind1) = skewness((gcorrtemp(interestregion &posgint(:,:,1,rind1))));

    myfinalmask(interestregion &posgint(:,:,1,rind1))...
      = myfinalmask(interestregion &posgint(:,:,1,rind1)) + ...
            (abs((gintendedalt(rind)) ...
            - (gcorrtemp(interestregion & posgint(:,:,1,rind1))))...
            /stdgcorr(rind,rind1)).^2;
  end;
end;
myfinalmask = myfinalmask/(nrr0 - 2*hwin/hstep+1)./posgintsum;



q1000gcorr = quantile((myfinalmask(interestregion)),0.999);
q100gcorr = quantile((myfinalmask(interestregion)),0.99);


meanrho=zeros(nrr1,1);
varrhoest=zeros(nrr1,1);
skewrho = zeros(nrr1,1);
rhointendedvec = zeros(nrr1,1);
myvarrhovec = zeros(nrr1,1);

myfinalmask=zeros(x_size,y_size);

for rind=1:nrr1
  rhointendedvec(rind) = (alphas*r_1(rind)^2)/(1-r_1(rind)^2);
  myvarrhovec(rind) = alphas^2*r_1(rind)^4/(2*pi*(1-r_1(rind)^2)^2)*...
        integral(@(t) intinvarrho(t, r_1(rind), alphas), -pi, pi);
  rhoesttemp = rhoest(:,:,rind);
  meanrho(rind) = mean(rhoesttemp(interestregion));
  varrhoest(rind) = var(rhoesttemp(interestregion));
  skewrho(rind) = skewness(rhoesttemp(interestregion));

  rhointended = rhointendedvec(rind);
  myvarrho = myvarrhovec(rind);
  myfinalmask(interestregion) = myfinalmask(interestregion) +  ...
          abs(rhointended - rhoesttemp(interestregion)).^2/myvarrho;
  
end;

myfinalmask = myfinalmask/nrr1;


q1000rhoest = quantile((myfinalmask(interestregion)),0.999);
q100rhoest = quantile((myfinalmask(interestregion)),0.99);


myfinalmask=zeros(x_size,y_size);
for rind = (2*hwin/hstep):nrr0
  for rind1 = 1:nrr1
    gcorrtemp=gcorr(:,:,rind,rind1);
    
    myfinalmask(interestregion &posgint(:,:,1,rind1))...
      = myfinalmask(interestregion &posgint(:,:,1,rind1)) + ...
            (abs((gintendedalt(rind)) ...
            - (gcorrtemp(interestregion & posgint(:,:,1,rind1))))...
            /stdgcorr(rind,rind1)).^2;
  end;

end;
myfinalmask = fact_corr*myfinalmask/(nrr0 - 2*hwin/hstep+1)./posgintsum;

for rind = 1:nrr1
  rhoesttemp = rhoest(:,:,rind);
  rhointended = rhointendedvec(rind);
  myvarrho = myvarrhovec(rind);
  myfinalmask(interestregion) = myfinalmask(interestregion) +  ...
          abs(rhointended - rhoesttemp(interestregion)).^2/myvarrho/nrr1;
  
end;
myfinalmask = myfinalmask/sqrt(2);

q1000comb = quantile((myfinalmask(interestregion)),0.999);
q100comb = quantile((myfinalmask(interestregion)),0.99);

