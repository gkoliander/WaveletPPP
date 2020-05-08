function matL = locminima(matA)
  
  matL = zeros(rows(matA), columns(matA));
  
  matright = matA(2:end-1,2:end-1)<matA(3:end,2:end-1);
  matleft = matA(2:end-1,2:end-1)<matA(1:end-2,2:end-1);
  mattop = matA(2:end-1, 2:end-1)<matA(2:end-1, 3:end);
  matbottom = matA(2:end-1, 2:end-1)<matA(2:end-1, 1:end-2);
  matL(2:end-1,2:end-1) = (matright + matleft + mattop + matbottom)>3;
  

endfunction
