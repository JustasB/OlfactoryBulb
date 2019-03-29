function val = corr_diff(M)
    tempM = uptri_1d(1-corrcoef(M));
    temp = atan(tempM(IX)./tempS(IX))/pi*4-1;
    val = mean(temp);
end