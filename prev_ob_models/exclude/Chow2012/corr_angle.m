function val = corr_angle(Scorr,Mcorr)
    tempS = 1-Scorr;
    tempM = 1-Mcorr;
    tempS(find(tempS<=0)) = NaN;
    val = atan(tempM./tempS)/pi*4-1;
end