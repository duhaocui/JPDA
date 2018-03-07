parfor i=1:10
   JPDA_singlesense(i);
end

ERRORS_all=0;
NN=10;
for i=1:NN
    i
    load(['JPDA_satsinglesense_',num2str(i),'.mat'])
    E=sqrt( sum(ERRORS.^2,3)/size(ERRORS,3) )
    
    ERRORS_all=ERRORS_all+E.^2;
end

sqrt(ERRORS_all/NN)
