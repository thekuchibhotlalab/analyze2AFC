function [bias,acc,dprime, criterion,acc1,acc2] = fn_getAccBiasByCount(mat)
    acc1 = mat(:,1)./sum(mat(:,1:2),2); acc2 = mat(:,4)./sum(mat(:,4:5),2);
    acc = (acc1+acc2)/2; bias = acc1-acc2;
    
    dprimeAccThreHigh = 0.95; dprimeAccThreLow = 0.05;
    acc1(acc1>dprimeAccThreHigh) = dprimeAccThreHigh; acc2(acc2>dprimeAccThreHigh) = dprimeAccThreHigh;
    acc1(acc1<dprimeAccThreLow) = dprimeAccThreLow; acc2(acc2<dprimeAccThreLow) = dprimeAccThreLow;

    dprime = norminv(acc1,0,1) - norminv(1-acc2,0,1);
    criterion = -0.5 * (norminv(acc1,0,1) + norminv(1-acc2,0,1));
end