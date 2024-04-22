function [bias,acc,dprime,criterion,acc_L,acc_R,ar_L,ar_R] = fn_getAccBias(stimulus,correct,miss)
    if nargin == 2;miss = false(size(correct)); end

    stimulus(stimulus==-1) = 2; % make all the L/R 1/-1 into 1/2

    acc_L = sum(correct & (stimulus==1 | stimulus==3)) / sum((stimulus==1 | stimulus==3) & miss~=1);
    acc_R = sum(correct & (stimulus==2 | stimulus==4)) / sum((stimulus==2 | stimulus==4) & miss~=1);
    bias = (acc_L - acc_R); % / (acc_s1 + acc_s2);
    acc  = (acc_L + acc_R)/2;

    %dprimeAccThre1 = 1-0.25/length(stimulus); dprimeAccThre2 = 0.25/length(stimulus);
    dprimeAccThre1 = 0.95; dprimeAccThre2 = 0.05;
    acc_L(acc_L>dprimeAccThre1) = dprimeAccThre1; acc_R(acc_R>dprimeAccThre1) = dprimeAccThre1;
    acc_L(acc_L<dprimeAccThre2) = dprimeAccThre2; acc_R(acc_R<dprimeAccThre2) = dprimeAccThre2;

    dprime = norminv(acc_L,0,1) - norminv(1-acc_R,0,1);
    criterion = -0.5 * (norminv(acc_L,0,1) + norminv(1-acc_R,0,1));
    ar_L = sum(~miss & (stimulus==1 | stimulus==3)) / sum((stimulus==1 | stimulus==3));
    ar_R = sum(~miss & (stimulus==2 | stimulus==4)) / sum((stimulus==2 | stimulus==4));

end