function histA = fn_hist_update(histA,currA,param)

histA = histA .* param.lambdaHist + currA .* (q-param.lambdaHist);

end