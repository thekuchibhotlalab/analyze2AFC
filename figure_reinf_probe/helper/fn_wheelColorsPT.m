function color = fn_wheelColorsPT(colorName,colorAlpha)
    
colorNames = {'Reinf','Probe','L','R','correct','incorrect','bias','ITI'};
colors = {[0.25 0.25 0.25],[0.8000    0.4549         0],[0.7843    0.3098    0.3176],[0.3725    0.3882    0.6471],...
     [62 108 35]/255,[240 25 64]/255, [0.5906    0.3400    0.6400],[0.6 0.6 0.6]};
color = colors{strcmpi(colorNames,colorName)};

if exist('colorAlpha')
    color = fn_colorAlpha(color,colorAlpha);   
end

end