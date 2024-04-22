clear; 
mice = {'zz109'};
task = 'puretone'; 
zz109 = fn_getObjPT_bin40(mice);zz109 = zz109{1};
%%
start = 101; finish = 500; biasStart1 = 201; biasEnd1 = 241;
biasStart2 = 244; biasEnd2 = 379;
biasStart3 = 380; biasEnd3 = 499;
figure; hold on;
fill([biasStart1 biasStart1:biasEnd1 biasEnd1],[0;zz109.behav.bias(biasStart1:biasEnd1);0],...
    fn_wheelColorsPT('R',0.2),'EdgeColor','none')
fill([biasStart2 biasStart2:biasEnd2 biasEnd2],[0;zz109.behav.bias(biasStart2:biasEnd2);0],...
    fn_wheelColorsPT('L',0.2),'EdgeColor','none')
fill([biasStart3 biasStart3:biasEnd3 biasEnd3],[0;zz109.behav.bias(biasStart3:biasEnd3);0],...
    fn_wheelColorsPT('R',0.2),'EdgeColor','none')
plot([start finish],[0 0],'Color',[0.8 0.8 0.8],'LineWidth',2)
plot(start:finish, zz109.behav.bias(start:finish),'LineWidth',2,'Color',fn_wheelColorsPT('bias'));
xlim([start finish]); ylim([-1 1]); yticks(-1:0.5:1); xticks([120 500])

%%
clear; 
mice = {'zz062'};
task = 'puretone'; 
zz062 = fn_getObjPT_bin40(mice);zz062 = zz062{1};

figure; plot( zz062.behav.bias)

%%
start = 15; finish = 500; 
biasStart1 = 15; biasEnd1 = 185;
biasStart2 = 218; biasEnd2 = 333;
biasStart3 = 334; biasEnd3 = 500;

figure; hold on;
fill([biasStart1 biasStart1:biasEnd1 biasEnd1],[0;zz062.behav.bias(biasStart1:biasEnd1);0],...
    fn_wheelColorsPT('R',0.2),'EdgeColor','none')
fill([biasStart2 biasStart2:biasEnd2 biasEnd2],[0;zz062.behav.bias(biasStart2:biasEnd2);0],...
    fn_wheelColorsPT('L',0.2),'EdgeColor','none')
fill([biasStart3 biasStart3:biasEnd3 biasEnd3],[0;zz062.behav.bias(biasStart3:biasEnd3);0],...
    fn_wheelColorsPT('R',0.2),'EdgeColor','none')
plot([start finish],[0 0],'Color',[0.8 0.8 0.8],'LineWidth',2)
plot(start:finish, zz062.behav.bias(start:finish),'LineWidth',2,'Color',fn_wheelColorsPT('bias'));
xlim([start finish]); ylim([-1 1]); yticks(-1:0.5:1); xticks([100 500])


%%
start = 1890; finish = 2633; 
biasStart1 = 1899; biasEnd1 = 2099;
biasStart2 = 2404; biasEnd2 = 2542;


figure; hold on;
fill([biasStart1 biasStart1:biasEnd1 biasEnd1],[0;zz062.behav.bias(biasStart1:biasEnd1);0],...
    fn_wheelColorsPT('R',0.2),'EdgeColor','none')
fill([biasStart2 biasStart2:biasEnd2 biasEnd2],[0;zz062.behav.bias(biasStart2:biasEnd2);0],...
    fn_wheelColorsPT('L',0.2),'EdgeColor','none')

plot([start finish],[0 0],'Color',[0.8 0.8 0.8],'LineWidth',2)
plot(start:finish, zz062.behav.bias(start:finish),'LineWidth',2,'Color',fn_wheelColorsPT('bias'));
xlim([start finish]); ylim([-1 1]); yticks(-1:0.5:1); xticks([2000 2500])

%% HIGH AND LOW BIAS BLOCKS
start = 15; finish = 500; 

figure; hold on;
fill([start start finish finish],[-0.2 0.2 0.2 -0.2],[0.8 0.8 0.8],'EdgeColor','none')

fill([start start finish finish],[0.2 0.4 0.4 0.2],fn_wheelColorsPT('L',0.2),'EdgeColor','none')
fill([start start finish finish],[0.4 1 1 0.4],fn_wheelColorsPT('L',0.4),'EdgeColor','none')

fill([start start finish finish],-[0.2 0.4 0.4 0.2],fn_wheelColorsPT('R',0.2),'EdgeColor','none')
fill([start start finish finish],-[0.4 1 1 0.4],fn_wheelColorsPT('R',0.4),'EdgeColor','none')

plot([start finish],[0 0],'Color',[0.8 0.8 0.8],'LineWidth',2)
plot(start:finish, zz062.behav.bias(start:finish),'LineWidth',2,'Color',fn_wheelColorsPT('bias'));
xlim([start finish]); ylim([-1 1]); yticks(-1:0.5:1); xticks([100 500])



%% HIGH AND LOW BIAS BLOCKS
start = 1890; finish = 2633; 

figure; hold on;
fill([start start finish finish],[-0.2 0.2 0.2 -0.2],[0.8 0.8 0.8],'EdgeColor','none')

fill([start start finish finish],[0.2 0.4 0.4 0.2],fn_wheelColorsPT('L',0.2),'EdgeColor','none')
fill([start start finish finish],[0.4 1 1 0.4],fn_wheelColorsPT('L',0.4),'EdgeColor','none')

fill([start start finish finish],-[0.2 0.4 0.4 0.2],fn_wheelColorsPT('R',0.2),'EdgeColor','none')
fill([start start finish finish],-[0.4 1 1 0.4],fn_wheelColorsPT('R',0.4),'EdgeColor','none')

plot([start finish],[0 0],'Color',[0.8 0.8 0.8],'LineWidth',2)
plot(start:finish, zz062.behav.bias(start:finish),'LineWidth',2,'Color',fn_wheelColorsPT('bias'));
xlim([start finish]); ylim([-1 1]); yticks(-1:0.5:1); xticks([2000 2500])



