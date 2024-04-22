clear; 
preTrial = 1500;
soundTime = 2000; 
afterSound = 2500;
totalTime = preTrial + soundTime + afterSound;
timeAxis = (-preTrial:soundTime + afterSound-1)/1000;

baseTemplate = zeros(totalTime,1);

soundOnL = baseTemplate; soundOnL(preTrial+1:preTrial+soundTime) = 1;
soundOnR = baseTemplate; soundOnR(preTrial+1:preTrial+soundTime) = 1;
wheelPos = baseTemplate; wheelPos(1:150) = linspace(0,0.5,150);
wheelPos(151:400) = linspace(0.5,-0.3,250); wheelPos(401:500) = linspace(-0.3,0,100);
wheelPos(preTrial+1801:preTrial+soundTime+100) = linspace(0,1.5,300);
wheelPos(preTrial+soundTime+101:end) = 1.5; 
reward = baseTemplate; reward(preTrial+soundTime+1:preTrial+soundTime+80) = 1;
rewardPort = baseTemplate; rewardPort(preTrial+soundTime+101:preTrial+soundTime+2000+100) = 1;
lick = baseTemplate; lickTime = preTrial+soundTime+300 + (0:150:900);
for i = 1:length(lickTime)
    lick(lickTime(i)+1:lickTime(i)+50) = 1; 
end
lick = lick(1:totalTime);

cellMat = {soundOnL,soundOnR,wheelPos,rewardPort,lick};


figure; subplot(1,2,1); hold on;
fill([0 2 2 0],[-18 -18 4 4],matlabColors(1,0.2),'LineStyle','None');
for i = 1:length(cellMat)
    h = plot(timeAxis,(-i+1)*4+cellMat{i});
    %if i == 1 || i == 2
    set(h,'Color',[0.2 0.2 0.2]); set(h,'LineWidth',1.5);
    %end
end
plot([min(timeAxis) max(timeAxis)],[-7 -7],'--','LineWidth',1,'Color',[0.2 0.2 0.2])
ylim([-18 4]); xlim([-1.5 totalTime/1000-1.5]); xticks(-1:1:5)

subplot(1,2,2); hold on;
fill([0 2 2 0],[-18 -18 4 4],matlabColors(1,0.2),'LineStyle','None');
for i = 1:length(cellMat)
    h = plot(timeAxis,(-i+1)*4+cellMat{i});
    %if i == 1 || i == 2
    set(h,'Color',[0.2 0.2 0.2]); set(h,'LineWidth',1.5);
    %end
end
plot([min(timeAxis) max(timeAxis)],[-7 -7],'--','LineWidth',1,'Color',[0.2 0.2 0.2])
ylim([-18 4]); xlim([2500/1000-1.5 totalTime/1000-1.5]); xticks([1 2 3 4 5])