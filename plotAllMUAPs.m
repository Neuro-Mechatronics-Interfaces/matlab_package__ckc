function plotAllMUAPs(MUAPs,MU_ID,Fsamp,BackColor,FrontColor)
Fontssize1 =10;
YLimMin=10000;
YLimMax=0;
col{1}= FrontColor; %[0.5,0.5,1.0];
col{2}= FrontColor; % [1,1,0];
YLimCorrect = 0.9;

h=figure;
%set(h,'Renderer','OpenGL');
set(h,'PaperPositionMode','auto');
Position = get(0,'ScreenSize');
set(h,'Position',Position);
set(h,'Color',BackColor);

% calculate the number of rows;
switch length(MUAPs)
    case {1,2,3}
        R=1; C=length(MUAPs);
    case {4,5,6}
        R=2; C=ceil(length(MUAPs)/2);
    case {7,8,9}
        R=3; C=3;
    case {10,11,12}
        R=3; C=4;
    case {13,14,15}
        R=3; C=5;
    case {16,17,18}
        R=3; C=6;
    case {19,20}
        R=4; C=5;
    case {21,22,23,24}
        R=4; C=6;
    case {25}
        R=5; C=5;
    case {26,27,28,29,30}
        R=5; C=6;
    case {31,32,33,34,35}
        R=5; C=7;
    case {36,37,38,39,40}
        R=5; C=8;
    otherwise
        R=8; C=8;
end

for MUid = 1:min(R*C,length(MUAPs))
    subplot(R,C,MUid); hold on;
    
    rY=size(MUAPs{MUid},1);
    cY=size(MUAPs{MUid},2);
    
    for r=1:rY
        for c=1:cY
            if ~isempty(MUAPs{MUid}{r,c})
                T=floor(size(MUAPs{MUid}{r,c},2)/2);
            end
        end
    end
    
    tmp = [];
    for r=1:rY
        for c=1:cY
            mMUAPs =[];
            if isempty(MUAPs{MUid}{r,c})
                tmp(end+1) = 0;
            else
                for k1=1:size(MUAPs{MUid}{r,c},1)
                    [mMUAPs(k1), iMUAPs]=max(abs(MUAPs{MUid}{r,c}(k1,:)));
                end
                tmp(end+1) = 1.0*max(mMUAPs)+eps;
            end
            tmp(end)=tmp(end)+eps;
        end
    end
    for r=1:rY
        dh(r) = 1.5*max(tmp);
    end
    yLineHight = round(max(tmp)*100000)/100;
    if yLineHight > 10
        yLineHight = floor(max(tmp)/10)*10;
    end
    xLineWidth = round(0.02*Fsamp);
    
    for r=1:rY
        for c=1:cY
            if ~isempty(MUAPs{MUid}{r,c})
               
                mMUAPs =[];
                if size(MUAPs{MUid}{r,c},1) > 1
                    for k1=1:size(MUAPs{MUid}{r,c},1)
                        plot((c-1)*T*1.1+[1:T],MUAPs{MUid}{r,c}(k1,:).'+sum(dh(1:r)),'Color',[0.6,0.6,0.6]);
                    end
                end
                plot((c-1)*(2*T+1)*1.1+[1:(2*T+1)],mean(MUAPs{MUid}{r,c},1).'+sum(dh(1:r)),'LineWidth',0.5,'Color',[0.1,0.1,0.7]);
                %plot3((c-1)*(2*T+1)*1.1+[1:(2*T+1)],sum(dh(1:r))+0*[1:(2*T+1)],mean(MUAPs{MUid}{r,c},1).','LineWidth',0.5,'Color',col{1}); view(-15,75);
                [mMUAPs(k1), iMUAPs]=max(abs(mean(MUAPs{MUid}{r,c},1)));
            end
            
            YLimMax=max([YLimMax,max(mMUAPs)+sum(dh(1:r))]);
            YLimMin=min([YLimMin,min(mMUAPs)+sum(dh(1:r))]);
        end
    end
    
    plot((1)*(2*T+1)*1.1+[0,0],-dh(1)+[0,yLineHight],'Color',FrontColor,'LineWidth',1);
    plot((1)*(2*T+1)*1.1-[-T/3,+T/3],-dh(1) + [0,0],'Color',FrontColor,'LineWidth',1);
    plot((1)*(2*T+1)*1.1-[-T/3,+T/3],-dh(1) + [yLineHight,yLineHight],'Color',FrontColor,'LineWidth',1);
    text((1)*(2*T+1)*1.15+T/3,-dh(1)+yLineHight/2,[num2str(yLineHight) ''],'VerticalAlignment','middle','HorizontalAlignment','left','Color',FrontColor,'FontSize',Fontssize1);
    % \muV
    axis tight;
    hAx=get(h,'CurrentAxes');
    set(hAx,'Color',BackColor);
    set(hAx,'XColor',BackColor);
    set(hAx,'YColor',BackColor);
    set(hAx,'ZColor',BackColor);
    set(hAx,'XTickLabel',[]);
    set(hAx,'YTickLabel',[]);
    title(['MU ' MU_ID{MUid}],'Color',FrontColor,'FontSize',Fontssize1)
end
