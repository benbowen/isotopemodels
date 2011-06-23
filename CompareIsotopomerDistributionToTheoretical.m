function F = CompareIsotopomerDistributionToTheoretical(ysim,ytheo,charge,params)
% ytheo is the theoretical distribution
% ysim is the measured distribution
%% this part is screwed up for molecules that are super labeled. Could be problematic to leave out
% xx=find(ytheo(:,1)>(min(ytheo(:,1))+params.massrange/charge));
% ytheo(xx,:)=[];
%% map the m/z values of real onto theoretical. shift relevant peaks in ysim to the m/z value in ytheo if they are within tolerance.
y=zeros(size(ysim,1),1);
y = ysim; %initial y to have all of the values of ysim.  The loop below will shift the m/z values of some peaks.
for i = 1:size(ytheo,1)         %for each peak in ytheo
    xx=find(abs(ysim(:,1)-ytheo(i,1)) < 0.1);        %find the index of peaks in ysim that match to within the tolerance of a ytheo peak
    if ~isempty(xx)    %if there are peaks that match %xx is an index that is between 1 and length(ysim);
        y(xx,1) = ytheo(i,1);   %set the m/z of the matched peaks to be equal to ytheo
        y(xx,2)=sum(ysim(xx,2));
    end
end
[u ua]=unique(y(:,1));
y=y(ua,:);
%% make a "full" ytheo that has all mass values in ysim and a "full" y that has all mass values in ytheo
all_mzs = union(y(:,1),ytheo(:,1));
ytheof = zeros(length(all_mzs),2);
yf = ytheof;
for i = 1:length(all_mzs)
    ytheof(i,1) = all_mzs(i);
    yf(i,1) = all_mzs(i);
    xx = ytheo(:,1) == all_mzs(i);
    zz = y(:,1) == all_mzs(i);
    if isempty(xx)
        ytheof(i,2) = 0;
    else
        ytheof(i,2) = sum(ytheo(xx,2));
    end
    if isempty(zz)
        yf(i,2) = 0;
    else
        yf(i,2) = sum(y(zz,2));
    end
end

%% take top n
[s sx]=sort(yf(:,2),1,'descend');
if length(s)>=params.minnum;
    yf=yf(sx(1:params.minnum),:);
    ytheof=ytheof(sx(1:params.minnum),:);
    all_mzs=all_mzs(sx(1:params.minnum));
    % compute the sse
yf(:,2)=yf(:,2)/sum(yf(:,2));
ytheof(:,2)=ytheof(:,2)/sum(ytheof(:,2));
F=sum((ytheof(:,2)-yf(:,2)).^2);


end
if ~isfinite(F)
    F=1000;
end
F
%% plot
if params.plot==1
    figure(1);
    subplot(2,1,1),
    s1=stem(all_mzs,yf(:,2),'ro-','linewidth',3,'markersize',8);
    hold on
    s2=stem(all_mzs,ytheof(:,2),'ko-','linewidth',2,'markersize',10);
    hold off
    set(s1,'markerfacecolor',[1 0 0])
    xlabel('m/z','fontsize',12,'fontweight','bold');
    ylabel('Intensity (au)','fontsize',12,'fontweight','bold');
    set(gca,'fontsize',12,'fontweight','bold','linewidth',2);
    drawnow
end