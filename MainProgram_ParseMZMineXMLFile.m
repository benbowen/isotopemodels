clear all
close all
clc
%%
addpath(genpath('/Users/bpb/matlabtoolbox/'))
%%
[pathname filename]=myfileselector('xml')
%%
for iii= 1:length(filename)
    % load the peaklist
    [features fname]=parse_mzmine_xmlfile([pathname{iii} filename{iii}]);
    % load the raw datafile
    d=Bmzxmlread([pathname{iii} fname]);
    for j= 1:length(features)
        features(j).data.rt=zeros(size(features(j).data.scan_id));
        for i = 1:length(features(j).data.rt)
            features(j).data.rt(i)=str2double(d.scan(features(j).data.scan_id(i)).retentionTime(3:end-1));
        end
        xx=find(features(j).data.intensity~=0);
        [m mx]=max(features(j).data.intensity(xx));
        features(j).centroidmz=sum(features(j).data.intensity(xx).*features(j).data.mz(xx))/sum(features(j).data.intensity(xx));
        features(j).centroidrt=sum(features(j).data.intensity(xx).*features(j).data.rt(xx))/sum(features(j).data.intensity(xx));
        features(j).maximamz=features(j).data.mz(xx(mx));
        features(j).maximart=features(j).data.rt(xx(mx));
        features(j).maxmz=max(features(j).data.mz(xx));
        features(j).maxrt=max(features(j).data.rt(xx));
        features(j).minmz=min(features(j).data.mz(xx));
        features(j).minrt=min(features(j).data.rt(xx));
        features(j).numpoints=length(xx);
        if ~isempty(strfind(fname,'pos'))
        features(j).polarity='+';
        else
            features(j).polarity='-';
        end
    end
    Data(iii).features=features;
    Data(iii).peaklist=filename{iii};
    Data(iii).rawfile=fname;
    Data(iii).pathname=pathname{iii};
    iii
end
save DarkCrustDataForAtlas Data