function mct = phadir2dt(phadir,hyp,jE)
% function mct = phadir2dt(phadir,hyp,jE)
%
% 2020-05-05
% This function makes a list of "catalog" differential times
% from a list of picks. The differential times (DT) are T1a-T1b, where 
% "1" is a station number and a/b are events. 
%
% Need to add in a way to limit number of "neighbors"...
%
% INPUTS
%
%     phadir = [EV,STA,PHA,T]
%        hyp = hypocenters [LAT,LON,DEP (km)]
%        jE  = list of event pairs to use [EV1,EV2]
%     %% maxD = maximum hypocentral distance, in km
%     %%  Nn = maximum number of neighbors
%     %%       only applies if more than 200 events are picked at the 
%     %%       given station
%
% OUTPUTS
%
%     mct = [EVa,EVb,STA,PHA,DT,(CC)]
%            note that the CC column is just zero in this case.

%if nargin < 4
%    if Ne > 500
%        Nn = 10;
%    else
%        Nn = Inf;
%   end
%end



% -- Initial cull of mct based inter-event distance
%Ne  = size(hyp,1);
%{
lat = repmat(hyp(:,1),1,Ne);
lon = repmat(hyp(:,2),1,Ne);
dep = repmat(hyp(:,3),1,Ne);
dEE = sqrt( distance(lat,lon,lat',lon',[6371,0]).^2 +(dep-dep').^2 );

% -- Pre-determine event combos within maximum distance
jj = find(dEE < maxD);
[j1,j2] = ind2sub([Ne Ne],jj);
jE = unique(sort([j1 j2],2),'rows');
jE(jE(:,1)==jE(:,2),:) = [];
%}
% -- Might need to further cull event combos...restrict maximum
% -- number of neighbors...


mct = zeros(1e7,6);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% -- First deal with P-waves
A = phadir(phadir(:,3)==1,:);

usta = unique(A(:,2));
jj   = 0;

for jss = usta'
    
    Asta = A(A(:,2)==jss,:);
    Ne   = size(Asta,1);
    
    if Ne < 2
        continue
    end
   
    [sw1,j1] = ismember(jE(:,1),Asta(:,1));
    [sw2,j2] = ismember(jE(:,2),Asta(:,1));
    jb = find(sw1.*sw2);
    N = size(jb,1);
    
    if N 
        dt = Asta(j1(jb),4)-Asta(j2(jb),4);
        mct(jj+[1:N]',:) = [Asta(j1(jb),1),Asta(j2(jb),1),repmat([jss,1],N,1),dt,ones(N,1)];
        jj = jj + N;
    end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% -- Repeat for S-waves
A = phadir(phadir(:,3)==2,:);

usta = unique(A(:,2));

for jss = usta'
    
    Asta = A(A(:,2)==jss,:);
    Ne   = size(Asta,1);
    
    if Ne < 2
        continue
    end
   
    [sw1,j1] = ismember(jE(:,1),Asta(:,1));
    [sw2,j2] = ismember(jE(:,2),Asta(:,1));
    jb = find(sw1.*sw2);
    N = size(jb,1);
    
    if N 
        dt = Asta(j1(jb),4)-Asta(j2(jb),4);
        mct(jj+[1:N]',:) = [Asta(j1(jb),1),Asta(j2(jb),1),repmat([jss,2],N,1),dt,ones(N,1)];
        jj = jj + N;
    end
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% -- Sort and cull

mct(jj+1:end,:) = [];
mct = sortrows(mct);



%mct(dEE(Ne*mct(:,1)+mct(:,2)) > maxD,:) = [];

% -- Should I also cull based on redundancy??? /self-consistency


%je   = nchoosek(1:Ne,2); 
%je   = je(dEE(Ne*(je(:,1)-1)+je(:,2)) < maxD,:);

 %dt = Asta(je(:,1),4)-Asta(je(:,2),4);
 %     mct(jj+[1:N]',:) = [Asta(je(:,1),1),Asta(je(:,2),1),repmat([jss,1],N,1),dt,ones(N,1)];
     