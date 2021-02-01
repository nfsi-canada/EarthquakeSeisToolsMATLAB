function mct = ph2dt(hyp,phadir,maxD,maxN)
% function mct = ph2dt(hyp,phadir,maxD,maxN)
%
% 2020-10-21
% This function mimics the program 'ph2dt' that forms double-differences
% as input for "hypoDD"
%
%  INPUTS
%
%       hyp == [LAT,LON,DEP]
%    phadir == [EV,STA,PHA,T]
%      maxD == maximum inter-event distance, in km
%              note that DEP must be in km
%      maxN == maximum number of neighbours an event can have
%
% OUTPUTS
%
%       mct == [EV1 EV2 STA PHA DT (CC)]
%              
%


Ne  = size(hyp,1);
jEt = cell(Ne-1,1);

NphM = zeros(Ne);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% PART 1
% -- First only consider inter-event distance
for j1 = 1:Ne-1
    
    j2  = (j1+1:Ne)'; 
    lat1 = hyp(j1,1);
    lon1 = hyp(j1,2);
    dep1 = hyp(j1,3);
    lat2 = hyp(j2,1);
    lon2 = hyp(j2,2);
    dep2 = hyp(j2,3);
    
    dE1 = sqrt(distance(lat1,lon1,lat2,lon2,[6371 0]).^2 +(dep1-dep2).^2);
    jj = find(dE1 < maxD);
    Nj = length(jj);
    
    jEt{j1} = [repmat(j1,Nj,1) j2(jj)];
    
end

jE0 = vertcat(jEt{:});

disp(['PH2DT: Inter-event distances complete.  ',datestr(datetime('now'))])


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% -- PART 2

% -- Select event pairs based on a number of comparisons
% -- Get list of [STA,PHA] picked for each event
phE = cell(Ne,1);
for ie = 1:Ne
    phE{ie} = phadir(phadir(:,1)==ie,2:3);
end



jE = zeros(size(jE0));
jj = 0;

for ie = 1:Ne
    jje = find( (jE0(:,1)==ie)+(jE0(:,2)==ie) );
    Nj = size(jje,1);
    if Nj <= maxN
        jE(jj+[1:Nj],:) = jE0(jje,:);
        jj = jj+Nj;
    else
        % -- If there are > maxN events within maxD, select
        % -- maxN events that have the most comparisons 
        
        % -- Indices of events 'ie' is compared to
        jcmp = jE0(jje,:)';
        jcmp = setdiff(jcmp(:)',ie,'stable');
        
        NphM = zeros(Nj,1);
        for ii = 1:Nj
            NphM(ii) = length(intersect(phE{ie},phE{jcmp(ii)},'rows'));
        end
        [~,jsrt] = sort(NphM,'descend');
        jE(jj+[1:maxN],:) = jE0(jje(jsrt(1:maxN)),:);
        jj = jj+maxN;
    end
end

jE(jj+1:end,:) = [];
jE = unique(jE,'rows');


disp(['PH2DT: Neighbour selection complete.  ',datestr(datetime('now'))])

mct = phadir2dt(phadir,jE);

disp(['PH2DT complete.  ',datestr(datetime('now'))])