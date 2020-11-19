function phadir = read_PhaseDAT(f,stadir)
% function phadir = read_PhaseDAT(f,stadir)
%
% 2020-10-09
% This function takes in a path to a hypoDD 'phase.dat' formatted file
% and returns results in a cell array "phadir". It requires a station
% directory (probably read in from a "station.dat" file
%
%   INPUTS
%
%        f == path to "phase.dat"-style file
%   stadir == {STA,LAT,LON,(ELV)}
%
%   OUTPUTS
%
%     phadir == [EV#,STA#,PHA,T,WGHT]
%
% 
fph = fopen(f);    

txtPH = textscan(fph,'%s %f %f %s %s %s %s %s %s %s %s %s %s %s %f');
nline = size(txtPH{2},1);
ievln = cellstrfind(txtPH{1},'#');
iphl  = setdiff([1:nline],ievln);
nevln = length(ievln);
nph   = nline-nevln;

% -- phase matrix all numeric [EV, STA, PHS, T, WGHT]
% -- EVs, STA have indices from stadir, evdir, for PHS P == 1, S == 2.
phadir = zeros(nph,5);

% -- T, WGHT can be set before loop
phadir(:,4) = txtPH{2}(iphl);
phadir(:,5) = txtPH{3}(iphl);

% -- Select P or S
indP = cellstrfind(txtPH{4}(iphl),'P');
indS = cellstrfind(txtPH{4}(iphl),'S');
phadir(indP,3) = 1;
phadir(indS,3) = 2;

% -- Get station indices, if station can't be found just remove line
for iph = nph:-1:1
    try
        phadir(iph,2) = cellstrfind(stadir(:,1),txtPH{1}{iphl(iph)});
        
    catch
        phadir(iph,:) = [];
        iphl(iph)    = [];
    end
end

% -- Write event indices to phase
% -- Add fake last index of new event pair so that find statement in loop works
ievln(end+1) = nline+1;
for ie = 1:nevln
    ev = txtPH{15}(ievln(ie));
    i1 = find( (iphl > ievln(ie)).*(iphl < ievln(ie+1)) );

    phadir(i1,1) = ev;

end
