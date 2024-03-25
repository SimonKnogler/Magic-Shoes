function ISres = ISanalysis(A,B,leveldef,KR)
 
 
% * INPUTS:
% - vector A and B of discrete time onsets
% - leveldef:  how far to go in the farey tree. 
% - KR = option to determine how many ratio are performed
%       > KR = 0; don't know how many ratio
%       > KR = x; there is x ratio performed (reduce the possibilities to x ratios performed)
 
% 
% * OUTPUT = ISres (structure)
% - ISres.CI:           Best circular variance found on 'leveldef' number of Farey level - phase coupling index
% - ISres.CIa:          Best circular variance found for onsets A only 
% - ISres.CImean:       Corresponding mean of the circular variance - phase shift 
% - ISres.CIlevel:      Level of the Farey tree where stronger phase coupling is (the highest if multiple choice)
% - ISres.CIstats:      (2,n) where n is the number of ratio used to compute the RP series leading to CI.
%                             First line: ratio from the most to the less occurring
%                             second line: occurrence (in percentage) of the ratio
% 
% 
% - ISres.RM:           Best Dwell time on the relative phase designed for each Farey level (defined from +- 40 degrees of relative phase differences)
% - ISres.RMa:          Level of the Farey tree where stronger DT is (the higher if multiple choice)
% - ISres.RMstats:      (2,n) where n is the number of ratio used to compute the RP leading to DT.
%                             First line: ratio from the most occurent to
%                             the less occurent
%                             second line: occurrence (in percentage) of the ratio
 
 
 
%% 
 
 
 NumVect = A;
 DenVect = B;
 
 % Put every vector in same dimension
 [l,c] = size(NumVect);
 if c > l;NumVect = NumVect';end;
 [l,c] = size(DenVect);
 if c > l;DenVect = DenVect';end;clear l c
 clear ref tar
 
 
 
 %% determining frequency ratio and phase at specific onset  
TimeOnset = sortrows (cat(1,NumVect,DenVect));
TimeOnset(diff(TimeOnset) == 0) = [];
 
NumPhase = - ones(numel(TimeOnset),1);
DenPhase = - ones(numel(TimeOnset),1);
 
for p = 1:numel(NumVect);
NumPhase(TimeOnset == NumVect(p)) = 0 + (p-1) * (2*pi);
end;
 
for p = 1:numel(DenVect);
DenPhase(TimeOnset == DenVect(p)) = 0 + (p-1) * (2*pi);
end;
 
%% option for knowing numphase rank - time onsets A only
optnum = (NumPhase ~= -1);
 
%% find the phase series for NumVect and DenVect
% fill the phase  == -1 by the real phase value (when possible <=> when a
% cycle is available)
 
DenFill = find(DenPhase == -1);
 
DenFill(DenFill > find(DenPhase > 0,1,'last')) = [];
DenFill(DenFill < find(DenPhase == 0)) = [];
 
    for p = 1:numel(DenFill);
        indrp = find(DenVect <= TimeOnset(DenFill(p)), 1, 'last');
        indra = find(DenVect > TimeOnset(DenFill(p)), 1, 'first');
    
        DenPhase(DenFill(p),:) =  DenPhase(TimeOnset == DenVect(indrp,:),:) ...
            + (2*pi) * (TimeOnset(DenFill(p))-DenVect(indrp))/(DenVect(indra)-DenVect(indrp));
    
        clear indrp indra
    end;
 
 
        % same for num 
        
            NumFill = find(NumPhase == -1);
 
            NumFill(NumFill > find(NumPhase > 0,1,'last')) = [];
            NumFill(NumFill < find(NumPhase == 0)) = [];
 
                for p = 1:numel(NumFill);
                      indrp = find(NumVect <= TimeOnset(NumFill(p)), 1, 'last');
                      indra = find(NumVect > TimeOnset(NumFill(p)), 1, 'first');
    
                      NumPhase(NumFill(p),:) =  NumPhase(TimeOnset == NumVect(indrp,:),:) ...
                         + (2*pi) * (TimeOnset(NumFill(p))-NumVect(indrp))/(NumVect(indra)-NumVect(indrp));
    
                       clear indrp indra
                end;
 
  
  
  DenPhase(DenPhase == -1) = NaN;  
        
  NumPhase(NumPhase == -1) = NaN;  
  
  clear NumFill DenFill
  
                
                
%% Determine the real quotient for each phase value
 
% First determine the interval where I can do that
 
realquotient = ones(numel(TimeOnset),1);
 
 
for p = 1:numel(TimeOnset);
    
    if (isnan(NumPhase(p)) == 1) || (isnan(DenPhase(p)) == 1); % means 
        
        realquotient(p,:) = NaN;
        
    else NumMin = find(NumVect <= TimeOnset(p,:),1,'last');
         DenMin = find(DenVect <= TimeOnset(p,:),1,'last');
         
         %Now possibility that the NumMin is the last Numvect available
         %(same reasonning for DenVect)
         
         if NumMin + 1 > numel(NumVect); NumMin = NumMin - 1;end;
         if DenMin + 1 > numel(DenVect); DenMin = DenMin - 1;end;
 
    
                    NumMax = NumMin + 1;
                    DenMax = DenMin + 1;
                    
         realquotient(p,:) = (DenVect(DenMax)-DenVect(DenMin))/(NumVect(NumMax)-NumVect(NumMin));
 
    
    end;
    
    clear NumMax DenMax NumMin DenMin 
    
end;
 
realquotient(realquotient == 0) = NaN;
realquotient(isinf(realquotient) == 1) = NaN;
 
%% filtering for NaN related to method
% here we have some NaN values that are related to the way we computed the
% different vectors. We can delete them so the following NaN values will only be
% related to the fareyquotient being 1/0 or 0/1 (impossible values)
 
 
 
TimeOnset(isnan(realquotient) == 1) = [];
NumPhase(isnan(realquotient) == 1) = [];
DenPhase(isnan(realquotient) == 1) = [];
 
optnum(isnan(realquotient) == 1) = [];
 
realquotient(isnan(realquotient) == 1) = [];
 
 
%% find Farey quotient per level
 
% Associate to each real quotient a farey ratio at a certain level n
% depending on ranges found in level n+1 
 
 
fareynum = ones(numel(realquotient),leveldef);
fareyden = ones(numel(realquotient),leveldef);
 
 
raj = 1;% 2 4 8 16 32 64 
num = [0 1];
den = [1 1];
 
for lev = 1:leveldef;
    
    % determining numerateur and denominateur of the next level; 
    num_nextlevel = ones(numel(num)+raj,1)';
    den_nextlevel = ones(numel(num)+raj,1)';
 
    num_nextlevel(1:2:end) = num;
    den_nextlevel(1:2:end) = den;
    
    num_nextlevel(2:2:end) = num(1:end-1)+num(2:end);
    den_nextlevel(2:2:end) = den(1:end-1)+den(2:end);
    
    % Do the operation you want on this level
    % percentage of occurrence for the ratio
    
    intervalle = num_nextlevel./den_nextlevel;
    intervalle(1:2:end) = [];
    intervalle  = cat(2,0,intervalle,1);
    
    for p = 1:numel(realquotient);
        
         if isnan(realquotient(p)) == 1;
             fareynum(p,lev) = NaN;
             fareyden(p,lev) = NaN;
             
     
        
         else if realquotient(p) <= 1;%case the actual interval work:
    % the reference frequency is lower than the target frequency
    % find the good interval and report value for farey ratio:
        ind = find(intervalle < realquotient(p),1,'last');
        fareynum(p,lev) = num(ind);
        fareyden(p,lev) = den(ind);
        
        else ind = find(intervalle < (1/realquotient(p)),1,'last');% in case target freq < ref freq 
        fareynum(p,lev) = den(ind);
        fareyden(p,lev) = num(ind);
              end;
              
        clear ind
         end;
    end;
    
    numfinal = num;
    denfinal = den;
    
    % preparing next level
    num = num_nextlevel;
    den = den_nextlevel;
    raj = 2*raj;
     
 end;
    
 
 
 
 
 
fareyratio = cat(2 , numfinal./denfinal , denfinal(end-1:-1:1)./numfinal(end-1:-1:1));
fareyrationum = cat(2, numfinal , denfinal(end-1:-1:1));
fareyratioden = cat(2, denfinal , numfinal(end-1:-1:1));
 
 
fareynum(fareynum == 0) = NaN;
fareyden(fareyden == 0) = NaN;
 
 
fareyquotient = fareynum./fareyden;
 
clear raj num den num_nextlevel den_nextlevel intervalle numfinal denfinal
 
 
 
%% 
%% Determining Individual Coordination Index
%% 
 
%% Definition of the relative phase per Farey level
 
% for each level, we can compute the relative phase following:
% RPab = NUM * DenPhase - DEN * NumPhase
% !!! this only works for non NaN values 
 
 
 
RP = ones(numel(TimeOnset),leveldef); %predefining RP
 
 
% determining the relative phase per level
    
 
for lev = 1:leveldef;
        
        for p = 1:numel(TimeOnset)
    
             if isnan(fareyquotient(p,lev)) == 1; 
        
                RP(p,lev) = NaN;
        
             else RP(p,lev) = fareyden(p,lev) * NumPhase(p,:) - fareynum(p,lev) * DenPhase(p,:);
             
             end;
         
        end;
        
end;
    
% Now we play with the option 'KR':
    
    
    % (1) if we don't know how many ratio (KR = 0) RP stays as it is   
    % (2) if there is x ratios performed (KR = x) we take into account RP for only these x ratios 
 
    
if KR > 0; % option (2)
 
    
            for lev = 1:leveldef;
 
    
    matRatio = fareyquotient(:,lev);
    matRatio(isnan(matRatio) == 1) = [];
    
    
    
    
    
 
    
    if isempty(matRatio) == 1; %fareyquotient vector has only NaN
        
        RP(:,lev) = NaN;
        
    else 
        
        % extract different ratio presents at the level
    RatioPres = sortrows(matRatio);
    RatioPres(diff(RatioPres) == 0) = [];
        
        % determine occurrence for each fareyquotient existing in the
        % fareyquotient vector
        
        ratioOcc = ones(numel(RatioPres),1);
        for p = 1:numel(RatioPres);
              ratioOcc(p,:) = numel(find(matRatio == RatioPres(p)))/numel(matRatio);
        end;
        
        % we rank the different fareyquotients per % of occurrence
 
        MatR = sortrows(cat(2,ratioOcc,RatioPres),1);
        MatR = MatR(end:-1:1,:); %ratios are ordered by occurrence    
        clear ratioOcc
        
        % now we have all fareyquotients existing at the level
        % that are ranked depending on their % of occurrence
         
        
        % we create indexGR = [0 1] matrix with 1 indicating a rank for the
        % correct fareyquotient (the 'hmr' most occurring) and 0 the ranks
        % of the other ones (typically the ones to work on)
        
        % remove the ratio the less occurring
        MatR(KR+1:end,:) = [];
 
        
        
            indexGR = zeros(numel(TimeOnset),1);
            
            
            for p = 1:size(MatR,1);
                indexGR = indexGR + (fareyquotient(:,lev) == MatR(p,2));
            end;
            
            
            % For the ranks where indexGR = 0, we change the RP values depending on how
            % close the realquotient is from the 'KR' fareyquotient
            % considered
            
            index0 = find(indexGR == 0);
            clear indexGR
            
            for p = 1:numel(index0);
                
                if isnan(realquotient(index0(p))) == 1;
                    
                     RP(index0(p),lev) = NaN;
                     fareyquotient(index0(p),lev) = NaN;
                     
                else
                    
                    % determine the closer ratio from the realquotient
                    % among the 'KR' most occurring fareyquotient on the
                    % level
 
                    Grat = find(abs(MatR(:,2) - realquotient(index0(p))) == ...
                                 max(abs(MatR(:,2) - realquotient(index0(p)))),1,'first');
                             
                    closerFR =  MatR(Grat,2);
                    clear Grat
                    
                   
                      fareynumlev = fareyrationum(fareyratio == closerFR);
                      fareydenlev = fareyratioden(fareyratio == closerFR);
                      
                      RP(index0(p),lev) = fareynumlev * DenPhase(index0(p),:) - fareydenlev * NumPhase(index0(p),:); 
                      fareyquotient(index0(p),lev) = closerFR;
                      
                      clear closerFR fareynumlev fareydenlev
 
                end;
                
                
            end;
            
            clear index0 MatR
            
     end;
    
    
   clear matRatio RatioPres
             end;
        
end;
 
 
          
         
 % transform RP within [0 2*pi]
 
         % transform negative into [0 2pi]
              RP(RP < 0) = RP(RP < 0) + (2*pi) * (floor(abs(RP(RP < 0))/(2*pi)) + 1);
 
         % transform RP > 2pi into [0 2pi]
              RP(RP >= 2*pi) = RP(RP >= 2*pi) - (2*pi) * (floor(RP(RP >= 2*pi)/(2*pi)));
 
 
%% Determining CI, CImean, CIlevel and CIstats             
       
CImeanlev = ones(1,leveldef);
CIlev = ones(1,leveldef);
 
 
for lev = 1:leveldef; 
    
    matRP = RP(:,lev);
    matRP(isnan(matRP) == 1) = [];
    
        [CImeanlev(:,lev),CIlev(:,lev)] = circmean( matRP );
 
        % security to make sure that the NaN values related to extreme ratios
        % such as 1/4 1/7 or 7/1 6/1 etc, do not affect the index
         
        CIlev(:,lev) = (numel(matRP)/numel(RP(:,lev))) * CIlev(:,lev);
        
        clear matRP
end;
 
ISres.CI = max(CIlev);
 
ISres.CIlev = CIlev;
 
ISres.CIlevel = find(CIlev == ISres.CI,1,'last');
ISres.CImean = CImeanlev(:,ISres.CIlevel);
 
clear CImeanlev CIlev
 
% We now determine CI Stats
 
    lev = ISres.CIlevel;
    
   
    matRatio = fareyquotient(:,lev);
    matRatio(isnan(matRatio) == 1) = [];
    
    RatioPres = sortrows(matRatio);
    RatioPres(diff(RatioPres) == 0) = [];
    
    
    if isempty(RatioPres) == 1;
        
        ISres.CIstats =  NaN;
        
    else
        
             ratioOcc = ones(numel(RatioPres),1);
             for p = 1:numel(RatioPres);
              ratioOcc(p,:) = numel(find(matRatio == RatioPres(p)))/numel(matRatio);
             end;
    
             matstats = sortrows(cat(2,RatioPres,ratioOcc),2);    
    
             ISres.CIstats =  matstats(end:-1:1,:)';
             
    end;
    
             clear matstats ratioOcc RatioPres matRatio lev
 
             
             
             
%% now we determine CI for A onsets only in the same level
  
    matRPa = RP(optnum == 1, ISres.CIlevel);
    matRPa(isnan(matRPa) == 1) = [];
    
        [~,ISres.CIa] = circmean( matRPa );
 
        % security to make sure that the NaN values related to extreme ratios
        % such as 1/4 1/7 or 7/1 6/1 etc, do not affect the index
         
        ISres.CIa = (numel(matRPa)/numel(RP(optnum == 1, ISres.CIlevel))) * ISres.CIa;
        
        clear matRPa
 
  
             
         clear matRP delta Wn Dn
 
             
 %% now we determine RM for A and B onsets with the improved Return Map method
            
             matRP = RP(:,ISres.CIlevel); 
             matRP(isnan(matRP) == 1) = [];
 
             delta = diff(matRP)/(2*pi)*360; %difference between two RP values 
                                                    % between[-360;+360]
                                                    
             Dn = delta/sqrt(2);% Euclidien Distance from identity line x=y for each points of the return map
             DnRefnormInf = 40/sqrt(2); %this is the Euclidian distance for a phase diff of 40 deg = 0.6981 rad = 0.1111 in [0 1] 
             DnRefnormSup = 320/sqrt(2);
             DnMax = 360/sqrt(2);
             
             %this is the euclidien distance for a phase diff of 40 deg = 0.6981 rad = 0.1111 in [0 1] 
 
             Wn = ones(numel(Dn),1);
             Wn(abs(Dn) <= DnRefnormInf) = 1 - (abs(Dn(abs(Dn) <= DnRefnormInf))/DnRefnormInf);
             Wn(abs(Dn) >= DnRefnormSup) = (abs(Dn(abs(Dn) >= DnRefnormSup)) - DnRefnormSup) / (DnMax - DnRefnormSup);
             Wn(abs(Dn)>DnRefnormInf & abs(Dn)<DnRefnormSup) = 0;
 
             % computing the final phase coupling index
     ISres.RM = sum(Wn)/numel(Wn);
     
     % Return map data
     ISres.ReturnMap.RP = matRP;
     ISres.ReturnMap.Dn = Dn;
     ISres.ReturnMap.Wn = Wn;
             
     
     clear matRP delta Wn Dn
     
     
   %% now we determine RM for A onsets only with the improved Return Map method
   
     
             matRP = RP(optnum == 1, ISres.CIlevel); 
             matRP(isnan(matRP) == 1) = [];
 
             delta = diff(matRP)/(2*pi)*360; %difference between two RP values 
                                                    % between[-360;+360]
                                                    
             Dn = delta/sqrt(2);% Euclidian Distance from identity line x=y for each points of the return map
             DnRefnormInf = 40/sqrt(2); %this is the Euclidian distance for a phase diff of 40 deg = 0.6981 rad = 0.1111 in [0 1] 
             DnRefnormSup = 320/sqrt(2);
             DnMax = 360/sqrt(2);
             
             %this is the euclidien distance for a phase diff of 40 deg = 0.6981 rad = 0.1111 in [0 1] 
 
             Wn = ones(numel(Dn),1);
             Wn(abs(Dn) <= DnRefnormInf) = 1 - (abs(Dn(abs(Dn) <= DnRefnormInf))/DnRefnormInf);
             Wn(abs(Dn) >= DnRefnormSup) = (abs(Dn(abs(Dn) >= DnRefnormSup)) - DnRefnormSup) / (DnMax - DnRefnormSup);
             Wn(abs(Dn)>DnRefnormInf & abs(Dn)<DnRefnormSup) = 0;
 
             % computing the final phase coupling index
     ISres.RMa = sum(Wn)/numel(Wn);
             
             
 
return
 
 
 
 
 
