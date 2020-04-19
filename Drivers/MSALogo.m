%% Multiple sequence alignment 
Query = 'FGFRVFGVVLIIVDIIVVIVDLAISEKKRGIREILEGVSLAIALFFLVDVLMRVFVEGFKNYFRSKLNTLDAVIVVGTLLINMTYSFSDLAATDQMPQMVTLLRVLRIVILIRIFRLASQ';
VSD = fastaread('D:\OneDrive - rice.edu\Francois\Paper\JEDI-2P\Figures\MSA\NCBI_blasp\91ZACGNW01R_seqdump.txt');
VSD(end+1)=struct('Header','ASAPVSD','Sequence',Query);
% % 
dist = seqpdist(VSD,'ScoringMatrix','BLOSUM62');
tree = seqlinkage(dist,'average',VSD)
ma1 = multialign(VSD,tree,'ScoringMatrix',...
                {'pam150','pam200','pam250'});
%%
% [W, ~]=seqlogo2(ma,'alphabet', 'aa')%, 'startAt', 1, 'endAt', 171)
 = 1:23;
S2 = 29:56;
S3 = 66:87;
S4 = 99:120;
StrandsSt = [57, 85, 122, 398];
MSA = cell2mat(extractfield(ma1,'Sequence')');
MSA_align2_ASAP = find(MSA(end,:) ~= '-');
seqlogo2(MSA(:,MSA_align2_ASAP), 'alphabet', 'aa', 'startAt', S1(1), 'endAt', S1(end),'ActualSt',StrandsSt(1))
seqlogo2(MSA(:,MSA_align2_ASAP), 'alphabet', 'aa', 'startAt', S2(1), 'endAt', S2(end),'ActualSt',StrandsSt(2))
seqlogo2(MSA(:,MSA_align2_ASAP), 'alphabet', 'aa', 'startAt', S3(1), 'endAt', S3(end),'ActualSt',StrandsSt(3))
[wm, h] = seqlogo2(MSA(:,MSA_align2_ASAP), 'alphabet', 'aa', 'startAt', S4(1), 'endAt', S4(end),'ActualSt',StrandsSt(4))
figure(1)
hold on
[npos, handle] = SeqLogoFig(MSA(:,MSA_align2_ASAP(S1)), 'alphabet', 'aa','CUTOFF',0);
set(gca,'xtick',[],'ytick',[])
[npos, handle] = SeqLogoFig(MSA(:,MSA_align2_ASAP), 'alphabet', 'aa','CUTOFF',0);
set(gca,'xtick',[],'ytick',[])
[npos, handle] = SeqLogoFig(MSA(:,MSA_align2_ASAP), 'alphabet', 'aa','CUTOFF',0);
set(gca,'xtick',[],'ytick',[])
[npos, handle] = SeqLogoFig(MSA(:,MSA_align2_ASAP), 'alphabet', 'aa','CUTOFF',0);
set(gca,'xtick',[],'ytick',[])
%-----