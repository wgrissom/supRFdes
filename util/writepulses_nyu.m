function writepulses_nyu(rf,g,fname,NominalFlipAngle,PulseName,Comment,InitialPhase,Asymmetry)

if ~exist('PulseName','var')
  PulseName = 'pTx';
end
if ~exist('Comment','var')
  Comment = 'pTx pulse';
end
Oversampling = '1 #BS'; % WAG: This seems to be how much finer the RF dt is than the grad dt
if ~exist('InitialPhase','var')
  InitialPhase = 0;
end
if ~exist('Asymmetry','var')
  Asymmetry = 0.5;
end


NUsedChannels = size(rf,2);
MaxAbsRF = max(abs(rf(:)));
Samples = size(rf,1);

fp = fopen(fname,'w');

% write header
fprintf(fp,'[pTXPulse]\n\n');
fprintf(fp,'NUsedChannels    = %d\n',NUsedChannels);
fprintf(fp,'MaxAbsRF         = %0.2f\n',MaxAbsRF);
fprintf(fp,'Samples          = %d\n',Samples);
fprintf(fp,'PulseName        = %s\n',PulseName);
fprintf(fp,'Comment          = %s\n',Comment);
fprintf(fp,'Oversampling     = %s\n',Oversampling);
fprintf(fp,'NominalFlipAngle = %d\n',NominalFlipAngle);
fprintf(fp,'InitialPhase     = %d\n',InitialPhase);
fprintf(fp,'Asymmetry        = %0.1f\n\n',Asymmetry);

% write gradient waveforms
fprintf(fp,'[Gradient]\n\n');
for ii = 1:Samples
  fprintf(fp,'G[%d]=\t%g\t%g\t%g\n',ii-1,...
          round(g(ii,1)*100)/100,round(g(ii,2)*100)/100,round(g(ii,3)*100)/100);
end

% write RF waveforms
rf = rf./max(abs(rf(:)));
for ii = 1:NUsedChannels
  fprintf(fp,'\n[pTXPulse_ch%d]\n\n',ii-1);
  for jj = 1:Samples
    fprintf(fp,'RF[%d]=\t%g\t%g\n',jj-1,...
            round(abs(rf(jj,ii))*1000)/1000,round((angle(rf(jj,ii))+pi)*1000)/1000);
  end
end

fclose(fp);