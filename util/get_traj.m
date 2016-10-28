disp('Determining k-space trajectory...')

kdimxy = floor(FOVsubsamp/deltax); %#k-space lines
if ~exist('dt')
  dt = 4e-6;     %sec; sampling period in pulse sequence
end

gambar = 4257;               % gamma/2pi in Hz/T
gam = gambar*2*pi;         % gamma in radians/g

if ~exist('gmax')
  gmax = 4;             %g/cm
end
if ~exist('dgdtmax')
  dgdtmax = 18000;      %g/cm/s; max slew rate. Use 150 for
                        %trapezoid
end

%dgdtmax = 150;

if strcmp(traj,'spiral')

kdimxy = floor(FOVsubsamp/deltax); %#k-space lines
%gmax = 4;             %g/cm
%dgdtmax = 18000;      %g/cm/s; max slew rate. Use 150 for trapezoid

if ~exist('densamp') | ~exist('dentrans')
  disp 'setting densamp = 75, dentrans = 75'
  densamp = 75;
  dentrans = 75;
end
if ~exist('nl')
  disp 'setting nl = 1'
  nl = 1;
end

[g,k,t,s,ddens,NN] = spiralgradlx6(FOVsubsamp,kdimxy,dt,dgdtmax/100,gmax,nl,densamp,dentrans);

g = [real(g(:)) imag(g(:))];

% reverse trajectory

% default to reverse if not specified
if ~exist('forwardspiral')
  forwardspiral = 1;
end

if ~forwardspiral
  g = flipdim(g,1);
end

gtemp = [];
for ii = 1:2
  gtemp(:,ii) = interp1([0:size(g,1)-1],g(:,ii),[0:1/8:size(g,1)-1/8],'spline',0);
end
ktemp = -flipud(cumsum(flipud(gtemp)))*dt/8*gam/2/pi;
k = ktemp([1:8:size(ktemp,1)],:);

end

if strcmp(traj,'flybackEPI')

% generate flyback EPI trajectory
gambar = 4257;               % gamma/2pi in Hz/T
gam = gambar*2*pi;         % gamma in radians/g
NN = [];
% move to starting position
g = -(1+1i)*dotrap(1/2/deltax/4257,gmax,dgdtmax/1,dt);
NN = [NN length(g)];
for ii = 1:kdimxy+1
   % readout
   g = [g dotrap(1/deltax/4257,gmax,dgdtmax/1,dt)]; 
   NN = [NN length(g) - sum(NN)];
   % rewind to beginning in kx, blip up in ky
   if ii ~= kdimxy+1
      g = [g (-1+1i/kdimxy)*dotrap(1/deltax/4257,gmax,dgdtmax/1,dt)];
      NN = [NN length(g) - sum(NN)];
   end
end
% blip to DC
g = [g -(1+1i)*dotrap(1/2/deltax/4257,gmax,dgdtmax/1,dt)];
NN = [NN length(g) - sum(NN)];
% forward integrate to get k
k = 4257*cumsum(g)*dt;
% reverse integrate to get k
krev = 4257*cumsum(fliplr(g))*dt;
kpointsb1 = [];
cumNN = cumsum(NN);
for ii = 2:2:length(NN)
     kpointsb1 = [kpointsb1 krev(cumNN(ii-1)+1:cumNN(ii))];
end

end

if strcmp(traj,'EPI')

% generate standard EPI trajectory
gambar = 4257;               % gamma/2pi in Hz/T
gam = gambar*2*pi;         % gamma in radians/g
% move to starting position
g = -(1+1i)*dotrap(1/2/deltax/4257,gmax,dgdtmax/1,dt);
NN = length(g)
sn = 1;
for ii = 1:kdimxy+1
   % readout
   g = [g sn*dotrap(1/deltax/4257,gmax,dgdtmax/1,dt)]; 
   % rewind to beginning in kx, blip up in ky
   if ii ~= kdimxy+1
      g = [g (1i/kdimxy)*dotrap(1/deltax/4257,gmax,dgdtmax/1,dt)];
   end
   sn = sn * -1; % alternate sign
end
NN = [length(g)-NN NN];
% blip to DC
g = [g -(1+1i)*dotrap(1/2/deltax/4257,gmax,dgdtmax/1,dt)];

kpointsb1 = -4257 * dt * (ones(size(g))*sum(g) - cumsum(g)) ;
kpointsb1 = kpointsb1(NN(2)+1:end-NN(2));
g = g(NN(2)+1:end);

end
