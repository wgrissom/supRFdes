function [g,k,t,s,dens,NN]=spiralgradlx6(opfov,opxres,gts,gslew,gamp,nl,densamp,dentrans);


% variable density paramters
if nargin == 5
densamp = 75;		% duration of full density sampling (# of samples)
dentrans = 75;		% duration of transition from higher to lower
                        % should be >= densamp/2
nl = 1.0;		% degree of undersampling outer part (any
                        % real #)
end
%     densamp = 200;
%     dentrans = 200;
%     nl = 3.0;

fsgcm = gamp;		% fullscale g/cm
risetime = gamp/gslew*10000;  	% us
ts = gts;		% sampling time
gts = gts;		% gradient sampling time
targetk = opxres/2;
A = 32766;		% output scaling of waveform (fullscale)

MAXDECRATIO = 32;
GAM = 4257.0;
S = (gts/1e-6)*A/risetime;
dr = ts/gts;
OMF = 2.0*pi * opfov/(1/(GAM*fsgcm*gts));
OM = 2.0*pi/nl * opfov/(1/(GAM*fsgcm*gts));
distance = 1.0 / (opfov*GAM*fsgcm*gts/A);
i = sqrt(-1);

ac = A;
loop = 1;
absk = 0;
decratio = 1;
S0 = gslew*100; 
ggx = [];
ggy = [];
dens = [];

while loop > 0,		% start over
	loop = 0;
	om = OM/decratio; 
	omf = OMF/decratio; 
	s = S/decratio;
	g0 = 0; 
	gx = g0; 
	gy = 0; 
	absg = abs(g0);
	oldkx = 0; 
	oldky = 0; 
	tkx = gx; 
	tky = gy; 
	kxt = tkx; 
	kyt = tky;
	thetan_1 = 0; 
	taun = 0; 
	n = 0;
	den1 = 0;
	while absk < targetk
	    realn = n/decratio;
	    if rem(realn,200) == 1
	       realn;
	    end
	    taun_1 = taun; 
	    taun = abs(tkx + i*tky)/A; 
	    tauhat = taun;
	    if realn > densamp
	      if den1 == 0 
	    kdenrad = abs(tkx + i*tky)/distance/decratio;
		den1 = 1;
	      end
	      if realn > densamp+dentrans
		scoffset = scthat;
		denoffset = taun_1;
	        scthat = scoffset + om*(tauhat - denoffset);
		fractrans = 1;
	      else
		scoffset = scthat; 
		denoffset = taun_1;
		fractrans = (realn - densamp)/dentrans;
		fractrans = 1 - ( (fractrans-1)*(fractrans-1));
	        scthat = scoffset + (omf + (om-omf)*fractrans)*(tauhat - denoffset);
	      end
	    else
 	      fractrans = 0;
	      scthat = omf*tauhat;
	    end
	    theta = atan2(scthat,1.0)+scthat;
	    if absg < ac
		deltheta = theta-thetan_1;
		B = 1.0/(1.0+tan(deltheta)*tan(deltheta));
		gtilde = absg;
		t1 = s*s;
		t2 = gtilde*gtilde*(1-B);
		if t2 > t1
		decratio = decratio * 2.0;
		    if decratio > MAXDECRATIO
			printf('k-space calculation failed.\n');
			return;
		    end
		    loop = 1;
		    break;
		end
		t3 = sqrt(t1-t2);
		absg = sqrt(B)*gtilde+t3;
		if (absg > ac)
		    absg = ac;
	        end
	    end 
	    tgx = absg*cos(theta);
	    tgy = absg*sin(theta);
	    tkx = tkx + tgx;
	    tky = tky + tgy;
	    thetan_1=theta;
	    if rem(n,decratio) == 0
		m = round(n/decratio);
		gx = round(((tkx-oldkx))/decratio);
		gx = gx - rem(gx,2);
		gy = round(((tky-oldky))/decratio);
		gy = gy - rem(gy,2);
		ggx(m+1) = gx;
		ggy(m+1) = gy;
		kxt = kxt + gx;
		kyt = kyt + gy;
		oldkx = tkx;
		oldky = tky;
	    	if rem(m,dr) == 0
		    m  = m / dr;
		    absk = abs(kxt + i*kyt)/distance;
		    dens(m+1) = omf/(omf + (om-omf)*fractrans);
		    if absk > targetk
			break;
		    end
		    kx(m+1) = kxt/distance;
		    ky(m+1) = kyt/distance;
		end
	    end
	    n = n+1;
	end %while
end
		npts = m+1;
		arraysize = (2*nl*om*taun/(2*pi));
res = opfov/arraysize*10;
gvec = (ggx + i.*ggy)./A.*fsgcm;
% svec = diff(gvec)/gts/1000;
% plot((1:npts-1).*gts.*1000,abs(svec),(1:npts).*gts.*1000,abs(gvec))
dt = gts*1000;
delk = 1/4.258/opfov; 	% (g ms)/cm

g = gvec;

% ramp down

l2 = length(g);
rsteps = ceil(abs(g(l2))/(S0*0.99)/gts);
ind3 = l2 + [1:rsteps];
g(ind3) = g(l2).*[rsteps-1:-1:0]./rsteps;
dens(ind3) = ones(size(ind3));
l4 = l2 + rsteps;

% rewinder

skx = real(sum(g));
sky = imag(sum(g));
if 0
  rstep2 = ceil(sqrt(abs(skx + i*sky)/(S0*0.99*gts)));
  tri = [(1:rstep2) (rstep2-1:-1:0)];
  stri = sum(tri);
  ind4 = l4 + [1:2*rstep2];
  g(ind4) = -(tri/stri*skx + i.*tri/stri*sky);
  dens(ind4) = ones(size(ind4));

  l5 = l4 + 2*rstep2;
end


rewx = dotrap(abs(real(sum(g)))*gts,gamp,gslew*50,gts);
rewy = dotrap(abs(imag(sum(g)))*gts,gamp,gslew*50,gts);
if length(rewx) > length(rewy)
  % use rewx and scale it to get appropriate rewy
  g(end+1:end+length(rewx)) = -sign(real(sum(g)))*rewx - ...
      sign(imag(sum(g))) * 1i * abs(imag(sum(g))/ ...
                                  real(sum(g)))*rewx;
else
  g(end+1:end+length(rewy)) = -sign(real(sum(g)))*abs(real(sum(g))/imag(sum(g)))*rewy - ...
                                1i*sign(imag(sum(g)))*rewy;
end

%k = cumsum(g).*dt/delk/opfov;
gtemp = interp(g,8);
ktemp = cumsum(gtemp).*dt/delk/opfov/8;
k = ktemp([1:8:length(ktemp)]);

t = [0:length(g)].*gts;
s = diff(g)./(gts*1000);  %slew rate vector
if 0
  NN = [l2 (rsteps+2*rstep2)];
end
NN = [l2 rsteps+max(length(rewx),length(rewy))];
