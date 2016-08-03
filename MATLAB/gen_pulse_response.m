function [pr,t,tstep] = gen_pulse_response(tf,f,maxfreq,pw)
% parameters
if nargin<3
    maxfreq=400; % desired bandwidth in gigahertz
end
if nargin<4
    pw=0.1; % pulse width in nanoseconds
end
% input check
f=f*1e-9;
if (length(tf)~=length(f))
    error('unequal input lengths');
end
if f(1)~=0
    error('no DC response')
end
fdiff=diff(f);
fstep=f(2)-f(1);
if max(fdiff)-min(fdiff)>1e-12
    warning('unequal frequency steps --> downsampling');
    tf=tf(rem(f,fstep)==0);
    f=f(rem(f,fstep)==0);
    fdiff = diff(f);
    if max(fdiff)-min(fdiff)>1e-12
        error('irregular frequency steps')
    end
end
%
fstop=maxfreq;
if rem(fstop/2,fstep)~=0
    error('half maximum frequency not evenly divisible by frequency step');
end
tstep=1/fstop;
% generate modified tf
tf1 = zeros(fstop/2/fstep,1);
tf1(1:length(tf)) = tf;
tf2 = [tf1;conj(flipud(tf1(2:end-1)))];
% convert to impulse response
ir = real(ifft(tf2));
% generate pulse
if rem(pw,tstep)~= 0 
    error('pulse width is not evenly divisible by the time step')
end
pu = ones(pw/tstep,1);
pr = conv(pu,ir);
t = 0:length(pr)-1; %linspace(0,(length(pr)-1)*tstep,length(pr));