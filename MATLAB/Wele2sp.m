function [Snp]=Wele2sp(rlgcfile, L, f);
%
% function [Snp]=Wele2sp(rlgcfile, L, f)
%
% Converts w-elements to scattering parameters (causal)
%
% Inputs:   rlgcfile - path to rlgc file representing a lossy transmissoin line( w-element)
%               L    - length of the transmission line
%               f    - freqeucny range of interest
%
% Outputs:      s    - the sparameter in 'touchstone format'
%
% William Peters (william.r.peters@intel.com) 3/8/2005
% w2s('z50.rlc', 0.1411, 10e6:10e6:1e9)
% w2s('z50.rlc', 0.1411, 25e6:25e6:5.12e9)

freq_ref = 1.0e9;       % Reference frequency is set to 1 GHz
[Lpl,Cpl,Rpl,Gpl] = W_conversion(rlgcfile,f,freq_ref);
f3d=zeros(size(Rpl,1),size(Rpl,2),length(f));
for nn = 1:size(Rpl,1)
    for mm = 1:size(Rpl,2)
        f3d(mm,nn,:)=f;
    end
end
zc = Rpl + i*2*pi*f3d.*Lpl;
yc = Gpl + i*2*pi*f3d.*Cpl;
y=zcyc2y(zc,yc,f,L);
s=y2s(y,50);
Snp= snpexport(s);
return;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    function [y]=zcyc2y(zc,yc,f,L)

[m,n,k]=size(zc);       % m,n should be equal
if(k~=length(f));       
    error('number of samples do not match impedance matrix data');
    return;
end
for i=1:k
    z0=zc(:,:,i);
    y0=yc(:,:,i);
    [s,t2]=eig(z0*y0);
    t=diag( sqrt(diag(t2)) );
    e1=diag( coth( diag(t)*L) );
    e2=diag( -1*csch( diag(t)*L) );
    y11=z0\s*t*e1/s;
    y12=z0\s*t*e2/s;
    y0=[y11,y12; y12, y11];
    y(i,:)=[f(i),reshape(y0,1,4*m*n)];
end;
return;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [s]=y2s(y,z0)

if( nargin < 2)
    z0=50;
end

np=floor(sqrt(length(y(1,:))));      
[m,n]= size(y);
if( n ~=(np*np+1) )
    error('wrong y paramter block');
end

z0=z0(:);
if( length(z0) ~=np )     
    z0=z0(1,1);                         
else 
    z0=z0*ones(1,np).*eye(np,np);     
end

ii=eye(np,np);
k=sqrt(z0)*ii;

for j=1:m
    yy=(reshape(y(j,2:end), np, np)).';  
     ss=2*ii/(ii+k*yy*k) - ii;
    s(j,:)=[y(j,1),reshape(ss.', 1, np*np)];
end
return;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [Snp] = snpexport( s )
if nargin ~= 1
    error('usage: snpexport(tsfile, s )');
end

np=floor(sqrt(length(s(1,:))));          % number of ports
[m,n]= size(s);
if( n ~=(np*np+1) )
    error('wrong s paramter block'); 
end

for i=1:m
   for j=1:np
      for k=1:np
         jdx=(j-1)*np+k+1;  
         Snp(j,k,i)=s(i,jdx);
      end
   end
end

return;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

 function [L_matrix, C_matrix, R_matrix, G_matrix] = W_conversion(rlgc_filename,freq,freq_ref);

% W_conversion function converts conventional W-element model to tabular
% W-element models for frequencies specified in the input
%   
%   [l,c,r,g] = W_conversion(rlgc_filename,freq,freq_ref)
%   Inputs:
%       * rlgc_filename - path to rlgc file representing a lossy transmissoin
%         line( w-element) The W-element file has units in meters
%       * freq - frequency points for W-element to be calculated
%       * freq_ref - reference frequency
%   Outputs:
%       * L_matrix, C_matrix, R_matrix and G_matrix are num_cond by
%       num_cond by N_freq matrices
%
% Author: Chi-te Chen (chi-te.chen@intel.com)

[lo,co,ro,go,rs,gd,num_cond] = readrlgc2(rlgc_filename);

% total number of frequency points
N = length(freq);

% ----------- Other fixed parameters in this version of the script ------- % 
% Surface roughness
sr_rms = 0;         % This version does not consider surface roughness

% copper conductivity
conductivity = 5.8e7;

% Speed of light 
cvel = 2.99e8;     % Speed of light in s/m

% Units of the RLGC parameters
unit_in = 'meter';        % This version does not include unit conversion
unit_out = 'meter';       % The units for in and out will be the same

% DC dielectric loss parameter, set to 0
sigma = 0;
% model valid from 10^m1
m1 = 1;
% model valid upto 10^m2
m2 = 16;

% ------------------ Start main calculation ----------------------- %
N1=(num_cond+1)*num_cond/2;

% Adjusting data containing zeros to 1e-15;
lo = nearzero(lo);
co = nearzero(co);
ro = nearzero(ro);
go = nearzero(go);
rs = nearzero(rs);
gd = nearzero(gd);

% Calculate RLGC matrices at the reference frequency
l_ref = lo';
c_ref = co';
r_ref = rs'*sqrt(freq_ref);
g_ref = gd'*freq_ref;
ro_ref = ro';

% Initialzation
%cc index = zeros(N,1);       % No need
omega = zeros(N,1);
factor = zeros(N,1);
er_double_prime = zeros(N,1);
er_prime = zeros(N,1);
tand = zeros(N,1);

l = zeros(N, N1);
c = zeros(N, N1);
r = zeros(N, N1);
g = zeros(N, N1);

% update 12/2004
% compute effective er and loss tangent from input model
% convert lower half matrices into full matrices
c_full = zeros(num_cond,num_cond);
l_full = zeros(num_cond,num_cond);
gd_full = zeros(num_cond,num_cond);
kn = 0;

for in=1:num_cond
    for jn=1:in
        kn=kn+1;
        c_full(in,jn) = c_ref(kn);
        c_full(jn,in) = c_ref(kn);
        l_full(in,jn) = l_ref(kn);
        l_full(jn,in) = l_ref(kn);
        gd_full(in,jn) = gd(kn);
        gd_full(jn,in) = gd(kn);
    end
end

% Compute c_fs
c_fs = zeros(num_cond,num_cond);
c0 = cvel;
% tao check unit input in inches
c_fs = 1/c0^2*inv(l_full);

% compute er_eff and tand_eff
er_effn = zeros(num_cond,num_cond);
tand_effn = zeros(num_cond,num_cond);
for in = 1:num_cond
    for jn=1:num_cond
        er_effn(in,jn) = c_full(in,jn)/c_fs(in,jn);
        tand_effn(in,jn) = gd_full(in,jn)/c_full(in,jn)/2/3.1416;
    end
end

er_eff = 0;
er_eff_sum = 0;
tand_eff = 0;
tand_eff_sum = 0;

for in=1:num_cond,
    er_eff_sum=er_eff_sum+er_effn(in,in);
    tand_eff_sum=tand_eff_sum+tand_effn(in,in);
end
er_eff=er_eff_sum/num_cond;
tand_eff=tand_eff_sum/num_cond;

% compute Debye model parameters from er_eff and tand_eff
er_slope=er_eff*tand_eff*log(10)*2/3.1416;
er_inf=er_eff-er_slope*log(10^m2/freq_ref)/log(10);

% calculate er and loss tan at discrete frequencies
% formula from paper "IEEE EMC Trans vol 43 no 4 Nov. 2001 pp662-667"
e0=9.98e-12;
omega1=6.28*10^m1;
omega2=6.28*10^m2;
%index_ref = M*log10(freq_ref);
omega = 6.28*freq;      % omega = 2*pi*freq;

for i=1:N,
%cc    index(i)=i/M;            % No need
%cc    freq(i)=10^(i/M);        % Frequency points are from function inputs
%cc    omega(i)=6.28*10^(i/M);
    factor(i)=er_slope*log(complex(omega2, omega(i))/complex(omega1, omega(i)))/log(10);
    er_prime(i)=er_inf+real(factor(i));
    er_double_prime(i)=-(imag(factor(i))-sigma/omega(i)/e0);
    tand(i)=er_double_prime(i)/er_prime(i);
end

surface_roughness_coeff=zeros(N,1);

% create table w-element model
for i=1:N;
    % surface_roughness_coeff is from Hammerstad
    surface_roughness_coeff(i) = 1 + 2/3.1416*atan(1.4*3.1416*4*3.1416*1e-7*conductivity*(sr_rms*1e-6)^2*freq(i));
    r(i,:) = r_ref*sqrt(freq(i)/freq_ref) * surface_roughness_coeff(i);
    % 2nd l(1,:) term for internal inductance
    l(i,:) = l_ref + r(i,:)/6.28/freq(i);
    % c scaled by er ratio
    c(i,:) = c_ref*er_prime(i)/er_eff;
    % g scaled with frequency and material properties
    g(i,:) = g_ref*freq(i)/freq_ref*er_prime(i)/er_eff*tand(i)/tand_eff;

    % add DC resistance
    temp = 1;
    temp2 = 2;
    if r(i,1) < ro_ref(1)
        while temp <= N1
            r(i,temp) = ro_ref(temp);
            temp = temp + temp2;
            temp2 = temp2 + 1;
        end
    end
end

% ------------ Reshape the LCRG to matrix format ---------------- %
L_matrix = zeros(num_cond,num_cond,N);
C_matrix = zeros(num_cond,num_cond,N);
R_matrix = zeros(num_cond,num_cond,N);
G_matrix = zeros(num_cond,num_cond,N);

for ii = 1:num_cond,
    for jj = 1:ii,
        L_matrix(ii,jj,:) = l(:,ii*(ii-1)/2+jj);
        L_matrix(jj,ii,:) = L_matrix(ii,jj,:);
        C_matrix(ii,jj,:) = c(:,ii*(ii-1)/2+jj);
        C_matrix(jj,ii,:) = C_matrix(ii,jj,:);
        R_matrix(ii,jj,:) = r(:,ii*(ii-1)/2+jj);
        R_matrix(jj,ii,:) = R_matrix(ii,jj,:);
        G_matrix(ii,jj,:) = g(:,ii*(ii-1)/2+jj);
        G_matrix(jj,ii,:) = G_matrix(ii,jj,:);
    end
end

return;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%
% Look for data function
%%%%%%%%%%%%%%%%%%%%%%%%%
% Application: 
%    while look_for_data(line)
%       line = fgetl(fid)
%    end
%
% assumption that data line containts numerical data
%
% returns 0 - this line contains data
% returns 1 - no data (continue loop)

function keeplooking = look_for_data(thisline)
      
   
    if length(regexp(thisline, '*')) > 0  
        %if regexp(thisline, '\.')
        keeplooking = 1;
    else
        if length(regexp(thisline, '[0-9]')) > 0 
            keeplooking = 0;
        else
            keeplooking = 1;
        end
    end  
    
%rewind to start of first data line

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function out = nearzero(in);

% function [out] = nearzero(in)
%   Modify incoming data containing zeros to 1e-15
%   Input
%       * in is a 1-dimensional array
%   Output:
%       * out is also a\ 1-dimensional array
%
% Author: Chi-te Chen (chi-te.chen@intel.com)

index = find(in == 0);
out = in;
out(index) = 1.0e-15;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [lo,co,ro,go,rs,gd,num_cond] = readrlgc2(input_model_filename)

fid_file = fopen(input_model_filename);

% Confirm that open worked
if fid_file<0
    error(['Could not open file: ' fid_file]);
end
done_file = 0;
nheaderlines = 0;

% -------------------- Look for Number of Conductors ---------------- %
% find Number of Conductors
while ~done_file
    nheaderlines = nheaderlines+1;
    line = fgetl(fid_file);
    if regexp(line,'Number of Conductors') > 0
        done_file = 1;
    end;
end


% look for next data line (past comments)
line = fgetl(fid_file);
while look_for_data(line)
    line = fgetl(fid_file);
end;
fseek(fid_file, -(length(line))-4, 0);

% Read Lo
data = fscanf(fid_file,'%e \n', 1);
num_cond = data;

N1=(num_cond+1)*num_cond/2;

% ----------------------------------------------------------------- %
% find where data (Lo) starts
while ~done_file
    nheaderlines = nheaderlines+1;
    line = fgetl(fid_file);
    if regexp(line,'Lo') > 0
        done_file = 1;
    end;
    %look for RLCG config name
    if length(regexp(line, 'config:')) > 0
        configname = line;
    end
end

% look for next data line (past comments)
line = fgetl(fid_file);
while look_for_data(line)
    line = fgetl(fid_file);
end;
fseek(fid_file, -(length(line))-4, 0);

% Read Lo
data = fscanf(fid_file,'%e \n', N1);
lo = data;

% look for next data line
line = fgetl(fid_file);
while look_for_data(line)
    line = fgetl(fid_file);
end
fseek(fid_file, -(length(line))-4, 0);

% read Co
data = fscanf(fid_file,'%e \n', N1);
co = data;

% look for next data line
line = fgetl(fid_file);
while look_for_data(line)
    line = fgetl(fid_file);
end
fseek(fid_file, -(length(line))-4, 0);

% read ro
data = fscanf(fid_file,'%e \n', N1);
ro = data;

% look for next data line
line = fgetl(fid_file);
while look_for_data(line)
    line = fgetl(fid_file);
end
fseek(fid_file, -(length(line))-4, 0);

% read go
data = fscanf(fid_file,'%e \n', N1);
go = data;

% look for next data line
line = fgetl(fid_file);
while look_for_data(line)
    line = fgetl(fid_file);
end
fseek(fid_file, -(length(line))-4, 0);

% read rs
data = fscanf(fid_file,'%e \n', N1);
rs = data;

% look for next data line
line = fgetl(fid_file);
while look_for_data(line)
    line = fgetl(fid_file);
end
fseek(fid_file, -(length(line))-4, 0);

% read gd
data = fscanf(fid_file,'%e \n', N1);
gd = data;

fclose(fid_file);;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function ecomp=get_ecomp(num_cond, lo, co, er_prime, cvel, m1, M)
% --- Solve for Er_eff & er compensation for correct er_prime_inf (using 1ghz)
er_prime=er_prime;
offset = 0;
lo_avg=0;
co_avg=0;
for i=1:num_cond                                 % avg lo & co using diag values
    offset = i + offset;
    lo_avg = lo_avg + lo( offset );
    co_avg = co_avg + co( offset );
end
lo_avg = lo_avg / num_cond;
co_avg = co_avg / num_cond;
Er_eff = (cvel*(sqrt(lo_avg*co_avg)))^2;        %forumula for er_eff
ecomp = Er_eff - er_prime((10-m1)*M);           %10-m1 * M is 1GHZ position in er_prime
%data(3) = data(3) + ecomp;                 %erprimeinf = erprimeinf + ecomp

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
