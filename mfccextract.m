function mfccfile = mfccextract(filename,sr)
% This function takes all 26 second FsReq (Hz) .wav files 
% and returns .mfcc file of MFCC feature vectors
% to be used with HTK
FsReq = sr;   % set the sampling frequency (Hz)
frameMs = 25; % set the frame size (ms)

N = (FsReq/1000)*frameMs; % # of samples in each frame

% ------ Mel-scale filterbank ------ %
hz2mel = @(hz)(1125*log(1+hz/700)); % Hertz to mel warping
mel2hz = @(mel)(700*exp(mel/1125)-700); % mel to Hertz warping
K = 257; % length of each filter vector
M = 26; % number of filters - SET HERE
[H1,~] = trifbank(M,K,[0 FsReq/2],FsReq,hz2mel,mel2hz);
% ---------------------------------- %

[y,Fs] = audioread(filename);      % get signal, sample rate of wav file
if Fs ~= FsReq                     % check wav file is correct sample rate    
    y = resampleaudio(y,Fs,FsReq); % resample to required sample rate 
end    
y = preemphasis(y);                % apply preemphasis to signal
y = buffer(y,Fs,Fs/2);             % divide into 50% overlapped frames
y =(  y' * diag(hamming(Fs))  )';  % apply hamming window
dim = size(y);                     % get dimensions of frame matrix
numFrames = dim(2);                % get number of frames
J = zeros(12,numFrames);           % matrix J for MFCC feature vectors
    
for frame = 1:numFrames   
    ps = magSpec(y(:,frame));      % get power spectrum vector       
    ev = filterbank(ps,H1,M);      % get energy vector from filterbank      
    ev = log(ev);                  % take log of each filter energy    
    fv = dcttrunc(ev);             % take dct and truncate feature vector
    %fve = loge(fv);               % add log energy component
    J(:,frame) = fv;               % store feature vector of each frame     
end 

mfccfile = writemfcc(J,filename,...
    frameMs);                      % write feature vectors to .mfcc file   
end

function ynew = resampleaudio(y,Fs,FsReq)
[P,Q] = rat(FsReq/Fs);
ynew = resample(y,P,Q);
end

function sigout = preemphasis(y)
% This function applies pre emphasis to
% a time-domain signal, returns signal
%sigout = sum(y,2)/size(y,2); % convert stereo signal to mono
b = [1 -0.95];
sigout = filter(b,1,y);       % high-pass filter
end

function ps = magSpec(frame)
% This function returns a 257 point 
% 1 sided estimated power spectral
% density of a column vector
dft = fft(frame,512); % 512 point DFT of frames
p = abs(dft).^2;      % periodogram of DFT
ps = p(1:257);        % 1 sided power spectrum is first 257 coefficients
end

function ev = filterbank(ps,H,M)
% H = filterbank matrix
% M = number of filters
ev = zeros(M,1);         % energy vector ev
for filter = 1:M
    fi = H(filter,:).';  % get filter vector
    fiev = zeros(257,1); % stores energy coefficients of ps for each filter
    for x = 1:257
        fiev(x) = (ps(x)*fi(x)); % calculate filter energy
    end
    fie = sum(fiev);     % take sum of coefficients
    ev(filter) = fie;    % store energy of each filter
end
end

function fv = dcttrunc(ev)
% This function takes the discrete cosine transform
% of an energy vector and removes the pitch-related
% quefrencies with truncation
ev = dct(ev);  % take DCT of ev to get cepstral domain
fv = ev(1:12); % take first 12 points to get MFCC feature vector
end

function fve =  loge(fv)
% This function returns the feature vector with
% log energy component
[N,~] = size(fv);
ec = (sum(fv))^2; % calculate energy component
lec = log(ec);  % take log of ec
fve = zeros(N+1,1);
for i = 1:N
    fve(i) = fv(i);
end
fve(N+1) = lec; % append energy components
end

function mfccfile = writemfcc(J,filename,ms)
% This function writes a matrix of MFCC feature vectors
% to the HTK file format .mfcc
mfccfile = strcat(filename(1:end-4),'.mfcc'); % create new filename
[vecDims,nSamples] = size(J);
sampPeriod = (ms/2)*10000;  % convert ms to 100ns
parmKind = 6;               % sample kind is MFCC
sampSize = 4*12;            % 12 4-byte float values 
% Open file for writing:
fid = fopen(mfccfile, 'w', 'ieee-be');
% Write the header information
fwrite(fid,nSamples,'int32');   % number of vectors in file (4 byteint)
fwrite(fid,sampPeriod,'int32'); % sample period in 100ns units (4 byte int)
fwrite(fid,sampSize,'int16');   % number of bytes per sample (2 byte int)
fwrite(fid,parmKind,'int16');   % code for the sample kind (2 byte int)
% Write the data: one coefficient at a time:
for i = 1:nSamples     
   for j = 1:vecDims    
     fwrite(fid, J(j, i), 'float32');    
   end
end
% close the file
fclose(fid);
end
    
function [ H, f, c ] = trifbank( M, K, R, fs, h2w, w2h )
% TRIFBANK Triangular filterbank.
%
%   [H,F,C]=TRIFBANK(M,K,R,FS,H2W,W2H) returns matrix of M triangular filters 
%   (one per row), each K coefficients long along with a K coefficient long 
%   frequency vector F and M+2 coefficient long cutoff frequency vector C. 
%   The triangular filters are between limits given in R (Hz) and are 
%   uniformly spaced on a warped scale defined by forward (H2W) and backward 
%   (W2H) warping functions.
%
%   Inputs
%           M is the number of filters, i.e., number of rows of H
%
%           K is the length of frequency response of each filter 
%             i.e., number of columns of H
%
%           R is a two element vector that specifies frequency limits (Hz), 
%             i.e., R = [ low_frequency high_frequency ];
%
%           FS is the sampling frequency (Hz)
%
%           H2W is a Hertz scale to warped scale function handle
%
%           W2H is a wared scale to Hertz scale function handle
%
%   Outputs
%           H is a M by K triangular filterbank matrix (one filter per row)
%
%           F is a frequency vector (Hz) of 1xK dimension
%
%           C is a vector of filter cutoff frequencies (Hz), 
%             note that C(2:end) also represents filter center frequencies,
%             and the dimension of C is 1x(M+2)
%
%   Example
%           fs = 16000;               % sampling frequency (Hz)
%           nfft = 2^12;              % fft size (number of frequency bins)
%           K = nfft/2+1;             % length of each filter
%           M = 23;                   % number of filters
%
%           hz2mel = @(hz)(1127*log(1+hz/700)); % Hertz to mel warping function
%           mel2hz = @(mel)(700*exp(mel/1127)-700); % mel to Hertz warping function
%
%           % Design mel filterbank of M filters each K coefficients long,
%           % filters are uniformly spaced on the mel scale between 0 and Fs/2 Hz
%           [ H1, freq ] = trifbank( M, K, [0 fs/2], fs, hz2mel, mel2hz );
%
%           % Design mel filterbank of M filters each K coefficients long,
%           % filters are uniformly spaced on the mel scale between 300 and 3750 Hz
%           [ H2, freq ] = trifbank( M, K, [300 3750], fs, hz2mel, mel2hz );
%
%           % Design mel filterbank of 18 filters each K coefficients long, 
%           % filters are uniformly spaced on the Hertz scale between 4 and 6 kHz
%           [ H3, freq ] = trifbank( 18, K, [4 6]*1E3, fs, @(h)(h), @(h)(h) );
%
%            hfig = figure('Position', [25 100 800 600], 'PaperPositionMode', ...
%                              'auto', 'Visible', 'on', 'color', 'w'); hold on; 
%           subplot( 3,1,1 ); 
%           plot( freq, H1 );
%           xlabel( 'Frequency (Hz)' ); ylabel( 'Weight' ); set( gca, 'box', 'off' ); 
%       
%           subplot( 3,1,2 );
%           plot( freq, H2 );
%           xlabel( 'Frequency (Hz)' ); ylabel( 'Weight' ); set( gca, 'box', 'off' ); 
%       
%           subplot( 3,1,3 ); 
%           plot( freq, H3 );
%           xlabel( 'Frequency (Hz)' ); ylabel( 'Weight' ); set( gca, 'box', 'off' ); 
%
%   Reference
%           [1] Huang, X., Acero, A., Hon, H., 2001. Spoken Language Processing: 
%               A guide to theory, algorithm, and system development. 
%               Prentice Hall, Upper Saddle River, NJ, USA (pp. 314-315).
%
%   Author
%           Kamil Wojcicki, June 2011

    if( nargin~= 6 ), help trifbank; return; end; % very lite input validation

    f_min = 0;          % filter coefficients start at this frequency (Hz)
    f_low = R(1);       % lower cutoff frequency (Hz) for the filterbank 
    f_high = R(2);      % upper cutoff frequency (Hz) for the filterbank 
    f_max = 0.5*fs;     % filter coefficients end at this frequency (Hz)
    f = linspace( f_min, f_max, K ); % frequency range (Hz), size 1xK

    % filter cutoff frequencies (Hz) for all filters, size 1x(M+2)
    c = w2h( h2w(f_low)+[0:M+1]*((h2w(f_high)-h2w(f_low))/(M+1)) );

    % implements Eq. (6.140) given in [1] (for a given set of warping functions)
    H = zeros( M, K );                  % zero otherwise
    for m = 1:M 
        k = f>=c(m)&f<=c(m+1);          % up-slope
        H(m,k) = 2*(f(k)-c(m)) / ((c(m+2)-c(m))*(c(m+1)-c(m)));
        k = f>=c(m+1)&f<=c(m+2);        % down-slope
        H(m,k) = 2*(c(m+2)-f(k)) / ((c(m+2)-c(m))*(c(m+2)-c(m+1)));
    end

    % H = H./repmat(max(H,[],2),1,K);  % normalize to unit height
    % H = H./repmat(trapz(f,H,2),1,K); % normalize to unit area (inherently done)
end

%%% EOF