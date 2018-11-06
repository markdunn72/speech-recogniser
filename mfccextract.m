function mfccextract(filename)
% parameters: filename of input audio file 
% returns:    MFCC file

[y,Fs] = audioread(filename); % get signal & sample rate of input file
ySize = length(y);            % get length of signal

% --------- Configuration ---------- %
disp('Analyzing signal...')
fDur = 25;                % frame duration (ms)
fSDur = 10;               % frame step     (ms)
fSize = (Fs/1000)*fDur;   % frame length
fStep = (Fs/1000)*fSDur;  % frame step length
fN = 1;                   % max # of frames in signal
i = 1;      
while i+fSize-1 <= ySize
    i = i+fStep;
    fN = fN+1;
end
fN = fN-1;
% ---------------------------------- %

% ------ Mel-scale filterbank ------ %
hz2mel = @(hz)(1125*log(1+hz/700));     % Hertz to mel warping
mel2hz = @(mel)(700*exp(mel/1125)-700); % mel to Hertz warping
K = 257;                                % length of each filter vector
M = 26;                                 % number of filters - SET HERE
[FB,~] = trifbank(M,K,[0 Fs/2],Fs,hz2mel,mel2hz);
% ---------------------------------- %

% ------- Signal processing -------- %
disp('Processing signal...')
y = preemphasis(y);                % apply preemphasis to signal

F = zeros(12,fN);                  % matrix F for MFCC feature vectors
i = 1;                             % signal start index

for f = 1:fN                       % for each frame
    iEnd = i+fSize-1;              % get end index
    ps = spectral(y(i:iEnd));      % get spectral domain  
    ev = filterbank(ps,FB,M);      % get energy vector from filterbank      
    ev = log(ev);                  % take log of each filter energy    
    fv = cepstral(ev);             % get cepstral coefficients
    %fve = loge(fv);               % add log energy component
    F(:,f) = fv;                   % store feature vector
    i = i+fStep;                   % step signal start index
end
% ---------------------------------- %

disp('Writing file...')

writehtk_lite(strcat(filename(1:end-4),'.mfcc'), F.', fSDur*1E-3, 9);
%mfccfile = writemfcc(F,filename, fSize); % write feature vectors to .mfcc file   

end

function sigout = preemphasis(y)
% This function applies pre emphasis to
% a time-domain signal, returns signal
b = [1 -0.95];
sigout = filter(b,1,y);       % high-pass filter
end

function ps = spectral(frame)
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

function fv = cepstral(ev)
% This function takes the discrete cosine transform
% of an energy vector and removes the pitch-related
% quefrencies with truncation
ev = fft(ev);  % take DCT of ev to get cepstral domain
fv = ev(1:12); % take first 12 points to get MFCC feature vector
end

function fve = loge(fv)
% This function returns the feature vector with
% log energy component
[N,~] = size(fv);
ec = (sum(fv))^2; % calculate energy component
lec = log(ec);    % take log of ec
fve = zeros(N+1,1);
for i = 1:N
    fve(i) = fv(i);
end
fve(N+1) = lec;   % append energy components
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


function mfccfile = writemfcc(J,filename,fSize)
% This function writes a matrix of MFCC feature vectors
% to the HTK file format .mfcc
mfccfile = strcat(filename(1:end-4),'.mfcc'); % create new filename
[vecDims,nSamples] = size(J);
sampPeriod = (fSize/2)*10000; % convert ms to 100ns
parmKind = 6;                 % sample kind is MFCC
sampSize = 4*vecDims;         % 12 4-byte float values 
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

function writehtk_lite( filename, features, sampPeriod, parmKind )
% WRITEHTK_LITE Simple routine for writing HTK feature files.
%
%   WRITEHTK_LITE( FILENAME, FEATURES, SAMPPERIOD, PARMKIND )
%   writes FEATURES to HTK [1] feature file specified by FILENAME,
%   with sample period (s) defined in SAMPPERIOD and parameter kind
%   in PARAMKIND. Note that this function provides a trivial 
%   implementation with limited functionality. For fully featured 
%   support of HTK I/O refer for example to the VOICEBOX toolbox [2].
%   
%   Inputs
%           FILENAME is a filename as string for a HTK feature file
%
%           FEATURES is a feature matrix with feature vectors 
%           as rows and feature dimensions as columns
%
%           SAMPPERIOD is a sample period (s)
%
%           PARMKIND is a code indicating a sample kind
%           (see Sec. 5.10.1 of [1], pp. 80-81)
%
%   Example
%           % write features to sp10_htk.mfc file with sample period 
%           % set to 10 ms and feature type specified as MFCC_0
%           readhtk_lite( 'sp10_htk.mfc', features, 10E-3, 6+8192 );
%
%   References
%
%           [1] Young, S., Evermann, G., Gales, M., Hain, T., Kershaw, D., 
%               Liu, X., Moore, G., Odell, J., Ollason, D., Povey, D., 
%               Valtchev, V., Woodland, P., 2006. The HTK Book (for HTK 
%               Version 3.4.1). Engineering Department, Cambridge University.
%               (see also: http://htk.eng.cam.ac.uk)
%
%           [2] VOICEBOX: MATLAB toolbox for speech processing by Mike Brookes
%               url: http://www.ee.ic.ac.uk/hp/staff/dmb/voicebox/voicebox.html
%   Author: Kamil Wojcicki, September 2011
    mfcfile = fopen( filename, 'w', 'b' );
    [ nSamples, sampSize ] = size( features );
    
    fwrite( mfcfile, nSamples, 'int32' );
    fwrite( mfcfile, sampPeriod*1E7, 'int32' );
    fwrite( mfcfile, 4*sampSize, 'int16' );
    fwrite( mfcfile, parmKind, 'int16' );
    
    count = fwrite( mfcfile, features.', 'float' );
    fclose( mfcfile );
    if count~=nSamples*sampSize
        error( sprintf('write_HTK_file: count!=nSamples*sampSize (%i!=%i), filename: %s', count, nSamples*sampSize, filename)); 
    end
end
% EOF
