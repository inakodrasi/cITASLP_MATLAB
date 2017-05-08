function coh_diff = comp_diffcoh(D_mic,nfft,fs,option)

if nargin<4
    option = 'spherical';
end


k = 0:nfft/2;
f = fs*k/nfft;

c = 340;
M = size(D_mic,1);
k = 0:nfft/2;
coh_diff = ones(M,M,length(k));
for m_out = 1:M-1
    for m_in = m_out+1:M
        if strcmp(option, 'spherical'),
            coh_diff(m_out,m_in,:) = sinc ( 2*f * D_mic(m_out,m_in)/c ); 
        else
            coh_diff(m_out,m_in,:) =  besselj(0, 2*pi*f * D_mic(m_out,m_in)/c);
        end
        coh_diff(m_in,m_out,:) = coh_diff(m_out,m_in,:);
    end
end