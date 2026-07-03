function [RPN,PSD,SPD]=evaluate1(imageyh,h)

if ~isreal(imageyh)
    imagen = angle(imageyh);
else
    imagen = imageyh;
end

[m,n]=size(imagen);
RPN = 0;
SPSD = 0;
SPD=0;
image2 = padarray(imagen,[h h],'symmetric');
image3 = padarray(imagen,[1 1],'symmetric');

% PSD
for i=1:m
    for j=1:n
        i1=i+h;
        j1=j+h;
        w=image2(i1-h:i1+h,j1-h:j1+h);
        % [mw,nw]=size(w);
        w1=mean2(w);
        % for i1=1:mw
        %     for j1=1:nw
        psdw = (w-w1).^2;
        psd = sqrt((sum(psdw(:)))/((2*h+1)^2-1));
        SPSD = SPSD+psd;
    end
end
PSD=SPSD/(m*n);

% SPD
for i=1:m
    for j=1:n
        ii=i+1;
        jj=j+1;
        APD = 1/8*(abs(image3(ii,jj)-image3(ii-1,jj-1))+abs(image3(ii,jj)-image3(ii,jj-1))+abs(image3(ii,jj)-image3(ii-1,jj))...
            +abs(image3(ii,jj)-image3(ii-1,jj+1))+abs(image3(ii,jj)-image3(ii+1,jj-1))+abs(image3(ii,jj)-image3(ii+1,jj))...
            +abs(image3(ii,jj)-image3(ii,jj+1))+abs(image3(ii,jj)-image3(ii+1,jj+1)));
        SPD = SPD+APD;
    end
end

% RPN
if ~isreal(imageyh)
    phase = angle(imageyh);
else
    phase = imageyh;
end

[M, N] = size(phase);
residue_map = zeros(M, N);

for i = 1:M-1
    for j = 1:N-1
        block = phase(i:i+1, j:j+1);
        delta1 = wrapToPi(block(1,2) - block(1,1));
        delta2 = wrapToPi(block(2,2) - block(1,2));
        delta3 = wrapToPi(block(2,1) - block(2,2));
        delta4 = wrapToPi(block(1,1) - block(2,1));
        total_delta = delta1 + delta2 + delta3 + delta4;
        
        if abs(total_delta - 2*pi) < 1e-6
            residue_map(i, j) = 1;
        elseif abs(total_delta + 2*pi) < 1e-6
            residue_map(i, j) = -1;
        end
    end
end
RPN = sum(abs(residue_map(:)));
end

function wrapped = wrapToPi(phase_diff)
    wrapped = mod(phase_diff + pi, 2*pi) - pi;
end