function [SSIM, RMSE, PSNR] = evaluate2(imageys, imageyh)

if ~isreal(imageys)
    image_ref = angle(imageys);
else
    image_ref = imageys;
end

if ~isreal(imageyh)
    image_test = angle(imageyh);
else
    image_test = imageyh;
end
[m,n] = size(image_test);

% --- RMSE ---
diff  = image_ref - image_test;
MSE   = sum(diff(:).^2) / (m*n);
RMSE  = sqrt(MSE);

% --- SSIM ---
SSIM = ssim(image_test, image_ref, 'DynamicRange', 2*pi);

% --- PSNR ---
MAX_I = pi;
PSNR = 10*log10((MAX_I^2) / MSE);
end




