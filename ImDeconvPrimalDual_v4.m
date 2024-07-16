function OutPut = ImDeconvPrimalDual_v4(g, H, Param)
%% [f, regpar, ISNR, iter_time] = ImDeconvPrimalDual(g, H, Param)
% Primal-dual method is applied to image deconvolution problem.
% Suppose the mathematical model is given by
%
% $$g = H*f + n$$
%
% where g, f and n are the observed image with size $n\times m$, the original image 
% and the noise respectively. H is the blur matrix. 
%
% To restore f, We solve the constrained minimization problem
% 
% $$\min_f TV(f)$$
% 
% subject to $||H*f-y||_2 \leq c$ .
% or the unconstrained minimization problem
%
% $$\min_f TV(f)+ \frac{\lambda}{2} ||H*f-y||_2^2$$
%
%%  Input:
%         Param.OrigIm:     original image              
%         Param.PrimalStep: primal step, default=1;     
%         Param.MaxIter:    default=100;                
%         Param.Sol:        initial solution, default=g 
%         Param.p:          initial dual variable
%         Param.Disp:       display the result or not. default=1;
%         Param.UpBound:    upper bound for the norm of the residual $c$ 
%         Param.GCV:        1 for using GCV method to estimate parameter
%         Param.Reglambda:  regularization parameter  $\lambda$
%         Param.SolRE:      stop criterion for the relative difference. default=1e-4;  
%%  OutPut:
%         OutPut.Sol:       restored image f
%         OutPut.p:         dual variable p;
%         OutPut.ISNR:      the evolution of ISNR against iteration number
%         OutPut.IterTime:  the evolution of CPU time against iteration number
%
%% Example 1: Constrainted Problem, estimate upbound
%  I = imread('cameraman.tif'); f = double(I);
%  psf = fspecial('average',9);
%  H = BlurMatrix(psf, size(f));
%  g = H * f + 2 * randn(size(f));
%
%  output = ImDeconvPrimalDual_v4(g, H);
%  figure; imshow(uint8(output.Sol));
%  figure; plot(output.Reglambda);
%
%% Example 2: Constrainted Problem, given upbound
%  I = imread('cameraman.tif'); f = double(I);
%  psf = fspecial('average',9);
%  H = BlurMatrix(psf, size(f));
%  g = H * f + 2 * randn(size(f));
%  Param.OrigIm     = f;
%  c = (size(g,1)*size(g,2))*ImageStdDev(g);
%  Param.UpBound    = c;
%  output = ImDeconvPrimalDual_v4(g, H, Param);
%  figure; imshow(uint8(output.Sol));
%  figure; plot(output.IterTime, output.ISNR);
%  figure; plot(output.Reglambda);
%
%% Example 3: Penalized Proble, given regularization parameter
%  I = imread('cameraman.tif'); f = double(I);
%  psf = fspecial('average',9);
%  H = BlurMatrix(psf, size(f));
%  g = H * f + 2 * randn(size(f));
%  Param.OrigIm     = f;
%  Param.Reglambda  = 7.78;
%  output = ImDeconvPrimalDual_v4(g, H, Param);
%  figure; imshow(uint8(output.Sol));
%
%% Example 4: Penalized Proble, given regularization parameter
%  I = imread('cameraman.tif'); f = double(I);
%  psf = fspecial('average',9);
%  H = BlurMatrix(psf, size(f));
%  g = H * f + 2 * randn(size(f));
%  Param.OrigIm     = f;
%  Param.GCV        = 1;
%  output = ImDeconvPrimalDual_v4(g, H, Param);
%  figure; imshow(uint8(output.Sol));
%  figure; plot(output.Reglambda);
%
%%=========================================================================
%  Copyright(c), April 2011, Dr.WEN You-wei(wenyouwei@gmail.com)
%  Copyright(c), Sep.  2012, Dr.WEN You-wei(wenyouwei@gmail.com)
%%=========================================================================


f       = g;               MaxIter = 100;        mu   = 0;
SolRE   = 1e-4;            t       = 2.;         flag = 1;
p       = OpGradient(f);   gcv     = 0;          

if nargin == 3
    if isfield(Param,'OrigIm'),     xtrue     = Param.OrigIm;     end;
    if isfield(Param,'PrimalStep'), t         = Param.PrimalStep; end
    if isfield(Param,'MaxIter'),    MaxIter   = Param.MaxIter;    end
    if isfield(Param,'Disp'),       flag      = Param.Disp;       end
    if isfield(Param,'SolRE'),      SolRE     = Param.SolRE;      end
    if isfield(Param,'Reglambda'),  lambda    = Param.Reglambda;  end
    if isfield(Param,'UpBound'),    c         = Param.UpBound;    
                                    RegMethod = 'Discrepancy';    end
    if isfield(Param,'GCV'),        gcv       = Param.GCV;        
                                    RegMethod = 'GCV';            end
    if isfield(Param,'Sol'),        f         = Param.Sol;        end
    if isfield(Param,'p');          p         = Param.p;          end
end

if ~isvar('c') && ~isvar('lambda') && ~gcv, 
    sigma = ImageStdDev(g);
    c = sigma * size(g,1) * size(g,2)*0.98;
    if sigma<1, c = 0.92 * c; end
    RegMethod = 'Discrepancy';
end
if isvar('lambda'), 
    mu = lambda * t; 
    RegMethod = 'RegFixed';
end


Param.Sol = f;
OutPut = Param;

HtH  = conj(H.eigblurmatrix) .*  H.eigblurmatrix;
Htg  = conj(H.eigblurmatrix) .*  H.diagop(g) ; %H' * y;    
gker = H.diagop(g) ;  

regpar = zeros(MaxIter,1); ISNR = regpar;  iter_time = regpar;
cont   = 1; k = 0; 
s      = 1/16/t;

start_time = cputime;
df         = OpGradient(f);

while cont
    k = k+1;  
    
    %% update the dual variable
    pold = p;
    p    = projection(pold - s * df);  

    %% update the primal variable 
    u    = f - t * div(p);
    uker = H.diagop(u); 
    bker = gker - H.eigblurmatrix .* uker; %    rsdl = H * fnew - y; 

    %% compute the parameter and estimate the image.
    switch RegMethod
        case 'GCV'
            [mu,fnew,pold] = RegSel_GCV2(H,HtH,bker,uker,gker, p,pold,Htg,t,s);
        case 'Discrepancy'
            if norm(bker,'fro') <= c,         
                mu   = 0; 
                fnew = u;
            else  %% projection of primal variable
                mu   = Projection_Regpara_Newtona(H.eigblurmatrix, bker, c, mu);
                fnew = H.diagopback((mu * Htg + uker)./(mu * HtH + 1)); 
            end 
        case 'RegFixed'
            fnew = H.diagopback((mu * Htg + uker)./(mu * HtH + 1));
    end
            
    %% update the dual variable
    df   = OpGradient(fnew);
    p    = projection(pold - s * df);
    
    %% relative error
    re   = norm(fnew-f,'fro')/norm(f,'fro');
    f    = fnew;
    cont = (k<MaxIter)&&(re>SolRE);% && (RE_eng>EngRE);    
    if isvar('xtrue'),     
        ISNR(k) = fun_ISNR(xtrue,g,f); 
        if mod(k,10)==0 && flag
            fprintf('%3d-th  isnr: %2.2f,   regpara: %1.2e  re: %1.2e\n', k,ISNR(k), mu, re);
        end
    end
    iter_time(k) = cputime - start_time;
    regpar(k)    = mu/t;         
end

OutPut.Sol      = f;
OutPut.p        = p;
OutPut.ISNR     = ISNR(1:k);
OutPut.IterTime = iter_time(1:k);
if isvar('c'), 
    tmp              = 1./(mu * HtH + 1);
    OutPut.AdjustTau = sum(tmp(:));
end
OutPut.Reglambda = regpar(1:k);


function sigma = ImageStdDev(y)
h = [0.03522629188571 0.08544127388203 -0.13501102001025 -0.45987750211849 0.80689150931109 -0.33267055295008];       
h = h(end:-1:1); %% flip is used only to give exactly same result as previous version of the code (which was using filter2)
        
z = conv2(y,h,'same');
z=conv2(z,h','same');
sigma = median(abs(z(:)))/.6745;


function [mu] = Projection_Regpara_Newtona(H, rsdl, c, mu)
%% This program is to find the parameter beta such that
%      || r ./(beta * |H_i|^2 +1)||<=c
% here
%      EigConvMat:  H_i
%      rsdl:        r
%      boundsig:
%%
if nargin <4, mu = 0; end
H2   = conj(H) .* H;
k = 0; 
absrsdl2 = conj(rsdl).*rsdl;
A = 1./(mu * H2 + 1);
rsdl2 = A .* rsdl;    %rsdl2 = inv(mu * H2 + 1) .* rsdl;

cont = 1;
while cont
    phi   = norm(rsdl2,'fro')^2;
    Htmp  = (A .* A .* A) .* H2;
    tmp   = Htmp .* absrsdl2;
    fdiff = -2 * real(sum(tmp(:)));
    dmu   = (c^2 - phi)/fdiff;
    munew = mu + dmu;
    RE    = abs(dmu)/mu;
    mu    = munew;
    if mu <0, mu = 1e-5;end
    A = 1./(mu * H2 + ones(size(H2)));
    rsdl2 = A .* rsdl;    %rsdl2 = inv(mu * H2 + 1) .* rsdl;
    
    k     = k + 1;
    cont  = (norm(rsdl2,'fro') > c)&&(k<50)&&(RE>5e-2);
end
return;


function ISNR = fun_ISNR(X,Y, Xest)
%  X:    original image
%  Y:    observed image
%  Xest: estimated image
%ISNR = 10*log10(sum(abs(X(:)-Y(:)).^2)/sum(abs(Xest(:)-X(:)).^2));
ISNR = 10*log10(norm(X-Y,'fro').^2/norm(Xest-X,'fro').^2);
