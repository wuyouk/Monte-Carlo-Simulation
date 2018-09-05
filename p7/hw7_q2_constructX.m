function xX = hw7_q2_constructX(x, typeOfPolynomial)
%
len1 = length(x);

switch typeOfPolynomial
    
    case 'Chebychev_firstKind'
        
        %% Chebychev typeOfPolynomials (first kind)
        u0 = ones(len1,1);
        u1 = x;
        u2 = 2*x.^2 - 1;
        u3 = 4*x.^3 - 3*x;
        u4 = 8*x.^4 - 8*x.^2 + 1;
        u5 = 16*x.^5 - 20*x.^3 + 5*x;
        xX = [u0 u1 u2 u3]; % u4 u5];
    
    case 'Chebychev_secondKind'
        
        %% Chebychev typeOfPolynomials (second kind)
        u0 = ones(len1,1);
        u1 = 2*x;
        u2 = 4*x.^2-1;
        u3 = 8*x.^3-4*x;
        u4 = 16*x.^4 - 12*x.^2 + 1;
        u5 = 32*x.^5 - 32*x.^3 + 6*x;
        %u6 = 64*x.^6 - 80*x.^4 + 24*x.^2 - 1;
        %u7 = 128*x.^7 - 192*x.^5 + 80*x.^3 - 8*x;
        %u8 = 256*x.^8 - 448*x.^6 + 240*x.^4 - 40*x.^2 + 1;
        %u9 = 512*x.^9 - 1024*x.^7 + 672*x.^5 - 160*x.^3 + 10*x;
        %
        %xX = [u0 u1 u2 u3 u4 u5 u6 u7 u8 u9];
        xX = [u0 u1 u2 u3]; % u4 u5];
        
    case 'Laguerre'
        %% Laguerre typeOfPolynomials
        u0 = ones(len1,1);
        u1 = -x+1;
        u2 = (x.^2-4*x+2)/2;
        u3 = (-x.^3+9*x.^2-18*x+6)/6;
        u4 = (x.^4-16*x.^3+72*x.^2-96*x+24)/24;
        xX = [u0 u1 u2 u3]; % u4];
end
        