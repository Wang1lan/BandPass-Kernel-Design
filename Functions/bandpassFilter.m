    function filter=bandpassFilter(fRadius,sigmaHighFreq,sigmaLowFreq,scaleLowFreq)
        if ~exist('fRadius','var')
            fRadius=3;
        end
        if ~exist('sigmaLowFreq','var')
            sigmaLowFreq=0.3;
        end
        if ~exist('sigmaHighFreq','var')
            sigmaHighFreq=1.0;
        end
        if ~exist('scaleLowFreq','var')
            scaleLowFreq=0.6;
        end
        sz=2*fRadius+1;
        scale=32;
        fx=linspace(-1,1, scale*sz);
        fy=fx;
        [x,y]=meshgrid(fx,fy);
        sqX=x.^2+y.^2;
        fL=exp(-sqX/(sigmaLowFreq^2));
        fH=exp(-sqX/(sigmaHighFreq^2));
        fBP=fH-scaleLowFreq*fL;
        filter=fsamp2(fx,fy,fBP,[sz,sz]);
        % filter=fsamp2(fBP);
        filter=filter/sum(filter(:));
        filter(fRadius+1,fRadius+1)=filter(fRadius+1,fRadius+1)-1;

        sumVal=sum(sum(filter(filter>0)));
        filter=filter/sumVal;
    end