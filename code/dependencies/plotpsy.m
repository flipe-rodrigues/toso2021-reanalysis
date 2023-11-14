function h = plotpsy(data, psy, plotopt)
    %UNTITLED Summary of this function goes here
    %   Detailed explanation goes here

    % default plot options
    if ~exist('plotopt','var')
        plotopt = struct();
    end
    if ~isfield(plotopt,'dataedgeclr')
        plotopt.dataedgeclr = [1,1,1] * 0;
    end
    if ~isfield(plotopt,'datafaceclr')
        plotopt.datafaceclr = [1,1,1] * .15;
    end
    if ~isfield(plotopt,'fitclr')
        plotopt.fitclr = plotopt.datafaceclr;
    end
    if ~isfield(plotopt,'gradeclrs')
        plotopt.gradeclrs = true;
    end
    if ~isfield(plotopt,'linestyle')
        plotopt.linestyle = '-';
    end
    if ~isfield(plotopt,'linewidth')
        plotopt.linewidth = 2;
    end
    if ~isfield(plotopt,'marker')
        plotopt.marker = 'o';
    end
    if ~isfield(plotopt,'markersize')
        plotopt.markersize = 7.5;
    end
    if ~isfield(plotopt,'normalizemarkersize')
        plotopt.normalizemarkersize = true;
    end
    if ~isfield(plotopt,'markersizenormalizer')
        plotopt.markersizenormalizer = mean(data.n(data.n>0));
    end
    if ~isfield(plotopt,'flipy')
        plotopt.flipy = false;
    end
    if ~isfield(plotopt,'plotdata')
        plotopt.plotdata = true;
    end
    if ~isfield(plotopt,'plotfit')
        plotopt.plotfit = true;
    end
    if ~isfield(plotopt,'overallvisibility')
        plotopt.overallvisibility = 'on';
    end
    if ~isfield(plotopt,'parent')
        plotopt.parent = gca;
    end

    % sigmoid function handle
    S = psy.options.sigmoidHandle;

    % psychometric function handle
    psi = @(x, m, w, g, l) ...
        g + (1 - l - g) * S(x,m,w);

    % data prep
    k = size(psy.data,1);
    x = data.x;
    y = data.y;
    if plotopt.flipy
        y = data.n - y;
    end
    n = data.n;
    err = data.err;

    mhat = psy.Fit(1);
    what = psy.Fit(2);
    lhat = psy.Fit(3);
    ghat = psy.Fit(4);

    xx = linspace(-.05,1.05,1e3);
    yy = psi(xx,mhat,what,ghat,lhat);
    if plotopt.flipy
        yy = 1 - yy;
    end

    % plot estimated psychometric curve
    if plotopt.plotfit
        p = gobjects(3,1);
        xrange = xx >= x(1) & xx <= x(end);
        p(1) = plot(xx(xrange),yy(xrange),...
            'color',plotopt.fitclr,...
            'linestyle',plotopt.linestyle,...
            'linewidth',plotopt.linewidth,...
            'handlevisibility','off',...
            'parent',plotopt.parent);
        p(2) = plot(xx(xx < x(1)),yy(xx < x(1)),...
            'color',plotopt.fitclr,...
            'linestyle',':',...
            'linewidth',plotopt.linewidth,...
            'handlevisibility','off',...
            'parent',plotopt.parent);
        p(3) = plot(xx(xx > x(end)),yy(xx > x(end)),...
            'color',plotopt.fitclr,...
            'linestyle',':',...
            'linewidth',plotopt.linewidth,...
            'handlevisibility','off',...
            'parent',plotopt.parent);
        uistack(p,'bottom');
    end
    
    % plot proportion of long choices
    if plotopt.plotdata
        if plotopt.gradeclrs
            datafaceclrs = colorlerp([plotopt.datafaceclr;[0,0,0]],k);
        else
            datafaceclrs = repmat(plotopt.datafaceclr,k,1);
        end
        p = gobjects(k,1);
        for ii = 1 : k
            if plotopt.normalizemarkersize
                markersize = max(n(ii) / ...
                    plotopt.markersizenormalizer * plotopt.markersize, 2.5);
            else
                markersize = plotopt.markersize;
            end
            p(ii) = errorbar(x(ii),y(ii)/n(ii),err(ii),...
                'color',plotopt.dataedgeclr,...
                'marker',plotopt.marker,...
                'markersize',markersize,...
                'markerfacecolor',datafaceclrs(ii,:),...
                'linewidth',plotopt.linewidth,...
                'capsize',0,...
                'handlevisibility','off',...
                'parent',plotopt.parent);
        end
        uistack(p,'top');
    end
    
    % dummy legend handle
    if plotopt.overallvisibility
        h = plot(max(xlim)*1.5,max(ylim)*1.5,...
            'color',plotopt.fitclr,...
            'marker',plotopt.marker,...
            'markerfacecolor',plotopt.datafaceclr,...
            'markeredgecolor',plotopt.dataedgeclr,...
            'linestyle','-',...
            'markersize',plotopt.markersize,...
            'linewidth',plotopt.linewidth);
    end
end