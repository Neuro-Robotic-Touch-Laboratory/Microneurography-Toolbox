function [a r phi zm] = circ_plot_ax(ax,alpha, format, formats, varargin)
%
% r = circ_plot(alpha, ...)
%   Plotting routines for circular data.
%
%   Input:
%     alpha     sample of angles in radians
%     [format		specifies style of plot
%                 pretty, histogram, density, []
%     [formats  standard matlab string for plot format (like '.r')]
%
%     The different plotting styles take optional arguments:
%         pretty:   fourth argument toggles between showing mean direction
%                     and not showing it
%         hist:     fourth argument determines number of bins/bin centers
%                   fifth argument determines whether normalized or count
%                     histogram is shown
%                   sixth argument toggles between showing mean direction
%                     and not showing it
%
%       All of these arguments can be left empty, i.e. set to [], so that
%       the default value will be used. If additional arguments are
%       supplied in the name-value style ('linewidth', 2, ...), these are
%       used to change the properties of the mean resultant vector plot.         
%
%   Output:
%     a         axis handle
%
%   Examples:
%     alpha = randn(60,1)*.4+pi/2;
%     figure
%     subplot(2,2,1)
%     circ_plot(alpha,'pretty','ro',true,'linewidth',2,'color','r'),
%     title('pretty plot style')
%     subplot(2,2,2)
%     circ_plot(alpha,'hist',[],20,true,true,'linewidth',2,'color','r')
%     title('hist plot style')
%     subplot(2,2,3)
%     circ_plot(alpha,[],'s')
%     title('non-fancy plot style')
%    
%
% Circular Statistics Toolbox for Matlab

% By Philipp Berens & Marc J. Velasco, 2009
% velasco@ccs.fau.edu, berens@tuebingen.mpg.de

if nargin < 2 || isempty(format)
    format = '';
end


switch format
  case 'pretty'
    % plot in 'pretty style'
    % draws unit circle and marks points around the circle
    % adds optionally the mean resultant vector
    
    if nargin < 3|| isempty(formats) 
      formats = 'o';
    end
    
    % convert angles to unit vectors
    z = exp(1i*alpha);

    % create unit circle
    zz = exp(1i*linspace(0, 2*pi, 101));

    plot(ax,real(z), imag(z), formats, real(zz), imag(zz), 'k', [-2 2], [0 0], 'k:', [0 0], [-2 2], 'k:');
    set(ax, 'XLim', [-1.1 1.1], 'YLim', [-1.1 1.1])

    % plot mean directions with an overlaid arrow if desired
    if nargin > 2 && ~isempty(varargin{1})
      s = varargin{1};
    else
      s = true;
    end
    
    if s
      r = circ_r(alpha);
      phi = circ_mean(alpha);
      hold (ax,'on');
      zm = r*exp(1i*phi);
      plot(ax,[0 real(zm)], [0, imag(zm)],varargin{2:end})
      hold(ax, 'off');
    end

    axis(ax,'square');
    set(ax,'box','off')
    set(ax,'xtick',[])
    set(ax,'ytick',[])
    text(ax,1.2, 0, '0'); text(ax,-.05, 1.2, '\pi/2');  text(ax,-1.35, 0, '�\pi');  text(ax,-.075, -1.2, '-\pi/2');

    
  case 'hist'
    % plot in  'hist style'
    % this is essentially a wrapper for the rose plot function of matlab
    % adds optionally the mean resultant vector
    
    if nargin < 3|| isempty(formats) 
      formats = '-';
    end
    
    if nargin > 3 && ~isempty(varargin{1})
      x = varargin{1};
    else
      x = 20;
    end
    
    [t,r] = rose(ax,alpha,x);%polarhistogram
    if nargin> 3 && varargin{2}
      polar(ax,t,2*r/sum(r),formats) % polarplot 
      mr = max(2*r/sum(r));
    else
      polar(ax,t,r,formats)
      mr = max(r);
    end
    
     % plot mean directions with an overlaid arrow if desired
    if nargin > 5 && ~isempty(varargin{3})
      s = varargin{3};
    else
      s = true;
    end
    
    if s
      r = circ_r(alpha) * mr;
      phi = circ_mean(alpha);
      hold(ax,'on');
      zm = r*exp(1i*phi);
      plot(ax,[0 real(zm)], [0, imag(zm)],varargin{4:end})
      hold(ax,'off');
    end
    
   
  otherwise
    if nargin < 3
      formats = 'o';
    end
    polar(ax,alpha, ones(size(alpha)), formats);
end

a = ax;
