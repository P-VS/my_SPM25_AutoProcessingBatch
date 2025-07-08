function out=save_rp_plot(sub,ses,run,task,datpath,params)

if sub<10
    substring = ['sub-0' num2str(sub)];
else
    substring = ['sub-' num2str(sub)];
end

subpath = fullfile(datpath,substring,['ses-00' num2str(ses)]);

preproc_func = fullfile(subpath,params.func_save_folder);

rpdir = dir([preproc_func filesep 'rp_*' task '_bold.txt']);
if ~isempty(rpdir)
    rpfile = fullfile(preproc_func,rpdir.name);
else
    out = 0;
    return
end

Params=load(rpfile);

fg = spm_figure('FindWin','Graphics');

spm_figure('Clear','Graphics');
ax = axes('Position',[0.1 0.65 0.8 0.2],'Parent',fg,'Visible','off');
set(get(ax,'Title'),'String','Image realignment',...
'FontSize',16,'FontWeight','Bold','Visible','on');
x     =  0.1;
y     =  0.9;
text(x,y,[substring ' Session ' num2str(ses) ' task: ' task],...
'FontSize',10,'Interpreter','none','Parent',ax);
y = y - 0.08;

ax = axes('Position',[0.1 0.35 0.8 0.2],'Parent',fg,'XGrid','on','YGrid','on',...
'NextPlot','replacechildren','ColorOrder',[0 0 1;0 0.5 0;1 0 0]);
plot(Params(:,1:3),'Parent',ax)
s  = {'x translation','y translation','z translation'};
legend(ax, s, 'Location','Best')
set(get(ax,'Title'),'String','translation','FontSize',16,'FontWeight','Bold');
set(get(ax,'Xlabel'),'String','image');
set(get(ax,'Ylabel'),'String','mm');

ax = axes('Position',[0.1 0.05 0.8 0.2],'Parent',fg,'XGrid','on','YGrid','on',...
'NextPlot','replacechildren','ColorOrder',[0 0 1;0 0.5 0;1 0 0]);
plot(Params(:,4:6)*180/pi,'Parent',ax)
s  = {'pitch','roll','yaw'};
%text([2 2 2], Params(2, 4:6)*180/pi, s, 'Fontsize',10,'Parent',ax)
legend(ax, s, 'Location','Best')
set(get(ax,'Title'),'String','rotation','FontSize',16,'FontWeight','Bold');
set(get(ax,'Xlabel'),'String','image');
set(get(ax,'Ylabel'),'String','degrees');

saveas(fg,fullfile(preproc_func,['rp_plots_' substring '_task-' task '_bold.png']));

out=1;

end