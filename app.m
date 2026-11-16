function fig=app(varargin)
p=inputParser;
addParameter(p,'Visible','on');
addParameter(p,'AutoStart',false);
parse(p,varargin{:});
windows=findall(groot,'Type','figure','Tag','FluctuationTeachingApp');
for k=1:numel(windows)
    if strcmp(getappdata(windows(k),'InterfaceVersion'),'3.2')
        fig=windows(k);
        if strcmp(p.Results.Visible,'on'), figure(fig); end
        return
    end
end
html=fullfile(fileparts(mfilename('fullpath')),'ui.html');

model=qiu(8,4);
pre=model;
started=false; paused=true; paused0=true; phase='settings';
speed=1; n0=8; left0=4;
count=0; avg=0; m2=0; ht=[]; hl=[]; hs=[];
tnext=0; tstat=0; dt0=.2; span=120;
ack=0; seq=0;
screen=get(groot,'ScreenSize');
width=min(1420,screen(3)-64); height=min(880,screen(4)-96);
fig=uifigure('Name','粒子数涨落 · 实验台','Tag','FluctuationTeachingApp', ...
    'Position',[max(20,(screen(3)-width)/2) max(30,(screen(4)-height)/2) width height], ...
    'Color',[.95 .96 .97],'Visible','off','CloseRequestFcn',@closeapp);
fig.Theme='light';
grid=uigridlayout(fig,[1 1],'Padding',[0 0 0 0],'RowSpacing',0,'ColumnSpacing',0);
view=uihtml(grid,'Tag','SimulationSurface');
tm=timer('ExecutionMode','fixedSpacing','Period',.04,'BusyMode','drop', ...
    'Name','FluctuationTeachingTimer','TimerFcn',@tick);
setappdata(fig,'InterfaceVersion','3.2');
setappdata(fig,'SimulationModel',model);
setappdata(fig,'SimulationTimer',tm);
setappdata(fig,'HTMLSurface',view);
setappdata(fig,'HTMLReady',false);
setappdata(fig,'DispatchAction',@dispatch);
sample(model.measure());
view.DataChangedFcn=@receive;
update();
view.HTMLSource=html;
fig.Visible=p.Results.Visible;
if p.Results.AutoStart, dispatch(struct('action','start','n',8,'left',4)); end

    function receive(~,event)
        packet=event.Data;
        if strcmp(packet.kind,'action')
            dispatch(packet);
        end
    end

    function dispatch(packet)
        if isfield(packet,'id')
            ack=packet.id;
        end
        switch char(packet.action)
            case 'ready'
                setappdata(fig,'HTMLReady',true);
            case 'draft'
                [n,left]=counts(packet);
                p1=qiu(n,left);
                n0=n; left0=left; pre=p1;
            case 'start'
                [n,left]=counts(packet);
                runexp(n,left,true);
            case 'restart'
                if started, runexp(model.N,model.left0,false); end
            case 'togglePause'
                if started&&strcmp(phase,'simulation')
                    paused=~paused;
                    if paused, stop(tm); else, start(tm); end
                end
            case 'settings'
                if started
                    paused0=paused; paused=true; stop(tm);
                    n0=model.N; left0=model.left0;
                    pre=qiu(n0,left0);
                    phase='settings';
                end
            case 'return'
                if started
                    phase='simulation'; paused=paused0;
                    if~paused, start(tm); end
                end
            case 'speed'
                speed=packet.value;
            case 'resetStats'
                count=0; avg=0; m2=0; tstat=model.Time;
        end
        update();
    end

    function [n,left]=counts(packet)
        n=packet.n; left=packet.left;
        left=double(left); n=double(n);
    end

    function runexp(n,left,usepre)
        if usepre&&strcmp(phase,'settings')&&pre.N==n&&pre.left0==left
            m1=pre;
        else
            m1=qiu(n,left);
        end
        stop(tm);
        model=m1; setappdata(fig,'SimulationModel',model);
        n0=n; left0=left;
        phase='simulation'; started=true; paused=false;
        ht=[]; hl=[]; hs=[];
        count=0; avg=0; m2=0; tstat=0; tnext=0;
        sample(model.measure());
        start(tm);
    end

    function sample(state)
        ht(end+1)=state.time; hl(end+1)=state.left; hs(end+1)=state.entropy;
        count=count+1;
        %左侧比例的均值和方差
        delta=state.fraction-avg;
        avg=avg+delta/count;
        m2=m2+delta*(state.fraction-avg);
        tnext=state.time+dt0;
        keep=ht>=state.time-span;
        hs=hs(keep); hl=hl(keep); ht=ht(keep);
    end

    function tick(~,~)
        if paused||~strcmp(phase,'simulation'), return; end
        remaining=.16*speed; tc=tic;
        while remaining>1e-10
            gap=max(0,tnext-model.Time);
            if gap<1e-10, sample(model.measure()); gap=dt0; end
            dt=model.advance(min(remaining,gap),max(.001,.030-toc(tc)));
            remaining=remaining-dt;
            if model.Time>=tnext-1e-10, sample(model.measure()); end
            if toc(tc)>=.030, break; end
        end
        update();
        drawnow limitrate nocallbacks
    end

    function update()
        if strcmp(phase,'settings'), shown=pre; else, shown=model; end
        stats=shown.measure();
        seq=seq+1;
        state=struct('kind','state','sequence',seq,'ackId',ack, ...
            'phase',phase,'started',started,'paused',paused, ...
            'draftN',n0,'draftLeft',left0,'draftRight',n0-left0, ...
            'n',shown.N,'initialLeft',shown.left0,'initialRight',shown.N-shown.left0, ...
            'positions',shown.Position,'box',shown.box,'diameter',shown.Diameter, ...
            'time',stats.time,'left',stats.left,'right',stats.right,'fraction',stats.fraction, ...
            'entropy',stats.entropy,'entropyMax',stats.entropyMaximum,'temperature',stats.temperature, ...
            'speed',speed,'statsCount',count,'sigma',0,'mean',avg, ...
            'statsFrom',tstat,'referenceSigma',1/(2*sqrt(shown.N)), ...
            'historyTime',ht,'historyLeft',hl,'historyEntropy',hs);
        if count>=2, state.sigma=sqrt(max(0,m2/count)); end
        setappdata(fig,'InterfaceState',state);
        setappdata(fig,'PreviewModel',pre);
        view.Data=state;
    end

    function closeapp(~,~)
        stop(tm); delete(tm);
        delete(fig);
    end
end
