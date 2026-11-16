classdef qiu<handle
    properties (Constant)
        Diameter=1
        Mass=1
        kt=1
        box=[54 32 24]
    end
    properties (SetAccess=private)
        N
        left0
        Seed
        Radius=0.5
        Position
        Velocity
        Time=0
    end
    properties (Access=private)
        pa
        pb
        ref
        Skin=1
        Low
        High
    end
    methods
        function obj=qiu(n,left,seed)
            if nargin<1, n=8; end
            if nargin<2, left=floor(n/2); end
            if nargin<3, seed=[]; end
            obj.reset(n,left,seed);
        end

        function reset(obj,n,left,seed)
            if nargin<4||isempty(seed), seed=obj.newseed(); end
            obj.Seed=double(seed); obj.N=double(n); obj.left0=double(left);
            obj.Time=0; obj.Low=-obj.box/2; obj.High=obj.box/2;
            stream=RandStream('mt19937ar','Seed',seed);
            hi=obj.High-obj.Radius; lo=obj.Low+obj.Radius;
            spacing=obj.Diameter;
            dims=floor((hi-lo)/spacing)+1;
            cells=cell(prod(dims),1);
            p=zeros(n,3);
            offsets=obj.offset(false);
            for k=1:n
                while true
                    candidate=lo+rand(stream,1,3).*(hi-lo);
                    if k<=left
                        candidate(1)=lo(1)+rand(stream)*(-lo(1));
                    else
                        candidate(1)=rand(stream)*hi(1);
                    end
                    cellpos=floor((candidate-lo)/spacing)+1;
                    neighbors=cellpos+offsets;
                    neighbors=neighbors(all(neighbors>=1&neighbors<=dims,2),:);
                    ids=neighbors(:,1)+(neighbors(:,2)-1)*dims(1)+...
                        (neighbors(:,3)-1)*dims(1)*dims(2);
                    nearby=[cells{ids}];
                    if isempty(nearby)||all(sum((p(nearby,:)-candidate).^2,2)>obj.Diameter^2)
                        p(k,:)=candidate;
                        id=cellpos(1)+(cellpos(2)-1)*dims(1)+(cellpos(3)-1)*dims(1)*dims(2);
                        cells{id}(end+1)=k;
                        break
                    end
                end
            end
            %按总动能缩放初速度
            v=randn(stream,n,3);
            v=v*sqrt(3*n/sum(v(:).^2));
            obj.Position=p; obj.Velocity=v;
            obj.pairs();
        end

        function elapsed=advance(obj,duration,budget)
            if nargin<3, budget=inf; end
            t0=obj.Time;
            tc=tic;
            t1=t0+duration;
            while obj.Time<t1-1e-12
                vmax=max(sqrt(sum(obj.Velocity.^2,2)));
                step=min([t1-obj.Time,0.25,obj.Skin/(4*max(vmax,eps))]);
                if obj.N>128
                    travel=max(sqrt(sum((obj.Position-obj.ref).^2,2)));
                    if travel+vmax*step>=obj.Skin/2
                        obj.pairs();
                    end
                end
                [tw,i,axis]=obj.wall();
                [tp,ip]=obj.pair();
                dt=min([step,tw,tp]);
                obj.Position=obj.Position+obj.Velocity*dt;
                obj.Time=obj.Time+dt;
                if tw<=step&&tw<=tp
                    if obj.Velocity(i,axis)>0
                        obj.Position(i,axis)=obj.High(axis)-obj.Radius;
                    else
                        obj.Position(i,axis)=obj.Low(axis)+obj.Radius;
                    end
                    obj.Velocity(i,axis)=-obj.Velocity(i,axis);
                elseif tp<=step
                    a=obj.pa(ip); b=obj.pb(ip);
                    dr=obj.Position(b,:)-obj.Position(a,:);
                    dist=norm(dr);
                    normal=dr/dist;
                    dv=dot(obj.Velocity(b,:)-obj.Velocity(a,:),normal);
                    obj.Velocity(a,:)=obj.Velocity(a,:)+dv*normal;
                    obj.Velocity(b,:)=obj.Velocity(b,:)-dv*normal;
                    shift=(obj.Diameter-dist)/2;
                    obj.Position(a,:)=obj.Position(a,:)-shift*normal;
                    obj.Position(b,:)=obj.Position(b,:)+shift*normal;
                end
                if toc(tc)>=budget, break; end
            end
            elapsed=obj.Time-t0;
        end

        function s=measure(obj)
            s.left=sum(obj.Position(:,1)<0);
            s.right=obj.N-s.left;
            s.fraction=s.left/obj.N;
            s.entropy=obj.entropy(s.left,obj.N);
            s.entropyMaximum=obj.entropy(floor(obj.N/2),obj.N);
            s.time=obj.Time;
            s.kineticEnergy=obj.Mass*sum(obj.Velocity(:).^2)/2;
            s.temperature=2*s.kineticEnergy/(3*obj.N*obj.kt);
        end
    end

    methods (Access=private)
        function [t,i,axis]=wall(obj)
            hi=obj.High-obj.Radius; lo=obj.Low+obj.Radius;
            tt=inf(obj.N,3);
            for dim=1:3
                v=obj.Velocity(:,dim);
                positive=v>0; negative=v<0;
                tt(positive,dim)=(hi(dim)-obj.Position(positive,dim))./v(positive);
                tt(negative,dim)=(lo(dim)-obj.Position(negative,dim))./v(negative);
            end
            [t,index]=min(tt(:));
            t=max(0,t);
            [i,axis]=ind2sub(size(tt),index);
        end

        function [t,index]=pair(obj)
            t=inf; index=0;
            if isempty(obj.pa), return; end
            r=obj.Position(obj.pb,:)-obj.Position(obj.pa,:);
            v=obj.Velocity(obj.pb,:)-obj.Velocity(obj.pa,:);
            bb=sum(r.*v,2); aa=sum(v.^2,2);
            cc=sum(r.^2,2)-obj.Diameter^2;
            dd=bb.^2-aa.*cc;
            valid=aa>eps&bb<-1e-12&dd>=0;
            validx=find(valid);
            if isempty(validx), return; end
            %两球碰撞时间
            tt=max(0,cc(validx))./(-bb(validx)+sqrt(dd(validx)));
            [t,index2]=min(tt); index=validx(index2);
        end

        function pairs(obj)
            obj.ref=obj.Position;
            if obj.N<=128
                [obj.pa,obj.pb]=find(triu(true(obj.N),1));
                return
            end
            rc=obj.Diameter+obj.Skin;
            dims=floor((obj.High-obj.Low)/rc)+1;
            xyz=floor((obj.Position-obj.Low)/rc)+1;
            xyz=min(max(xyz,1),dims);
            ids=xyz(:,1)+(xyz(:,2)-1)*dims(1)+(xyz(:,3)-1)*dims(1)*dims(2);
            [~,order]=sort(ids);
            counts=accumarray(ids,1,[prod(dims),1]);
            first=cumsum([1;counts(1:end-1)]);
            offsets=obj.offset(true);
            as=cell(size(offsets,1),1); bs=as;
            for k=1:size(offsets,1)
                q=xyz+offsets(k,:);
                inside=all(q>=1&q<=dims,2);
                src=find(inside); q=q(inside,:);
                dst=q(:,1)+(q(:,2)-1)*dims(1)+(q(:,3)-1)*dims(1)*dims(2);
                count=counts(dst); occupied=count>0;
                src=src(occupied); dst=dst(occupied); count=count(occupied);
                if isempty(src), continue; end
                a=repelem(src,count); a=a(:);
                base=repelem(cumsum(count)-count,count); base=base(:);
                start=repelem(first(dst),count); start=start(:);
                local=(0:sum(count)-1)'-base;
                b=order(start+local); b=b(:);
                if all(offsets(k,:)==0)
                    keep=b>a; a=a(keep); b=b(keep);
                end
                keep=sum((obj.Position(a,:)-obj.Position(b,:)).^2,2)<=rc^2;
                as{k}=a(keep); bs{k}=b(keep);
            end
            obj.pa=vertcat(as{:}); obj.pb=vertcat(bs{:});
        end
    end

    methods (Static)
        function s1=entropy(left,n)
            s1=(gammaln(n+1)-gammaln(left+1)-gammaln(n-left+1))/(n*log(2));
        end
    end
    methods (Static,Access=private)
        function seed=newseed()
            persistent rng0 seed0
            if isempty(rng0)
                rng0=RandStream('mt19937ar','Seed','shuffle');
                seed0=-1;
            end
            seed=floor(rand(rng0)*2^32);
            while seed==seed0, seed=floor(rand(rng0)*2^32); end
            seed0=seed;
        end

        function offsets=offset(half)
            [x,y,z]=ndgrid(-1:1,-1:1,-1:1);
            offsets=[x(:),y(:),z(:)];
            if half
                offsets=offsets(z(:)>0|(z(:)==0&y(:)>0)|...
                    (z(:)==0&y(:)==0&x(:)>=0),:);
            end
        end
    end
end
