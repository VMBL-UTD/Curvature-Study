% F : geometry faces
% V : geometry vertices
function [surfs] = separateMeshSurfaces(F,V,ecnorm_xyz,ecnorm_uvw,offset)

% Create storage for output struct
surfs = struct();
surfs.caps{1,1}.V = [];
surfs.caps{1,1}.F = [];
surfs.lumen.V = [];
surfs.lumen.F   = [];
surfs.outer.V = [];
surfs.outer.F   = [];

% Get face connectivity
C = patchConnectivity(F,V);

% Get the number of end caps
num_caps = size(ecnorm_xyz,1);

% Offset end cap vectors
ecnorm_xyz = ecnorm_xyz - ecnorm_uvw.*offset;

% Create blank masks for nodes
ec_nodes = zeros(size(V,1),1);
ec_nodesg = zeros(size(V,1),1);
for i=1:num_caps

    % Get potential end cap nodes
    [~,s] = projectPointOntoPlane(V,ecnorm_uvw(i,:),ecnorm_xyz(i,:));
    ppos_v_inds = find(s==1);
    
    % There may be multiple groups of faces so group them based on connectivity
    [vgroups] = connVertsInSet(ppos_v_inds,C.vertex.vertex);
    
    % Check if there are multiple groups
    if length(unique(vgroups)) > 1
        
        % Find the end cap node nearest to the end cap outlet vector
        [~,min_ind] = min(sqrt(sum((V(ppos_v_inds,:) - ecnorm_xyz(i,:)).^2,2)));
        gid = vgroups(min_ind);
    
        % Find all nodes corresponding to this group
        ppos_v_inds = ppos_v_inds(vgroups==gid);
    
    end
    
    % Update the nodes mask
    ec_nodes(ppos_v_inds) = 1;
    ec_nodesg(ppos_v_inds) = i;
    
end

% Store the end cap nodes
capsv = find(ec_nodes==1);
capsvg = ec_nodesg(capsv);

% Get faces in end caps
capsf = zeros(size(F,1),1);
capsfg = zeros(size(F,1),1);
for i=1:num_caps
    fec_cg = find(all(ismember(F,capsv(capsvg==i)),2));
    capsf(fec_cg) = 1;
    capsfg(fec_cg) = i;
end
capsf = find(capsf);
capsfg = capsfg(capsf);

% Create copy of faces and remove faces connected to end cap nodes
Ff = F;
Ff(any(ismember(F,capsv),2),:) = [];

% Group based on connectivity
optstr.outputType='label';
g = tesgroup(Ff,optstr);
gu = unique(g);

% Should only have 2 surfaces
if length(gu)==1
    warning("Surfaces may not have been separated correctly");
end
    
sa = zeros(length(gu),1);
for i=1:length(gu)
    sa(i) = sum(patchArea(Ff(g==i,:),V));
end
[~,sa_maxi] = max(sa); outeri = gu(sa_maxi);
[~,sa_mini] = min(sa); lumeni = gu(sa_mini);

% sa1 is the outer surface, sa2 is the lumen surface
outerf = find(any(ismember(F,unique(Ff(g==outeri,:))),2));
outere = patchBoundary(F(outerf,:));
lumenf = find(any(ismember(F,unique(Ff(g==lumeni,:))),2));
lumene = patchBoundary(F(lumenf,:));

% Store faces
surfs.outer.F = outerf;
surfs.lumen.F = lumenf;
surfs.caps{1,1}.F = capsf;
surfs.caps{1,1}.Fg = capsfg;

% Store vertices
surfs.outer.V = unique(F(outerf,:));
surfs.lumen.V = unique(F(lumenf,:));
surfs.caps{1,1}.V = capsv;
surfs.caps{1,1}.Vg = capsvg;

% Store edge nodes
surfs.outer.E = outere;
surfs.lumen.E = lumene;

% Store edge node groups
optstr.outputType='label';
surfs.outer.Eg = tesgroup(outere,optstr);
surfs.lumen.Eg = tesgroup(lumene,optstr);

end


% Nice litte function adapted from: https://math.stackexchange.com/a/100766
function [point, sideofplane] = projectPointOntoPlane(P, uvw, xyz)
t = (uvw(1)*xyz(1) - uvw(1).*P(:,1) + uvw(2)*xyz(2) - uvw(2).*P(:,2) + uvw(3)*xyz(3) - uvw(3).*P(:,3))...
    ./ (uvw(1)^2 + uvw(2)^2 + uvw(3)^2);

point = P+t.*uvw;

sideofplane = sideOfPlane(P, uvw, xyz);
end


% output:
%   0: on plane
%   1: on positive side of plane
%   -1: on negative side of plane
function [side] = sideOfPlane(P, plane_uvw, plane_xyz)
dots = zeros(size(P,1),1);
P = P - plane_xyz;
for i=1:size(P,1)
    dots(i) = dot(P(i,:)./norm(P(i,:)),plane_uvw);
end
side = ones(size(P,1),1);
side(dots<0) = -1;
side(dots==0) = 0;
end


function [vgroups] = connVertsInSet(vin,C)
parsed_v = zeros(length(vin),1);
vgroups = zeros(length(vin),1);
group_count=1;
while sum(parsed_v)~=length(parsed_v)
    f = vin(find(parsed_v==0,1));
    conn_v=f;
    len_conn=0;
    while(len_conn<length(conn_v))
        len_conn=length(conn_v);
        vlist = unique([reshape(C(conn_v,:).',[],1); f]);
        conn_v = vlist(ismember(vlist,vin));
    end
    fi = ismember(vin,conn_v);
    parsed_v(fi) = 1;
    vgroups(fi) = group_count;
    group_count = group_count+1;
end
end