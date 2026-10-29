function testPowerBasedSpectrum()
% Compare both implementations with the wave generator's JONSWAP spectrum.
root = fileparts(fileparts(mfilename('fullpath')));
mssRoot = fullfile(fileparts(root),'MSS');
helpers = tempname;
mkdir(helpers);
originalPath = path;
cleanup = onCleanup(@() cleanupHelpers(helpers,originalPath)); %#ok<NASGU>
copyfile(fullfile(root,'matlab','private','computeVesselDamping.m'),helpers);
addpath(genpath(mssRoot));
addpath(helpers);

v.freqs = [0 linspace(.05,4,80) 10];
v.MRB = diag([100 120 140 160 180 200]);
v.main.submerged = false;
nf = numel(v.freqs);
v.A = zeros(6,6,nf);
v.B = zeros(6,6,nf);
v.C = repmat(diag([0 0 40 30 80 0]),1,1,nf);
for i = 1:6
    v.A(i,i,:) = 5 + i*v.freqs + .2*i*v.freqs.^2;
    v.B(i,i,:) = .5 + .1*i*v.freqs.^2;
end
v.A(2,4,:) = .01*v.freqs;
v.A(4,2,:) = v.A(2,4,:);
v.B(2,4,:) = .001*v.freqs;
v.B(4,2,:) = v.B(2,4,:);
v.powerBased = struct('T_126',[1 1 1],'delta_zeta_345',[0 .1 0]);

for wp = [.3 .7 1 2]
    standalone = computeVesselDamping(v,wp);
    mss = computeManeuveringModel(v,wp,[1 1 1],[0 .1 0],0);
    compareImplementations(standalone,mss);
    assert(isequal(standalone.A,v.A) && isequal(standalone.B,v.B));
    assert(isequal(standalone.C,v.C));

    % Independent dense integration using the actual wave-generation function.
    w = linspace(1e-6,4,20001)';
    spectrum = waveSpectrum(7,[1 wp 3.3],w,0);
    weight = spectrum / trapz(w,spectrum);
    finite = v.freqs < 10;
    for i = 1:6
        a = interp1(v.freqs(finite),reshape(v.A(i,i,finite),[],1),w,'pchip');
        b = interp1(v.freqs(finite),reshape(v.B(i,i,finite),[],1),w,'pchip');
        expectedA = trapz(w,a.*weight);
        expectedB = trapz(w,b.*weight);
        assert(abs(mss.powerBased.A_eq(i,i)-expectedA) < 1e-5*expectedA);
        assert(abs(mss.powerBased.B_eq(i,i)-expectedB) < 1e-5*expectedB);
    end
    [~,peak] = max(spectrum);
    assert(abs(w(peak)-wp) < 3e-4);
    assert(mss.powerBased.Bv(3,3)==0 && mss.powerBased.Bv(5,5)==0);
    expectedRoll = .2*sqrt(mss.M(4,4)*mss.G(4,4));
    assert(abs(mss.powerBased.Bv(4,4)-expectedRoll) < 1e-12*expectedRoll);
end

% Missing/empty peaks use the same fixed reference, even with an old cache.
v.powerBased.omega_p = .7;
standaloneDefault = computeVesselDamping(v);
mssDefault = computeManeuveringModel(v);
compareImplementations(standaloneDefault,mssDefault);
assert(mssDefault.powerBased.omega_p==1);
explicitDefault = computeManeuveringModel(v,1,[1 1 1],[0 .1 0],0);
assert(isequal(mssDefault.powerBased.A_eq,explicitDefault.powerBased.A_eq));
emptyDefault = computeManeuveringModel(v,[],[1 1 1],[0 .1 0],0);
assert(isequal(emptyDefault.powerBased.B_eq,explicitDefault.powerBased.B_eq));

% Infinite-frequency markers must not affect the spectrum average.
changed = v;
changed.A(:,:,end) = 1e30;
changed.B(:,:,end) = 1e30;
marker = computeManeuveringModel(changed,.7,[1 1 1],[0 .1 0],0);
reference = computeManeuveringModel(v,.7,[1 1 1],[0 .1 0],0);
assert(isequal(marker.powerBased.A_eq,reference.powerBased.A_eq));
assert(isequal(marker.powerBased.B_eq,reference.powerBased.B_eq));

% Keep the submerged four-time-constant formulation aligned as well.
v.main.submerged = true;
v.C(3,3,:) = 0;
v.powerBased = struct('T_1236',[1 1 1 1],'delta_zeta_45',[.1 0]);
standalone = computeVesselDamping(v,.7);
mss = computeManeuveringModel(v,.7,[1 1 1 1],[.1 0],0);
compareImplementations(standalone,mss);
fprintf('JONSWAP baseline, independent spectrum integration, shared damping and MSS alignment passed.\n');
end

function compareImplementations(standalone,mss)
assert(strcmp(mss.powerBased.spectrumType,'JONSWAP'));
assert(mss.powerBased.gamma==3.3);
fields = {'A_eq','B_eq','Bv'};
for k = 1:numel(fields)
    name = fields{k};
    a = standalone.powerBased.(name);
    b = mss.powerBased.(name);
    assert(norm(a-b,'fro') < 1e-12*max(1,norm(a,'fro')));
end
assert(norm(standalone.D-mss.D,'fro') < 1e-12*norm(mss.D,'fro'));
end

function cleanupHelpers(folder,originalPath)
path(originalPath);
rmdir(folder,'s');
end
