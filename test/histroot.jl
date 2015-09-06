using ROOT, Histograms, ROOTHistograms
using Base.Test

h = ErrorHistogram([-Inf, 0,1,2,3,4, Inf])
for i=1:100
    push!(h, 0.3)
end
th = to_root(h, "hist1")
_h = from_root(th)
@test all(edges(h)[1] .== edges(_h)[1])
@test all(contents(h) .== contents(_h))
@test all(errors(h) .== errors(_h))

h = ErrorHistogram([-Inf, 0,1,2,3,4, Inf])
for i=1:100
    push!(h, 0.3, 0.1)
end
th = to_root(h, "hist2")
_h = from_root(th)
@test all(edges(h)[1] .== edges(_h)[1])
@test all(contents(h) .== contents(_h))
@test all(errors(h) .== errors(_h))

nh = ErrorHistogram([-Inf, 0,1,2, Inf], [-Inf, 0,1,2, Inf])
for i=1:100
    push!(nh, (0.3, 0.3))
end
th = to_root(nh, "nhist")
nh2 = from_root(th)
@test nbins(nh) == nbins(nh2)
@test all(edges(nh)[1] .== edges(nh2)[1])
@test all(edges(nh)[2] .== edges(nh2)[2])

@test all(errors(nh) .== errors(nh2))
@test all(contents(nh) .== contents(nh2))

tf = ROOT.TFile("tfile.root", "RECREATE")
h = ErrorHistogram([-Inf, linspace(-1, 1, 30)..., Inf])

for i=1:1000000
    push!(h, rand(-1:0.0001:1))
end
th = to_root(h, "hist3")
h2 = from_root(th)
Write(th)
Close(tf)

hd = load_hists_from_file("tfile.root")
h = hd["hist3"]

ROOT.gROOT.process_line("TFile* tf = new TFile(\"tfile.root\");")
ROOT.gROOT.process_line("TH1D* h = (TH1D*)(tf->Get(\"hist3\"));")
for i=1:nbins(h)[1]
    @test contents(h)[i] == ROOT.gROOT.process_line("h->GetBinContent($(i-1));")
end
