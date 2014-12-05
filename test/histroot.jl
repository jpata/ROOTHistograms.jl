using ROOT, Histograms, ROOTHistograms
using Base.Test

h = Histogram([-Inf, 0,1,2,3,4, Inf])
for i=1:100
    hfill!(h, 0.3)
end
th = to_root(h, "hist1")
@test from_root(th) == h

h = Histogram([-Inf, 0,1,2,3,4, Inf])
for i=1:100
    hfill!(h, 0.3, 0.1)
end
th = to_root(h, "hist2")
@test from_root(th) == h


nh = NHistogram({[-Inf, 0,1,2, Inf], [-Inf, 0,1,2, Inf]})
for i=1:100
    hfill!(nh, {0.3, 0.3})
end
th = to_root(nh, "nhist")
nh2 = from_root(th)
@test size(nh) == size(nh2)
@test all(nh2.baseh.bin_entries .== nh.baseh.bin_entries)
@test all(nh2.baseh.bin_contents .== nh.baseh.bin_contents)
@test all(contents(nh) .== contents(nh2))
@test all(entries(nh) .== entries(nh2))

tf = ROOT.TFile("tfile.root", "RECREATE")
h = Histogram([-Inf, linspace(-1, 1, 30)..., Inf])

for i=1:1000000
    hfill!(h, rand(-1:0.0001:1))
end
th = to_root(h, "hist3")
h2 = from_root(th)
Write(th)
Close(tf)

hd = load_hists_from_file("tfile.root")
h = hd["hist3"]

ROOT.gROOT.process_line("TFile* tf = new TFile(\"tfile.root\");")
ROOT.gROOT.process_line("TH1D* h = (TH1D*)(tf->Get(\"hist3\"));")
for i=1:nbins(h)
    @test contents(h)[i] == ROOT.gROOT.process_line("h->GetBinContent($(i-1));")
end
