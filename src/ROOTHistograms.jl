module ROOTHistograms

using ROOT, Histograms

import ROOT.to_root
function to_root(h::NHistogram, name="hist")
    length(size(h))==2 || error("to_root not defined for nd=$(length(size(h)))")

    nbins_x = nbins(h, 1)
    nbins_y = nbins(h, 2)

    conts, errs = contents(h), errors(h)

    #remove overflow bins
    e1 = convert(Vector{Float64}, h.edges[1][2:end-1])
    e2 = convert(Vector{Float64}, h.edges[2][2:end-1])

    hi = TH2D(name, name, int32(nbins_x - 3), pointer(e1), int32(nbins_y - 3), pointer(e2))
    for i=0:nbins_x-1
        for j=0:nbins_y-1
            SetBinContent(hi, int32(i), int32(j), conts[i + 1, j + 1])
            SetBinError(hi, int32(i), int32(j), errs[i + 1, j + 1])
        end
    end
    SetEntries(hi, sum(entries(h)))
    return hi
end

function to_root(h::Histogram, name="hist")
    edges = copy(h.bin_edges)

    #remove underflow low and overflow high
    edges = h.bin_edges[2:end-1]
    hi = TH1D(name, name, int32(nbins(h) - 3), pointer(edges))
    for i=0:GetNbinsX(hi)+1
        SetBinContent(hi, int32(i), contents(h)[i + 1])
        SetBinError(hi, int32(i), errors(h)[i + 1])
    end
    SetEntries(hi, sum(entries(h)))

    # println(h)
    # Print(hi, "ALL")
    conts, errs, ents, edges = get_hist_bins(hi)

    if !all(conts .==  contents(h))
        warn("mismatch in contents: \n", hcat(conts, contents(h))|>string, "\n")
    end
    if !all(errs .==  errors(h))
        warn("mismatch in errors: \n", hcat(errs, errors(h))|>string, "\n")
    end
    # if !all(ents .==  entries(h))
    #     warn("mismatch in entries: \n", hcat(ents, entries(h))|>string, "\n")
    # end
    if !all(edges .==  h.bin_edges)
        warn("mismatch in edges: \n", hcat(edges, h.bin_edges)|>string, "\n")
    end
    return hi
end

function get_hist_bins(h::Union(TH1D, TH1A, TH1); error_type=:contents)
    nb = int64(GetNbinsX(h))
    nb>0 || error("nbins = $nb")

    #+3 = underflow, overflow, superfluous extra bin
    conts = zeros(Float64, nb+3)
    errs = zeros(Float64, nb+3)
    ents = zeros(Float64, nb+3)

    #underflow low, overflow low, overflow high
    edges = zeros(Float64, nb+3)
    for n=0:nb+1
        conts[n+1] = GetBinContent(h, int32(n))
        errs[n+1] = GetBinError(h, int32(n))
        #entries[n+1] = GetEntries(h) * conts[n+1]
        edges[n+1] = GetBinLowEdge(h, int32(n))
    end

    #this work for histograms for which the bin errors have been manually set
    #to non-Poisson, GetEntries is meaningless
    if error_type == :errors
        ents = (conts ./ errs).^2
        ents[isnan(ents)] = 0
        ents[ents .== Inf] = 1
        #println(hcat(conts, errs, ents, conts./sqrt(ents)))
    end

    #this works for Poisson bin errors
    if error_type == :contents
        ents = conts ./ sum(conts) .* GetEntries(h)
        ents[isnan(ents)] = 0
    end
    #ents = int(round(ents))

    edges[1] = -Inf
    edges[nb+2] = edges[nb+1] + GetBinWidth(h, int32(nb))
    edges[nb+3] = Inf

    if error_type != :errors
        if GetEntries(h)>0
            @assert abs(sum(ents) - GetEntries(h))/sum(ents)<0.0001 string("entries unequal ", sum(ents), "!=", GetEntries(h))
        end
    end
    if Integral(h)>0
        @assert abs(sum(conts) - Integral(h, int32(0), int32(nb+1)))/sum(conts)<0.0001 string("contents unequal ", sum(conts), "!=", Integral(h, int32(0), int32(nb+1)))
    end
    # if (abs(sum(conts) - Integral(h)) > 100000 * eps(Float64))
    #     warn(
    #         GetName(h)|>bytestring,
    #         " integral mismatch: $(sum(conts)) != $(Integral(h))"
    #     )
    # end
    return conts, errs, ents, edges
end

function load_with_errors(f::TDirectoryA, k::ASCIIString; kwargs...)
    th = root_cast(TH1, Get(f, k))
    #println(th)
    conts, errs, ents, edgs = get_hist_bins(th; kwargs...)
    h = Histogram(ents, conts, edgs)
    return h
end

function from_root(o::TH1A)
    conts, errs, ents, edgs = get_hist_bins(o)
    h = Histogram(ents, conts, edgs)
    # for i=1:nbins(h)
    #     @assert abs(errors(h)[i] - errs[i])<0.0001 string("errors unequal $(hcat(errors(h), errs))")
    # end
    return h
end

function from_root(o::TH2A)
    nx, ny = size(o)
    #println("nx=$nx ny=$ny")
    arr = zeros(nx+3, ny+3)
    ents = zeros(nx+3, ny+3)
    for x=0:nx+1
        for y=0:ny+1
            arr[x+1, y+1] =
                GetBinContent(root_cast(TH1, o), int32(x), int32(y))
        end
    end

    ents = arr .* GetEntries(o) ./ sum(arr)
    ents[isnan(ents)] = 0
    ents = round(ents)

    return NHistogram(ents, arr, {[0:nx+2], [0:ny+2]})
end


import Base.size
Base.size(o::TH2D) = (GetNbinsX(root_cast(TH1, o)), GetNbinsY(root_cast(TH1, o)))

function load_hists_from_file(fn, hfilter=(name->true))
    tf = TFile(fn)
    @assert tf.p != C_NULL

    kl = GetListOfKeys(tf)
    @assert kl.p != C_NULL

    key_iterator = TListIter(kl.p)
    #kl = GetListOfKeys(tf)
    #objs = GetList(tf)

    ret = Dict()

    tic()
    for i=1:length(kl)
        #(i % 10000 == 0) && (println("$(i)/$(length(kl)) ", toq());tic())
        _k = Next(key_iterator)
        @assert _k != C_NULL
        const k = TKey(_k)
        const n = GetName(k) |> bytestring
        hfilter(n) || continue
        const o = ReadObj(k) |> to_root
        ret[n] = from_root(o)
    end

    Close(tf)
    return ret
end

function write_hists_to_file(hd::Associative, fn)
    tf = TFile(convert(ASCIIString, fn), "RECREATE")
    Cd(tf, "")
    for (k, v) in hd
        v::Union(Histogram, NHistogram)
        println("$k N=$(sum(entries(v)))")
        hi = to_root(v, string(k))
    end
    Write(tf)
    Close(tf)
end

export to_root, get_hist_bins, load_hists_from_file, from_root, write_hists_to_file
export load_with_errors

end
