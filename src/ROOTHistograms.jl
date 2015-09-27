module ROOTHistograms

using ROOT, Histograms

import ROOT.to_root
function to_root{T, E}(h::ErrorHistogram{T, 2, E}, name="hist")
    nbins_x = nbins(h, 1)
    nbins_y = nbins(h, 2)

    conts, errs = contents(h), errors(h)

    #remove overflow bins
    e1 = convert(Vector{Float64}, h.edges[1][2:end-1])
    e2 = convert(Vector{Float64}, h.edges[2][2:end-1])

    hi = TH2D(name, name, Int32(nbins_x - 3), pointer(e1), Int32(nbins_y - 3), pointer(e2))
    for i=0:nbins_x-2
        for j=0:nbins_y-2
            SetBinContent(hi, Int32(i), Int32(j), conts[i + 1, j + 1])
            SetBinError(hi, Int32(i), Int32(j), errs[i + 1, j + 1])
        end
    end
    SetEntries(hi, sum(entries(h)))
    return hi
end

function to_root{T, E}(h::ErrorHistogram{T, 1, E}, name="hist")

    #remove underflow low and overflow high
    _edges = collect(edges(h)[1])[2:end-1]
    hi = TH1D(name, name, Int32(nbins(h)[1] - 2), pointer(_edges))
    for i=0:GetNbinsX(hi)+1
        SetBinContent(hi, Int32(i), contents(h)[i + 1])
        SetBinError(hi, Int32(i), errors(h)[i + 1])
    end
    SetEntries(hi, sum(entries(h)))

    # # println(h)
    # # Print(hi, "ALL")
    # conts, errs, ents, _edges2 = get_hist_bins(hi)
    # 
    # if !all(conts .==  contents(h))
    #     warn("mismatch in contents: \n", hcat(conts, contents(h))|>string, "\n")
    # end
    # if !all(errs .==  errors(h))
    #     warn("mismatch in errors: \n", hcat(errs, errors(h))|>string, "\n")
    # end
    # # if !all(ents .==  entries(h))
    # #     warn("mismatch in entries: \n", hcat(ents, entries(h))|>string, "\n")
    # # end
    # if !all(_edges2 .==  _edges)
    #     warn("mismatch in edges: \n", hcat(_edges, edges(h)[1])|>string, "\n")
    # end
    return hi
end

function get_hist_bins(h::Union{TH1D, TH1A, TH1}; error_type=:contents)
    nb = Int64(GetNbinsX(h))
    nb>0 || error("nbins = $nb")

    #+3 = underflow, overflow, superfluous extra bin
    conts = zeros(Float64, nb+2)
    errs = zeros(Float64, nb+2)
    ents = zeros(Float64, nb+2)

    #underflow low, overflow low, overflow high
    edges = zeros(Float64, nb+3)
    for n=0:nb+1
        #println("$n ", GetBinLowEdge(h, Int32(n)), " ", GetBinContent(h, Int32(n)), " ", GetBinError(h, Int32(n)))
        conts[n+1] = GetBinContent(h, Int32(n))
        errs[n+1] = GetBinError(h, Int32(n))
        #entries[n+1] = GetEntries(h) * conts[n+1]
        edges[n+1] = GetBinLowEdge(h, Int32(n))
    end

    #this works for histograms for which the bin errors have been manually set
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
    edges[nb+2] = edges[nb+1] + GetBinWidth(h, Int32(nb))
    edges[nb+3] = Inf

    if error_type != :errors
        if GetEntries(h)>0
            @assert abs(sum(ents) - GetEntries(h))/sum(ents)<0.0001 string("entries unequal ", sum(ents), "!=", GetEntries(h))
        end
    end
    if Integral(h)>0
        @assert abs(sum(conts) - Integral(h, Int32(0), Int32(nb+1)))/sum(conts)<0.0001 string("contents unequal ", sum(conts), "!=", Integral(h, Int32(0), Int32(nb+1)))
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
    h = ErrorHistogram{Float64, 1, Tuple{Vector{Float64}}}(
        (edgs, ), conts, errs.^2, :left
    )
    # for i=1:nbins(h)
    #     @assert abs(errors(h)[i] - errs[i])<0.0001 string("errors unequal $(hcat(errors(h), errs))")
    # end
    return h
end

function from_root(o::TH2A)
    nx, ny = size(o)
    #println("nx=$nx ny=$ny")
    conts = zeros(nx+2, ny+2)
    errs = zeros(nx+2, ny+2)
    edgs = (zeros(Float64, nx+3), zeros(Float64, ny+3))
    for x=0:nx+1
        edgs[1][x+1] = GetBinLowEdgeX(o, Int32(x))
    end
    for y=0:ny+1
        edgs[2][y+1] = GetBinLowEdgeY(o, Int32(y))
    end
    edgs[1][1] = -Inf
    edgs[1][end] = Inf
    edgs[2][1] = -Inf
    edgs[2][end] = Inf
    for x=0:nx+1
        for y=0:ny+1
            conts[x+1, y+1] =
                GetBinContent(root_cast(TH1, o), Int32(x), Int32(y))
            errs[x+1, y+1] =
                GetBinError(root_cast(TH1, o), Int32(x), Int32(y))
        end
    end

    return ErrorHistogram{Float64, 2, Tuple{Vector{Float64}, Vector{Float64}}}(
        edgs, conts, errs.^2, :left
    )
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

# function recurse_write(tf, hd::Associative, pref="")
#     if pref != ""
#         println("making $pref")
#         d = ROOT.mkdir(tf, pref)
#         @show d
#         if ROOT.root_pointer(d) == C_NULL
#             d2 =  root_cast(TDirectory, Get(tf, pref))
#             @show d2
#         end
#         d = root_cast(TDirectory, d)
#         assert(d.p != C_NULL)
#     else
#         d = tf
#     end
#     Cd(tf, pref)
#     for (k, v) in hd
#         s = length(pref)>0 ? "$pref/$k" : string(k)
#         recurse_write(d, v, s)
#     end
#     pref!= "" && Cd(tf, "..")
# end

# function recurse_write(tf::TDirectoryA, v::ErrorHistogram, pref="")
#     th = to_root(v)
#     @show tf
#     SetDirectory(th, ROOT.root_pointer(tf))
#     Write(th)
# end

function write_hists_to_file(fn, hd::Associative; verbose=false)
    tf = TFile(convert(ASCIIString, fn), "RECREATE")
    Cd(tf, "")
    for (objpath, obj) in hd
        spl = rsearch(objpath, "/")
        dirname = objpath[1:spl.start-1]
        objname = objpath[spl.start+1:end]
        if length(dirname) > 0
            verbose && println("dir=$dirname obj=$objname")
            d = Get(tf, dirname)
            if d.p == C_NULL
                verbose && println("making $dirname")
                d = ROOT.mkdir(tf, dirname)
            end
            d = root_cast(TDirectory, Get(tf, dirname))
            println(d) 
            assert(d.p != C_NULL)
        else
            d = tf
        end
        verbose && println("TDirectory=$d")
        robj = to_root(obj, objname)
        SetDirectory(robj, convert(Ptr{TDirectory}, d.p))
    end
    nb = Write(tf)
    verbose && println("Wrote $nb bytes")
    Close(tf)
end

export to_root, get_hist_bins, load_hists_from_file, from_root, write_hists_to_file
export load_with_errors

end
