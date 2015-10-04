using PyPlot

#we can overload the default plot method
import PyPlot.bar
function bar{T <: Any,M <: Any}(h::ErrorHistogram{T,1,M}; kwargs...)
    es = edges(h)[1]
    cs = contents(h)
    PyPlot.bar(es[1:end-1], cs, width=diff(es); kwargs...)
end

import PyPlot.plot
function plot{T <: Any,M <: Any}(h::ErrorHistogram{T,1,M}; kwargs...)
    es = edges(h)[1]
    cs = contents(h)
    PyPlot.plot(Base.midpoints(es), cs; kwargs...)
end

import PyPlot.errorbar
function errorbar{T <: Any,M <: Any}(h::ErrorHistogram{T,1,M}; kwargs...)
    es = edges(h)[1]
    cs = contents(h)
    errs = errors(h)
    PyPlot.errorbar(Base.midpoints(es), cs, errs; kwargs...)
end
