# Consume the output of the -S option of nplotnn.

using ArgMacros

using CSV
using Plots
using Dierckx

# This supports the --dist option.
function crossing_point(values::Array, cutoff::Real)::Int
  for i in eachindex(values)
    values[i] <= cutoff && return i
  end

  return 1;
end

# This supports the --points option
# Assumes that `values` is an output from spread.
# Return the indices at which `values` crosses the values implied by `delta`.
# Delta is a distance, like 0.1. In that case we return the indices at which
# `values` crosses [0.9, 0.8, 0.7...]
function crossing_points(values::Array, delta::Real)::Array{Int}
  result_size = convert(Int, round(1.0 / delta)) + 1
  result = zeros(Int32, result_size)

  dists = collect(range(1.0, stop=0.0, length=convert(Int, round(1.0 / delta) + 1)))
  # An index into the `dists` array. The last value will be 0, so skip that.
  ndx = length(dists) - 1
  # Scan in reverse order
  for i in length(values):-1:1
    v = values[i]
    v == 0.0 && continue
    if v > dists[ndx]
      result[ndx] = i
      ndx -= 1
    end
  end

  # get rid of zero's
  return filter(e -> e > 0, result)

end

function main()
  @inlinearguments begin
    @helpusage "nplotnn -S <file> file.spr > file.spr.smi; spreadplot.jl -S spread <file>    generates spread.png"
    @helpdescription """Consumes the output of the -S option of nplotnn, generating a plot of distance vs number selected.
    -S : specifies the name of the output file - .png will be added.
    --title : option is the title of the plot, --title 'Title'
    --label : legend of the curve, --legend 'foo'
    --lwd : controls the width of the line, --lwd 8
    --color : the color of the line, --color red
    --dist : draw a horizontal line at the point where a distance is crossed, --dist 0.15
    --spline : use spline interpolation. Useful for large datasets and svg output, --spline <npoints>
    --point : draw dots on the curve every specified interval, --point 0.1
    The 
    """
    @argumentrequired String stem "-S"
    @argumentoptional Float32 threshold "--dist"
    @argumentoptional String label "--label"
    @argumentoptional String title "--title"
    @argumentoptional Int32 spline "--spline"
    @argumentoptional Float32 point "--point"
    @argumentdefault String "png" format "--format"
    @argumentdefault Int32 4 lwd "--lwd"
    @argumentdefault String "green" colour "--color"

    @positionalrequired String input_file "input_file"
  end

  data = CSV.File(input_file, header=true, delim= ",")
  distance = data.distance
  if isnothing(label)
    label = "Selected"
  end
  
  if isnothing(spline)
    plot(range(1,length(distance)), distance, ylim=(0.0, distance[2]), xlabel="Number Selected",
         ylabel="Distance", lwd=lwd, color=Symbol(colour), label=label)
  else
    x = range(1, stop=length(distance), length=spline)
    spl = Spline1D(range(1,length(distance)), distance)
    tmp = [spl(q) for q in x]
    plot(x, tmp, ylim=(0.0, distance[2]), xlabel="Number Selected",
         ylabel="Distance", lwd=lwd, color=Symbol(colour), label=label)
  end

   if ! isnothing(title)
     plot!(title=title)
   end

   if ! isnothing(point)
     x = crossing_points(distance, point)
     y = [distance[i] for i in x]
     scatter!(x, y, label=false)
   end

  if ! isnothing(threshold)
    ndx = crossing_point(distance, threshold)
    y = distance[ndx]
    plot!([0, length(distance)], [y, y], linestyle=:dash, color=:grey, label="Dist $(threshold) count $(ndx)")
  end

  savefig("$(stem).$(format)")
end

main()

