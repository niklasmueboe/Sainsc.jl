#=
# Example analysis

This is an example workflow demonstrating all the different functionalities the 
`StereoSSAM.jl` package offers.

## General workflow

We start by loading the package and defining the path of some example data.
=#

using StereoSSAM

using ImageTransformations
using ImageShow
using Printf

data_path = "./data";

stereoseq_file = joinpath(data_path, "Mouse_brain_Adult_GEM_bin1.tsv.gz");
signature_file = joinpath(data_path, "signatures.csv");
color_file = joinpath(data_path, "celltype_colors.json");

#=
### Loading data

We read the data from file using [`readstereoseq`](@ref). The input will 
usually be a GEM file of a StereoSeq experiment. The output is a 
[`GridCounts`](@ref) of gene counts stored as 
[`SparseArrays.SparseMatrixCSC`](@extref). The input data stems from the
original [Stereo-seq publication](https://doi.org/10.1016/j.cell.2022.04.003) and can
be downloaded [here](https://db.cngb.org/stomics/mosta/download/).
=#

counts = readstereoseq(stereoseq_file)

#=
We also load pre-defined cell-type signatures
=#

using CSV
using DataFrames

## read pre-determined cell-type signatures
signatures = CSV.read(signature_file, DataFrame);
celltypes = signatures.Celltype;
select!(signatures, Not(:Celltype));

## remove genes that are not detected in Stereo-seq experiment
signatures = signatures[:, names(signatures) .∈ [keys(counts)]];

print(signatures[1:5, 1:8])

#=
We can have a look at the overall structure by calculating the [`totalrna`](@ref) 
which is basically the sum across all genes per pixel.
=#

total_counts = totalrna(counts);
simshow(imresize(total_counts; ratio=1 / 45))

#=
### Kernel

Next, we define a [`gaussiankernel`](@ref) to use for smoothing and thus integrating 
the gene expression locally across pixels. The relevant parameters are the abndwidth of 
the kernel and the radius i.e. how large the kernel will be measured in bandwidths.
For example, setting the bandwidth ``\sigma=8`` and the radius ``r=2`` will result in a 
kernel of size ``2r*\sigma+1=33``.

We can also convert the kernel to `Float32` to reduce memory usage later on.
=#

kernel = gaussiankernel(8, 2);
kernel = convert.(Float32, kernel);

#=
We can evaluate the effects by e.g. applying the kernel to the totalRNA.
=#

total_rna = kde(Matrix(totalrna(counts)), kernel);
simshow(imresize(total_rna; ratio=1 / 45))

#=
### Local Maxima

The smoothed totalRNA is used to find local maxima of that can be used as cell proxies 
to identify gene signatures or for further analysis.

[`findlocalmaxima`](@ref) allows you to specify a minimum distance between to 
neighboring local maxima.
=#

lm = findlocalmaxima(total_rna, 6);
println("$(length(lm)) local maxima detected")

#=
Once we identified the local maxima we can load them for a defined set of genes. 
[`getlocalmaxima`](@ref) allows you to load the maxima in a variety of output types.
Have a look at the existing `StereoSSAM.jl` extensions. Here we are going to load them
as `Muon.AnnData` object.
=#

using Muon

adata = getlocalmaxima(AnnData, counts, lm, kernel; genes=names(signatures))

#=
### Cell-type map

Once we have identified gene signatures from literature, previous experiments, 
or by analysing local maxima we can assign a celltype to each pixel. This is usually the 
most time-intensive processing step. [`assigncelltype`](@ref) returns 2 matrices;
a map of assigned celltypes of each pixel as 
`CategoricalArrays.CategoricalMatrix` and the cosine similarity of the 
assigned cell type for each pixel as matrix.
=#

celltype_map, cosine_sim = assigncelltype(counts, signatures, kernel; celltypes=celltypes);

#=
In the final step we can visualize the cell-type map.
=#

using Colors
using JSON

## generate the lookup-table for cell-type colors
lut = Dict(Iterators.map(((k, v),) -> k => RGB(v...), pairs(JSON.parsefile(color_file))));

background_color = RGB(0, 0, 0);

#-

celltype_img = map(x -> get(lut, x, background_color), celltype_map);
simshow(imresize(celltype_img; ratio=1 / 45))

#=
When visualizing the cell-type map, it might be beneficial to remove background noise.
A histogram of the totalRNA can be helpful to define a threshold.
=#

using SparseArrays
using Plots

histogram(nonzeros(sparse(total_rna)); title="KDE of totalRNA", bins=200)

#-

celltype_img[total_rna .< 1.5] .= background_color;

simshow(imresize(celltype_img; ratio=1 / 45))

#=
## Utils

Some additional useful functions are included in `StereoSSAM.jl`.

### Cropping & Masking

[`crop!`](@ref) allows you to slice the counts across all genes in x-y dimension.
=#
crop!(counts, (4_001:7_000, 6_001:9_000));

simshow(imresize(kde(totalrna(counts), kernel); ratio=1 / 10))

#=
[`mask!`](@ref) allows you to filter the field of view by an arbitrary binary mask.
=#

## create a mask
roi_mask = zeros(Bool, 3_000, 3_000)
roi_mask[1:1_500, 1_501:3_000] .= true
roi_mask[1_501:3_000, 1:1_500] .= true

mask!(counts, roi_mask);

simshow(imresize(kde(totalrna(counts), kernel); ratio=1 / 10))

#=
### Binning

The same extensions that are available to load local maxima are also available for 
reading the StereoSeq data into bins using [`readstereoseqbinned`](@ref). This can be 
useful for e.g. identifying gene signatures from binned data, or detecting 
highly variable genes that can be used for analysis. We again here demonstrate it using
`Muon.AnnData`.
=#

adata = readstereoseqbinned(AnnData, stereoseq_file, 50)