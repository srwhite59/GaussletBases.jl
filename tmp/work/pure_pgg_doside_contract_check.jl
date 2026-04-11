push!(LOAD_PATH, "/Users/srw/Dropbox/GaussletModules")

using Test
using Gausslets
using PureGaussianGausslet

function build_case(doside::Int)
    atom = PureGaussianGausslet.Atom(1, 0.0, 0.0, 0.0, 0.5)
    gaulet = gausslet(cf1092)
    bas1dx, bas1dy, bas1dz, bas3dd, shells = mktempdir() do dir
        cd(dir) do
            PureGaussianGausslet.getNGgaussletonly(
                gaulet,
                [atom],
                doside,
                0.5,
                2.0,
                8.0;
                dwidth = 8.0,
                doInvsqrt = true,
            )
        end
    end
    return (
        basis_dim = length(bas3dd),
        n1x = size(bas1dx.all1dbmat, 2),
        n1y = size(bas1dy.all1dbmat, 2),
        n1z = size(bas1dz.all1dbmat, 2),
        shell_count = length(shells),
    )
end

low = build_case(3)
high = build_case(5)

@test low.basis_dim < high.basis_dim
@test low != high

println("pure-pgg-doside-contract-ok")
println((doside3 = low, doside5 = high))
