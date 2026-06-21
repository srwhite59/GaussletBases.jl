#!/usr/bin/env julia

using GaussletBases
using LinearAlgebra

fixture = joinpath(
    dirname(@__DIR__),
    "test/driver_inputs/h2_pqs_q5_independent_source_box_r4_base_hamiltonian.jl",
)
hamfile = tempname() * ".jld2"

empty!(ARGS)
append!(ARGS, [
    fixture,
    "save_artifact=false",
    "save_tsv=false",
    "save_basis_artifact=false",
    "save_ham_artifact=true",
    "materialize_route=true",
    "hamiltonian_output=:cartesian_ida_hamiltonian",
    "hamfile=$(repr(hamfile))",
])

harness_elapsed = @elapsed include(joinpath(@__DIR__, "cartesian_driver_harness.jl"))

const H2_H1_LOWEST = -0.79460371733658908
const H2_SELF_COULOMB = 0.4569117646737212
const ATOL = 1.0e-10

function _check_le(label, observed, bound)
    observed <= bound || error("$label expected <= $bound, got $observed")
    println("  ", label, " = ", observed)
    return nothing
end

function _self_coulomb(V, orbital)
    density = orbital * transpose(orbital)
    rho = 0.5 .* (density .+ transpose(density))
    occupations = vec(diag(rho))
    sym = 0.5 .* (V .+ transpose(V))
    return 2.0 * dot(occupations, sym * occupations) -
        dot(vec(rho), vec(sym .* rho))
end

function _check_hamiltonian(ham)
    typeof(ham) == CartesianIDAHamiltonian{Float64} ||
        error("materialization did not return CartesianIDAHamiltonian{Float64}")
    size(ham.kinetic) == (471, 471) || error("H2 final dimension mismatch")
    (ham.nup, ham.ndn) == (1, 1) || error("H2 electron count mismatch")
    ham.nuclear_charges == [1.0, 1.0] || error("H2 nuclear charge mismatch")
    ham.nuclear_positions == [0.0 0.0 -2.0; 0.0 0.0 2.0] ||
        error("H2 nuclear position mismatch")

    _check_le("K_symmetry", norm(ham.kinetic - transpose(ham.kinetic), Inf), ATOL)
    _check_le("V_symmetry",
        norm(ham.electron_electron_ida - transpose(ham.electron_electron_ida), Inf), ATOL)
    for (index, matrix) in enumerate(ham.nuclear_attraction_unit_by_center)
        all(isfinite, matrix) || error("U_$index contains nonfinite values")
        _check_le("U_$(index)_symmetry", norm(matrix - transpose(matrix), Inf), ATOL)
    end
    all(isfinite, ham.kinetic) || error("K contains nonfinite values")
    all(isfinite, ham.electron_electron_ida) || error("V contains nonfinite values")

    H = one_body_hamiltonian(ham)
    eig = eigen(Symmetric((H + transpose(H)) ./ 2))
    lowest_index = argmin(eig.values)
    lowest = eig.values[lowest_index]
    orbital = eig.vectors[:, lowest_index]
    self_coulomb = _self_coulomb(ham.electron_electron_ida, orbital)
    _check_le("h1_lowest_delta", abs(lowest - H2_H1_LOWEST), ATOL)
    _check_le("self_coulomb_delta", abs(self_coulomb - H2_SELF_COULOMB), ATOL)
    return H
end

println("h2_pqs_base_hamiltonian_smoke_checks")
direct_H = _check_hamiltonian(materialization)
isfile(hamfile) || error("Hamiltonian artifact was not written")
reloaded = read_cartesian_ida_hamiltonian(hamfile)
reloaded_H = _check_hamiltonian(reloaded)
_check_le("readback_one_body_delta", norm(direct_H - reloaded_H, Inf), 0.0)
println("h2_pqs_base_hamiltonian_smoke_elapsed_s=", harness_elapsed)
