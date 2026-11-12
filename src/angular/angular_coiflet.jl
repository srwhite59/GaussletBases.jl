function _coiflet_weights_digest(weights::AbstractVector{<:Real}, expected_len::Int, whiten_target::Bool)
    length(weights) == expected_len ||
        throw(ArgumentError("cubature_weights length ($(length(weights))) must match cardinality ($expected_len)"))
    all(w -> isfinite(w) && w > 0.0, weights) ||
        throw(ArgumentError("cubature_weights must be finite and strictly positive"))
    w64 = Float64[Float64(w) for w in weights]
    all(w -> isfinite(w) && w > 0.0, w64) ||
        throw(ArgumentError("converted cubature_weights must be finite and strictly positive"))
    abs(sum(w64) - 4.0 * pi) <= (whiten_target ? 1.0e-5 : 1.0e-11) ||
        throw(ArgumentError("cubature_weights must sum to 4*pi (got $(sum(w64)))"))
    return bytes2hex(SHA.sha1(reinterpret(UInt8, w64)))
end

function _polar_coiflet_isometry(
    parent_moments::AbstractMatrix{<:Real},
    target_moments::AbstractMatrix{<:Real},
    ;
    whiten_target::Bool = false,
)
    nf, m_parent = size(parent_moments)
    nf_tar, n_coarse = size(target_moments)
    (nf == nf_tar && m_parent == n_coarse && n_coarse >= nf) || throw(DimensionMismatch("invalid coiflet moment dimensions"))
    raw_orth_err = opnorm(target_moments * transpose(target_moments) - Matrix{Float64}(I, nf, nf), Inf)
    if whiten_target
        raw_orth_err <= 1.0e-5 || throw(ArgumentError("target moments exceed whiten_target tol ($raw_orth_err)"))
        ft = svd(Matrix{Float64}(target_moments); full = false)
        target_eff = ft.U * ft.Vt
    else
        raw_orth_err <= 1.0e-11 || throw(ArgumentError("target coiflet moments not row-orthonormal ($raw_orth_err)"))
        target_eff = Matrix{Float64}(target_moments)
    end
    f_parent = svd(Matrix{Float64}(parent_moments); full = false)
    minimum(f_parent.S) >= 1.0 - 1.0e-11 ||
        throw(ArgumentError("parent basis lacks exact Y_lm through l_flat (min_sigma=$(minimum(f_parent.S)))"))
    s_parent = Float64[w >= 0.0 ? 1.0 : -1.0 for w in view(parent_moments, 1, :)]
    proto = Matrix{Float64}(Diagonal(s_parent))
    m0_iso = f_parent.U * f_parent.Vt
    q0, r0 = nullspace(m0_iso), nullspace(target_eff)
    f_comp = svd(transpose(q0) * proto * r0; full = false)
    u_iso = transpose(m0_iso) * target_eff + q0 * (f_comp.U * f_comp.Vt) * transpose(r0)
    maximum(abs, parent_moments * u_iso - target_eff) <= 1.0e-11 ||
        throw(ArgumentError("constructed coiflet rotation failed exact moment check"))
    min_align = minimum(s_parent .* diag(u_iso))
    min_align >= 0.90 ||
        throw(ArgumentError("coiflet rotation lost parent site alignment (min_align=$min_align)"))
    return u_iso, min_align, raw_orth_err
end

function _build_coiflet_profile_from_parent(
    parent_profile::ShellLocalAngularProfile,
    coarse_coords::AbstractMatrix{<:Real},
    cubature_weights::AbstractVector{<:Real},
    l_flat::Int,
    key::ShellLocalAngularProfileKey;
    whiten_target::Bool = false,
)
    pb = parent_profile.basis
    pb.l_inject >= l_flat >= 0 || throw(ArgumentError("coiflet shell requires l_inject >= l_flat >= 0"))
    n_coarse = size(coarse_coords, 1)
    flat_channels = ylm_channels(l_flat)
    parent_moments = reduce(vcat, [_shell_local_ylm_moment_block(pb, l) for l in 0:l_flat])
    w_sqrt = sqrt.(Float64.(cubature_weights))
    target_moments = [w_sqrt[a] * _real_spherical_harmonic(flat_channels[c], view(coarse_coords, a, :))
                      for c in eachindex(flat_channels.channel_data), a in 1:n_coarse]
    u_iso, min_align, raw_orth_err =
        @timeg "coiflet polar rotation" _polar_coiflet_isometry(parent_moments, target_moments; whiten_target = whiten_target)
    ov_r = transpose(u_iso) * pb.final_overlap * u_iso
    kin_r = transpose(u_iso) * pb.final_kinetic * u_iso
    final_kin = 0.5 .* (kin_r .+ transpose(kin_r))
    inj_ov = pb.injected_overlap * u_iso
    inj_kin_r = inj_ov * final_kin * transpose(inj_ov)
    basis = ShellLocalInjectedAngularBasis(
        pb.point_set, pb.beta, copy(pb.shell_weights), pb.theta_nn, pb.kappa,
        pb.l_inject, pb.injected_channels, pb.prototype_count, pb.injected_count,
        pb.whitened_complement_count, n_coarse, pb.prototype_overlap, pb.prototype_kinetic,
        pb.prototype_orthonormalizer, pb.ylm_prototype_coupling, pb.raw_to_grand,
        pb.grand_coefficients * u_iso, pb.ylm_coefficients * u_iso, pb.prototype_coefficients * u_iso,
        0.5 .* (ov_r .+ transpose(ov_r)), final_kin, inj_ov, 0.5 .* (inj_kin_r .+ transpose(inj_kin_r)),
    )
    all(>(0.0), _shell_local_ylm_moment_block(basis, 0)) ||
        throw(ArgumentError("coiflet profile must have positive final monopoles"))
    meta = merge(parent_profile.gauge_metadata, (coiflet_gauge = key.gauge_version, l_flat = l_flat,
        cubature_weights = Float64.(cubature_weights), min_site_alignment = min_align, target_moment_orth_err = raw_orth_err))
    return _wrap_shell_local_angular_profile(key, basis, meta)
end

"""
    coiflet_angular_profile(point_set, cubature_weights; l_flat, beta=2.0,
        l_inject=:auto, tau=1e-12, whiten=:svd, whiten_target=false)
    coiflet_angular_profile(cs_or_name; beta, l_inject=:auto, l_flat,
        tau=1e-12, whiten=:svd)

Return a cached same-size rotation of an exact-injected angular parent.
Catalog calls use stored beta/injection/flat defaults and declared degree.
Strict targets are `sqrt(weight[a])*Y_lm(n_a)` with row-Gram and weight-sum
tolerances of 1e-11. Arbitrary sets may opt into `whiten_target=true` for
rounded inputs (raw row-Gram and weight-sum errors at most 1e-5): the matched
target is then the row-orthonormal polar target, NOT the original raw moments.
`whiten` remains the separate parent-complement control. The rotation preserves
the one-body span, not general basis-dependent IDA. Treat cached arrays as read-only.
"""
function coiflet_angular_profile(point_set::SpherePointSet, cubature_weights::AbstractVector{<:Real};
    l_flat::Int, beta::Real = 2.0, l_inject::Union{Int,Symbol} = :auto,
    tau::Real = 1.0e-12, whiten::Symbol = :svd, whiten_target::Bool = false)
    resolved_l = l_inject === :auto ? _choose_l_inject(point_set.cardinality) : Int(l_inject)
    resolved_l >= l_flat >= 0 || throw(ArgumentError("coiflet_angular_profile requires l_inject >= l_flat >= 0"))
    w_tag = _coiflet_weights_digest(cubature_weights, point_set.cardinality, whiten_target)
    mode_tag = Symbol("polar_coiflet_lflat_", l_flat, "_wt", Int(whiten_target), "_w_", w_tag)
    key = ShellLocalAngularProfileKey(point_set.provenance.source_tag, _shell_local_point_set_checksum(point_set),
        point_set.order, Float64(beta), resolved_l, Float64(tau), whiten, mode_tag)
    return get!(_shell_local_angular_profile_cache(), key) do
        parent = @timeg "coiflet parent profile" shell_local_angular_profile(point_set; beta = beta, l_inject = resolved_l, tau = tau, whiten = whiten)
        _build_coiflet_profile_from_parent(parent, point_set.coordinates, cubature_weights, l_flat, key; whiten_target = whiten_target)
    end
end

function coiflet_angular_profile(cs::CubatureSpherePointSet;
    beta::Real = cs.default_beta, l_inject::Union{Int,Symbol} = :auto,
    l_flat::Int = cs.default_l_flat, tau::Real = 1.0e-12, whiten::Symbol = :svd)
    resolved_l = l_inject === :auto ? cs.default_l_inject : Int(l_inject)
    2 * l_flat <= cs.cubature_degree || throw(ArgumentError("cubature_degree=$(cs.cubature_degree) < 2*l_flat=$(2 * l_flat)"))
    return coiflet_angular_profile(cs.point_set, cs.cubature_weights;
        l_flat = l_flat, beta = beta, l_inject = resolved_l, tau = tau, whiten = whiten, whiten_target = false)
end

coiflet_angular_profile(name::Union{Symbol,AbstractString}; kwargs...) =
    coiflet_angular_profile(cubature_sphere_point_set(name); kwargs...)
