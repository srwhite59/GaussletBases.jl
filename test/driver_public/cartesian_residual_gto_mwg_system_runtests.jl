using GaussletBases
using LinearAlgebra
using SHA
using Test
using TOML

@testset "Finite collinear residual-GTO connection" begin
    GB=GaussletBases; C=GB.CartesianFinalBasisRealization; R=GB.CartesianResidualGaussians
    function candidates(z)
        orbitals = [CartesianGaussianShellOrbitalRepresentation3D("$(i)_$(p)", p,
            (0.,0.,v), [.8], [1.], :axiswise_normalized_cartesian_gaussian)
            for (i,v) in pairs(z) for p in ((0,0,0),(1,0,0),(0,1,0),(0,0,1))]
        CartesianGaussianShellSupplementRepresentation3D(:scratch_sp, orbitals, (;))
    end
     # Per-center reference is independent of streamed production accumulation.
    function reference(w,z,Z,A,e)
        t,bundles=w.terminal_basis,w.parent_axis_bundles
        loc=[(0.,0.,v) for v in z]
        result=cartesian_residual_gto_mwg_system(w,z,Z;supplement=A,expansion=e)
        base=cartesian_collinear_operators(w,z,Z;expansion=e)
        raw=C._r3a_qw_blocks(t,bundles,A,loc,e)
        products=C.pqs_terminal_residual_gto_augmented_products(t,bundles,nothing,A,
            result.residual,loc,Z;expansion=e,supplement_blocks=raw)
        Hga=raw.mixed.kinetic+sum(Z[i]*raw.mixed.nuclear[i] for i in eachindex(Z))
        Haa=raw.self.kinetic+sum(Z[i]*raw.self.nuclear[i] for i in eachindex(Z))
        (;result,base,raw,products,Hga,Haa)
    end
    primitive(o,i)=CartesianGaussianShellOrbitalRepresentation3D(o.label,o.angular_powers,
        o.center,[o.exponents[i]],[1.],o.primitive_normalization)
    aaoracle(a,b,z,Z,e)=map(k->sum(
        a.coefficients[i]*b.coefficients[j]*aaprimitive(primitive(a,i),primitive(b,j),z,Z,e)[k] for i in eachindex(a.exponents),j in eachindex(b.exponents)),(1,2))
    gaoracle(w,a,z,Z,e)=map(k->sum(
        a.coefficients[i]*gaprimitive(w,primitive(a,i),z,Z,e)[k] for i in eachindex(a.exponents)),(1,2))
     # Five-point GH integrates these polynomial Gaussian products exactly.
    GH = eigen(SymTridiagonal(zeros(5), sqrt.((1:4)./2)))
    hermite_nodes = GH.values
    hermite_weights = sqrt(pi).*abs2.(GH.vectors[1,:])
    norm1(a,l) = (2a/pi)^.25 * (l==0 ? 1. : sqrt(4a))
    function axisint(a,A,l,na,b,D,m,nb; t=0.,E=0.,kind=:overlap)
        g=a+b+t; mu=(a*A+b*D+t*E)/g
        damp=exp(-(a*(A-mu)^2+b*(D-mu)^2+t*(E-mu)^2))
        val=0.
        for (u,w) in zip(hermite_nodes,hermite_weights)
            x=mu+u/sqrt(g); y=x-D
            poly = if kind==:kinetic
                -.5*((m>=2 ? m*(m-1)*y^(m-2) : 0.) - 2b*(2m+1)*y^m + 4b*b*y^(m+2))
            else
                y^m
            end
            val += w*(x-A)^l*poly
        end
        return na*nb*damp*val/sqrt(g)
    end
    function aaprimitive(a,b,z,Z,e)
        s=ntuple(k->axisint(a.exponents[1],a.center[k],a.angular_powers[k],norm1(a.exponents[1],a.angular_powers[k]),
            b.exponents[1],b.center[k],b.angular_powers[k],norm1(b.exponents[1],b.angular_powers[k])),3)
        K=sum(axisint(a.exponents[1],a.center[k],a.angular_powers[k],norm1(a.exponents[1],a.angular_powers[k]),
            b.exponents[1],b.center[k],b.angular_powers[k],norm1(b.exponents[1],b.angular_powers[k]);kind=:kinetic)*
            prod(s[j] for j in 1:3 if j!=k) for k in 1:3)
        U=0.
        for (v,q) in zip(z,Z), (c,t) in zip(e.coefficients,e.exponents)
            U-=q*c*prod(axisint(a.exponents[1],a.center[k],a.angular_powers[k],norm1(a.exponents[1],a.angular_powers[k]),
                b.exponents[1],b.center[k],b.angular_powers[k],norm1(b.exponents[1],b.angular_powers[k]);t,E=k==3 ? v : 0.) for k in 1:3)
        end
        return prod(s),K+U
    end
    function gaprimitive(w, orbital, z, Z, e)
        proxy=C._r3a_qw_proxy_layers(w.parent_axis_bundles)
        function axisvalues(k; t=0.,E=0.,kind=:overlap)
            layer=proxy[(:x,:y,:z)[k]]
            gs=primitives(primitive_set(layer))
            v=[axisint(inv(2g.width^2),g.center_value,0,1.,orbital.exponents[1],orbital.center[k],
                orbital.angular_powers[k],norm1(orbital.exponents[1],orbital.angular_powers[k]);t,E,kind) for g in gs]
            return transpose(stencil_matrix(layer))*v
        end
        ov=ntuple(k->axisvalues(k),3)
        ki=ntuple(k->axisvalues(k;kind=:kinetic),3)
        so=zeros(w.terminal_basis.final_dimension); ho=copy(so)
        function addproduct!(out,vs,scale)
            for block in w.terminal_basis.blocks
                localv=[prod(vs[k][s[k]] for k in 1:3) for s in block.support_states]
                value=isnothing(block.coefficients) ? localv : block.coefficients'*localv
                out[block.column_range] .+= scale.*value
            end
        end
        addproduct!(so,ov,1.)
        for k in 1:3
            addproduct!(ho,ntuple(j->j==k ? ki[j] : ov[j],3),1.)
        end
        for (v,q) in zip(z,Z),(c,t) in zip(e.coefficients,e.exponents)
            addproduct!(ho,ntuple(k->axisvalues(k;t,E=k==3 ? v : 0.),3),-q*c)
        end
        return so,ho
    end
    function mixed_mwg_oracle(w,centers,widths,e,rindex)
        proxy=C._r3a_qw_proxy_layers(w.parent_axis_bundles)
        terms=ntuple(3) do k
            layer=proxy[(:x,:y,:z)[k]]; gs=primitives(primitive_set(layer)); M=stencil_matrix(layer)
            iw=M'*[sqrt(2pi)*g.width for g in gs]
            [begin
                v=[sqrt(2pi)*g.width*exp(-t*(g.center_value-centers[rindex,k])^2/
                    (1+2t*(g.width^2+widths[rindex,k]^2)))/sqrt(1+2t*(g.width^2+widths[rindex,k]^2)) for g in gs]
                (M'*v)./iw
            end for t in e.exponents]
        end
        out=zeros(w.terminal_basis.final_dimension)
        pg=ntuple(k->GB._nested_axis_pgdg(w.parent_axis_bundles,(:x,:y,:z)[k]),3)
        for b in w.terminal_basis.blocks
            values=[sum(c*prod(terms[k][j][s[k]] for k in 1:3) for (j,c) in pairs(e.coefficients)) for s in b.support_states]
            if isnothing(b.coefficients)
                out[b.column_range].=values
            else
                weights=[prod(pg[k].weights[s[k]] for k in 1:3) for s in b.support_states]
                out[b.column_range].=(b.coefficients'*(weights.*values))./(b.coefficients'*weights)
            end
        end
        out
    end
    function validate(w,data,z,Z,e; primitive_case=true)
        s=data.result; r=s.residual; n=r.base_dimension; nr=r.residual_dimension
        @test 0 < nr <= length(s.supplement.orbitals)
        @test size(s.hamiltonian.one_body,1)<600
        @test Base.summarysize(s)<512*1024^2
        X=data.raw.mixed.overlap; S_AA=data.raw.self.overlap
        T=[Matrix{Float64}(I,n,n) r.T_G; zeros(size(S_AA,1),n) r.T_A]
        pg=ntuple(k->GB._nested_axis_pgdg(w.parent_axis_bundles,(:x,:y,:z)[k]),3)
        Sgg=zeros(n,n)
        C._assemble_terminal_product_operator!(Sgg,w.terminal_basis,pg[1].overlap,
            pg[2].overlap,pg[3].overlap,C._terminal_operator_buffers(w.terminal_basis)...)
        metric=[Sgg X; X' S_AA]
        S=T'*metric*T
        ortho=norm(S-I,Inf)
        @test ortho<5e-8
        full=[data.base.one_body data.Hga; data.Hga' data.Haa]
        herr=norm(s.hamiltonian.one_body-T'*full*T,Inf)
        @test herr<1e-10
        @test s.hamiltonian.electron_electron_ida[1:n,1:n]==data.base.electron_electron_ida
        @test norm(s.hamiltonian.one_body-s.hamiltonian.one_body',Inf)<1e-10
        @test norm(s.hamiltonian.electron_electron_ida-s.hamiltonian.electron_electron_ida',Inf)<1e-10
        errS=0.; errH=0.
        for i in eachindex(s.supplement.orbitals), j in eachindex(s.supplement.orbitals)
            so,ho=aaoracle(s.supplement.orbitals[i],s.supplement.orbitals[j],z,Z,e)
            errS=max(errS,abs(S_AA[i,j]-so)); errH=max(errH,abs(data.Haa[i,j]-ho))
        end
        @test errS<1e-12
        @test errH<1e-10
        gaS=0.; gaH=0.
        for j in unique([1,min(6,size(S_AA,1)),min(7,size(S_AA,1)),size(S_AA,1)])
            so,ho=gaoracle(w,s.supplement.orbitals[j],z,Z,e)
            gaS=max(gaS,norm(so-X[:,j],Inf)); gaH=max(gaH,norm(ho-data.Hga[:,j],Inf))
        end
        @test gaS<1e-12
        @test gaH<1e-10
        centers,widths=R.moment_matched_gaussians(data.products,r)
        rr=zeros(nr,nr)
        for i in 1:nr,j in 1:nr,(c,t) in zip(e.coefficients,e.exponents)
            rr[i,j]+=c*prod(exp(-t*(centers[i,k]-centers[j,k])^2/
                (1+2t*(widths[i,k]^2+widths[j,k]^2)))/sqrt(1+2t*(widths[i,k]^2+widths[j,k]^2)) for k in 1:3)
        end
        rerr=norm(rr-s.hamiltonian.electron_electron_ida[n+1:end,n+1:end],Inf)
        @test rerr<1e-10
        gmerr=maximum(norm(mixed_mwg_oracle(w,centers,widths,e,j)-
            s.hamiltonian.electron_electron_ida[1:n,n+j],Inf) for j in (1,div(nr,2),nr))
        @test gmerr<1e-10
        if primitive_case
        probes=CartesianGaussianShellSupplementRepresentation3D(:pxpy,s.supplement.orbitals[6:7],(;))
        packet=ExternalGTOOrbitalPacket(probes,Matrix{Float64}(I,2,2),
            ExternalGTOOrbitalSpinBlock(:restricted,Matrix{Float64}(I,2,2),[1.,1.]))
        cross=gto_overlap_matrix(s,probes)
        imported=import_external_gto_orbitals(s,packet)
        @test imported.alpha.imported_coefficients==cross
        @test norm(cross-T'*metric[:,n+6:n+7],Inf)<1e-10
        @test abs(imported.alpha.orbital_captures[1]-1)<5e-8
        @test abs(imported.alpha.orbital_captures[2]-1)<5e-8
        end
    end
     # Contracted candidates use a different potential on the unchanged working basis.
    e=coulomb_gaussian_expansion(doacc=false)
    build(z,Z)=cartesian_collinear_working_basis(z,Z;core_spacing=.6,transverse_spacing=.6,
        padding_parallel=3.,padding_transverse=3.,core_side=3,angular_reference_count=3,
        outer_face_count=3,tail_spacing=2.8,angular_resolution_scale=1.4,expansion=e)
    for (z,Z,n,repulsion) in (([-1.2,0.,1.2],ones(3),231,25/12),([-2.4,0.,2.4],[1.,2.,1.],223,1.875))
        w=build(z,Z); A=candidates(z); data=reference(w,z,Z,A,e)
        @test size(data.result.hamiltonian.one_body)==(n+12,n+12)
        @test data.result.residual.owner_retained_counts==[4,4,4]
        @test abs(data.result.hamiltonian.nuclear_repulsion-repulsion)<=1e-12
        validate(w,data,z,Z,e)
    end
    z=[-1.2,0.,1.2]; Z=ones(3); w=build(z,Z)
    z2=[-1.3,.1,1.4]; Z2=[1.,2.,1.]
    A=CartesianGaussianShellSupplementRepresentation3D(:contracted,
        [CartesianGaussianShellOrbitalRepresentation3D(string(i),(0,0,0),(0.,0.,v),
            [.8,1.6],[.7,.2],:axiswise_normalized_cartesian_gaussian) for (i,v) in pairs(z2)],(;))
    data=reference(w,z2,Z2,A,e); validate(w,data,z2,Z2,e;primitive_case=false)
    packet=ExternalGTOOrbitalPacket(A,data.raw.self.overlap,ExternalGTOOrbitalSpinBlock(:restricted,inv(cholesky(Symmetric(data.raw.self.overlap)).U),ones(3)))
    @test import_external_gto_orbitals(data.result,packet).alpha.imported_coefficients==gto_overlap_matrix(data.result,A)*packet.alpha.coefficients
    explicit=cartesian_residual_gto_mwg_system(w,z2,Z2;supplement=A,expansion=e,residual_occupation_cutoff=1e-8)
    @test data.result.residual.occupation_cutoff==1e-8
    @test explicit.hamiltonian==data.result.hamiltonian
    @test explicit.residual.T_G==data.result.residual.T_G
    @test explicit.residual.T_A==data.result.residual.T_A
    for cutoff in (1e-10,1e-6)
        selected=cartesian_residual_gto_mwg_system(w,z2,Z2;supplement=A,expansion=e,residual_occupation_cutoff=cutoff)
        expected=C.pqs_terminal_residual_gto_augmentation(w.terminal_basis,w.parent_axis_bundles,A,
            [(0.,0.,v) for v in z2];residual_occupation_cutoff=cutoff)
        @test selected.residual.occupation_cutoff==cutoff
        @test selected.residual.T_G==expected.T_G
        @test selected.residual.T_A==expected.T_A
        cross=gto_overlap_matrix(selected,A)
        @test cross==gto_overlap_matrix(data.result,A)
        @test import_external_gto_orbitals(selected,packet).alpha.imported_coefficients==cross*packet.alpha.coefficients
    end
    @test_throws ArgumentError cartesian_residual_gto_mwg_system(w,z2,Z2;
        supplement=A,expansion=e,residual_occupation_cutoff=1e6)
    for cutoff in (-1.,NaN,Inf,big"1e1000")
        @test_throws ArgumentError cartesian_residual_gto_mwg_system(w,z2,Z2;
            supplement=A,expansion=e,residual_occupation_cutoff=cutoff)
    end
    A=candidates(z); call(zv,Zv,a=A,exp=e)=cartesian_residual_gto_mwg_system(w,zv,Zv;supplement=a,expansion=exp)
    for (zv,Zv) in ((reverse(z),Z),(z,[-1.,1.,1.]),(z,[1.,NaN,1.]),([0.,0.,1.],Z))
        @test_throws ArgumentError call(zv,Zv)
    end
    @test_throws TypeError call(z,Z,nothing)
    @test_throws ArgumentError call(z2,Z)
    @test_throws ArgumentError call(z,Z,A,coulomb_gaussian_expansion(doacc=true))
    orb=A.orbitals[1]
    for (exps,coefs,center,powers,normalization) in (
        (Float64[],Float64[],orb.center,(0,0,0),orb.primitive_normalization),
        ([.8],[1.,2.],orb.center,(0,0,0),orb.primitive_normalization),
        ([-1.],[1.],orb.center,(0,0,0),orb.primitive_normalization),
        ([.8],[NaN],orb.center,(0,0,0),orb.primitive_normalization),
        ([.8],[1.],(NaN,0.,0.),(0,0,0),orb.primitive_normalization),
        ([.8],[1.],orb.center,(-1,0,0),orb.primitive_normalization),
        ([.8],[1.],orb.center,(0,0,0),:unsupported),
        ([.8],[0.],orb.center,(0,0,0),orb.primitive_normalization),
    )
        bad=CartesianGaussianShellSupplementRepresentation3D(:invalid,[CartesianGaussianShellOrbitalRepresentation3D("bad",powers,center,exps,coefs,normalization)],(;))
        @test_throws ArgumentError call(z,Z,bad)
    end
    r=data.result
    for ham in ((;one_body=r.hamiltonian.one_body),merge(r.hamiltonian,(;extra=1)),merge(r.hamiltonian,(;nuclear_repulsion=1)))
        @test_throws ArgumentError GB._validate_cartesian_residual_gto_mwg_system(ham,r.terminal_basis,r.supplement,r.residual)
    end
    @test_throws DimensionMismatch GB._validate_cartesian_residual_gto_mwg_system(merge(r.hamiltonian,(;electron_electron_ida=zeros(1,1))),r.terminal_basis,r.supplement,r.residual)
    @test_throws ArgumentError GB._validate_cartesian_residual_gto_mwg_system(merge(r.hamiltonian,(;nuclear_repulsion=NaN)),r.terminal_basis,r.supplement,r.residual)
end


@testset "Public residual-GTO/MWG system" begin
    nuclei = NTuple{3,Float64}[(0.0, 0.0, -2.0), (0.0, 0.0, 2.0)]
    system = (;
        atom_symbols = ["H", "H"], nuclear_charges = [1.0, 1.0],
        atom_locations = nuclei, nup = 1, ndn = 1)
    basis = (;
        q = 5, core_spacing = 0.5,
        xmax_parallel = 6.0, xmax_transverse = 4.0)
    supplement = (;
        basis_by_center = ["cc-pVTZ", "cc-pVTZ"], lmax = 1,
        uncontracted = false, width_filtering = nothing)

    result = cartesian_residual_gto_mwg_system(
        system; basis, supplement)
    ham = result.hamiltonian
    nfinal = 487 + 18
    matrices = [ham.kinetic, ham.electron_electron_ida,
        ham.nuclear_attraction_unit_by_center...]
    @test ham isa CartesianIDAHamiltonian{Float64}
    @test size(ham.kinetic) == size(ham.electron_electron_ida) == (nfinal, nfinal)
    @test all(matrix -> all(isfinite, matrix), matrices)
    @test all(matrix -> norm(matrix - transpose(matrix), Inf) <= 1.0e-10, matrices)
    @test ham.nup == ham.ndn == 1
    @test ham.nuclear_charges == [1.0, 1.0]
    @test ham.nuclear_positions == [nucleus[axis] for nucleus in nuclei, axis in 1:3]

    h1 = one_body_hamiltonian(ham)
    orbital = eigen(Symmetric(h1)).vectors[:, 1]
    density = abs2.(orbital)
    @test dot(density, ham.electron_electron_ida * density) ≈
        0.4574161883692301 atol = 1.0e-10

    probe_source = legacy_bond_aligned_diatomic_gaussian_supplement(
        "H", "cc-pVTZ", nuclei; lmax = 1)
    source_representation = basis_representation(probe_source)
    probe_orbitals = source_representation.orbitals[[4, 5]]
    @test getproperty.(probe_orbitals, :angular_powers) == [(1, 0, 0), (0, 1, 0)]
    @test all(orbital -> orbital.center == first(nuclei), probe_orbitals)
    @test probe_orbitals[1].exponents == probe_orbitals[2].exponents
    @test all(orbital -> length(orbital.exponents) == 1 &&
        orbital.coefficients == [1.0] &&
        orbital.primitive_normalization == :axiswise_normalized_cartesian_gaussian,
        probe_orbitals)
    probes = CartesianGaussianShellSupplementRepresentation3D(
        :public_ci_probe, probe_orbitals,
        basis_metadata(source_representation))
    overlap = gto_overlap_matrix(result, probes)
    rows = [nfinal, 1, 488, 1]
    @test size(overlap) == (nfinal, 2)
    @test all(isfinite, overlap)
    @test gto_overlap_matrix(result, probes; block_indices = rows) == overlap[rows, :]
    @test_throws BoundsError gto_overlap_matrix(result, probes; block_indices = [0])

    S_GG = Matrix{Float64}(I, 2, 2)
    coefficients = Matrix{Float64}(I, 2, 2)
    source_block = ExternalGTOOrbitalSpinBlock(:restricted, coefficients, [1.0, 1.0])
    packet = ExternalGTOOrbitalPacket(probes, S_GG, source_block)
    imported = import_external_gto_orbitals(result, packet)
    @test imported.cross_overlap_size == size(overlap)
    @test imported.cross_overlap_finite
    @test imported.ordering_fingerprint_valid
    @test imported.S_GG_fingerprint_valid
    @test imported.S_GG_symmetry_error <= 1.0e-12
    @test imported.S_GG_expected_error <= 1.0e-12
    @test imported.alpha.imported_coefficients == overlap
    @test imported.alpha.source_orthogonality_error <= 1.0e-12
    @test imported.alpha.capture_matrix ≈ transpose(overlap) * overlap
    @test imported.alpha.density_trace_source == 2.0
    @test imported.alpha.density_trace_capture ≈
        dot(source_block.occupations, imported.alpha.orbital_captures)
    @test abs(imported.alpha.worst_orbital_capture - 1.0) <= 1.0e-8

    theta = 0.37
    rotation = [cos(theta) -sin(theta); sin(theta) cos(theta)]
    rotated_block = ExternalGTOOrbitalSpinBlock(
        :restricted, coefficients * rotation, [1.0, 1.0])
    rotated_packet = ExternalGTOOrbitalPacket(probes, S_GG, rotated_block)
    rotated = import_external_gto_orbitals(result, rotated_packet)
    @test rotated.alpha.imported_coefficients ≈ overlap * rotation
    @test rotated.alpha.source_orthogonality_error <= 1.0e-12
    @test rotated.alpha.density_trace_capture ≈
        imported.alpha.density_trace_capture atol = 1.0e-12 rtol = 1.0e-12

    alpha = ExternalGTOOrbitalSpinBlock(:alpha, coefficients[:, 1:1], [1.0])
    beta = ExternalGTOOrbitalSpinBlock(:beta, coefficients[:, 2:2], [1.0])
    spin_packet = ExternalGTOOrbitalPacket(probes, S_GG, alpha; beta)
    spin_import = import_external_gto_orbitals(result, spin_packet)
    @test spin_import.alpha.spin == :alpha
    @test spin_import.beta !== nothing
    @test spin_import.beta.spin == :beta
    @test spin_import.alpha.imported_coefficients == overlap[:, 1:1]
    @test spin_import.beta.imported_coefficients == overlap[:, 2:2]

    fixture_manifest = joinpath(@__DIR__, "external_cartesian_gto_h2_ccpvtz_v1.toml")
    fixture_payload = splitext(fixture_manifest)[1] * ".f64"
    fixture_packet = read_external_cartesian_gto_packet(fixture_manifest)
    @test fixture_packet.provenance.external_cartesian_gto.producer.exporter_version == "1.0.0"
    @test length(fixture_packet.ao_labels) == 30
    @test any(sum(powers) == 2 for powers in fixture_packet.angular_powers)
    @test fixture_packet.alpha.spin == :restricted
    @test fixture_packet.alpha.occupations == [2.0]
    @test fixture_packet.beta === nothing
    @test norm(transpose(fixture_packet.alpha.coefficients) * fixture_packet.S_GG *
        fixture_packet.alpha.coefficients - I, Inf) <= 1.0e-12

    fixture_overlap = gto_overlap_matrix(result, fixture_packet.probes)
    fixture_expected = fixture_overlap * fixture_packet.alpha.coefficients
    fixture_import = import_external_gto_orbitals(result, fixture_packet)
    @test fixture_import.alpha.imported_coefficients == fixture_expected
    @test fixture_import.ordering_fingerprint_valid
    @test fixture_import.S_GG_fingerprint_valid
    closest = closest_external_gto_determinant(
        result, fixture_packet; minimum_gram_eigenvalue = 0.99)
    @test closest.imported.alpha.imported_coefficients == fixture_expected
    @test closest.alpha.minimum_gram_eigenvalue >= 0.99
    @test closest.alpha.gram_eigenvalues ≈
        eigvals(Symmetric(transpose(fixture_expected) * fixture_expected))
    @test closest.alpha.principal_angles_radians ≈
        acos.(sqrt.(clamp.(closest.alpha.gram_eigenvalues, 0.0, 1.0)))
    @test closest.alpha.maximum_principal_angle_radians ==
        maximum(closest.alpha.principal_angles_radians)
    @test closest.alpha.orthonormality_error <= 1.0e-10
    @test norm(transpose(closest.alpha.coefficients) * closest.alpha.coefficients - I, Inf) <= 1.0e-10
    @test_throws ArgumentError closest_external_gto_determinant(
        result, fixture_packet; minimum_gram_eigenvalue = 0.0)
    @test_throws ArgumentError closest_external_gto_determinant(
        result, fixture_packet; minimum_gram_eigenvalue = 1.0)

    fixture_identity = fixture_packet.provenance.external_cartesian_gto
    packet_with = function (identity, alpha = fixture_packet.alpha)
        return ExternalGTOOrbitalPacket(fixture_packet.probes, fixture_packet.S_GG, alpha;
            ao_labels = fixture_packet.ao_labels,
            provenance = (; external_cartesian_gto = identity))
    end
    wrong_positions = copy(fixture_identity.nuclear_positions)
    wrong_positions[1, 3] += 0.01
    @test_throws ArgumentError closest_external_gto_determinant(
        result, packet_with(merge(fixture_identity, (; nuclear_positions = wrong_positions)));
        minimum_gram_eigenvalue = 0.99)
    @test_throws ArgumentError closest_external_gto_determinant(
        result, packet_with(merge(fixture_identity, (; nalpha = 2)));
        minimum_gram_eigenvalue = 0.99)
    fractional = ExternalGTOOrbitalSpinBlock(
        :restricted, fixture_packet.alpha.coefficients, [1.5])
    @test_throws ArgumentError closest_external_gto_determinant(
        result, packet_with(fixture_identity, fractional); minimum_gram_eigenvalue = 0.99)

    with_fixture = function (edit)
        mktempdir() do directory
            manifest = joinpath(directory, basename(fixture_manifest))
            payload = joinpath(directory, basename(fixture_payload))
            cp(fixture_manifest, manifest)
            cp(fixture_payload, payload)
            edit(manifest, payload)
            @test_throws ArgumentError read_external_cartesian_gto_packet(manifest)
        end
    end
    replace_manifest = function (old, new)
        with_fixture() do manifest, _
            contents = read(manifest, String)
            @test occursin(old, contents)
            write(manifest, replace(contents, old => new; count = 1))
        end
    end
    replace_manifest("exponent = \"bohr_inverse_square\"", "exponent = \"angstrom_inverse_square\"")
    replace_manifest("shape = [30, 30]", "shape = [29, 30]")
    replace_manifest("de579e9188f46d7e60cb6c7a0ad86263ad66edc22d5a905a64a4129d2669b0a1", repeat("0", 64))
    replace_manifest("ao_index_1based = 1", "ao_index_1based = 2")
    replace_manifest("0.025494863234680084", "0.12549486323468008")
    with_fixture() do manifest, payload
        bytes = read(payload)
        values = collect(reinterpret(Float64, bytes))
        values[2] += 1.0e-4
        values[31] += 1.0e-4
        changed = collect(reinterpret(UInt8, values))
        write(payload, changed)
        document = TOML.parsefile(manifest)
        document["payload"]["sha256"] = bytes2hex(sha256(changed))
        document["arrays"][1]["sha256"] = bytes2hex(sha256(changed[1:(30 * 30 * 8)]))
        open(manifest, "w") do io
            TOML.print(io, document; sorted = true)
        end
    end

    stale_order = ExternalGTOOrbitalPacket(
        probes, S_GG, source_block; ordering_fingerprint = "stale")
    @test_throws ArgumentError import_external_gto_orbitals(result, stale_order)
    stale_overlap = ExternalGTOOrbitalPacket(
        probes, S_GG, source_block; S_GG_fingerprint = "stale")
    @test_throws ArgumentError import_external_gto_orbitals(result, stale_overlap)

    wrong_metric = copy(S_GG)
    wrong_metric[1, 1] += 0.05
    wrong_metric_packet = ExternalGTOOrbitalPacket(
        probes, wrong_metric, source_block;
        S_GG_fingerprint = external_gto_overlap_fingerprint(wrong_metric))
    @test_throws ArgumentError import_external_gto_orbitals(result, wrong_metric_packet)
    nonsymmetric_metric = copy(S_GG)
    nonsymmetric_metric[1, 2] += 1.0e-3
    nonsymmetric_packet = ExternalGTOOrbitalPacket(
        probes, nonsymmetric_metric, source_block;
        S_GG_fingerprint = external_gto_overlap_fingerprint(nonsymmetric_metric))
    @test_throws ArgumentError import_external_gto_orbitals(result, nonsymmetric_packet)

    beta_only = ExternalGTOOrbitalSpinBlock(:beta, coefficients[:, 2:2], [1.0])
    @test_throws ArgumentError ExternalGTOOrbitalPacket(probes, S_GG, beta_only)
    @test_throws ArgumentError ExternalGTOOrbitalPacket(
        probes, S_GG, source_block; beta)
    wrong_beta = ExternalGTOOrbitalSpinBlock(:alpha, coefficients[:, 2:2], [1.0])
    @test_throws ArgumentError ExternalGTOOrbitalPacket(probes, S_GG, alpha; beta = wrong_beta)

    nonorthogonal = ExternalGTOOrbitalSpinBlock(
        :restricted, 2.0 .* coefficients, [1.0, 1.0])
    nonorthogonal_packet = ExternalGTOOrbitalPacket(probes, S_GG, nonorthogonal)
    @test_throws ArgumentError import_external_gto_orbitals(result, nonorthogonal_packet)

    @test_throws ArgumentError cartesian_residual_gto_mwg_system(
        merge(system, (; extra = true)); basis, supplement)
    @test_throws ArgumentError cartesian_residual_gto_mwg_system(
        system; basis = merge(basis, (; radius = 4.0)), supplement)
    @test_throws ArgumentError cartesian_residual_gto_mwg_system(
        system; basis, supplement = merge(supplement, (; extra = true)))
end
