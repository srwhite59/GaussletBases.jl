# Compare the self-contained candidate with the still-authoritative prose during transition.
module CartesianAuthorityTransition

using SHA
using TOML

include(joinpath(@__DIR__, "check_cartesian_authority_shadow.jl"))
using .CartesianAuthorityShadow
include(joinpath(@__DIR__, "check_cartesian_authority_candidate.jl"))
using .CartesianAuthorityCandidate

const ROOT = CartesianAuthorityCandidate.ROOT
const TRANSITION_PATH = joinpath(
    CartesianAuthorityCandidate.PRODUCER_DOCS,
    "authority_transition_snapshot.toml",
)
const ARTIFACT_KIND = "cartesian_authority_transition_snapshot"
const RECONCILIATION_SOURCE =
    "docs/src/developer/designs/cartesian_hamiltonian_producer/reviews/authority_transition_rehearsal_pass394_2026-07-12.md"
const TOP_LEVEL_KEYS = Set([
    "schema_version",
    "artifact_kind",
    "authoritative",
    "authorization_complete",
    "candidate_sha256",
    "registry_sha256",
    "registry_file_sha256",
    "agents_file_sha256",
    "agents_whitelist_sha256",
    "agents_whitelist_block_sha256",
    "record_count",
    "execution_whitelist_count",
    "semantic_parity",
    "semantic_review_status",
    "reconciliation_source",
])

digest(text::AbstractString) = bytes2hex(sha256(codeunits(String(text))))

function serialized(data)
    io = IOBuffer()
    TOML.print(io, data; sorted = true)
    text = String(take!(io))
    return endswith(text, '\n') ? text : text * "\n"
end

function _legacy_document_paths(record)
    paths = String[]
    for link in record["owner_document_links"]
        relative = first(split(link, '#'; limit = 2))
        absolute = normpath(joinpath(dirname(CartesianAuthorityShadow.REGISTRY_PATH), relative))
        push!(paths, relpath(absolute, ROOT))
    end
    return sort!(unique!(paths))
end

function parity_errors(candidate_data, legacy_records, legacy_whitelist)
    errors = String[]
    candidate_records = candidate_data["records"]
    candidate_by_id = Dict(record["id"] => record for record in candidate_records)
    legacy_by_id = Dict(record["id"] => record for record in legacy_records)
    candidate_ids = sort!(collect(keys(candidate_by_id)))
    legacy_ids = sort!(collect(keys(legacy_by_id)))
    candidate_ids == legacy_ids || push!(errors, "candidate/registry ID coverage mismatch")

    for id in intersect(Set(candidate_ids), Set(legacy_ids))
        candidate = candidate_by_id[id]
        legacy = legacy_by_id[id]
        candidate["title"] == legacy["heading_title"] || push!(errors, "$id title mismatch")
        candidate_paths = sort!(unique!(String[ref["path"] for ref in candidate["document_refs"]]))
        candidate_paths == _legacy_document_paths(legacy) ||
            push!(errors, "$id document-owner path mismatch")
    end

    derived = CartesianAuthorityCandidate.execution_ids(candidate_data)
    derived == sort!(String.(legacy_whitelist)) ||
        push!(errors, "candidate-derived execution whitelist mismatch")
    return sort!(unique!(errors))
end

function expected_snapshot(candidate, registry_text, agents_text)
    CartesianAuthorityCandidate.check_candidate(candidate)
    legacy_records, registry_digest = CartesianAuthorityShadow._registry_records(registry_text)
    whitelist, whitelist_digest, whitelist_block_digest =
        CartesianAuthorityShadow._parse_agents_whitelist(agents_text)
    errors = parity_errors(candidate.data, legacy_records, whitelist)
    candidate_block_digest = CartesianAuthorityCandidate.digest(
        CartesianAuthorityCandidate.render_whitelist_block(candidate),
    )
    candidate_block_digest == whitelist_block_digest ||
        push!(errors, "candidate-derived whitelist block differs from marked AGENTS block")
    isempty(errors) || error("Cartesian authority transition parity errors:\n" * join(errors, "\n"))
    return Dict{String,Any}(
        "schema_version" => 2,
        "artifact_kind" => ARTIFACT_KIND,
        "authoritative" => false,
        "authorization_complete" => false,
        "candidate_sha256" => candidate.sha256,
        "registry_sha256" => registry_digest,
        "registry_file_sha256" => digest(registry_text),
        "agents_file_sha256" => digest(agents_text),
        "agents_whitelist_sha256" => whitelist_digest,
        "agents_whitelist_block_sha256" => whitelist_block_digest,
        "record_count" => length(candidate.data["records"]),
        "execution_whitelist_count" => length(whitelist),
        "semantic_parity" => "manual_review_complete",
        "semantic_review_status" => "independent_rehearsal_and_focused_reconciliation_passed",
        "reconciliation_source" => RECONCILIATION_SOURCE,
    )
end

expected_snapshot() = expected_snapshot(
    CartesianAuthorityCandidate.load_snapshot(),
    read(CartesianAuthorityShadow.REGISTRY_PATH, String),
    read(CartesianAuthorityShadow.AGENTS_PATH, String),
)

function _git_head()
    head = strip(read(`git -C $ROOT rev-parse HEAD`, String))
    occursin(r"^[0-9a-f]{40}$", head) || error("invalid Git HEAD: $head")
    return head
end

function _transition_bindings(
    candidate,
    transition,
    transition_text;
    git_head,
    checker_sha256,
    shadow_checker_sha256,
)
    transition["candidate_sha256"] == candidate.sha256 ||
        error("transition candidate SHA-256 does not match captured candidate")
    return Dict{String,Any}(
        "git_head" => git_head,
        "transition_snapshot_sha256" => digest(transition_text),
        "transition_checker_sha256" => checker_sha256,
        "shadow_checker_sha256" => shadow_checker_sha256,
        "registry_sha256" => transition["registry_sha256"],
        "registry_file_sha256" => transition["registry_file_sha256"],
        "agents_file_sha256" => transition["agents_file_sha256"],
        "agents_whitelist_sha256" => transition["agents_whitelist_sha256"],
        "agents_whitelist_block_sha256" => transition["agents_whitelist_block_sha256"],
        "semantic_parity" => transition["semantic_parity"],
        "semantic_review_status" => transition["semantic_review_status"],
        "reconciliation_source" => transition["reconciliation_source"],
        "candidate_sha256" => candidate.sha256,
    )
end

function _validated_transition(candidate, transition_text, registry_text, agents_text)
    actual = TOML.parse(transition_text)
    Set(String.(keys(actual))) == TOP_LEVEL_KEYS || error("transition snapshot has unexpected keys")
    actual["candidate_sha256"] == candidate.sha256 ||
        error("transition candidate SHA-256 does not match captured candidate")
    expected = expected_snapshot(candidate, registry_text, agents_text)
    actual == expected || error("Cartesian authority transition snapshot drift")
    transition_text == serialized(actual) || error("transition snapshot serialization is not canonical")
    return actual
end

function write_transition_rehearsal(output_dir)
    candidate = CartesianAuthorityCandidate.load_snapshot()
    transition_text = read(TRANSITION_PATH, String)
    registry_text = read(CartesianAuthorityShadow.REGISTRY_PATH, String)
    agents_text = read(CartesianAuthorityShadow.AGENTS_PATH, String)
    checker_text = read(@__FILE__, String)
    shadow_checker_path = joinpath(@__DIR__, "check_cartesian_authority_shadow.jl")
    shadow_checker_text = read(shadow_checker_path, String)
    git_head = _git_head()
    transition = _validated_transition(candidate, transition_text, registry_text, agents_text)
    bindings = _transition_bindings(
        candidate,
        transition,
        transition_text;
        git_head = git_head,
        checker_sha256 = digest(checker_text),
        shadow_checker_sha256 = digest(shadow_checker_text),
    )

    CartesianAuthorityCandidate.write_rehearsal(
        output_dir,
        candidate;
        transition_bindings = bindings,
    )

    read(candidate.path, String) == candidate.text || error("candidate changed during rehearsal")
    read(TRANSITION_PATH, String) == transition_text || error("transition snapshot changed during rehearsal")
    read(CartesianAuthorityShadow.REGISTRY_PATH, String) == registry_text ||
        error("registry changed during rehearsal")
    read(CartesianAuthorityShadow.AGENTS_PATH, String) == agents_text ||
        error("AGENTS changed during rehearsal")
    read(@__FILE__, String) == checker_text || error("transition checker changed during rehearsal")
    read(shadow_checker_path, String) == shadow_checker_text ||
        error("shadow checker changed during rehearsal")
    _git_head() == git_head || error("Git HEAD changed during rehearsal")
    manifest = TOML.parsefile(joinpath(abspath(output_dir), "manifest.toml"))
    manifest["transition_bindings"] == bindings || error("rehearsal manifest binding mismatch")
    return abspath(output_dir)
end

function check_transition()
    isfile(TRANSITION_PATH) || error("missing transition snapshot: $TRANSITION_PATH")
    candidate = CartesianAuthorityCandidate.load_snapshot()
    transition_text = read(TRANSITION_PATH, String)
    registry_text = read(CartesianAuthorityShadow.REGISTRY_PATH, String)
    agents_text = read(CartesianAuthorityShadow.AGENTS_PATH, String)
    _validated_transition(candidate, transition_text, registry_text, agents_text)
    read(candidate.path, String) == candidate.text || error("candidate changed during transition check")
    read(TRANSITION_PATH, String) == transition_text || error("transition snapshot changed during transition check")
    read(CartesianAuthorityShadow.REGISTRY_PATH, String) == registry_text ||
        error("registry changed during transition check")
    read(CartesianAuthorityShadow.AGENTS_PATH, String) == agents_text ||
        error("AGENTS changed during transition check")
    return nothing
end

function write_transition()
    content = serialized(expected_snapshot())
    temporary, io = mktemp(dirname(TRANSITION_PATH))
    try
        write(io, content)
        flush(io)
        close(io)
        mv(temporary, TRANSITION_PATH; force = true)
    catch
        isopen(io) && close(io)
        ispath(temporary) && rm(temporary; force = true)
        rethrow()
    end
    return TRANSITION_PATH
end

function _expect_failure(candidate, records, whitelist, needle)
    errors = parity_errors(candidate, records, whitelist)
    any(error -> occursin(needle, error), errors) ||
        error("negative transition check did not produce $(repr(needle)); got $(join(errors, "; "))")
end

function self_test()
    candidate_snapshot = CartesianAuthorityCandidate.load_snapshot()
    candidate = candidate_snapshot.data
    registry_text = read(CartesianAuthorityShadow.REGISTRY_PATH, String)
    agents_text = read(CartesianAuthorityShadow.AGENTS_PATH, String)
    records, _ = CartesianAuthorityShadow._registry_records(registry_text)
    whitelist, _, _ = CartesianAuthorityShadow._parse_agents_whitelist(agents_text)
    isempty(parity_errors(candidate, records, whitelist)) || error("transition parity must pass before self-test")

    mismatched = expected_snapshot(candidate_snapshot, registry_text, agents_text)
    mismatched["candidate_sha256"] = repeat("0", 64)
    try
        _validated_transition(candidate_snapshot, serialized(mismatched), registry_text, agents_text)
        error("candidate/transition mismatch unexpectedly passed")
    catch error_value
        occursin("candidate SHA-256", sprint(showerror, error_value)) || rethrow()
    end

    broken = deepcopy(candidate)
    broken["records"][1]["title"] *= " changed"
    _expect_failure(broken, records, whitelist, "title mismatch")

    broken = deepcopy(candidate)
    original_path = broken["records"][1]["document_refs"][1]["path"]
    replacement = first(
        record for record in broken["records"]
        if record["document_refs"][1]["path"] != original_path
    )
    broken["records"][1]["document_refs"][1]["path"] =
        replacement["document_refs"][1]["path"]
    _expect_failure(broken, records, whitelist, "document-owner path mismatch")

    shortened = whitelist[2:end]
    _expect_failure(candidate, records, shortened, "execution whitelist mismatch")
    return nothing
end

function main(args = ARGS)
    if isempty(args) || args == ["--check"]
        check_transition()
    elseif args == ["--write"]
        println(write_transition())
    elseif args == ["--self-test"]
        self_test()
    elseif length(args) == 2 && args[1] == "--write-rehearsal"
        println(write_transition_rehearsal(args[2]))
    else
        error("usage: julia --project=docs docs/check_cartesian_authority_transition.jl [--check|--write|--self-test|--write-rehearsal DIR]")
    end
end

end

if abspath(PROGRAM_FILE) == @__FILE__
    CartesianAuthorityTransition.main()
end
