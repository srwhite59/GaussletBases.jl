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
const TOP_LEVEL_KEYS = Set([
    "schema_version",
    "artifact_kind",
    "authoritative",
    "authorization_complete",
    "candidate_sha256",
    "registry_sha256",
    "agents_whitelist_sha256",
    "agents_whitelist_block_sha256",
    "record_count",
    "execution_whitelist_count",
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

function expected_snapshot()
    candidate = CartesianAuthorityCandidate.load_snapshot()
    CartesianAuthorityCandidate.check_candidate(candidate)
    legacy_records, registry_digest = CartesianAuthorityShadow._registry_records()
    whitelist, whitelist_digest, whitelist_block_digest =
        CartesianAuthorityShadow._agents_whitelist()
    errors = parity_errors(candidate.data, legacy_records, whitelist)
    candidate_block_digest = CartesianAuthorityCandidate.digest(
        CartesianAuthorityCandidate.render_whitelist_block(candidate),
    )
    candidate_block_digest == whitelist_block_digest ||
        push!(errors, "candidate-derived whitelist block differs from marked AGENTS block")
    isempty(errors) || error("Cartesian authority transition parity errors:\n" * join(errors, "\n"))
    return Dict{String,Any}(
        "schema_version" => 1,
        "artifact_kind" => ARTIFACT_KIND,
        "authoritative" => false,
        "authorization_complete" => false,
        "candidate_sha256" => candidate.sha256,
        "registry_sha256" => registry_digest,
        "agents_whitelist_sha256" => whitelist_digest,
        "agents_whitelist_block_sha256" => whitelist_block_digest,
        "record_count" => length(candidate.data["records"]),
        "execution_whitelist_count" => length(whitelist),
    )
end

function check_transition()
    isfile(TRANSITION_PATH) || error("missing transition snapshot: $TRANSITION_PATH")
    expected = expected_snapshot()
    actual_text = read(TRANSITION_PATH, String)
    actual = TOML.parse(actual_text)
    Set(String.(keys(actual))) == TOP_LEVEL_KEYS || error("transition snapshot has unexpected keys")
    actual == expected || error("Cartesian authority transition snapshot drift")
    actual_text == serialized(actual) || error("transition snapshot serialization is not canonical")
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
    candidate = CartesianAuthorityCandidate.load_snapshot().data
    records, _ = CartesianAuthorityShadow._registry_records()
    whitelist, _, _ = CartesianAuthorityShadow._agents_whitelist()
    isempty(parity_errors(candidate, records, whitelist)) || error("transition parity must pass before self-test")

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
    else
        error("usage: julia --project=docs docs/check_cartesian_authority_transition.jl [--check|--write|--self-test]")
    end
end

end

if abspath(PROGRAM_FILE) == @__FILE__
    CartesianAuthorityTransition.main()
end
