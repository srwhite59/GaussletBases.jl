# Validate explicit semantic metadata without making it execution authority.
module CartesianAuthorityCandidate

using SHA
using TOML

include(joinpath(@__DIR__, "check_cartesian_authority_shadow.jl"))
using .CartesianAuthorityShadow

const ROOT = CartesianAuthorityShadow.ROOT
const ROOT_REAL = realpath(ROOT)
const CANDIDATE_PATH = joinpath(
    CartesianAuthorityShadow.PRODUCER_DOCS,
    "authority_candidate.toml",
)
const ARTIFACT_KIND = "cartesian_authority_candidate"
const SCOPE = "explicit_non_authoritative_candidate_metadata_only"
const DEPENDENCY_SEMANTICS = "non_topological_references_cycles_allowed"
const REVIEWED_SEMANTICS_SHA256 = "506c5dcfcc15db9726969702c87d9c26835e4c40ae9f4f643c1a2ebd28629418"
const LIFECYCLES = Set([
    "proposed",
    "approved",
    "implemented",
    "completed",
    "deferred",
    "suspended",
    "superseded",
    "retired",
    "rejected",
])
const GRANTS = Set([
    "none",
    "design",
    "measurement",
    "implementation",
    "maintenance",
    "preservation",
    "retirement",
])
const SURFACES = Set([
    "source",
    "tests",
    "docs",
    "tools",
    "artifacts",
    "driver",
])
const EXECUTION_GRANTS = Set(["implementation", "maintenance", "preservation", "retirement"])
const EXECUTION_SURFACES = Set(["source", "tests", "tools", "artifacts", "driver"])
const CLOSED_LIFECYCLES = Set(["superseded", "retired", "rejected"])
const TOP_LEVEL_KEYS = Set([
    "schema_version",
    "artifact_kind",
    "authoritative",
    "authorization_complete",
    "generated",
    "semantic_review_complete",
    "scope",
    "dependency_semantics",
    "registry_source",
    "agents_source",
    "registry_sha256",
    "agents_whitelist_sha256",
    "registry_record_count",
    "agents_whitelist_count",
    "canonical_documents",
    "records",
])
const RECORD_KEYS = Set([
    "id",
    "heading_title",
    "lifecycle",
    "grant",
    "surfaces",
    "agents_whitelisted",
    "canonical",
    "source_paths",
    "test_paths",
    "planned_test_paths",
    "dependency_ids",
    "scope",
    "registry_record_sha256",
])
const CANONICAL_KEYS = Set(["path", "heading"])
const CANONICAL_DOCUMENT_KEYS = Set(["path", "sha256"])
const SEMANTIC_RECORD_FIELDS = (
    "id",
    "lifecycle",
    "grant",
    "surfaces",
    "agents_whitelisted",
    "canonical",
    "source_paths",
    "test_paths",
    "planned_test_paths",
    "dependency_ids",
    "scope",
)

function _serialized(data)
    io = IOBuffer()
    TOML.print(io, data; sorted = true)
    text = String(take!(io))
    return endswith(text, '\n') ? text : text * "\n"
end

function _sorted_unique_strings(values)
    values isa AbstractVector || return false
    all(value -> value isa AbstractString, values) || return false
    return values == sort!(unique!(String.(values)))
end

function _safe_repo_path(
    path;
    allow_missing = false,
    root = ROOT,
    root_real = ROOT_REAL,
)
    path isa AbstractString || return nothing
    isabspath(path) && return nothing
    normalized = normpath(String(path))
    normalized == path || return nothing
    (normalized == ".." || startswith(normalized, "../")) && return nothing
    absolute = joinpath(root, normalized)
    islink(absolute) && return nothing
    probe = if ispath(absolute)
        realpath(absolute)
    elseif allow_missing && isdir(dirname(absolute))
        realpath(dirname(absolute))
    else
        return nothing
    end
    relative = relpath(probe, root_real)
    (!isabspath(relative) && relative != ".." && !startswith(relative, "../")) ||
        return nothing
    return absolute
end

function _semantic_digest(data)
    semantic_records = Dict{String,Any}[]
    for record in get(data, "records", Any[])
        record isa AbstractDict || continue
        push!(semantic_records, Dict(
            field => deepcopy(get(record, field, nothing))
            for field in SEMANTIC_RECORD_FIELDS
        ))
    end
    canonical_paths = String[]
    documents = get(data, "canonical_documents", Any[])
    if documents isa AbstractVector
        for document in documents
            document isa AbstractDict || continue
            path = get(document, "path", nothing)
            path isa AbstractString && push!(canonical_paths, String(path))
        end
    end
    projection = Dict{String,Any}(
        "canonical_document_paths" => sort!(unique!(canonical_paths)),
        "dependency_semantics" => get(data, "dependency_semantics", nothing),
        "records" => semantic_records,
    )
    return bytes2hex(sha256(codeunits(_serialized(projection))))
end

function _markdown_headings(path)
    lines = CartesianAuthorityShadow._normalized_lines(read(path, String))
    headings = String[]
    for line in CartesianAuthorityShadow._visible_markdown_lines(lines)
        isnothing(line) && continue
        matched = match(r"^ {0,3}#{1,6}\s+(.+?)\s*$", line)
        isnothing(matched) || push!(headings, strip(matched.captures[1]))
    end
    return headings
end

function _expected_canonical(record)
    sections = Dict{String,Any}[]
    for link in record["owner_document_links"]
        relative = first(split(link, '#'; limit = 2))
        absolute = normpath(joinpath(dirname(CartesianAuthorityShadow.REGISTRY_PATH), relative))
        push!(sections, Dict{String,Any}(
            "path" => relpath(absolute, ROOT),
            "heading" => first(_markdown_headings(absolute)),
        ))
    end
    unique_sections = Dict{String,Any}[]
    seen = Set{Tuple{String,String}}()
    for section in sort!(sections; by = item -> (item["path"], item["heading"]))
        key = (section["path"], section["heading"])
        key in seen && continue
        push!(seen, key)
        push!(unique_sections, section)
    end
    return unique_sections
end

function _canonical_document_errors(data)
    errors = String[]
    documents = get(data, "canonical_documents", nothing)
    documents isa AbstractVector || return ["canonical_documents must be a list"], Set{String}()
    paths = String[]
    for (index, document) in enumerate(documents)
        document isa AbstractDict || begin
            push!(errors, "canonical_documents[$index] must be a table")
            continue
        end
        Set(String.(keys(document))) == CANONICAL_DOCUMENT_KEYS ||
            push!(errors, "canonical_documents[$index] has unexpected keys")
        path = get(document, "path", nothing)
        digest = get(document, "sha256", nothing)
        path isa AbstractString && digest isa AbstractString || begin
            push!(errors, "canonical_documents[$index] path/sha256 must be strings")
            continue
        end
        push!(paths, String(path))
        absolute = _safe_repo_path(path)
        isnothing(absolute) && begin
            push!(errors, "canonical_documents[$index] path is invalid: $path")
            continue
        end
        (!startswith(path, "docs/") || !endswith(path, ".md")) &&
            push!(errors, "canonical_documents[$index] must name docs Markdown: $path")
        occursin(r"^[0-9a-f]{64}$", digest) ||
            push!(errors, "canonical_documents[$index] has invalid sha256")
        bytes2hex(sha256(read(absolute))) == digest ||
            push!(errors, "canonical document content drift: $path")
    end
    paths == sort!(unique!(paths)) ||
        push!(errors, "canonical_documents must be sorted and unique")
    return errors, Set(paths)
end

function _path_errors(record, id)
    errors = String[]
    for (field, prefix, must_exist) in (
        ("source_paths", r"^(?:src|bin|tools)/", true),
        ("test_paths", "test/", true),
        ("planned_test_paths", "test/", false),
    )
        values = get(record, field, nothing)
        _sorted_unique_strings(values) || begin
            push!(errors, "$id.$field must be a sorted unique string list")
            continue
        end
        for value in values
            absolute = _safe_repo_path(value; allow_missing = !must_exist)
            isnothing(absolute) && begin
                push!(errors, "$id.$field is noncanonical or escapes the repository: $value")
                continue
            end
            if prefix isa Regex
                isnothing(match(prefix, value)) &&
                    push!(errors, "$id.$field has invalid prefix: $value")
            elseif !isnothing(prefix) && !startswith(value, prefix)
                push!(errors, "$id.$field has invalid prefix: $value")
            end
            if must_exist
                isfile(absolute) || push!(errors, "$id.$field path does not exist: $value")
            else
                ispath(absolute) && push!(errors, "$id.$field is already present: $value")
            end
        end
    end
    return errors
end

function _canonical_errors(record, id, document_paths)
    errors = String[]
    canonical = get(record, "canonical", nothing)
    canonical isa AbstractVector || return ["$id.canonical must be a list"]
    isempty(canonical) && push!(errors, "$id.canonical must not be empty")
    keys_seen = Tuple{String,String}[]
    for (index, section) in enumerate(canonical)
        section isa AbstractDict || begin
            push!(errors, "$id.canonical[$index] must be a table")
            continue
        end
        Set(String.(keys(section))) == CANONICAL_KEYS ||
            push!(errors, "$id.canonical[$index] has unexpected keys")
        path = get(section, "path", nothing)
        heading = get(section, "heading", nothing)
        path isa AbstractString && heading isa AbstractString || begin
            push!(errors, "$id.canonical[$index] path/heading must be strings")
            continue
        end
        (!startswith(path, "docs/") || !endswith(path, ".md")) &&
            push!(errors, "$id.canonical[$index] must name a docs Markdown file: $path")
        path in document_paths ||
            push!(errors, "$id.canonical[$index] is absent from canonical_documents: $path")
        absolute = _safe_repo_path(path)
        isnothing(absolute) && begin
            push!(errors, "$id.canonical[$index] escapes the repository: $path")
            continue
        end
        isfile(absolute) || begin
            push!(errors, "$id.canonical[$index] path does not exist: $path")
            continue
        end
        heading in _markdown_headings(absolute) ||
            push!(errors, "$id.canonical[$index] heading does not exist: $heading")
        push!(keys_seen, (String(path), String(heading)))
    end
    keys_seen == sort!(unique!(keys_seen)) ||
        push!(errors, "$id.canonical must be sorted and unique")
    return errors
end

function validation_errors(data)
    errors = String[]
    expected_records, registry_hash = CartesianAuthorityShadow._registry_records()
    whitelist, whitelist_hash = CartesianAuthorityShadow._agents_whitelist()
    expected_by_id = Dict(record["id"] => record for record in expected_records)
    whitelist_set = Set(whitelist)
    canonical_errors, canonical_document_paths = _canonical_document_errors(data)
    append!(errors, canonical_errors)

    Set(String.(keys(data))) == TOP_LEVEL_KEYS ||
        push!(errors, "candidate has unexpected or missing top-level keys")
    schema_version = get(data, "schema_version", nothing)
    schema_version isa Integer && schema_version == 2 ||
        push!(errors, "schema_version must be integer 2")
    get(data, "artifact_kind", nothing) == ARTIFACT_KIND ||
        push!(errors, "artifact_kind mismatch")
    get(data, "authoritative", nothing) === false ||
        push!(errors, "candidate must remain authoritative=false")
    get(data, "authorization_complete", nothing) === false ||
        push!(errors, "candidate must remain authorization_complete=false")
    get(data, "generated", nothing) === false ||
        push!(errors, "candidate must remain generated=false")
    get(data, "semantic_review_complete", nothing) === true ||
        push!(errors, "semantic_review_complete must be true after review")
    get(data, "scope", nothing) == SCOPE || push!(errors, "scope mismatch")
    get(data, "dependency_semantics", nothing) == DEPENDENCY_SEMANTICS ||
        push!(errors, "dependency_semantics mismatch")
    get(data, "registry_source", nothing) == CartesianAuthorityShadow.REGISTRY_SOURCE ||
        push!(errors, "registry_source mismatch")
    get(data, "agents_source", nothing) == CartesianAuthorityShadow.AGENTS_SOURCE ||
        push!(errors, "agents_source mismatch")
    get(data, "registry_sha256", nothing) == registry_hash ||
        push!(errors, "registry_sha256 drift")
    get(data, "agents_whitelist_sha256", nothing) == whitelist_hash ||
        push!(errors, "agents_whitelist_sha256 drift")
    registry_record_count = get(data, "registry_record_count", nothing)
    registry_record_count isa Integer && registry_record_count == length(expected_records) ||
        push!(errors, "registry_record_count mismatch")
    agents_whitelist_count = get(data, "agents_whitelist_count", nothing)
    agents_whitelist_count isa Integer && agents_whitelist_count == length(whitelist) ||
        push!(errors, "agents_whitelist_count mismatch")
    _semantic_digest(data) == REVIEWED_SEMANTICS_SHA256 ||
        push!(errors, "reviewed semantic metadata drift")

    records = get(data, "records", nothing)
    records isa AbstractVector || return vcat(errors, ["records must be a list"])
    ids = String[]
    referenced_canonical_paths = Set{String}()
    for (index, record) in enumerate(records)
        record isa AbstractDict || begin
            push!(errors, "records[$index] must be a table")
            continue
        end
        Set(String.(keys(record))) == RECORD_KEYS ||
            push!(errors, "records[$index] has unexpected or missing keys")
        id = get(record, "id", nothing)
        id isa AbstractString || begin
            push!(errors, "records[$index].id must be a string")
            continue
        end
        id = String(id)
        push!(ids, id)
        expected = get(expected_by_id, id, nothing)
        isnothing(expected) && begin
            push!(errors, "$id is absent from registry")
            continue
        end
        get(record, "heading_title", nothing) == expected["heading_title"] ||
            push!(errors, "$id heading_title drift")
        get(record, "registry_record_sha256", nothing) == expected["record_sha256"] ||
            push!(errors, "$id registry_record_sha256 drift")
        get(record, "canonical", nothing) == _expected_canonical(expected) ||
            push!(errors, "$id canonical metadata does not match registry ownership")
        canonical_value = get(record, "canonical", Any[])
        if canonical_value isa AbstractVector
            for section in canonical_value
                section isa AbstractDict || continue
                path = get(section, "path", nothing)
                path isa AbstractString && push!(referenced_canonical_paths, String(path))
            end
        end
        whitelisted = get(record, "agents_whitelisted", nothing)
        whitelisted isa Bool || push!(errors, "$id agents_whitelisted must be Boolean")
        whitelisted == (id in whitelist_set) || push!(errors, "$id whitelist mismatch")

        lifecycle = get(record, "lifecycle", nothing)
        grant = get(record, "grant", nothing)
        lifecycle in LIFECYCLES || push!(errors, "$id has unknown lifecycle: $lifecycle")
        grant in GRANTS || push!(errors, "$id has unknown grant: $grant")
        surfaces = get(record, "surfaces", nothing)
        _sorted_unique_strings(surfaces) ||
            push!(errors, "$id.surfaces must be a sorted unique string list")
        if surfaces isa AbstractVector
            for surface in surfaces
                surface in SURFACES || push!(errors, "$id has unknown surface: $surface")
            end
        else
            surfaces = String[]
        end

        scope = get(record, "scope", nothing)
        scope isa AbstractString && !isempty(strip(scope)) ||
            push!(errors, "$id.scope must be nonempty")
        dependencies = get(record, "dependency_ids", nothing)
        _sorted_unique_strings(dependencies) ||
            push!(errors, "$id.dependency_ids must be a sorted unique string list")
        if dependencies isa AbstractVector
            id in dependencies && push!(errors, "$id depends on itself")
            for dependency in dependencies
                haskey(expected_by_id, dependency) ||
                    push!(errors, "$id has unknown dependency: $dependency")
            end
        end

        grant == "none" && !isempty(surfaces) &&
            push!(errors, "$id grant none requires no surfaces")
        lifecycle in CLOSED_LIFECYCLES && whitelisted === true &&
            push!(errors, "$id closed lifecycle cannot be whitelisted")
        lifecycle in ("retired", "rejected") && grant != "none" &&
            push!(errors, "$id retired/rejected lifecycle requires grant none")
        lifecycle == "superseded" && !(grant in ("none", "preservation")) &&
            push!(errors, "$id superseded lifecycle permits only none/preservation")
        lifecycle in ("proposed", "deferred", "suspended") && whitelisted === true &&
            push!(errors, "$id non-executable lifecycle cannot be whitelisted")
        if whitelisted === true
            grant in EXECUTION_GRANTS ||
                push!(errors, "$id whitelist requires an execution grant")
            isempty(intersect(Set(surfaces), EXECUTION_SURFACES)) &&
                push!(errors, "$id whitelist requires an execution surface")
        elseif grant in ("implementation", "preservation", "retirement")
            push!(errors, "$id $grant grant requires whitelist membership")
        end
        whitelisted === false && grant == "maintenance" && surfaces != ["docs"] &&
            push!(errors, "$id non-whitelisted maintenance must be docs-only")
        grant == "measurement" && !isempty(surfaces) &&
            push!(errors, "$id measurement grant must not claim repository surfaces")
        grant == "design" && surfaces != ["docs"] &&
            push!(errors, "$id design grant must be docs-only")
        lifecycle == "approved" && whitelisted === true && grant != "implementation" &&
            push!(errors, "$id approved whitelist record requires implementation grant")
        lifecycle == "deferred" && whitelisted === true &&
            push!(errors, "$id deferred record cannot be whitelisted")
        grant == "implementation" && !(lifecycle in ("approved", "implemented")) &&
            push!(errors, "$id implementation grant has incompatible lifecycle")
        grant == "measurement" && !(lifecycle in ("approved", "completed", "deferred")) &&
            push!(errors, "$id measurement grant has incompatible lifecycle")
        grant == "design" && !(lifecycle in ("proposed", "approved", "completed")) &&
            push!(errors, "$id design grant has incompatible lifecycle")

        source_value = get(record, "source_paths", String[])
        test_value = get(record, "test_paths", String[])
        planned_value = get(record, "planned_test_paths", String[])
        source_paths = _sorted_unique_strings(source_value) ? String.(source_value) : String[]
        test_paths = _sorted_unique_strings(test_value) ? String.(test_value) : String[]
        planned_test_paths = _sorted_unique_strings(planned_value) ? String.(planned_value) : String[]
        "source" in surfaces && !any(startswith(path, "src/") for path in source_paths) &&
            push!(errors, "$id source surface requires an src path")
        "driver" in surfaces && !any(startswith(path, "bin/") for path in source_paths) &&
            push!(errors, "$id driver surface requires a bin path")
        "tools" in surfaces && !any(startswith(path, "tools/") for path in source_paths) &&
            push!(errors, "$id tools surface requires a tools path")
        "artifacts" in surfaces && isempty(source_paths) &&
            push!(errors, "$id artifacts surface requires an implementation path")
        "tests" in surfaces && isempty(test_paths) && isempty(planned_test_paths) &&
            push!(errors, "$id tests surface requires a present or planned test path")
        whitelisted === true && occursin("-TEST-", id) && surfaces != ["tests"] &&
            push!(errors, "$id whitelisted TEST record must be tests-only")
        !isempty(planned_test_paths) &&
            !(lifecycle == "approved" && grant == "implementation" && "tests" in surfaces) &&
            push!(errors, "$id planned tests require approved test implementation")

        append!(errors, _path_errors(record, id))
        append!(errors, _canonical_errors(record, id, canonical_document_paths))
    end

    ids == sort!(unique!(ids)) || push!(errors, "records must be sorted and unique by ID")
    Set(ids) == Set(keys(expected_by_id)) || push!(errors, "candidate registry ID coverage mismatch")
    length(records) == length(expected_records) || push!(errors, "candidate record count mismatch")
    canonical_document_paths == referenced_canonical_paths ||
        push!(errors, "canonical_documents must equal the referenced canonical path set")
    return sort!(unique!(errors))
end

function check_candidate()
    isfile(CANDIDATE_PATH) || error("missing authority candidate: $CANDIDATE_PATH")
    data = TOML.parsefile(CANDIDATE_PATH)
    errors = validation_errors(data)
    read(CANDIDATE_PATH, String) == _serialized(data) ||
        push!(errors, "candidate TOML serialization is not canonical")
    isempty(errors) || error("Cartesian authority candidate errors:\n" * join(first(errors, 100), "\n"))
    return nothing
end

function _expect_failure(data, needle)
    errors = validation_errors(data)
    any(error -> occursin(needle, error), errors) ||
        error("negative check did not produce $(repr(needle)); got $(join(errors, "; "))")
end

function self_test()
    data = TOML.parsefile(CANDIDATE_PATH)
    isempty(validation_errors(data)) || error("candidate must pass before self-test")

    broken = deepcopy(data)
    broken["records"][1]["lifecycle"] = "historical"
    _expect_failure(broken, "unknown lifecycle")

    broken = deepcopy(data)
    none_record = first(record for record in broken["records"] if record["grant"] == "none" && isempty(record["surfaces"]))
    none_record["surfaces"] = ["source"]
    _expect_failure(broken, "grant none requires no surfaces")

    broken = deepcopy(data)
    source_record = first(record for record in broken["records"] if "source" in record["surfaces"])
    empty!(source_record["source_paths"])
    _expect_failure(broken, "source surface requires an src path")

    broken = deepcopy(data)
    broken["records"][1]["canonical"][1]["heading"] = "Missing Candidate Heading"
    _expect_failure(broken, "heading does not exist")

    broken = deepcopy(data)
    whitelisted_record = first(record for record in broken["records"] if record["agents_whitelisted"])
    whitelisted_record["agents_whitelisted"] = false
    _expect_failure(broken, "whitelist mismatch")

    broken = deepcopy(data)
    traversal_record = first(record for record in broken["records"] if !isempty(record["source_paths"]))
    traversal_record["source_paths"][1] = "src/../Project.toml"
    _expect_failure(broken, "noncanonical or escapes the repository")

    broken = deepcopy(data)
    original_canonical = broken["records"][1]["canonical"]
    replacement = first(
        record for record in broken["records"]
        if record["canonical"] != original_canonical
    )
    broken["records"][1]["canonical"] = deepcopy(replacement["canonical"])
    _expect_failure(broken, "canonical metadata does not match registry ownership")

    broken = deepcopy(data)
    broken["schema_version"] = 2.0
    _expect_failure(broken, "schema_version must be integer 2")

    broken = deepcopy(data)
    broken["canonical_documents"][1]["sha256"] = repeat("0", 64)
    _expect_failure(broken, "canonical document content drift")

    broken = deepcopy(data)
    broken["records"][1]["scope"] *= " changed"
    _expect_failure(broken, "reviewed semantic metadata drift")

    broken = deepcopy(data)
    malformed_record = first(record for record in broken["records"] if !isempty(record["source_paths"]))
    malformed_record["source_paths"] = Any[1]
    _expect_failure(broken, "source_paths must be a sorted unique string list")

    broken = deepcopy(data)
    broken["dependency_semantics"] = "topological_build_order"
    _expect_failure(broken, "dependency_semantics mismatch")

    broken = deepcopy(data)
    extra_path = "docs/src/index.md"
    push!(broken["canonical_documents"], Dict(
        "path" => extra_path,
        "sha256" => bytes2hex(sha256(read(joinpath(ROOT, extra_path)))),
    ))
    sort!(broken["canonical_documents"]; by = document -> document["path"])
    _expect_failure(broken, "must equal the referenced canonical path set")

    mktempdir() do root
        mkpath(joinpath(root, "test"))
        link = joinpath(root, "test", "planned.jl")
        symlink(joinpath(dirname(root), "outside-missing.jl"), link)
        isnothing(_safe_repo_path(
            "test/planned.jl";
            allow_missing = true,
            root = root,
            root_real = realpath(root),
        )) || error("dangling planned-test symlink was accepted")
    end

    return nothing
end

function main(args = ARGS)
    if isempty(args) || args == ["--check"]
        check_candidate()
    elseif args == ["--self-test"]
        self_test()
    elseif args == ["--digest"]
        println(_semantic_digest(TOML.parsefile(CANDIDATE_PATH)))
    else
        error("usage: julia --project=docs docs/check_cartesian_authority_candidate.jl [--check|--self-test|--digest]")
    end
end

end

if abspath(PROGRAM_FILE) == @__FILE__
    CartesianAuthorityCandidate.main()
end
