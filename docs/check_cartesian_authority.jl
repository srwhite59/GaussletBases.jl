# Validate the authoritative Cartesian producer metadata and its generated views.
module CartesianAuthority

using Markdown
using SHA
using TOML

const ROOT = normpath(joinpath(@__DIR__, ".."))
const ROOT_REAL = realpath(ROOT)
const PRODUCER_DOCS = joinpath(
    ROOT,
    "docs",
    "src",
    "developer",
    "designs",
    "cartesian_hamiltonian_producer",
)
const AUTHORITY_PATH = joinpath(PRODUCER_DOCS, "authority.toml")
const REGISTRY_PATH = joinpath(PRODUCER_DOCS, "registry.md")
const AGENTS_PATH = joinpath(ROOT, "AGENTS.md")
const ARTIFACT_KIND = "cartesian_authority"
const WHITELIST_BEGIN = "<!-- BEGIN CARTESIAN HAMILTONIAN PRODUCER EXECUTION WHITELIST -->"
const WHITELIST_END = "<!-- END CARTESIAN HAMILTONIAN PRODUCER EXECUTION WHITELIST -->"
const LEGACY_PATHS = [
    joinpath(PRODUCER_DOCS, "authority_candidate.toml"),
    joinpath(PRODUCER_DOCS, "authority_transition_snapshot.toml"),
    joinpath(PRODUCER_DOCS, "registry_whitelist_shadow.toml"),
    joinpath(@__DIR__, "check_cartesian_authority_candidate.jl"),
    joinpath(@__DIR__, "check_cartesian_authority_shadow.jl"),
    joinpath(@__DIR__, "check_cartesian_authority_transition.jl"),
]

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
    "measurement",
])
const DOCUMENT_KINDS = Set(["canonical", "history", "evidence"])
const PATH_KINDS = Set(["source", "test", "tool", "driver", "docs", "measurement"])
const PATH_STATES = Set(["existing", "planned", "optional_local"])
const EVIDENCE_KINDS = Set(["git_commit", "manager_pass", "repo_path", "external_path"])
const SAFE_REPO_PATH = r"^[A-Za-z0-9._/-]+$"
const EXECUTION_GRANTS = Set(["implementation", "maintenance", "preservation", "retirement"])
const EXECUTION_SURFACES = Set(["source", "tests", "tools", "artifacts", "driver"])
const CLOSED_LIFECYCLES = Set(["superseded", "retired", "rejected"])

const TOP_LEVEL_KEYS = Set([
    "schema_version",
    "artifact_kind",
    "authoritative",
    "authorization_complete",
    "documents",
    "records",
])
const DOCUMENT_KEYS = Set(["path", "sha256"])
const RECORD_KEYS = Set([
    "id",
    "title",
    "lifecycle",
    "grant",
    "surfaces",
    "scope",
    "dependency_ids",
    "document_refs",
    "paths",
    "evidence_refs",
])
const DOCUMENT_REF_KEYS = Set(["kind", "path", "heading"])
const PATH_KEYS = Set(["kind", "state", "path"])
const EVIDENCE_REF_KEYS = Set(["kind", "value"])

struct AuthoritySnapshot
    path::String
    text::String
    data::Dict{String,Any}
    sha256::String
end

digest(bytes) = bytes2hex(sha256(bytes))
digest(text::AbstractString) = digest(codeunits(String(text)))

function serialized(data)
    io = IOBuffer()
    TOML.print(io, data; sorted = true)
    text = String(take!(io))
    return endswith(text, '\n') ? text : text * "\n"
end

function load_snapshot(path = AUTHORITY_PATH)
    text = read(path, String)
    data = TOML.parse(text)
    return AuthoritySnapshot(abspath(path), text, data, digest(text))
end

function _sorted_unique_strings(values)
    values isa AbstractVector || return false
    all(value -> value isa AbstractString, values) || return false
    return values == sort!(unique!(String.(values)))
end

function _sorted_unique_tables(values, key)
    values isa AbstractVector || return false
    all(value -> value isa AbstractDict, values) || return false
    keys = try
        [key(value) for value in values]
    catch
        return false
    end
    return keys == sort!(unique!(keys))
end

function _has_glob_or_ellipsis(path)
    return occursin(r"[*?\[\]{}]", path) || occursin("...", path)
end

function _single_line_text(value)
    value isa AbstractString || return false
    isempty(strip(value)) && return false
    occursin('\n', value) && return false
    occursin('\r', value) && return false
    return all(character -> character >= ' ' && character != Char(0x7f), value)
end

function _contains_symlink(path, root = ROOT)
    relative = relpath(path, root)
    current = root
    for component in splitpath(relative)
        current = joinpath(current, component)
        islink(current) && return true
        ispath(current) || break
    end
    return false
end

function safe_repo_path(path; allow_missing = false, root = ROOT, root_real = ROOT_REAL)
    path isa AbstractString || return nothing
    isabspath(path) && return nothing
    normalized = normpath(String(path))
    normalized == path || return nothing
    (normalized == ".." || startswith(normalized, "../")) && return nothing
    occursin(SAFE_REPO_PATH, normalized) || return nothing
    _has_glob_or_ellipsis(normalized) && return nothing
    absolute = joinpath(root, normalized)
    _contains_symlink(absolute, root) && return nothing
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

function _tracked_paths(root = ROOT)
    output = read(`git -C $root ls-files -z`, String)
    return Set(filter(!isempty, split(output, '\0')))
end

function _git_ignored(path, root = ROOT)
    return success(pipeline(
        `git -C $root check-ignore --quiet -- $path`;
        stdout = devnull,
        stderr = devnull,
    ))
end

function _normalized_lines(text)
    return split(replace(String(text), "\r\n" => "\n", '\r' => '\n'), '\n'; keepempty = true)
end

function _without_html_comments(line, in_comment)
    parts = String[]
    rest = String(line)
    while true
        if in_comment
            pieces = split(rest, "-->"; limit = 2)
            length(pieces) == 1 && return join(parts), true
            rest = String(pieces[2])
            in_comment = false
        else
            pieces = split(rest, "<!--"; limit = 2)
            push!(parts, String(pieces[1]))
            length(pieces) == 1 && return join(parts), false
            rest = String(pieces[2])
            in_comment = true
        end
    end
end

function markdown_headings(text)
    headings = String[]
    fence = nothing
    in_comment = false
    for line in _normalized_lines(text)
        stripped = lstrip(line)
        marker = startswith(stripped, "```") ? "```" : startswith(stripped, "~~~") ? "~~~" : nothing
        if !isnothing(fence)
            marker == fence && (fence = nothing)
            continue
        end
        uncommented, in_comment = _without_html_comments(line, in_comment)
        stripped = lstrip(uncommented)
        marker = startswith(stripped, "```") ? "```" : startswith(stripped, "~~~") ? "~~~" : nothing
        if !isnothing(marker)
            fence = marker
            continue
        end
        matched = match(r"^ {0,3}#{1,6}\s+(.+?)\s*$", uncommented)
        isnothing(matched) || push!(headings, strip(matched.captures[1]))
    end
    isnothing(fence) || error("unterminated Markdown fence")
    in_comment && error("unterminated Markdown HTML comment")
    return headings
end

function execution_ids(data)
    ids = String[]
    for record in data["records"]
        grant = get(record, "grant", nothing)
        surfaces = Set(String.(get(record, "surfaces", String[])))
        grant in EXECUTION_GRANTS && !isempty(intersect(surfaces, EXECUTION_SURFACES)) &&
            push!(ids, String(record["id"]))
    end
    return sort!(ids)
end

function _document_errors(data, tracked)
    errors = String[]
    documents = get(data, "documents", nothing)
    documents isa AbstractVector || return ["documents must be a list"], Dict{String,String}(), Dict{String,Vector{String}}()
    paths = String[]
    bytes_by_path = Dict{String,String}()
    headings_by_path = Dict{String,Vector{String}}()
    for (index, document) in enumerate(documents)
        document isa AbstractDict || begin
            push!(errors, "documents[$index] must be a table")
            continue
        end
        Set(String.(keys(document))) == DOCUMENT_KEYS ||
            push!(errors, "documents[$index] has unexpected keys")
        path = get(document, "path", nothing)
        expected_digest = get(document, "sha256", nothing)
        _single_line_text(path) && expected_digest isa AbstractString || begin
            push!(errors, "documents[$index] path/sha256 must be strings")
            continue
        end
        push!(paths, String(path))
        absolute = safe_repo_path(path)
        isnothing(absolute) && begin
            push!(errors, "documents[$index] path is invalid: $path")
            continue
        end
        startswith(path, "docs/") && endswith(path, ".md") ||
            push!(errors, "documents[$index] must name docs Markdown: $path")
        path in tracked || push!(errors, "documents[$index] is not tracked: $path")
        text = read(absolute, String)
        bytes_by_path[String(path)] = text
        headings_by_path[String(path)] = markdown_headings(text)
        occursin(r"^[0-9a-f]{64}$", expected_digest) ||
            push!(errors, "documents[$index] has invalid sha256")
        digest(text) == expected_digest || push!(errors, "document content drift: $path")
    end
    paths == sort!(unique!(paths)) || push!(errors, "documents must be sorted and unique")
    return errors, bytes_by_path, headings_by_path
end

function _path_errors(item, id, index, tracked)
    errors = String[]
    Set(String.(keys(item))) == PATH_KEYS || push!(errors, "$id.paths[$index] has unexpected keys")
    kind = get(item, "kind", nothing)
    state = get(item, "state", nothing)
    path = get(item, "path", nothing)
    kind in PATH_KINDS || push!(errors, "$id.paths[$index] has unknown kind: $kind")
    state in PATH_STATES || push!(errors, "$id.paths[$index] has unknown state: $state")
    _single_line_text(path) || return vcat(errors, ["$id.paths[$index].path must be nonempty single-line text"])
    kind in PATH_KINDS || return errors
    state in PATH_STATES || return errors
    allow_missing = state != "existing"
    absolute = safe_repo_path(path; allow_missing = allow_missing)
    isnothing(absolute) && return vcat(errors, ["$id.paths[$index] is invalid: $path"])
    ispath(absolute) && !isfile(absolute) &&
        push!(errors, "$id.paths[$index] existing filesystem object must be a file: $path")

    prefixes = Dict(
        "source" => "src/",
        "test" => "test/",
        "tool" => "tools/",
        "driver" => "bin/",
        "docs" => "docs/",
        "measurement" => "tmp/work/",
    )
    startswith(path, prefixes[kind]) || push!(errors, "$id.paths[$index] has wrong prefix: $path")
    state == "existing" && (!isfile(absolute) || !(path in tracked)) &&
        push!(errors, "$id.paths[$index] existing path must be a tracked file: $path")
    state == "planned" && path in tracked && push!(errors, "$id.paths[$index] planned path is tracked: $path")
    state == "planned" && kind != "test" && push!(errors, "$id.paths[$index] planned state is test-only")
    state == "optional_local" && kind != "measurement" &&
        push!(errors, "$id.paths[$index] optional_local state is measurement-only")
    state == "optional_local" && path in tracked &&
        push!(errors, "$id.paths[$index] optional_local path is tracked: $path")
    state == "optional_local" && !_git_ignored(path) &&
        push!(errors, "$id.paths[$index] optional_local path is not ignored: $path")
    kind == "measurement" && state != "optional_local" &&
        push!(errors, "$id.paths[$index] measurement path must be optional_local")
    kind == "docs" && !endswith(path, ".md") && push!(errors, "$id.paths[$index] docs path must be Markdown")
    return errors
end

function _dependency_errors(records)
    errors = String[]
    by_id = Dict{String,Any}()
    for record in records
        record isa AbstractDict || continue
        id = get(record, "id", nothing)
        id isa AbstractString || continue
        by_id[String(id)] = record
    end

    dependencies(record) = begin
        values = get(record, "dependency_ids", nothing)
        values isa AbstractVector || return String[]
        return String[value for value in values if value isa AbstractString]
    end

    for (id, record) in by_id
        grant = get(record, "grant", nothing)
        surface_values = get(record, "surfaces", nothing)
        surfaces = surface_values isa AbstractVector ?
            Set(String[value for value in surface_values if value isa AbstractString]) :
            Set{String}()
        executable = grant in EXECUTION_GRANTS &&
            !isempty(intersect(surfaces, EXECUTION_SURFACES))
        executable || continue
        for dependency in dependencies(record)
            haskey(by_id, dependency) || continue
            prerequisite = by_id[dependency]
            prerequisite_lifecycle = get(prerequisite, "lifecycle", nothing)
            prerequisite_grant = get(prerequisite, "grant", nothing)
            (prerequisite_lifecycle in CLOSED_LIFECYCLES || prerequisite_grant == "none") &&
                push!(errors,
                    "$id execution record has inactive dependency: $dependency")
        end
    end

    state = Dict(id => 0 for id in keys(by_id))
    stack = String[]
    function visit(id)
        state[id] == 2 && return
        if state[id] == 1
            start = findfirst(==(id), stack)
            cycle = vcat(stack[start:end], id)
            push!(errors, "dependency cycle: $(join(cycle, " -> "))")
            return
        end
        state[id] = 1
        push!(stack, id)
        for dependency in dependencies(by_id[id])
            haskey(by_id, dependency) && visit(dependency)
        end
        pop!(stack)
        state[id] = 2
        return
    end
    foreach(visit, sort!(collect(keys(by_id))))
    return errors
end

function validation_errors(data; tracked = _tracked_paths())
    errors = String[]
    Set(String.(keys(data))) == TOP_LEVEL_KEYS ||
        push!(errors, "authority has unexpected or missing top-level keys")
    get(data, "schema_version", nothing) === 3 || push!(errors, "schema_version must be integer 3")
    get(data, "artifact_kind", nothing) == ARTIFACT_KIND || push!(errors, "artifact_kind mismatch")
    get(data, "authoritative", nothing) === true || push!(errors, "authority must set authoritative=true")
    get(data, "authorization_complete", nothing) === true ||
        push!(errors, "authority must set authorization_complete=true")

    document_errors, document_bytes, headings_by_path = _document_errors(data, tracked)
    append!(errors, document_errors)
    document_paths = Set(keys(document_bytes))
    referenced_documents = Set{String}()

    records = get(data, "records", nothing)
    records isa AbstractVector || return vcat(errors, ["records must be a list"])
    ids = String[]
    for (index, record) in enumerate(records)
        record isa AbstractDict || begin
            push!(errors, "records[$index] must be a table")
            continue
        end
        Set(String.(keys(record))) == RECORD_KEYS || push!(errors, "records[$index] has unexpected keys")
        id = get(record, "id", nothing)
        id isa AbstractString || begin
            push!(errors, "records[$index].id must be a string")
            continue
        end
        id = String(id)
        push!(ids, id)
        occursin(r"^HP-[A-Z0-9]+(?:-[A-Z0-9]+)*$", id) || push!(errors, "$id is malformed")
        title = get(record, "title", nothing)
        _single_line_text(title) || push!(errors, "$id.title must be nonempty single-line text")
        lifecycle = get(record, "lifecycle", nothing)
        grant = get(record, "grant", nothing)
        lifecycle in LIFECYCLES || push!(errors, "$id has unknown lifecycle: $lifecycle")
        grant in GRANTS || push!(errors, "$id has unknown grant: $grant")
        surfaces = get(record, "surfaces", nothing)
        _sorted_unique_strings(surfaces) || push!(errors, "$id.surfaces must be sorted unique strings")
        surfaces = surfaces isa AbstractVector ? String.(surfaces) : String[]
        all(surface -> surface in SURFACES, surfaces) || push!(errors, "$id has unknown surface")
        scope = get(record, "scope", nothing)
        _single_line_text(scope) || push!(errors, "$id.scope must be nonempty single-line text")
        dependencies = get(record, "dependency_ids", nothing)
        _sorted_unique_strings(dependencies) || push!(errors, "$id.dependency_ids must be sorted unique strings")
        dependencies = dependencies isa AbstractVector ? String.(dependencies) : String[]
        id in dependencies && push!(errors, "$id depends on itself")

        refs = get(record, "document_refs", nothing)
        _sorted_unique_tables(refs, item -> (item["kind"], item["path"], item["heading"])) ||
            push!(errors, "$id.document_refs must be sorted unique tables")
        refs = refs isa AbstractVector ? refs : Any[]
        canonical_count = 0
        for (ref_index, ref) in enumerate(refs)
            ref isa AbstractDict || continue
            Set(String.(keys(ref))) == DOCUMENT_REF_KEYS ||
                push!(errors, "$id.document_refs[$ref_index] has unexpected keys")
            kind = get(ref, "kind", nothing)
            path = get(ref, "path", nothing)
            heading = get(ref, "heading", nothing)
            kind in DOCUMENT_KINDS || push!(errors, "$id.document_refs[$ref_index] has unknown kind")
            _single_line_text(path) && _single_line_text(heading) || begin
                push!(errors, "$id.document_refs[$ref_index] path/heading must be nonempty single-line text")
                continue
            end
            push!(referenced_documents, String(path))
            path in document_paths || push!(errors, "$id.document_refs[$ref_index] has unknown document: $path")
            available = get(headings_by_path, String(path), String[])
            count(==(heading), available) == 1 ||
                push!(errors, "$id.document_refs[$ref_index] heading is absent or nonunique: $heading")
            kind == "canonical" && (canonical_count += 1)
        end
        isempty(refs) && push!(errors, "$id.document_refs must not be empty")
        grant in EXECUTION_GRANTS && canonical_count == 0 && push!(errors, "$id execution grant needs a canonical document")
        grant == "design" && canonical_count == 0 && push!(errors, "$id design grant needs a canonical document")
        grant == "measurement" && lifecycle in ("approved", "deferred") && canonical_count == 0 &&
            push!(errors, "$id active measurement needs a canonical document")

        paths = get(record, "paths", nothing)
        _sorted_unique_tables(paths, item -> (item["kind"], item["state"], item["path"])) ||
            push!(errors, "$id.paths must be sorted unique tables")
        paths = paths isa AbstractVector ? paths : Any[]
        for (path_index, item) in enumerate(paths)
            item isa AbstractDict || continue
            append!(errors, _path_errors(item, id, path_index, tracked))
        end
        path_kinds = Set(String(item["kind"]) for item in paths if item isa AbstractDict && haskey(item, "kind"))

        evidence = get(record, "evidence_refs", nothing)
        _sorted_unique_tables(evidence, item -> (item["kind"], item["value"])) ||
            push!(errors, "$id.evidence_refs must be sorted unique tables")
        evidence = evidence isa AbstractVector ? evidence : Any[]
        for (evidence_index, item) in enumerate(evidence)
            item isa AbstractDict || continue
            Set(String.(keys(item))) == EVIDENCE_REF_KEYS ||
                push!(errors, "$id.evidence_refs[$evidence_index] has unexpected keys")
            kind = get(item, "kind", nothing)
            value = get(item, "value", nothing)
            kind in EVIDENCE_KINDS || push!(errors, "$id.evidence_refs[$evidence_index] has unknown kind")
            _single_line_text(value) || begin
                push!(errors, "$id.evidence_refs[$evidence_index].value must be nonempty single-line text")
                continue
            end
            kind == "git_commit" && !occursin(r"^[0-9a-f]{7,40}$", value) &&
                push!(errors, "$id has malformed git_commit evidence")
            kind == "manager_pass" && !occursin(r"^[1-9][0-9]*$", value) &&
                push!(errors, "$id has malformed manager_pass evidence")
            if kind == "repo_path"
                absolute = safe_repo_path(value)
                isnothing(absolute) && push!(errors, "$id has invalid repo_path evidence: $value")
                value in tracked || push!(errors, "$id repo_path evidence is not tracked: $value")
            elseif kind == "external_path"
                isabspath(value) || push!(errors, "$id external_path evidence must be absolute")
            end
        end

        grant == "none" && (!isempty(surfaces) || !isempty(paths)) &&
            push!(errors, "$id grant none requires no surfaces or owned paths")
        lifecycle in CLOSED_LIFECYCLES && grant != "none" &&
            push!(errors, "$id closed lifecycle requires grant none")
        lifecycle in ("proposed", "deferred", "suspended") && grant in EXECUTION_GRANTS &&
            push!(errors, "$id non-executable lifecycle has execution grant")
        lifecycle in ("retired", "rejected") && grant != "none" &&
            push!(errors, "$id retired/rejected lifecycle requires grant none")
        grant == "design" && !(lifecycle in ("proposed", "approved", "completed")) &&
            push!(errors, "$id design grant has incompatible lifecycle")
        grant == "measurement" && !(lifecycle in ("approved", "completed", "deferred")) &&
            push!(errors, "$id measurement grant has incompatible lifecycle")
        grant == "implementation" && !(lifecycle in ("approved", "implemented")) &&
            push!(errors, "$id implementation grant has incompatible lifecycle")
        grant == "maintenance" && !(lifecycle in ("implemented", "completed")) &&
            push!(errors, "$id maintenance grant has incompatible lifecycle")
        grant == "preservation" && lifecycle != "implemented" &&
            push!(errors, "$id preservation grant requires implemented lifecycle")
        grant == "retirement" && lifecycle != "approved" &&
            push!(errors, "$id retirement grant requires approved lifecycle")
        grant == "design" && surfaces != ["docs"] && push!(errors, "$id design grant must be docs-only")
        grant == "measurement" && !(Set(surfaces) <= Set(["docs", "measurement"])) &&
            push!(errors, "$id measurement grant may own only docs/measurement surfaces")
        grant == "measurement" && !("measurement" in surfaces) &&
            push!(errors, "$id measurement grant requires measurement surface")
        grant in EXECUTION_GRANTS && isempty(intersect(Set(surfaces), EXECUTION_SURFACES)) &&
            !(grant == "maintenance" && surfaces == ["docs"]) &&
            push!(errors, "$id execution grant lacks execution surface")
        grant in ("implementation", "preservation", "retirement") &&
            isempty(intersect(Set(surfaces), EXECUTION_SURFACES)) &&
            push!(errors, "$id $grant grant lacks execution surface")

        required_path_kind = Dict(
            "source" => "source",
            "tests" => "test",
            "tools" => "tool",
            "driver" => "driver",
            "docs" => "docs",
            "measurement" => "measurement",
        )
        for surface in surfaces
            surface == "artifacts" && continue
            required_path_kind[surface] in path_kinds ||
                push!(errors, "$id surface $surface lacks matching owned path")
        end
        "artifacts" in surfaces && isempty(intersect(path_kinds, Set(["source", "driver"]))) &&
            push!(errors, "$id artifacts surface lacks source/driver implementation path")
        for kind in path_kinds
            surface = Dict("source" => "source", "test" => "tests", "tool" => "tools", "driver" => "driver", "docs" => "docs", "measurement" => "measurement")[kind]
            surface in surfaces || push!(errors, "$id owned path kind $kind lacks surface $surface")
        end
    end

    ids == sort!(unique!(ids)) || push!(errors, "records must be sorted and unique by ID")
    id_set = Set(ids)
    for record in records, dependency in get(record, "dependency_ids", String[])
        dependency in id_set || push!(errors, "$(record["id"]) has unknown dependency: $dependency")
    end
    append!(errors, _dependency_errors(records))
    document_paths == referenced_documents ||
        push!(errors, "documents must exactly equal the referenced document set")
    return sort!(unique!(errors))
end

function assert_snapshot_unchanged(snapshot, document_paths)
    read(snapshot.path, String) == snapshot.text || error("authority changed during validation/rendering")
    for document in snapshot.data["documents"]
        document["path"] in document_paths || continue
        digest(read(joinpath(ROOT, document["path"]))) == document["sha256"] ||
            error("canonical document changed during validation/rendering: $(document["path"])")
    end
    return nothing
end

function legacy_live_errors(paths = LEGACY_PATHS; tracked = _tracked_paths())
    errors = String[]
    for path in paths
        relative = relpath(path, ROOT)
        (ispath(path) || islink(path) || relative in tracked) &&
            push!(errors, "legacy authority artifact remains live: $path")
    end
    return errors
end

function check_authority(snapshot = load_snapshot())
    checker_text = read(@__FILE__, String)
    registry_text = read(REGISTRY_PATH, String)
    agents_text = read(AGENTS_PATH, String)
    errors = validation_errors(snapshot.data)
    snapshot.text == serialized(snapshot.data) || push!(errors, "authority TOML serialization is not canonical")
    append!(errors, legacy_live_errors())
    isempty(errors) || error("Cartesian authority errors:\n" * join(first(errors, 100), "\n"))
    assert_snapshot_unchanged(snapshot, Set(document["path"] for document in snapshot.data["documents"]))
    registry = render_registry(snapshot)
    whitelist = render_whitelist_block(snapshot)
    registry == render_registry(snapshot) || error("registry render is nondeterministic")
    whitelist == render_whitelist_block(snapshot) || error("whitelist render is nondeterministic")
    assert_render_structure(snapshot, registry, whitelist)
    validate_committed_views(snapshot, registry, whitelist; registry_text, agents_text)
    assert_snapshot_unchanged(snapshot, Set(document["path"] for document in snapshot.data["documents"]))
    read(REGISTRY_PATH, String) == registry_text || error("registry changed during authority check")
    read(AGENTS_PATH, String) == agents_text || error("AGENTS changed during authority check")
    read(@__FILE__, String) == checker_text || error("authority checker changed during authority check")
    return nothing
end

function code_span(value)
    text = String(value)
    runs = [length(m.match) for m in eachmatch(r"`+", text)]
    fence = repeat("`", isempty(runs) ? 1 : maximum(runs) + 1)
    padding = startswith(text, '`') || endswith(text, '`') || startswith(text, ' ') || endswith(text, ' ') ? " " : ""
    return fence * padding * text * padding * fence
end

function markdown_text(value)
    text = String(value)
    for character in ('\\', '`', '*', '_', '[', ']', '<', '>', '#', '!', '|')
        text = replace(text, string(character) => "\\" * string(character))
    end
    return replace(text, '\n' => ' ')
end

function _list(values; empty = "none")
    isempty(values) && return empty
    return join((code_span(value) for value in values), ", ")
end

function render_whitelist_block(snapshot = load_snapshot())
    io = IOBuffer()
    println(io, WHITELIST_BEGIN)
    println(io, "> **Generated authority view. Do not edit this block.**")
    println(io, "> Source: `docs/src/developer/designs/cartesian_hamiltonian_producer/authority.toml`.")
    println(io, "> Authority SHA-256: `$(snapshot.sha256)`.")
    println(io)
    println(io, "Cartesian Hamiltonian producer source work is currently authorized only for")
    println(io, "these approved design IDs:")
    println(io)
    for id in execution_ids(snapshot.data)
        println(io, "- `", id, "`")
    end
    println(io)
    println(io, WHITELIST_END)
    return chomp(String(take!(io))) * "\n"
end

function render_registry(snapshot = load_snapshot())
    io = IOBuffer()
    println(io, "# Cartesian Hamiltonian Producer Authority Registry")
    println(io)
    println(io, "> **Generated authority view. Do not edit.** The record-level source is")
    println(io, "> [authority.toml](authority.toml), SHA-256 `$(snapshot.sha256)`.")
    println(io)
    println(io, "Tracked producer work is authorized only when a unique record has an")
    println(io, "execution grant and surface, and the requested change stays within its exact")
    println(io, "owned paths, scope, `current.md`, `invariants.md`, and canonical contract.")
    println(io, "Lifecycle never grants work by itself. Any missing or conflicting fact fails closed.")
    println(io)
    println(io, "## Records")
    println(io)
    for record in snapshot.data["records"]
        println(io, "### $(record["id"]) - $(markdown_text(record["title"]))")
        println(io)
        println(io, "- **Lifecycle:** $(code_span(record["lifecycle"]))")
        println(io, "- **Grant:** $(code_span(record["grant"]))")
        println(io, "- **Surfaces:** $(_list(record["surfaces"]))")
        executable = record["id"] in execution_ids(snapshot.data)
        println(io, "- **Execution whitelist:** $(code_span(string(executable)))")
        println(io, "- **Documents:**")
        for ref in record["document_refs"]
            relative = relpath(joinpath(ROOT, ref["path"]), PRODUCER_DOCS)
            label = markdown_text(basename(ref["path"]))
            println(io, "  - $(code_span(ref["kind"])) [$(label)]($(relative)); heading $(code_span(ref["heading"]))")
        end
        paths = record["paths"]
        if isempty(paths)
            println(io, "- **Owned paths:** none")
        else
            println(io, "- **Owned paths:**")
            for item in paths
                println(io, "  - $(code_span(item["kind"])) / $(code_span(item["state"])): $(code_span(item["path"]))")
            end
        end
        evidence = record["evidence_refs"]
        if isempty(evidence)
            println(io, "- **Evidence:** none")
        else
            println(io, "- **Evidence:**")
            for item in evidence
                println(io, "  - $(code_span(item["kind"])): $(code_span(item["value"]))")
            end
        end
        println(io, "- **Dependencies:** $(_list(record["dependency_ids"]))")
        println(io, "- **Scope:** $(markdown_text(record["scope"]))")
        println(io)
    end
    return rstrip(String(take!(io))) * "\n"
end

function assert_render_structure(snapshot, registry, whitelist)
    registry_lines = _normalized_lines(registry)
    record_headings = count(line -> startswith(line, "### HP-"), registry_lines)
    record_headings == length(snapshot.data["records"]) ||
        error("registry record-heading count mismatch")
    count(==(WHITELIST_BEGIN), _normalized_lines(whitelist)) == 1 ||
        error("whitelist begin-marker count mismatch")
    count(==(WHITELIST_END), _normalized_lines(whitelist)) == 1 ||
        error("whitelist end-marker count mismatch")
    occursin(WHITELIST_BEGIN, registry) && error("registry contains whitelist begin marker")
    occursin(WHITELIST_END, registry) && error("registry contains whitelist end marker")
    Markdown.parse(registry)
    Markdown.parse(whitelist)
    return nothing
end

function _standalone_marker(text, range)
    before = first(range) == firstindex(text) || text[prevind(text, first(range))] == '\n'
    after_index = nextind(text, last(range))
    after = after_index > lastindex(text) || text[after_index] in ('\r', '\n')
    return before && after
end

function marked_whitelist_block(text)
    starts = collect(findall(WHITELIST_BEGIN, text))
    stops = collect(findall(WHITELIST_END, text))
    length(starts) == 1 || error("expected one AGENTS whitelist begin marker")
    length(stops) == 1 || error("expected one AGENTS whitelist end marker")
    start, stop = only(starts), only(stops)
    first(start) < first(stop) || error("AGENTS whitelist markers are reversed")
    _standalone_marker(text, start) || error("AGENTS whitelist begin marker must occupy one line")
    _standalone_marker(text, stop) || error("AGENTS whitelist end marker must occupy one line")

    block_end = last(stop)
    following = nextind(text, block_end)
    if following <= lastindex(text) && text[following] == '\r'
        block_end = following
        following = nextind(text, following)
        following <= lastindex(text) && text[following] == '\n' && (block_end = following)
    elseif following <= lastindex(text) && text[following] == '\n'
        block_end = following
    end
    return String(SubString(text, first(start), block_end))
end

function validate_committed_views(
    snapshot,
    registry = render_registry(snapshot),
    whitelist = render_whitelist_block(snapshot);
    registry_text = read(REGISTRY_PATH, String),
    agents_text = read(AGENTS_PATH, String),
)
    registry_text == registry || error("generated registry is stale or partially edited")
    marked_whitelist_block(agents_text) == whitelist ||
        error("generated AGENTS whitelist block is stale or partially edited")
    return nothing
end

function _inside(path, root)
    relative = relpath(path, root)
    return !isabspath(relative) && relative != ".." && !startswith(relative, "../")
end

function _nearest_existing_ancestor(path)
    current = normpath(path)
    while !ispath(current)
        parent = dirname(current)
        parent == current && return nothing
        current = parent
    end
    return current
end

function _prepare_external_output_dir(path)
    output = normpath(abspath(path))
    ispath(output) && error("render output must be a new directory: $output")
    ancestor = _nearest_existing_ancestor(output)
    isnothing(ancestor) && error("render output has no existing ancestor")
    isdir(ancestor) || error("render output ancestor is not a directory: $ancestor")
    realpath(ancestor) == normpath(ancestor) ||
        error("render output path contains a symlink: $ancestor")
    _inside(realpath(ancestor), ROOT_REAL) &&
        error("render output must be outside the repository")
    mkpath(output)
    isdir(output) || error("render output is not a directory: $output")
    real = realpath(output)
    real == output || error("render output path contains a symlink: $output")
    _inside(real, ROOT_REAL) && error("render output must be outside the repository")
    return real
end

function _atomic_write(path, content; parent_real)
    dirname(path) == parent_real || error("atomic-write destination has unexpected parent")
    realpath(dirname(path)) == parent_real || error("atomic-write parent changed")
    islink(path) && error("atomic-write destination is a symlink: $path")
    ispath(path) && !isfile(path) && error("atomic-write destination is a non-file: $path")
    temporary, io = mktemp(dirname(path))
    try
        write(io, content)
        flush(io)
        ccall(:fsync, Cint, (Cint,), fd(io)) == 0 || error("fsync failed for $path")
        close(io)
        realpath(dirname(path)) == parent_real || error("atomic-write parent changed")
        islink(path) && error("atomic-write destination became a symlink: $path")
        ispath(path) && !isfile(path) && error("atomic-write destination became a non-file: $path")
        mv(temporary, path; force = true)
        realpath(dirname(path)) == parent_real || error("atomic-write parent changed")
    catch
        isopen(io) && close(io)
        ispath(temporary) && rm(temporary; force = true)
        rethrow()
    end
    return path
end

function render_external(output_dir, snapshot = load_snapshot())
    checker_text = read(@__FILE__, String)
    errors = validation_errors(snapshot.data)
    snapshot.text == serialized(snapshot.data) || push!(errors, "authority TOML serialization is not canonical")
    append!(errors, legacy_live_errors())
    isempty(errors) || error("Cartesian authority errors:\n" * join(first(errors, 100), "\n"))
    registry = render_registry(snapshot)
    whitelist = render_whitelist_block(snapshot)
    assert_render_structure(snapshot, registry, whitelist)
    assert_snapshot_unchanged(snapshot, Set(document["path"] for document in snapshot.data["documents"]))
    read(@__FILE__, String) == checker_text || error("authority checker changed during rendering")
    output_dir = _prepare_external_output_dir(output_dir)
    outputs = Dict(
        "agents_whitelist_block.md" => whitelist,
        "registry.md" => registry,
    )
    for (name, content) in outputs
        _atomic_write(joinpath(output_dir, name), content; parent_real = output_dir)
    end
    assert_snapshot_unchanged(snapshot, Set(document["path"] for document in snapshot.data["documents"]))
    read(@__FILE__, String) == checker_text || error("authority checker changed during rendering")
    return output_dir
end

function _expect_failure(data, needle)
    errors = validation_errors(data)
    any(error -> occursin(needle, error), errors) ||
        error("negative check did not produce $(repr(needle)); got $(join(errors, "; "))")
end

function _expect_error(operation, needle)
    error_value = try
        operation()
        nothing
    catch caught
        caught
    end
    isnothing(error_value) &&
        error("negative check unexpectedly passed; expected $(repr(needle))")
    message = sprint(showerror, error_value)
    occursin(needle, message) || throw(error_value)
    return nothing
end

function self_test()
    snapshot = load_snapshot()
    isempty(validation_errors(snapshot.data)) || error("authority must pass before self-test")
    registry = render_registry(snapshot)
    whitelist = render_whitelist_block(snapshot)
    validate_committed_views(snapshot, registry, whitelist)

    broken = deepcopy(snapshot.data)
    broken["authoritative"] = false
    _expect_failure(broken, "authoritative=true")

    broken = deepcopy(snapshot.data)
    broken["authorization_complete"] = false
    _expect_failure(broken, "authorization_complete=true")

    broken = deepcopy(snapshot.data)
    broken["schema_version"] = 4
    _expect_failure(broken, "schema_version")

    broken = deepcopy(snapshot.data)
    broken["artifact_kind"] = "cartesian_authority_candidate"
    _expect_failure(broken, "artifact_kind")

    broken = deepcopy(snapshot.data)
    broken["records"][1]["lifecycle"] = "historical"
    _expect_failure(broken, "unknown lifecycle")

    broken = deepcopy(snapshot.data)
    broken["records"][1]["grant"] = "caller-maintenance"
    _expect_failure(broken, "unknown grant")

    broken = deepcopy(snapshot.data)
    broken["records"][1]["surfaces"] = 1
    _expect_failure(broken, "surfaces must be sorted unique strings")

    broken = deepcopy(snapshot.data)
    source_record = first(record for record in broken["records"] if "source" in record["surfaces"])
    empty!(source_record["paths"])
    _expect_failure(broken, "lacks matching owned path")

    broken = deepcopy(snapshot.data)
    measurement_record = first(record for record in broken["records"] if record["grant"] == "measurement")
    measurement_record["paths"][end]["path"] = "tmp/work/*.jl"
    _expect_failure(broken, "is invalid")

    broken = deepcopy(snapshot.data)
    measurement_record = first(record for record in broken["records"] if record["grant"] == "measurement")
    measurement_record["paths"][end]["path"] = "tmp/work/unsafe).jl"
    _expect_failure(broken, "is invalid")

    broken = deepcopy(snapshot.data)
    broken["records"][1]["document_refs"][1]["heading"] = "Missing Heading"
    _expect_failure(broken, "heading is absent or nonunique")

    broken = deepcopy(snapshot.data)
    broken["documents"][1]["sha256"] = repeat("0", 64)
    _expect_failure(broken, "document content drift")

    broken = deepcopy(snapshot.data)
    broken["records"][1]["dependency_ids"] = ["HP-MISSING-01"]
    _expect_failure(broken, "unknown dependency")

    broken = deepcopy(snapshot.data)
    broken["records"][1]["dependency_ids"] = [broken["records"][1]["id"]]
    _expect_failure(broken, "depends on itself")

    broken = deepcopy(snapshot.data)
    first_record, second_record = broken["records"][1:2]
    first_record["dependency_ids"] = [second_record["id"]]
    second_record["dependency_ids"] = [first_record["id"]]
    _expect_failure(broken, "dependency cycle")

    broken = deepcopy(snapshot.data)
    execution_record = first(record for record in broken["records"] if
        record["grant"] in EXECUTION_GRANTS &&
        !isempty(intersect(Set(record["surfaces"]), EXECUTION_SURFACES)))
    inactive_record = first(record for record in broken["records"] if
        record["grant"] == "none")
    execution_record["dependency_ids"] = [inactive_record["id"]]
    _expect_failure(broken, "inactive dependency")

    broken = deepcopy(snapshot.data)
    broken["records"][1]["scope"] = "valid prefix\n## injected heading"
    _expect_failure(broken, "single-line text")

    broken = deepcopy(snapshot.data)
    planned_record = first(record for record in broken["records"] if any(
        item -> item["state"] == "planned",
        record["paths"],
    ))
    planned_path = first(item for item in planned_record["paths"] if item["state"] == "planned")
    planned_path["path"] = "test/nested"
    _expect_failure(broken, "must be a file")

    broken = deepcopy(snapshot.data)
    evidence_record = first(record for record in broken["records"] if !isempty(record["evidence_refs"]))
    evidence_record["evidence_refs"][1]["kind"] = "authority"
    _expect_failure(broken, "unknown kind")

    registry == render_registry(snapshot) || error("registry render changed")
    whitelist == render_whitelist_block(snapshot) || error("whitelist render changed")
    assert_render_structure(snapshot, registry, whitelist)

    agents_text = read(AGENTS_PATH, String)
    _expect_error("registry is stale") do
        validate_committed_views(
            snapshot,
            registry,
            whitelist;
            registry_text = registry * "\n",
            agents_text,
        )
    end
    _expect_error("AGENTS whitelist block is stale") do
        changed = replace(agents_text, first(execution_ids(snapshot.data)) => "HP-BROKEN-ID"; count = 1)
        validate_committed_views(snapshot, registry, whitelist; registry_text = registry, agents_text = changed)
    end
    _expect_error("begin marker") do
        changed = replace(agents_text, WHITELIST_BEGIN => ""; count = 1)
        validate_committed_views(snapshot, registry, whitelist; registry_text = registry, agents_text = changed)
    end

    mktempdir() do directory
        missing = joinpath(directory, "missing.toml")
        _expect_error("No such file") do
            load_snapshot(missing)
        end
        corrupt = joinpath(directory, "corrupt.toml")
        write(corrupt, "not = [valid")
        _expect_error("TOML") do
            load_snapshot(corrupt)
        end
        legacy = joinpath(directory, "legacy.toml")
        write(legacy, "legacy")
        isempty(legacy_live_errors([legacy])) &&
            error("legacy live artifact was not rejected")
        dangling = joinpath(directory, "dangling.toml")
        symlink(joinpath(directory, "missing-target"), dangling)
        isempty(legacy_live_errors([dangling])) &&
            error("dangling legacy symlink was not rejected")
    end

    mktempdir() do directory
        symlink(ROOT, joinpath(directory, "repo-link"))
        _expect_error("symlink") do
            _prepare_external_output_dir(joinpath(directory, "repo-link", "render"))
        end
        isnothing(safe_repo_path("../escape"; root = directory, root_real = realpath(directory))) ||
            error("repository traversal path was not rejected")

        output = joinpath(directory, "render")
        render_external(output, snapshot)
        read(joinpath(output, "registry.md"), String) == registry ||
            error("external registry render mismatch")
        read(joinpath(output, "agents_whitelist_block.md"), String) == whitelist ||
            error("external whitelist render mismatch")
        _expect_error("new directory") do
            _prepare_external_output_dir(output)
        end

        parent = realpath(directory)
        destination = joinpath(parent, "existing-directory")
        mkpath(destination)
        sentinel = joinpath(destination, "sentinel")
        write(sentinel, "preserve")
        try
            _atomic_write(destination, "replacement"; parent_real = parent)
            error("directory-valued atomic-write destination was not rejected")
        catch error_value
            occursin("non-file", sprint(showerror, error_value)) || rethrow()
        end
        read(sentinel, String) == "preserve" || error("atomic-write directory rejection lost data")
    end

    _expect_error("outside the repository") do
        _prepare_external_output_dir(tempname(joinpath(ROOT, "tmp")))
    end

    mktempdir() do directory
        copy = joinpath(directory, "authority.toml")
        write(copy, snapshot.text)
        local_snapshot = load_snapshot(copy)
        write(copy, snapshot.text * "\n")
        _expect_error("authority changed") do
            assert_snapshot_unchanged(local_snapshot, Set{String}())
        end
    end
    return nothing
end

function main(args = ARGS)
    if isempty(args) || args == ["--check"]
        check_authority()
    elseif args == ["--self-test"]
        self_test()
    elseif length(args) == 2 && args[1] == "--render"
        snapshot = load_snapshot()
        println(render_external(args[2], snapshot))
    else
        error("usage: julia --project=docs docs/check_cartesian_authority.jl [--check|--self-test|--render DIR]")
    end
end

end

if abspath(PROGRAM_FILE) == @__FILE__
    CartesianAuthority.main()
end
