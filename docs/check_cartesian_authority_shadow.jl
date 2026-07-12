# This mirrors exact registry text and raw whitelist membership only. It must
# never infer or grant effective source authorization.
module CartesianAuthorityShadow

using SHA
using TOML

const ROOT = normpath(joinpath(@__DIR__, ".."))
const PRODUCER_DOCS = joinpath(
    ROOT,
    "docs",
    "src",
    "developer",
    "designs",
    "cartesian_hamiltonian_producer",
)
const REGISTRY_PATH = joinpath(PRODUCER_DOCS, "registry.md")
const AGENTS_PATH = joinpath(ROOT, "AGENTS.md")
const SHADOW_PATH = joinpath(PRODUCER_DOCS, "registry_whitelist_shadow.toml")
const REGISTRY_SOURCE = "docs/src/developer/designs/cartesian_hamiltonian_producer/registry.md"
const AGENTS_SOURCE = "AGENTS.md"
const ARTIFACT_KIND = "cartesian_registry_whitelist_shadow"
const SCOPE = "registry_metadata_and_agents_whitelist_membership_only"
const WHITELIST_BEGIN = "<!-- BEGIN CARTESIAN HAMILTONIAN PRODUCER EXECUTION WHITELIST -->"
const WHITELIST_END = "<!-- END CARTESIAN HAMILTONIAN PRODUCER EXECUTION WHITELIST -->"
const EM_DASH_SEPARATOR = string(' ', Char(0x2014), ' ')
const ASCII_DASH_SEPARATOR = " - "
const OWNERSHIP_PREFIXES = (
    "Owner/canonical:",
    "Owner/history:",
    "Owner/source:",
    "Owner and canonical contract:",
    "Owner and canonical outcome:",
    "Owner and evidence:",
    "Owner:",
    "Canonical contract:",
    "Canonical contracts:",
    "Canonical boundary:",
    "Canonical record:",
    "Canonical design:",
    "Canonical lifecycle:",
    "Canonical owner:",
)

_normalized_lines(text::AbstractString) =
    split(replace(String(text), "\r\n" => "\n", '\r' => '\n'), '\n'; keepempty = true)

function _normalized_block(lines)
    cleaned = String[rstrip(String(line)) for line in lines]
    while !isempty(cleaned) && isempty(cleaned[end])
        pop!(cleaned)
    end
    return join(cleaned, "\n") * "\n"
end

_digest(text::AbstractString) = bytes2hex(sha256(codeunits(String(text))))

function _fence_marker(line::AbstractString)
    stripped = lstrip(String(line))
    startswith(stripped, "```") && return "```"
    startswith(stripped, "~~~") && return "~~~"
    return nothing
end

function _without_html_comments(line::AbstractString, in_comment::Bool)
    parts = String[]
    rest = String(line)
    while true
        if in_comment
            split_rest = split(rest, "-->"; limit = 2)
            length(split_rest) == 1 && return join(parts), true
            rest = String(split_rest[2])
            in_comment = false
        else
            split_rest = split(rest, "<!--"; limit = 2)
            push!(parts, String(split_rest[1]))
            length(split_rest) == 1 && return join(parts), false
            rest = String(split_rest[2])
            in_comment = true
        end
    end
end

function _visible_markdown_lines(lines)
    visible = Vector{Union{Nothing,String}}(undef, length(lines))
    fill!(visible, nothing)
    fence = nothing
    in_comment = false
    for (index, line) in pairs(lines)
        if !isnothing(fence)
            marker = _fence_marker(line)
            marker == fence && (fence = nothing)
            continue
        end
        uncommented, in_comment = _without_html_comments(line, in_comment)
        marker = _fence_marker(uncommented)
        if !isnothing(marker)
            fence = marker
            continue
        end
        visible[index] = uncommented
    end
    isnothing(fence) || error("unterminated Markdown fence")
    in_comment && error("unterminated Markdown HTML comment")
    return visible
end

function _visible_heading_levels(lines)
    levels = Dict{Int,Int}()
    for (index, line) in pairs(_visible_markdown_lines(lines))
        isnothing(line) && continue
        matched = match(r"^ {0,3}(#{1,3})\s+", line)
        isnothing(matched) || (levels[index] = length(matched.captures[1]))
    end
    return levels
end

function _parse_heading(line::AbstractString)
    stripped = lstrip(String(line), ' ')
    startswith(stripped, "### ") ||
        error("expected level-3 registry heading: $(repr(line))")
    body = String(stripped[5:end])
    parts = if occursin(EM_DASH_SEPARATOR, body)
        split(body, EM_DASH_SEPARATOR; limit = 2)
    elseif occursin(ASCII_DASH_SEPARATOR, body)
        split(body, ASCII_DASH_SEPARATOR; limit = 2)
    else
        error("unsupported registry heading separator: $(repr(line))")
    end
    length(parts) == 2 || error("malformed registry heading: $(repr(line))")
    id, title = String.(strip.(parts))
    occursin(r"^HP-[A-Z0-9]+(?:-[A-Z0-9]+)*$", id) ||
        error("malformed registry ID $(repr(id))")
    isempty(title) && error("empty registry title for $id")
    return id, title
end

function _paragraphs(lines)
    paragraphs = String[]
    current = String[]
    function finish!()
        isempty(current) && return
        push!(paragraphs, replace(join(strip.(current), " "), r"\s+" => " "))
        empty!(current)
    end
    for line in _visible_markdown_lines(lines)
        if isnothing(line) || isempty(strip(line))
            finish!()
        else
            push!(current, String(line))
        end
    end
    finish!()
    return paragraphs
end

_ownership_paragraph(text::AbstractString) =
    any(prefix -> startswith(text, prefix), OWNERSHIP_PREFIXES)

_ownership_followup(text::AbstractString) =
    startswith(text, "- ") || startswith(text, "[")

function _document_links(texts)
    links = String[]
    pattern = r"\]\(<?([^)>]+\.md(?:#[^)>]+)?)>?\)"
    for text in texts, matched in eachmatch(pattern, text)
        push!(links, matched.captures[1])
    end
    return sort!(unique!(links))
end

function _ownership_fields(paragraphs, id)
    blocks = String[]
    index = 1
    while index <= length(paragraphs)
        paragraph = paragraphs[index]
        if _ownership_paragraph(paragraph)
            parts = String[paragraph]
            follower = index + 1
            while follower <= length(paragraphs) &&
                    _ownership_followup(paragraphs[follower])
                push!(parts, paragraphs[follower])
                follower += 1
            end
            push!(blocks, join(parts, "\n"))
            index = follower
        else
            index += 1
        end
    end
    isempty(blocks) && error("$id has no explicit Owner/Canonical block")
    links = _document_links(blocks)
    isempty(links) && error("$id has no Owner/Canonical Markdown document link")
    for link in links
        relative = first(split(link, '#'; limit = 2))
        startswith(relative, "http://") && error("$id has nonlocal owner document $link")
        startswith(relative, "https://") && error("$id has nonlocal owner document $link")
        isfile(normpath(joinpath(dirname(REGISTRY_PATH), relative))) ||
            error("$id owner document does not exist: $link")
    end
    return blocks, links
end

function _state_fields(paragraphs, id)
    states = Tuple{String,String}[]
    for paragraph in paragraphs
        for label in ("Lifecycle", "Status")
            prefix = label * ":"
            startswith(paragraph, prefix) || continue
            value = strip(paragraph[(length(prefix) + 1):end])
            isempty(value) && error("$id has an empty $label field")
            push!(states, (label, value))
        end
    end
    length(states) == 1 || error("$id must have exactly one Lifecycle/Status paragraph")
    return only(states)
end

function _permission_fields(paragraphs, id)
    values = String[]
    for paragraph in paragraphs
        eligible = startswith(paragraph, "Permission:") ||
            startswith(paragraph, "Lifecycle:") || startswith(paragraph, "Status:")
        eligible || continue
        ranges = findall("Permission:", paragraph)
        length(ranges) <= 1 || error("$id has multiple Permission fields in one paragraph")
        isempty(ranges) && continue
        range = only(ranges)
        value = strip(paragraph[(last(range) + 1):end])
        isempty(value) && error("$id has an empty Permission field")
        push!(values, value)
    end
    isempty(values) && error("$id has no explicit Permission field")
    return values
end

function _registry_records(text::AbstractString)
    lines = _normalized_lines(text)
    levels = _visible_heading_levels(lines)
    boundaries = sort!(collect(keys(levels)))
    starts = Int[]
    for index in boundaries
        levels[index] == 3 || continue
        startswith(lstrip(lines[index], ' '), "### HP-") && push!(starts, index)
    end
    records = Dict{String,Any}[]
    seen = Set{String}()
    for start in starts
        id, title = _parse_heading(lines[start])
        id in seen && error("duplicate registry ID $id")
        push!(seen, id)
        next_boundary = findfirst(index -> index > start, boundaries)
        stop = isnothing(next_boundary) ? length(lines) : boundaries[next_boundary] - 1
        section_lines = lines[start:stop]
        paragraphs = _paragraphs(section_lines[2:end])
        state_label, state_text = _state_fields(paragraphs, id)
        permission_texts = _permission_fields(paragraphs, id)
        ownership_texts, owner_document_links = _ownership_fields(paragraphs, id)
        push!(records, Dict{String,Any}(
            "id" => id,
            "heading_title" => title,
            "state_label" => state_label,
            "state_text" => state_text,
            "permission_texts" => permission_texts,
            "ownership_texts" => ownership_texts,
            "owner_document_links" => owner_document_links,
            "record_sha256" => _digest(_normalized_block(section_lines)),
        ))
    end
    sort!(records; by = record -> record["id"])
    return records, _digest(_normalized_block(lines))
end

_registry_records() = _registry_records(read(REGISTRY_PATH, String))

function _standalone_marker(text, range)
    before = first(range) == firstindex(text) || text[prevind(text, first(range))] == '\n'
    after_index = nextind(text, last(range))
    after = after_index > lastindex(text) || text[after_index] in ('\r', '\n')
    return before && after
end

function _raw_marked_whitelist_block(text)
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

function _active_whitelist_marker_indices(lines)
    starts = Int[]
    stops = Int[]
    fence = nothing
    in_comment = false
    for (index, line) in pairs(lines)
        stripped = strip(line)
        if !isnothing(fence)
            _fence_marker(line) == fence && (fence = nothing)
            continue
        end
        if !in_comment && stripped in (WHITELIST_BEGIN, WHITELIST_END)
            push!(stripped == WHITELIST_BEGIN ? starts : stops, index)
            continue
        end
        uncommented, in_comment = _without_html_comments(line, in_comment)
        marker = _fence_marker(uncommented)
        isnothing(marker) || (fence = marker)
    end
    isnothing(fence) || error("unterminated Markdown fence")
    in_comment && error("unterminated Markdown HTML comment")
    length(starts) == 1 || error("expected one visible AGENTS whitelist begin marker")
    length(stops) == 1 || error("expected one visible AGENTS whitelist end marker")
    only(starts) < only(stops) || error("visible AGENTS whitelist markers are reversed")
    return only(starts), only(stops)
end

function _marked_whitelist_lines(text)
    _active_whitelist_marker_indices(_normalized_lines(text))

    raw_block = _raw_marked_whitelist_block(text)
    lines = _normalized_lines(raw_block)
    starts = findall(line -> line == WHITELIST_BEGIN, lines)
    stops = findall(line -> line == WHITELIST_END, lines)
    start, stop = only(starts), only(stops)
    return lines[(start + 1):(stop - 1)], _digest(raw_block)
end

function _parse_agents_whitelist(text)
    block_lines, block_digest = _marked_whitelist_lines(text)
    visible = _visible_markdown_lines(block_lines)
    start_phrase = "Cartesian Hamiltonian producer source work is currently authorized only for"
    starts = findall(line -> !isnothing(line) && strip(line) == start_phrase, visible)
    length(starts) == 1 || error("expected one AGENTS whitelist start marker")
    start = only(starts)
    ids = String[]
    for line in visible[(start + 1):end]
        isnothing(line) && continue
        stripped = strip(line)
        if startswith(stripped, "- ")
            matched = match(r"^- `(HP-[A-Z0-9]+(?:-[A-Z0-9]+)*)`$", stripped)
            isnothing(matched) &&
                error("malformed AGENTS whitelist item: $(repr(line))")
            push!(ids, matched.captures[1])
        elseif startswith(stripped, "* ") || startswith(stripped, "+ ")
            error("unsupported AGENTS whitelist bullet marker: $(repr(line))")
        elseif occursin("HP-", stripped)
            error("unsupported AGENTS whitelist syntax: $(repr(line))")
        end
    end
    length(ids) == length(unique(ids)) || error("duplicate AGENTS whitelist ID")
    sorted = sort(ids)
    digest = _digest(join(sorted, "\n") * "\n")
    return sorted, digest, block_digest
end

_agents_whitelist() =
    _parse_agents_whitelist(read(AGENTS_PATH, String))

function expected_shadow()
    records, registry_digest = _registry_records()
    whitelist, whitelist_digest, whitelist_block_digest = _agents_whitelist()
    registry_ids = Set(record["id"] for record in records)
    missing = sort!(collect(setdiff(Set(whitelist), registry_ids)))
    isempty(missing) || error("AGENTS whitelist IDs missing from registry: $(join(missing, ", "))")
    whitelist_set = Set(whitelist)
    for record in records
        record["agents_whitelisted"] = record["id"] in whitelist_set
    end
    return Dict{String,Any}(
        "schema_version" => 3,
        "artifact_kind" => ARTIFACT_KIND,
        "authoritative" => false,
        "authorization_complete" => false,
        "generated" => true,
        "scope" => SCOPE,
        "registry_source" => REGISTRY_SOURCE,
        "agents_source" => AGENTS_SOURCE,
        "registry_sha256" => registry_digest,
        "registry_file_sha256" => _digest(read(REGISTRY_PATH, String)),
        "agents_file_sha256" => _digest(read(AGENTS_PATH, String)),
        "agents_whitelist_sha256" => whitelist_digest,
        "agents_whitelist_block_sha256" => whitelist_block_digest,
        "registry_record_count" => length(records),
        "agents_whitelist_count" => length(whitelist),
        "records" => records,
    )
end

function _serialized(data)
    io = IOBuffer()
    TOML.print(io, data; sorted = true)
    text = String(take!(io))
    return endswith(text, '\n') ? text : text * "\n"
end

function _differences(expected, actual; path = "shadow")
    errors = String[]
    if expected isa AbstractDict && actual isa AbstractDict
        expected_keys = Set(string.(keys(expected)))
        actual_keys = Set(string.(keys(actual)))
        for key in sort!(collect(setdiff(expected_keys, actual_keys)))
            push!(errors, "$path missing key $key")
        end
        for key in sort!(collect(setdiff(actual_keys, expected_keys)))
            push!(errors, "$path has extra key $key")
        end
        for key in sort!(collect(intersect(expected_keys, actual_keys)))
            append!(errors, _differences(expected[key], actual[key]; path = "$path.$key"))
        end
    elseif expected isa AbstractVector && actual isa AbstractVector
        length(expected) == length(actual) ||
            push!(errors, "$path length $(length(actual)) != $(length(expected))")
        for index in 1:min(length(expected), length(actual))
            append!(errors, _differences(expected[index], actual[index]; path = "$path[$index]"))
        end
    elseif !(typeof(expected) === typeof(actual) && isequal(expected, actual))
        push!(errors, "$path mismatch: $(repr(actual)) != $(repr(expected))")
    end
    return errors
end

function check_shadow()
    isfile(SHADOW_PATH) || error("missing generated shadow file: $SHADOW_PATH")
    expected = expected_shadow()
    actual = TOML.parsefile(SHADOW_PATH)
    errors = _differences(expected, actual)
    committed = read(SHADOW_PATH, String)
    canonical = _serialized(actual)
    committed == canonical || push!(errors, "shadow TOML serialization is not canonical")
    isempty(errors) || error("Cartesian authority shadow drift:\n" * join(first(errors, 50), "\n"))
    return nothing
end

function write_shadow()
    temporary, io = mktemp(dirname(SHADOW_PATH))
    try
        write(io, _serialized(expected_shadow()))
        flush(io)
        close(io)
        mv(temporary, SHADOW_PATH; force = true)
    catch
        isopen(io) && close(io)
        ispath(temporary) && rm(temporary; force = true)
        rethrow()
    end
    return SHADOW_PATH
end

function _expect_parse_failure(text, needle)
    try
        _parse_agents_whitelist(text)
        error("negative whitelist check unexpectedly passed")
    catch error_value
        message = sprint(showerror, error_value)
        occursin(needle, message) || rethrow()
    end
end

function self_test()
    text = read(AGENTS_PATH, String)
    ids, _, block_digest = _parse_agents_whitelist(text)
    _expect_parse_failure(replace(text, WHITELIST_BEGIN => ""), "begin marker")
    _expect_parse_failure(text * "\n" * WHITELIST_BEGIN, "begin marker")
    reversed = replace(
        replace(text, WHITELIST_BEGIN => "TEMP-WHITELIST-MARKER"),
        WHITELIST_END => WHITELIST_BEGIN,
    )
    reversed = replace(reversed, "TEMP-WHITELIST-MARKER" => WHITELIST_END)
    _expect_parse_failure(reversed, "reversed")
    malformed = replace(text, "- `HP-OBJ-01`" => "- HP-OBJ-01")
    _expect_parse_failure(malformed, "malformed AGENTS whitelist item")
    whitespace_changed = replace(text, "- `HP-OBJ-01`" => "- `HP-OBJ-01`  ")
    changed_ids, _, changed_block_digest = _parse_agents_whitelist(whitespace_changed)
    ids == changed_ids || error("whitelist whitespace mutation changed parsed IDs")
    block_digest != changed_block_digest || error("raw whitelist block digest ignored byte drift")
    fenced = replace(text, WHITELIST_BEGIN => "```text\n" * WHITELIST_BEGIN)
    fenced = replace(fenced, WHITELIST_END => WHITELIST_END * "\n```")
    _expect_parse_failure(fenced, "visible AGENTS whitelist begin marker")
    commented = replace(text, WHITELIST_BEGIN => "<!--\n" * WHITELIST_BEGIN)
    commented = replace(commented, WHITELIST_END => WHITELIST_END * "\n-->")
    _expect_parse_failure(commented, "visible AGENTS whitelist begin marker")
    return nothing
end

function main(args = ARGS)
    if isempty(args) || args == ["--check"]
        check_shadow()
    elseif args == ["--write"]
        println(write_shadow())
    elseif args == ["--self-test"]
        self_test()
    else
        error("usage: julia --project=docs docs/check_cartesian_authority_shadow.jl [--check|--write|--self-test]")
    end
end

end

if abspath(PROGRAM_FILE) == @__FILE__
    CartesianAuthorityShadow.main()
end
