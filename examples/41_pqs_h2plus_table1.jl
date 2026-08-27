using GaussletBases

comparison = pqs_h2plus_comparison(
    independent_reference_total_energy_hartree = -0.6026342144949465,
)
fields = fieldnames(PQSH2PlusRow)
output_path = isempty(ARGS) ? "/tmp/pqs_h2plus_table1.tsv" : abspath(only(ARGS))
open(output_path, "w") do io
    println(io, join(string.(fields), '\t'))
    for row in comparison.rows
        println(io, join((repr(getfield(row, field)) for field in fields), '\t'))
    end
end

for row in comparison.rows
    println(row.basis, ": dimension=", row.dimension,
        ", capture=", row.capture,
        ", total energy=", row.total_energy_hartree,
        ", total error=", row.total_error_hartree)
end
println("TSV: ", output_path)

comparison
