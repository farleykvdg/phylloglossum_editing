using FASTX
using BioSequences

function main(protein_alignment_file::String, nucleotide_file::String, outfile::String)

    reader = open(FASTA.Reader, protein_alignment_file)
    proteins = Vector{FASTA.Record}(undef, 0)
    for record in reader
        push!(proteins, record)
    end
    close(reader)

    reader = open(FASTA.Reader, nucleotide_file)
    nts = Vector{FASTA.Record}(undef, 0)
    for record in reader
        push!(nts, record)
    end
    close(reader)

    aligned_nts = Vector{FASTA.Record}(undef, 0)
    gap_codon = LongDNA{4}([DNA_Gap, DNA_Gap, DNA_Gap])
    for (p, nt) in zip(proteins, nts)
        @assert FASTA.identifier(p) == FASTA.identifier(nt)
        aligned_ntseq = LongDNA{4}()
        ntseq = FASTA.sequence(LongDNA{4}, nt)
        pointer = 1
        for aa in FASTA.sequence(p)
            aa == AA_Term && break
            if aa == AA_Gap
                append!(aligned_ntseq, gap_codon)
            else
                append!(aligned_ntseq, ntseq[pointer:pointer+2])
                pointer += 3
            end
        end
        push!(aligned_nts, FASTA.Record(FASTA.identifier(nt), aligned_ntseq))
    end

    writer = open(FASTA.Writer, outfile)
    for record in aligned_nts
        write(writer, record)
    end
    close(writer)
end

main(ARGS[1], ARGS[2], ARGS[3])
