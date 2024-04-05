using ArgParse
using BioSequences
using BioAlignments
using XAM
using FASTX

include("countbases.jl")

function parse_commandline()
    s = ArgParseSettings()

    @add_arg_table! s begin
        "--mapQ"
            help = "mapping quality threshold"
            arg_type = Int
            default = 0
        "--baseQ"
            help = "base quality threshold"
            arg_type = Int
            default = 30
        "--mismatches", "-m"
            help = "mismatch threshold; the number of non-editing mismatches allowed before a read is rejected"
            arg_type = Int
            default = 0
        "--contextwindow", "-w"
            help = "window before and after the base to be counted that should not contain mismatches or indels"
            arg_type = Int
            default = 6
        "--forwardorientation", "-f"
            help = "flag indicating reads (or first read in a pair) are the forward RNA strand, unlike the usual Illumina libraries"
            action = :store_true
        "--UtoC", "-u"
            help = "flag indicating whether to count U to C mismatches as editing events"
            action = :store_true
        "--outfile","-o"
            help = "path to results file"
            arg_type = String
        "reference"
            help = "reference sequence(s) in fasta format"
            required = true
        "bam"
            help = "bam file; paired reads should be merged before mapping"
            required = true
    end

    return parse_args(s)
end

function main()
    parsed_args = parse_commandline()

    references = FASTA.Record[]
    open(FASTA.Reader, parsed_args["reference"]) do reader
        for record in reader
            push!(references,record)
        end
    end
    refseqs = Array{LongDNA}(undef,length(references))
    for (index,record) in enumerate(references)
        refseqs[index] = FASTA.sequence(LongDNA{4},record)
    end

    mapfile = parsed_args["bam"]
    reader = open(BAM.Reader, mapfile)

    mapQ_threshold = parsed_args["mapQ"]
    baseQ_threshold = parsed_args["baseQ"]
    mismatch_threshold = parsed_args["mismatches"]
    contextwindow = parsed_args["contextwindow"]
    utoc = parsed_args["UtoC"]

    out = parsed_args["outfile"]
    if isnothing(out)
        io = Base.stdout
    else
        io = open(out, "w")
    end

    fwd_base_counts, rev_base_counts = scanbam(refseqs,reader,mapQ_threshold,baseQ_threshold,mismatch_threshold,contextwindow,parsed_args["forwardorientation"], utoc,out)

    write(io,join(["refID","pos","ref","fA","fC","fG","fU","rA","rC","rG","rU","total"],"\t"),"\n")
    
    for (index,reference) in enumerate(references)
        for nt in 1:length(refseqs[index])
            write(io, join([FASTA.identifier(reference),nt,refseqs[index][nt],fwd_base_counts[index][nt,1],fwd_base_counts[index][nt,2],fwd_base_counts[index][nt,3],fwd_base_counts[index][nt,4],rev_base_counts[index][nt,1],rev_base_counts[index][nt,2],rev_base_counts[index][nt,3],rev_base_counts[index][nt,4],sum(fwd_base_counts[index][nt,:])+sum(rev_base_counts[index][nt,:])],"\t"),"\n")
        end
    end

    close(io)
    
    close(reader)

end

main()
