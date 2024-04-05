function isforwardstrand(record::BAM.Record)
    if BAM.flag(record) & SAM.FLAG_READ1 ≠ 0               # if it's read1
        if BAM.flag(record) & SAM.FLAG_REVERSE ≠ 0         # and if it's reverse complement
            return true
        end
    elseif BAM.flag(record) & SAM.FLAG_READ2 ≠ 0           # if it's read2
        if BAM.flag(record) & SAM.FLAG_REVERSE == 0        # and it's forward
            return true
        end
    elseif BAM.flag(record) & SAM.FLAG_REVERSE ≠ 0        # if it's single-end & reverse complement
            return true
    end
    return false
end

struct SequenceMatch
    ref_index::Int
    strand::Bool
    refpos::Int
    readpos::Int
    match::DNA
end

const blank = SequenceMatch(0,true,0,0,DNA_Gap)

function countbases(refseqs, next_record::BAM.Record, baseQ_threshold::Int, mismatch_threshold::Int, contextwindow::Int, forwardorientation::Bool, utoc::Bool, mismatches::IO)
    refindex = BAM.refid(next_record)
    aln = BAM.alignment(next_record)
    seq = BAM.sequence(next_record)
    qual = BAM.quality(next_record)
    forward = (isforwardstrand(next_record) && !forwardorientation) || (!isforwardstrand(next_record) && forwardorientation)
    potential_matches = fill(blank, length(seq))
    new_mismatches = SequenceMatch[]
    mismatch_mask = trues(length(seq))

    for pos in eachindex(seq)
        qual[pos] < baseQ_threshold && continue
        t1nt = seq[pos]
        t1nt == DNA_N && continue
        local refpos, op
        try
            refpos, op = seq2ref(aln, pos)
        catch
            return empty!(potential_matches)
        end
        if refpos < 1 || refpos > length(refseqs[refindex])
            #println(next_record)
            return empty!(potential_matches)
        end
        refnt = refseqs[refindex][refpos]
        if ismatchop(op) && (refnt == t1nt || (refnt == DNA_C && forward && t1nt == DNA_T)
                            || (refnt == DNA_G && !forward && t1nt == DNA_A)
                            || (refnt == DNA_T && forward && utoc && t1nt == DNA_C)
                            || (refnt == DNA_A && !forward && utoc && t1nt == DNA_G))
            # don't add potential editing sites at the end of reads
            if refnt ≠ t1nt && refpos == BAM.rightposition(next_record)
                # do nothing
            else
                potential_matches[pos] = SequenceMatch(refindex, forward, refpos, pos, t1nt)
            end
        elseif ismetaop(op)
            continue
        else
            push!(new_mismatches, SequenceMatch(refindex, forward, refpos, pos, t1nt))
            length(new_mismatches) > mismatch_threshold && return empty!(potential_matches)
        end
    end
    #filter matches around mismatches
    for mm in new_mismatches
        mismatch_mask[max(1, mm.readpos - contextwindow):min(length(mismatch_mask), mm.readpos + contextwindow)] .= false
        write(mismatches, join([refindex,mm.strand,string(mm.refpos),mm.match], '\t'), '\n')
    end
    
    return filter(m->m.readpos > 0 && mismatch_mask[m.readpos], potential_matches)
end



function scanbam(refseqs::Array{LongDNA},reader::BAM.Reader,mapQ_threshold::Int,baseQ_threshold::Int,mismatch_threshold::Int,contextwindow::Int,forwardorientation::Bool,utoc::Bool,out::String)

    fwd_base_counts = Array{Array{Int}}(undef,length(refseqs))
    rev_base_counts = Array{Array{Int}}(undef,length(refseqs))

    for (index,ref) in enumerate(refseqs)
        fwd_base_counts[index] = zeros(Int64, length(ref), 4)
        rev_base_counts[index] = zeros(Int64, length(ref), 4)
    end

    fwd_bases = Dict(DNA_A=>1, DNA_C=>2, DNA_G=>3, DNA_T=>4)
    rev_bases = Dict(DNA_A=>4, DNA_C=>3, DNA_G=>2, DNA_T=>1)

    mismatches = open(out * ".mismatches", "w")

    threads = 1
    records = Vector{BAM.Record}(undef, threads)

    countslock = ReentrantLock()

    while !eof(reader)
        n = 0
        for t in 1:threads
            records[t] = first(reader)
            n += 1
            eof(reader) && break
        end
        Threads.@threads for r in 1:n
            next_record = records[r]
            !BAM.isfilled(next_record) && continue
            # ignore bad read pairs (not properly mapped, duplicates etc)
            if !passesmuster(BAM.flag(next_record)) || BAM.mappingquality(next_record) < mapQ_threshold
                continue
            end
            matches = countbases(refseqs, next_record, baseQ_threshold, mismatch_threshold,contextwindow, forwardorientation, utoc, mismatches)
            lock(countslock) do
                for m in matches
                    counts_array = m.strand ? fwd_base_counts : rev_base_counts
                    base_index = get(m.strand ? fwd_bases : rev_bases, m.match, 0)
                    counts_array[m.ref_index][m.refpos, base_index] += 1
                end
            end
        end
    end
    close(mismatches)
    return fwd_base_counts, rev_base_counts
end

function passesmuster(flags)
    flags & SAM.FLAG_SECONDARY ≠ 0 && return false
    flags & SAM.FLAG_QCFAIL ≠ 0 && return false
    flags & SAM.FLAG_DUP ≠ 0 && return false
    flags & SAM.FLAG_SUPPLEMENTARY ≠ 0 && return false
    return true
end
