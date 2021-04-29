PRIMER_BED = params.refdir+"/nCoV-2019.primer.bed"
genome_fasta = params.refdir+"/nCoV-2019.reference.fasta"
GFF = params.refdir+"/MN908947.3.gff"
MIN_DEPTH=10
MIN_FREQ=0.75
MIN_QUAL_VAR=20

Channel
    .fromPath(params.csv).splitCsv(header:true)
    .map{ row-> tuple(row.id, file(row.fq1), file(row.fq1)) }
    .set { input_data }

process index_references {
    cpus 1
    memory '1 GB'
    time '1h'

    when:
        params.build_indexes
    
    input:
        file(genome_fasta)

    """
    bwa index $genome_fasta
    """
    
}

process subsample_fastq {
    cpus 2
    memory '1 GB'
    time '1h'

    input:
    set id, fq1, fq2 from input_data

    output:
    set id, file("${id}_subsample_R1_001.fastq.gz"), file("${id}_subsample_R2_001.fastq.gz") into subsampled_data

    """
    seqtk sample -s 1314 $fq1 $params.max_readpairs | gzip -c > ${id}_subsample_R1_001.fastq.gz &
    seqtk sample -s 1314 $fq2 $params.max_readpairs | gzip -c > ${id}_subsample_R2_001.fastq.gz &
    wait
    """
}

process map_and_trim {
    cpus 4
    memory '2 GB'
    time '1h'

    input:
        set id, fq1, fq2 from subsampled_data

    output:
        set id, file("${id}.trim.sort.bam"), file("${id}.trim.sort.bam.bai") into bam_depth, bam_consensus, bam_variantcalling, bam_qc
    
    """
    bwa mem -t 4 $genome_fasta $fq1 $fq2 | samtools sort | tee ${id}.qc.bam | samtools view -F 4 -o ${id}.sort.bam
    samtools flagstat ${id}.qc.bam > ${id}.flagstat
    samtools index ${id}.sort.bam
    ivar trim -e -i ${id}.sort.bam -b $PRIMER_BED -p ${id}.trim -m 30 -q 20
    samtools sort ${id}.trim.bam -o ${id}.trim.sort.bam
    samtools index ${id}.trim.sort.bam
    rm ${id}.sort.bam ${id}.sort.bam.bai ${id}.trim.bam ${id}.qc.bam
    """
}


process depth {
    cpus 1
    memory '1 GB'
    time '1h'

    input:
        set id, file(bam), file(bai) from bam_depth

    output:
        set id, file("${id}.depth") into depth_data

    """
    sambamba depth base -c0 $bam -o ${id}.depth
    """
}


process consensus {
    cpus 1
    memory '1 GB'
    time '1h'

    input:
        set id, file(bam), file(bai) from bam_consensus

    output:
        set id, file("${id}.consensus.fa"), file("${id}.consensus.qual.txt") into consensus_qc, consensus_nextclade, consensus_pangolin

    """
    samtools mpileup -aa -A -B -d 6000000 -Q 0 --reference $genome_fasta $bam | ivar consensus -p ${id}.consensus -n N -m $MIN_DEPTH -t $MIN_FREQ 
    """
}

process variant_calling {
    cpus 1
    memory '1 GB'
    time '1h'

    input:
        set id, file(bam), file(bai) from bam_variantcalling

    output:
        set id, file("${id}.freebayes.vep.vcf") into vcf

    """
    freebayes -p 1 --min-coverage $MIN_DEPTH --min-base-quality $MIN_QUAL_VAR -f $genome_fasta -F $MIN_FREQ -m 60 $bam > ${id}.freebayes.raw.vcf
    bcftools norm ${id}.freebayes.raw.vcf -f $genome_fasta -o ${id}.freebayes.vcf
    rm ${id}.freebayes.raw.vcf

    bgzip -i -f -c ${id}.freebayes.vcf > ${id}.freebayes.vcf.gz
    source activate vep
    vep -i ${id}.freebayes.vcf.gz --format vcf --gff ${GFF}.gz --fasta $genome_fasta -o ${id}.freebayes.vep.vcf --vcf --force_overwrite --no_stats --distance 10 --hgvs
    rm ${id}.freebayes.vcf.gz*
    """
}


process qc {
    cpus 1
    memory '1 GB'

    input:
        set id, file(fasta), file(qual), file(bam), file(bai) from consensus_qc.join(bam_qc)
    
    output:
        set id, file("${id}.qc.csv")

    """
    #python qc.py --illumina --outfile ${id}.qc.csv --sample ${id} --ref $genome_fasta --bam $bam --fasta $fasta
    echo "dsa" > ${id}.qc.csv
    """
}


process nextclade {
    cpus 1
    memory '1 GB'

    input:
	set id, file(fasta), file(qual) from consensus_nextclade

    output:
	set id, file("${id}.nextclade.tsv"), file("${id}.auspice.json")

   """
   nextclade --input-fasta $fasta --output-tsv ${id}.nextclade.tsv --output-tree ${id}.auspice.json
   """
}    


process pangolin {
    cpus 1
    memory '1 GB'

    input:
	file(fastas) from consensus_pangolin.collect {it[1]}

    output:
	file("pangolin_all.csv") 

    """
    source activate pangolin
    cat ${fastas.join(" ")} > all.fasta
    pangolin all.fasta -o all_pangolin
    cp all_pangolin/lineage_report.csv pangolin_all.csv
    """
}    
