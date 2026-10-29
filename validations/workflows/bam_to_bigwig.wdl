version 1.0

workflow bam_to_bigwig {
    input {
        File  bam
        File? bai
        File  chrom_sizes
    }

    if (!defined(bai)) {
        call index_bam { input: bam = bam }
    }

    File bai_file = select_first([bai, index_bam.bai])

    call count_mapped_reads { input: bam = bam, bai = bai_file }

    call bam_to_bedgraph {
        input:
            bam          = bam,
            bai          = bai_file,
            mapped_reads = count_mapped_reads.mapped_reads
    }

    call bedgraph_to_bigwig {
        input:
            bedgraph    = bam_to_bedgraph.bedgraph,
            chrom_sizes = chrom_sizes
    }

    output {
        File bigwig = bedgraph_to_bigwig.bigwig
        File bai_out = bai_file
    }
}

task index_bam {
    input {
        File bam
        Int  cpu    = 4
        Int  memory = 8
    }

    command <<<
        set -euo pipefail
        samtools index -@ ~{cpu} ~{bam}
        # samtools writes <bam>.bai next to the input; copy to cwd for delocalization
        cp ~{bam}.bai ~{basename(bam)}.bai
    >>>

    output {
        File bai = "~{basename(bam)}.bai"
    }

    runtime {
        docker: "quay.io/biocontainers/samtools:1.17--h00cdaf9_0"
        cpu:    cpu
        memory: "~{memory} GiB"
        disks:  "local-disk 100 HDD"
    }
}

task count_mapped_reads {
    input {
        File bam
        File bai
        Int  cpu    = 4
        Int  memory = 8
    }

    command <<<
        samtools view -c -F 4 -@ ~{cpu} ~{bam}
    >>>

    output {
        Int mapped_reads = read_int(stdout())
    }

    runtime {
        docker: "quay.io/biocontainers/samtools:1.17--h00cdaf9_0"
        cpu:    cpu
        memory: "~{memory} GiB"
        disks:  "local-disk 50 HDD"
    }
}

task bam_to_bedgraph {
    input {
        File bam
        File bai
        Int  mapped_reads
        Int  cpu    = 4
        Int  memory = 16
    }

    String sample_name  = basename(bam, ".bam")
    Float  scale_factor = 1000000.0 / mapped_reads

    command <<<
        set -euo pipefail

        bedtools genomecov \
            -ibam ~{bam} \
            -bg \
            -scale ~{scale_factor} \
            -pc \
        | sort -k1,1 -k2,2n \
        > ~{sample_name}.bedGraph
    >>>

    output {
        File bedgraph = "~{sample_name}.bedGraph"
    }

    runtime {
        docker: "quay.io/biocontainers/bedtools:2.31.0--hf5e1c6e_2"
        cpu:    cpu
        memory: "~{memory} GiB"
        disks:  "local-disk 200 HDD"
    }
}

task bedgraph_to_bigwig {
    input {
        File bedgraph
        File chrom_sizes
        Int  cpu    = 2
        Int  memory = 8
    }

    String sample_name = basename(bedgraph, ".bedGraph")

    command <<<
        set -euo pipefail

        bedGraphToBigWig \
            ~{bedgraph} \
            ~{chrom_sizes} \
            ~{sample_name}.bw
    >>>

    output {
        File bigwig = "~{sample_name}.bw"
    }

    runtime {
        docker: "quay.io/biocontainers/ucsc-bedgraphtobigwig:377--h0b8a92a_2"
        cpu:    cpu
        memory: "~{memory} GiB"
        disks:  "local-disk 100 HDD"
    }
}
