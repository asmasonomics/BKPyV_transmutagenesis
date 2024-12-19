configfile: "config.yaml"

rule all:
    input:
        qc="03_qc/multiqc_report.html",
        analysis=expand("10_analysis/{sample}_{duplex_type}-vs-{undiluted_type}/", sample=config["SAMPLES"], duplex_type=config["DUPLEX_TYPES"], undiluted_type=config["UNDILUTED_TYPES"]),
        contamination=expand("11_contamination_check/{sample}_{type}.selfSM", sample=config["SAMPLES"], type=config["TYPES"]),
        efficiency=expand("12_efficiency_estimate/RBs/{sample}_{type}.RBs", sample=config["SAMPLES"], type=config["TYPES"]),
        wg_metrics="12_efficiency_estimate/efficiency-metrics-whole-genome.tsv",
        vcfs=expand("13_vcf/{sample}_{duplex_type}-vs-{undiluted_type}.SNV.vcf.gz", sample=config["SAMPLES"], duplex_type=config["DUPLEX_TYPES"], undiluted_type=config["UNDILUTED_TYPES"])

rule fastqc:
    input:
        "00_raw/{sample}_{type}_{readnum}.fastq.gz",
    output:
        directory("03_qc/{sample}_{type}_{readnum}_fastqc")
    params:
        threads=8
    resources:
        runtime=60,
        mem_mb=2000,
        cpus_per_task=9
    shell:
        r"""
        module purge
        module load FastQC/0.11.9-Java-11

        mkdir {output}
        fastqc -t {params.threads} -o {output} {input}
        """

rule multiqc:
    input:
        expand("03_qc/{sample}_{type}_{readnum}_fastqc", sample=config["SAMPLES"], type=config["TYPES"], readnum=["read1", "read2"])
    output:
        "03_qc/multiqc_report.html"
    resources:
        runtime=30,
        cpus_per_task=9
    shell:
        r"""
        module purge
        module load MultiQC/1.13-intel-2021b

        multiqc 03_qc/*_fastqc/ -o 03_qc/
        """

rule extract_tags:
    input:
        fq1="00_raw/{sample}_{type}_read1.fastq.gz",
        fq2="00_raw/{sample}_{type}_read2.fastq.gz"
    output:
        fq1=temp("04_fastq_extract_tags/{sample}_{type}_read1.extract_tags.fastq.gz"),
        fq2=temp("04_fastq_extract_tags/{sample}_{type}_read2.extract_tags.fastq.gz")
    params:
        read_length=config["READ_LENGTH"],
        library_type=config["LIBRARY_TYPE"]
    resources:
        runtime=300,
        mem_mb=1000,
        cpus_per_task=1
    run:
        if params.library_type == "HpyCH4V":
            shell("module purge; module load NanoSeq/3.2.1-foss-2020b-R-4.0.3; extract_tags.py -a {input.fq1} -b {input.fq2} -c {output.fq1} -d {output.fq2} -m 3 -s 4 -l {params.read_length}")
        elif params.library_type == "sonicated":
            shell("module purge; module load NanoSeq/3.2.1-foss-2020b-R-4.0.3; extract_tags.py -a {input.fq1} -b {input.fq2} -c {output.fq1} -d {output.fq2} -m 3 -s 2 -l {params.read_length}")
        else:
            print("Specify library type as sonicated or HpyCH4V in config.yaml")

rule bwa_index:
    input:
        "01_ref/{genome}"
    output:
        amb="01_ref/{genome}.amb",
        ann="01_ref/{genome}.ann",
        bwt="01_ref/{genome}.bwt",
        pac="01_ref/{genome}.pac",
        sa="01_ref/{genome}.sa",
    resources:
        runtime=360,
        mem_mb=32000,
        cpus_per_task=8
    shell:
        r"""
        module purge
        module load BWA/0.7.17-GCCcore-11.2.0

        bwa index {input}
        """

rule bwa_mem:
    input:
        fq1="04_fastq_extract_tags/{sample}_{type}_read1.extract_tags.fastq.gz",
        fq2="04_fastq_extract_tags/{sample}_{type}_read2.extract_tags.fastq.gz",
        fa=expand("01_ref/{genome}", genome=config["GENOME"]),
        amb=expand("01_ref/{genome}.amb", genome=config["GENOME"]),
        ann=expand("01_ref/{genome}.ann", genome=config["GENOME"]),
        bwt=expand("01_ref/{genome}.bwt", genome=config["GENOME"]),
        pac=expand("01_ref/{genome}.pac", genome=config["GENOME"]),
        sa=expand("01_ref/{genome}.sa", genome=config["GENOME"])
    output:
        temp("05_sam/{sample}_{type}.sam")
    params:
        threads=8
    resources:
        runtime=2880, # 48h
        mem_mb=48000,
        cpus_per_task=8
    shell:
        r"""
        module purge
        module load BWA/0.7.17-GCCcore-11.2.0

        bwa mem -C -t {params.threads} \
          {input.fa} \
          {input.fq1} \
          {input.fq2} > {output}
        """

rule sort_and_add_rc_mc_tags:
    input:
        "05_sam/{sample}_{type}.sam"
    output:
        temp("06_bam_rc_mc_tags/{sample}_{type}.rc_mc_tags.bam")
    params:
        threads=4
    resources:
        runtime=540,
        mem_mb=24000,
        cpus_per_task=4
    shell:
        r"""
        module purge
        module load NanoSeq/3.2.1-foss-2020b-R-4.0.3

        bamsormadup inputformat=sam rcsupport=1 threads={params.threads} < {input} > {output}
        """

rule mark_dups:
    input:
        "06_bam_rc_mc_tags/{sample}_{type}.rc_mc_tags.bam"
    output:
        bam=temp("07_bam_mark_dups/{sample}_{type}.mark_dups.bam"),
        bai=temp("07_bam_mark_dups/{sample}_{type}.mark_dups.bam.bai")
    params:
        threads=4
    resources:
        runtime=540,
        mem_mb=3000,
        cpus_per_task=4
    shell:
        r"""
        module purge
        module load NanoSeq/3.2.1-foss-2020b-R-4.0.3

        bammarkduplicatesopt inputthreads={params.threads} optminpixeldif=2500 I={input} O={output.bam}

        module purge
        module load SAMtools/1.16.1-GCC-11.3.0

        samtools index -o {output.bai} {output.bam}
        """

rule filter_and_append_rb_tags:
    input:
        "07_bam_mark_dups/{sample}_{type}.mark_dups.bam"
    output:
        bam="08_bam_rb_tags/{sample}_{type}.rb_tags.bam",
        bai="08_bam_rb_tags/{sample}_{type}.rb_tags.bam.bai"
    resources:
        runtime=480,
        mem_mb=1000,
        cpus_per_task=1
    shell:
        r"""
        module purge
        module load NanoSeq/3.2.1-foss-2020b-R-4.0.3

        bamaddreadbundles -I {input} -O {output.bam}

        module purge
        module load SAMtools/1.16.1-GCC-11.3.0

        samtools index -o {output.bai} {output.bam}
        """

rule keep_random_read:
    input:
        "08_bam_rb_tags/{sample}_{type}.rb_tags.bam"
    output:
        bam="09_bam_random_read/{sample}_{type}.random_read.bam",
        bai="09_bam_random_read/{sample}_{type}.random_read.bam.bai"
    resources:
        runtime=420,
        mem_mb=4000,
        cpus_per_task=1
    shell:
        r"""
        module purge
        module load NanoSeq/3.2.1-foss-2020b-R-4.0.3

        randomreadinbundle -I {input} -O {output.bam}

        module purge
        module load SAMtools/1.16.1-GCC-11.3.0

        samtools index -o {output.bai} {output.bam}
        """

rule run_analysis:
    input:
        bam_undiluted="09_bam_random_read/{sample}_{undiluted_type}.random_read.bam",
        bam_duplex="08_bam_rb_tags/{sample}_{duplex_type}.rb_tags.bam",
        genome=expand("01_ref/{genome}", genome=config["GENOME"]),
        mask_SNP=expand("01_ref/SNP.sorted.{build}.bed.gz", build=config["GENOME_BUILD"]),
        mask_NOISE=expand("01_ref/NOISE.sorted.{build}.bed.gz", build=config["GENOME_BUILD"]),
    output:
        dir=directory("10_analysis/{sample}_{duplex_type}-vs-{undiluted_type}"),
        muts="10_analysis/{sample}_{duplex_type}-vs-{undiluted_type}/tmpNanoSeq/post/results.muts.vcf.gz",
        indel="10_analysis/{sample}_{duplex_type}-vs-{undiluted_type}/tmpNanoSeq/post/results.indel.vcf.gz"
    resources:
        runtime=840,
        mem_mb=96000,
        cpus_per_task=21
    shell:
        r"""
        module purge
        module load NanoSeq/3.2.1-foss-2020b-R-4.0.3

        cd {output.dir}

        # cov
        runNanoSeq.py -t 10 \
          -A ../../{input.bam_undiluted} \
          -B ../../{input.bam_duplex} \
          -R ../../{input.genome} \
          cov \
          -Q 0
        # part
        runNanoSeq.py -t 1 \
          -A ../../{input.bam_undiluted} \
          -B ../../{input.bam_duplex} \
          -R ../../{input.genome} \
          part \
          -n 20
        # dsa
        runNanoSeq.py -t 20 \
          -A ../../{input.bam_undiluted} \
          -B ../../{input.bam_duplex} \
          -R ../../{input.genome} \
          dsa \
          -C ../../{input.mask_SNP} \
          -D ../../{input.mask_NOISE} \
          -d 2 \
          -q 30
        # var
        runNanoSeq.py -t 20 \
          -A ../../{input.bam_undiluted} \
          -B ../../{input.bam_duplex} \
          -R ../../{input.genome} \
          var \
          -a 50 \
          -b 0 \
          -c 0.02 \
          -d 2 \
          -f 0.9 \
          -i 1 \
          -m 8 \
          -n 3 \
          -p 0 \
          -q 60 \
          -r 144 \
          -v 0.01 \
          -x 8 \
          -z 15
        # indel
        runNanoSeq.py -t 20 \
          -A ../../{input.bam_undiluted} \
          -B ../../{input.bam_duplex} \
          -R ../../{input.genome} \
          indel \
          -s {wildcards.sample}_{wildcards.duplex_type}-vs-{wildcards.undiluted_type} \
          --rb 2 \
          --t3 136 \
          --t5 8 \
          -a 50 \
          -c 0.02 \
          -z 15 \
          -v 0.01
        # post
        runNanoSeq.py -t 2 \
          -A ../../{input.bam_undiluted} \
          -B ../../{input.bam_duplex} \
          -R ../../{input.genome} \
          post
        """

rule check_contamination:
    input:
        bam="09_bam_random_read/{sample}_{type}.random_read.bam",
        hla_fa=expand("01_ref/{ref}", ref=config["1000G_REFERENCE_GENOME"])
    output:
        selfsm="11_contamination_check/{sample}_{type}.selfSM"
    resources:
        runtime=60,
        mem_mb=2000,
        cpus_per_task=4
    shell:
        r"""
        module purge
        module load VerifyBamID2/2.0.1-GCC-11.3.0

        mkdir -p 11_contamination_check
        VerifyBamID --SVDPrefix /opt/apps/eb/software/VerifyBamID2/2.0.1-GCC-11.3.0/resource/1000g.phase3.100k.b38.vcf.gz.dat --Reference {input.hla_fa} --BamFile {input.bam} --Output 11_contamination_check/{wildcards.sample}_{wildcards.type}
        """

rule estimate_efficiency:
    input:
        bam_rbtags="08_bam_rb_tags/{sample}_{type}.rb_tags.bam",
        bam_randomread="09_bam_random_read/{sample}_{type}.random_read.bam",
        genome=expand("01_ref/{genome}", genome=config["GENOME"])
    output:
        rbs="12_efficiency_estimate/RBs/{sample}_{type}.RBs",
        gc_inserts="12_efficiency_estimate/RBs.GC_inserts/{sample}_{type}.RBs.GC_inserts.tsv",
        pdf="12_efficiency_estimate/RBs.pdf/{sample}_{type}.RBs.pdf",
        tsv="12_efficiency_estimate/{sample}_{type}.tsv"
    params:
        threads=2
    resources:
        runtime=60,
        mem_mb=20000,
        cpus_per_task=4
    shell:
        r"""
        module purge
        module load NanoSeq/3.2.1-foss-2020b-R-4.0.3

        mkdir -p 12_efficiency_estimate
        cd 12_efficiency_estimate
        efficiency_nanoseq.pl -t {params.threads} -d ../{input.bam_randomread} -x ../{input.bam_rbtags} -r ../{input.genome} -o {wildcards.sample}_{wildcards.type}

        mkdir -p RBs/ RBs.GC_inserts/ RBs.pdf/
        mv {wildcards.sample}_{wildcards.type}.RBs RBs/
        mv {wildcards.sample}_{wildcards.type}.RBs.GC_inserts.tsv RBs.GC_inserts/
        mv {wildcards.sample}_{wildcards.type}.RBs.pdf RBs.pdf/
        """

rule collate_efficiency:
    input:
        expand("12_efficiency_estimate/{sample}_{type}.tsv", sample=config["SAMPLES"], type=config["TYPES"])
    output:
        wg_metrics="12_efficiency_estimate/efficiency-metrics-whole-genome.tsv",
        rb_metrics="12_efficiency_estimate/read-bundle-metrics-chr1-only.tsv"
    resources:
        runtime=15,
        mem_mb=500,
        cpus_per_task=1
    shell:
        r"""
        module purge
        module load NanoSeq/3.2.1-foss-2020b-R-4.0.3

        Rscript collate-efficiency-metrics.R
        """

rule organise_outputs:
    input:
        vcf_snv="10_analysis/{sample}_{duplex_type}-vs-{undiluted_type}/tmpNanoSeq/post/results.muts.vcf.gz",
        vcf_indel="10_analysis/{sample}_{duplex_type}-vs-{undiluted_type}/tmpNanoSeq/post/results.indel.vcf.gz"
    output:
        vcf_snv="13_vcf/{sample}_{duplex_type}-vs-{undiluted_type}.SNV.vcf.gz",
        vcf_indel="13_vcf/{sample}_{duplex_type}-vs-{undiluted_type}.indel.vcf.gz"
    resources:
        runtime=15,
        mem_mb=500,
        cpus_per_task=1
    shell:
        r"""
        mkdir -p 13_vcf
        cp {input.vcf_snv} {output.vcf_snv}
        cp {input.vcf_indel} {output.vcf_indel}
        """
