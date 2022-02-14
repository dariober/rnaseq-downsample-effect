# SAMPLES = {
#     "ap2-o2_O_2-1": ["SRR3437941_1.fastq.gz", "SRR3437941_2.fastq.gz"],
#     "ap2-o2_O_2-2": ["SRR3437942_1.fastq.gz", "SRR3437942_2.fastq.gz"],
#     "ap2-o2_O_3-1": ["SRR3437957_1.fastq.gz", "SRR3437957_2.fastq.gz"],
#     "ap2-o2_O_3-2": ["SRR3437958_1.fastq.gz", "SRR3437958_2.fastq.gz"],
#     "ap2-o_O_1-1": ["SRR3437929_1.fastq.gz", "SRR3437929_2.fastq.gz"],
#     "ap2-o_O_1-2": ["SRR3437930_1.fastq.gz", "SRR3437930_2.fastq.gz"],
#     "ap2-o_O_2-1": ["SRR3437945_1.fastq.gz", "SRR3437945_2.fastq.gz"],
#     "ap2-o_O_2-2": ["SRR3437946_1.fastq.gz", "SRR3437946_2.fastq.gz"],
#     "ap2-o_O_3-1": ["SRR3437961_1.fastq.gz", "SRR3437961_2.fastq.gz"],
#     "ap2-o_O_3-2": ["SRR3437962_1.fastq.gz", "SRR3437962_2.fastq.gz"],
#     "wt_O_1-1": ["SRR3437923_1.fastq.gz", "SRR3437923_2.fastq.gz"],
#     "wt_O_1-2": ["SRR3437924_1.fastq.gz", "SRR3437924_2.fastq.gz"],
#     "wt_O_2-1": ["SRR3437937_1.fastq.gz", "SRR3437937_2.fastq.gz"],
#     "wt_O_2-2": ["SRR3437938_1.fastq.gz", "SRR3437938_2.fastq.gz"],
#     "wt_O_3-1": ["SRR3437953_1.fastq.gz", "SRR3437953_2.fastq.gz"],
#     "wt_O_3-2": ["SRR3437954_1.fastq.gz", "SRR3437954_2.fastq.gz"],
# }

import pandas

ss = pandas.read_csv(config['sample_sheet'], sep= '\t', comment= '#')

UPTO = [15000000, 5000000, 1000000, 500000, 50000]

wildcard_constraints:
    upto = '|'.join([re.escape(str(x)) for x in UPTO]),
    sample_id = '|'.join([re.escape(str(x)) for x in ss.sample_id]),


rule all:
    input:
        os.path.join(workflow.basedir, 'results/dge_by_nreads.pdf'),


rule download_fasta_reference:
    output:
        genome="ref/PlasmoDB-53_PbergheiANKA_Genome.fasta",
    shell:
        r"""
        curl -s -L https://plasmodb.org/common/downloads/release-53/PbergheiANKA/fasta/data/PlasmoDB-53_PbergheiANKA_Genome.fasta > {output.genome}
        """


rule download_gff_reference:
    output:
        gff="ref/PlasmoDB-53_PbergheiANKA.gff",
    shell:
        r"""
        curl -s -L https://plasmodb.org/common/downloads/release-53/PbergheiANKA/gff/data/PlasmoDB-53_PbergheiANKA.gff > {output.gff}
        """


rule hisat2_index:
    input:
        genome="ref/PlasmoDB-53_PbergheiANKA_Genome.fasta",
    output:
        idx="ref/PlasmoDB-53_PbergheiANKA_Genome.fasta.8.ht2",
    shell:
        r"""
        hisat2-build -p 4 --seed 1234 -f {input.genome} {input.genome}
        """


rule downsample_fastq:
    input:
        fastq=(
            lambda wc: ss[ss.sample_id == wc.sample_id].R1.iloc[0]
            if wc.r == "1"
            else ss[ss.sample_id == wc.sample_id].R2.iloc[0]
        ),
    output:
        fastq= temp("{upto}/downsample_fastq/{sample_id}_{r}.fastq.gz"),
    shell:
        r"""
        zcat {input.fastq} \
        | paste - - - - \
        | awk -v OFS='\t' -v FS='\t' '{{print "chr1", 0, 1, $0}}' \
        | bedtools sample -seed 1234 -n {wildcards.upto} -i - \
        | cut -f 4- \
        | tr '\t' '\n' \
        | pigz > {output.fastq}
        """


rule cutadapt:
    input:
        read1="{upto}/downsample_fastq/{sample_id}_1.fastq.gz",
        read2="{upto}/downsample_fastq/{sample_id}_2.fastq.gz",
    output:
        read1=temp("{upto}/cutadapt/{sample_id}_1.fastq.gz"),
        read2=temp("{upto}/cutadapt/{sample_id}_2.fastq.gz"),
    shell:
        r"""
        cutadapt --quality-cutoff 15 --minimum-length 10 --cores 2 -a AGATCGGAAGAGC \
                -o {output.read1} -p {output.read2} {input.read1} {input.read2}
        """


rule hisat2_align:
    input:
        idx="ref/PlasmoDB-53_PbergheiANKA_Genome.fasta.8.ht2",
        genome="ref/PlasmoDB-53_PbergheiANKA_Genome.fasta",
        read1="{upto}/cutadapt/{sample_id}_1.fastq.gz",
        read2="{upto}/cutadapt/{sample_id}_2.fastq.gz",
    output:
        bam="{upto}/hisat2/{sample_id}.bam",
        hlog="{upto}/hisat2/{sample_id}.log",
    shell:
        r"""
        hisat2 --summary-file {output.hlog} --new-summary --fr \
           --rna-strandness RF --max-intronlen 5000 --threads 8 -x {input.genome} \
           -1 {input.read1} -2 {input.read2} \
        | samtools sort > {output.bam}
        """


rule index_bam:
    input:
        bam="{upto}/hisat2/{sample_id}.bam",
    output:
        bai="{upto}/hisat2/{sample_id}.bam.bai",
    shell:
        r"""
        samtools index {output.bam}
        """


rule featureCounts:
    input:
        bam="{upto}/hisat2/{sample_id}.bam",
        gff="ref/PlasmoDB-53_PbergheiANKA.gff",
    output:
        cnt="{upto}/featureCounts/{sample_id}.tsv",
    params:
        strand= lambda wc: ss[ss.sample_id == wc.sample_id].strand.iloc[0],
    shell:
        r"""
        featureCounts -p -T 4 -Q 10 -s {params.strand} -t exon -g gene_id -a {input.gff} \
            -o {output.cnt} {input.bam}
        """


rule count_matrix:
    input:
        cnt=expand("{{upto}}/featureCounts/{sample_id}.tsv", sample_id=ss.sample_id),
    output:
        cnt_mat="{upto}/featureCounts/counts.tsv",
    run:
        import pandas
        import os
        import re

        counts = []
        for tsv in input.cnt:
            cnt = pandas.read_csv(tsv, sep="\t", comment="#")
            cnt.rename(columns={cnt.columns[cnt.shape[1] - 1]: "count"}, inplace=True)
            sample_id = re.sub("\.tsv$", "", os.path.basename(tsv))
            cnt["sample_id"] = sample_id
            counts.append(cnt[["Geneid", "sample_id", "count"]])
        counts = pandas.concat(counts)
        mat = counts.pivot_table(
            index="Geneid", columns="sample_id", values="count"
        ).reset_index()
        mat.to_csv(output.cnt_mat, sep="\t", index=False)


rule edger:
    input:
        cnt_mat="{upto}/featureCounts/counts.tsv",
        ss = config['sample_sheet'],
    output:
        dge="{upto}/edger/dge.tsv",
    shell:
        r"""
cat <<'EOF'> {rule}.$$.tmp.R
library(edgeR)
library(data.table)

ss <- fread(cmd= 'grep -v "^#" {input.ss}')

mat <- fread('{input.cnt_mat}')
mat <- as.matrix(mat, rownames= 'Geneid')

mat <- mat[,match(ss$sample_id, colnames(mat))]
stopifnot(identical(ss$sample_id, colnames(mat)))

group <- ss$group
design <- model.matrix(~0 + group)
colnames(design) <- sub('group', '', colnames(design), fixed= TRUE)

y <- DGEList(counts= mat, group= group)
keep <- filterByExpr(y)
y <- y[keep, , keep.lib.sizes=FALSE]
y <- calcNormFactors(y)
y <- estimateDisp(y, design)

contrasts <- makeContrasts(
    `16h vs 22h` = h16 - h22,
    `Ook vs Gam` = Ook - Gam,
    levels= design)

# fit <- glmQLFit(y, design, prior.count= 1)
# dge <- list()
# for(cntr in colnames(contrasts)){{
#     lfc <- glmTreat(fit, contrast= contrasts[, cntr], lfc= log2(1.5))
#     detable <- topTags(lfc, n= nrow(y))$table
#     detable$gene_id <- row.names(detable)
#     detable <- data.table(detable)
#     detable[, contrast := cntr]
#     dge[[length(dge)+1]] <- detable
# }}
# dge <- rbindlist(dge)

fit <- glmFit(y, design, prior.count= 1)
dge <- list()
for(cntr in colnames(contrasts)){{
    lfc <- glmTreat(fit, contrast= contrasts[, cntr], lfc= log2(1.5))
    detable <- topTags(lfc, n= nrow(y))$table
    detable$gene_id <- row.names(detable)
    detable <- data.table(detable)
    detable[, contrast := cntr]
    dge[[length(dge)+1]] <- detable
}}
dge <- rbindlist(dge)

options(scipen= 9)
dge[, downsample_size := {wildcards.upto}]
write.table(dge, '{output.dge}', row.names= FALSE, sep= '\t', quote= FALSE)

EOF
Rscript {rule}.$$.tmp.R
rm {rule}.$$.tmp.R
        """


rule summarize_dge:
    input:
        dge= expand('{upto}/edger/dge.tsv', upto= UPTO),
    output:
        dge= os.path.join(workflow.basedir, 'results/dge_by_nreads.pdf'),
    shell:
        r"""
cat <<'EOF'> {rule}.$$.tmp.R
library(data.table)
library(ggplot2)

ff <- strsplit('{input.dge}', ' ')[[1]]

dge <- list()
for(x in ff) {{
    dge[[length(dge) + 1]] <- fread(x)
}}
dge <- rbindlist(dge)

smry <- dge[, list(.N, N_FDR= sum(FDR < 0.01)), list(downsample_size, contrast)][order(contrast, downsample_size)]

lab_pos <- log10(sort(unique(smry$downsample_size/1e6)))
lab_text <- sprintf('%.2f', sort(unique(smry$downsample_size/1e6)))

gg <- ggplot(data= smry, aes(x= log10(downsample_size/1e6), y= N_FDR)) +
    geom_col() +
    facet_wrap(~contrast, scales= 'free_y') +
    scale_x_continuous(breaks= lab_pos, labels= lab_text) +
    xlab('Millions of reads\ndownsampled from original fastq files') +
    ylab('N genes with\nFDR < 0.01') +
    theme_light() +
    theme(strip.text= element_text(colour= 'black'), axis.text.x= element_text(angle= 90, vjust = 0.5, hjust=1))
ggsave('{output.dge}', width= 16, height= 8, units= 'cm')

EOF
Rscript {rule}.$$.tmp.R
rm {rule}.$$.tmp.R
        """
        
