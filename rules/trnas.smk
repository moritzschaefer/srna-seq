from snakemake.remote import FTP, HTTP
HTTP = HTTP.RemoteProvider()

postfixes = [
    # from sample/fastq directory
    '',
    # from working directory
    '.mapped',
    '.mapped.covgPre.txt',
    '.mapped.covg.txt',
    '.mapped.speciesInfo.txt',
    '.mapped.top50covgPre.txt',
    # '.mapped.top50covgPre.txt.top50Pre_tdr.pdf',  # not always generated.. (if no pre_tdr reads found..)
    '.mapped.top50covg.txt',
    '.mapped.top50covg.txt.top50_tdr.pdf']

rule tdrmapper:
    input:
        trimmed_reads="trimmed/{sample}.fastq.gz",  # trimmed reads
    # params:
        reference_trnas=HTTP.remote('https://raw.githubusercontent.com/sararselitsky/tDRmapper/master/mm10_mature_pre_for_tdrMapper.fa', keep_local=True)
    output:
        expand('tdrmapper/{{sample}}.fastq.gz.hq_cs{postfix}', postfix=postfixes)
    conda:
        "../envs/tdrmapper.yaml"
    log:
        out="logs/tdrmapper/{sample}.log",
        err="logs/tdrmapper/{sample}.err"
    shell:
        "TdrMappingScripts.pl {input.reference_trnas} {input.trimmed_reads} 2> {log.err} > {log.out} && rm trimmed/{wildcards.sample}.fastq.gz.hq_cs.*not* && mv trimmed/{wildcards.sample}.fastq.gz.hq_cs {wildcards.sample}.fastq.gz.hq_cs.mapped* tdrmapper/"

rule combine_counts:
    '''
    We need a dedicated normalization function
    '''
    input:
        expand("tdrmapper/{unit.sample}-{unit.unit}.fastq.gz.hq_cs.mapped.speciesInfo.txt", unit=units.itertuples())  # 'unit' is known from Snakefile
    output:
        all='counts/tdrs.tsv',
        all_cpm='counts/tdrs.cpm.tsv'
    script:
        '../scripts/normalize_cpm.py'

rule tdr_pca:
    '''
    Run a PCA of the analyzed samples
    '''
    input:
        'counts/tdrs.cpm.tsv'
    output:
        'plot/tdrs_pca.png'
    conda:
        '../envs/pandas.yaml'
    script:
        '../scripts/plot_tdr_pca.py'

rule tdr_correlation:
    '''
    Correlation plot for all samples
    '''
    input:
        'counts/tdrs.cpm.tsv'
    output:
        'plot/tdrs_correlation.png'
    conda:
        '../envs/pandas.yaml'
    script:
        '../scripts/plot_tdr_correlation.py'
