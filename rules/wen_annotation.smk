from snakemake.remote import FTP, HTTP
HTTP = HTTP.RemoteProvider()

rule wen15_supp_s3:
    input:
        HTTP.remote("doi.org/10.1371/journal.pcbi.1004441.s011", static=True, keep_local=True, allow_redirects=True)
    conda:
        "../envs/pandas.yaml"
    output:
        "ref/wen15_mirtrons.gff3"
    script:
        "../scripts/convert_wen15.py"
