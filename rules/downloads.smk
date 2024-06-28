ncbi_api_key = ""
try:
    with open("config/ncbi_api_key.txt", "r") as f:
        ncbi_api_key = f.read().strip()
except:
    print("No NCBI API key found. Save your key in config/ncbi_api_key.txt")


rule download_gbk:
    output:
        "data/gbk/{acc}.gbk",
    conda:
        "../conda_env/bioinfo.yml"
    params:
        api_key=f"--api-key {ncbi_api_key}" if len(ncbi_api_key) > 0 else "",
    shell:
        """
        ncbi-acc-download {wildcards.acc} \
            -e all \
            -F genbank \
            {params.api_key}

        mv {wildcards.acc}.gbk {output}
        """


rule gbk_to_fa:
    input:
        gbk=rules.download_gbk.output,
    output:
        fa="data/fa/{acc}.fa",
    conda:
        "../conda_env/bioinfo.yml"
    shell:
        """
        python3 scripts/utils/gbk_to_fa.py --gbk {input.gbk} --fa {output.fa}
        """
