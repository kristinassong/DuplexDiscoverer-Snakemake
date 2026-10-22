###############################################################
# Download FASTQ files
###############################################################

rule wget_fq:
    input:
        config["sra_ids"]
    output:
        "resources/data/download_status.txt"
    message:
        "Downloading FASTQ files via wget."
    shell:
        """
        echo 'Downloading fastq files...' > {output} 
        while IFS=$'\t' read protocol url; do 
            srr=$(basename $url | sed 's/.fastq.gz//')
            outfile="resources/data/${{protocol}}_${{srr}}.fastq.gz"
            echo "Downloading $outfile" >> {output}
            wget -q -O $outfile $url && echo "Done: $outfile" >> {output}
        done < {input}
        echo 'Download complete.' >> {output}
        """