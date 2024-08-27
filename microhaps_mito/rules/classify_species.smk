# get the vivaxGEN ngs-pipeline base directory
import os

NGS_PIPELINE_BASE = os.environ["NGS_PIPELINE_BASE"]
NGSENV_BASEDIR = os.environ['NGSENV_BASEDIR']

include: f"{NGS_PIPELINE_BASE}/rules/utilities.smk"

infiles = config['infiles']
outdir = config['outdir']
underscore = config['underscore']

bed_file = NGSENV_BASEDIR + '/' + config['bedfile']
Reference = NGSENV_BASEDIR + '/' + config["reference"]
prim_fw = NGSENV_BASEDIR + '/' + config["prim_fw"]
prim_rv = NGSENV_BASEDIR + '/' + config["prim_rv"]
model_to_run = config.get("model", ["cnb", "nc", "pa", "ensemble"])
fasta_spec = NGSENV_BASEDIR + '/' + config["similarity_comparison"]
cnb_model = NGSENV_BASEDIR + '/' + "refs/model/cnb.pickle"
nc_model = NGSENV_BASEDIR + '/' + "refs/model/nc.pickle"
pa_model = NGSENV_BASEDIR + '/' + "refs/model/pa.pickle"
ensemble_model = NGSENV_BASEDIR + '/' + "refs/model/ensemble.pickle"

models_file_loc = {
    "cnb": cnb_model,
    "nc": nc_model,
    "pa": pa_model,
    "ensemble": ensemble_model
}

from ngs_pipeline import fileutils

read_files = fileutils.ReadFileDict(infiles, underscore=underscore)

def get_read_file(w):
    return {"R1": read_files[w.sample][0][0], "R2": read_files[w.sample][0][1]}

IDs = read_files.keys()

rule all:
    input:
        *[f"{outdir}/{sample}/paired_fasta/{sample}.fasta" for sample in IDs],
        *[f"{outdir}/{sample}/{sample}.{model}.tsv" for sample in IDs for model in model_to_run],
        *[f"{outdir}/{sample}/Report_{sample}_{model}.{strictness}.tsv" for sample in IDs for model in model_to_run for strictness in config['strictness']],
        *[f"{outdir}/combined_{model}_report.tsv" for model in model_to_run],


rule combine_sample_report:
    input:
        report_files = expand(f"{outdir}/{{sample}}/Report_{{sample}}_{{model}}.{{strictness}}.tsv",
            sample = IDs, model = model_to_run,
            strictness = config['strictness']),
    output:
        f"{outdir}/combined_{{model}}_report.tsv"
    run:
        import pandas as pd

        dfs = []
        for file in input.report_files:
            df = pd.read_csv(file, sep="\t")
            dfs.append(df)
        
        combined_df = pd.concat(dfs)
        columns = combined_df.columns
        sorted_columns = ["sample", "model", "reads(non-filtered)", "strictness", "predicted_species"]
        sorted_columns = sorted_columns + sorted([col for col in columns if col not in sorted_columns])
        combined_df = combined_df.round(3) 
        combined_df[sorted_columns].to_csv(output[0], sep="\t", index = False)

rule process_per_sample_result:
    input:
        sample_results = [f"{outdir}/{{sample}}/{{sample}}.{model}.tsv" for 
            model in model_to_run],
    output:
        f"{outdir}/{{sample}}/Report_{{sample}}_{{model}}.{{strictness}}.tsv"
    params:
        filter_thresholds = lambda w: config.get('filtering_criteria').get(w.strictness),
        
        min_reads_consideration = lambda w: config.get('min_reads_consideration').get(w.strictness),

        strictness = lambda w: w.strictness,
        models = model_to_run,
    run:
        import pandas as pd
        import numpy as np

        dfs = []
        for file, model_name in zip(input.sample_results, params.models):
            full_model_name = {"cnb": "Categorical_NaiveBayes", "nc": "Nearest_Centroid",
                "pa": "Pairwise_consensus_aligner_score", 'ensemble':"Ensemble"}.get(model_name)

            df = pd.read_csv(file, sep="\t")
            total_read = df.read_id.dropna().count()

            filter_threshold = params.filter_thresholds[model_name]
            ordered_cols = ["sample", "read_id", "species", "f_pair_diff"]
            class_cols = [a for a in sorted(df.columns) if not a in ordered_cols]
            ordered_cols.extend(class_cols)
            df = df[ordered_cols]
            df = df.query(f"f_pair_diff {filter_threshold}")
            class_name = df.columns[4:]
            prop_dict = {"prop_"+col: len(df[df["species"] == col])/len(df) if len(df) else 0 for col in class_name}
            n_read_dict = {col: len(df[df["species"] == col]) for col in class_name}
            total_filtered_read = sum([nr for nr in n_read_dict.values() if nr >= params.min_reads_consideration])
            
            n_read_str = f"{total_filtered_read} ({total_read})"
            sp_filtered = [f"{sp}: {np.round(nr/total_filtered_read, 2)}" for sp, nr in n_read_dict.items() if nr >= params.min_reads_consideration]
            predicted_sp_str = "; ".join(sp_filtered)
            model_result = pd.DataFrame({
                "sample": [wildcards.sample],
                "model": [full_model_name],
                "reads(non-filtered)": [n_read_str],
                "strictness": [params.strictness],
                "predicted_species": [predicted_sp_str],
                **n_read_dict,
                **prop_dict}
            )
            dfs.append(model_result)
        pd.concat(dfs).to_csv(output[0], sep="\t", index = False)


rule infer_species_with_model:
    input:
        fasta = f"{outdir}/{{sample}}/paired_fasta/{{sample}}.fasta",
    output:
        result = f"{outdir}/{{sample}}/{{sample}}.{{model}}.tsv",
    params:
        sample = lambda w: w.sample,
        model_file = lambda w: models_file_loc.get(w.model),
    shell:
        "ngs-pl run-model --output {output.result} --model {params.model_file} --sample {params.sample} {input.fasta}"

rule bam_to_mated_reads_fasta:
    input:
        bam = f"{outdir}/{{sample}}/nomerge/{{sample}}.bam",
        bai = f"{outdir}/{{sample}}/nomerge/{{sample}}.bam.bai",
        bed = bed_file,
        # ref = Reference,
    output:
        fasta = f"{outdir}/{{sample}}/paired_fasta/{{sample}}.fasta"
    shell:
        "ngs-pl bam2matedfasta --fasta {output.fasta} --bed {input.bed} --bam {input.bam}"


# https://jbloomlab.github.io/dms_tools2/index.html
# https://github.com/lh3/minimap2/issues/521 
rule map_reads:
    input:
        R1 = f"{outdir}/{{sample}}/trimmed/{{sample}}_R1.fastq.gz",
        R2 = f"{outdir}/{{sample}}/trimmed/{{sample}}_R2.fastq.gz",
        ref = Reference
    params:
        minimap2_extra_param = "--sam-hit-only --MD -A2 -B4 -O4,24 -E2,1 --end-bonus 10" # -A2 -B4 -O4,24 -E2,1 --end-bonus 10 
    output:
        bam = f"{outdir}/{{sample}}/nomerge/{{sample}}.bam"
    shell:
        "minimap2 -ax sr {params.minimap2_extra_param} {input[2]} {input[0]} {input[1]} | samtools view -bSh -f 3 | samtools sort -o {output.bam}"


rule trim_reads:
    input:
        R1 = f"{outdir}/{{sample}}/reads/raw_R1.fastq.gz",
        R2 = f"{outdir}/{{sample}}/reads/raw_R2.fastq.gz",
        prim_fw = prim_fw,
        prim_rv = prim_rv,
    output:
        R1 = f"{outdir}/{{sample}}/trimmed/{{sample}}_R1.fastq.gz",
        R2 = f"{outdir}/{{sample}}/trimmed/{{sample}}_R2.fastq.gz",
    shell: 
        "cutadapt -g file:{input.prim_fw} -G file:{input.prim_rv} -o {output.R1} -p {output.R2} --pair-adapters --discard-untrimmed --action=trim {input.R1} {input.R2}"


rule link_reads:
    input:
        unpack(get_read_file)
    output:
        R1 = f"{outdir}/{{sample}}/reads/raw_R1.fastq.gz",
        R2 = f"{outdir}/{{sample}}/reads/raw_R2.fastq.gz",
    run:
        import pathlib

        dest_file_R1 = pathlib.Path(output.R1)
        dest_file_R2 = pathlib.Path(output.R2)
        src_file_R1 = pathlib.Path(input.R1).resolve()
        src_file_R2 = pathlib.Path(input.R2).resolve()
        dest_file_R1.symlink_to(src_file_R1)
        dest_file_R2.symlink_to(src_file_R2)
