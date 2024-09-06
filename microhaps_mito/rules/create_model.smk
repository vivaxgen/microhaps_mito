import os

NGSENV_BASEDIR = os.environ['NGSENV_BASEDIR']
REFS_DIR = os.path.join(NGSENV_BASEDIR, 'refs')
REFS_MODEL_DIR = os.path.join(REFS_DIR, 'model')
BEDFILE = os.path.join(NGSENV_BASEDIR, config.get("bedfile", "refs/model/mit_cons.bed"))
if not os.path.exists(REFS_MODEL_DIR):
    os.makedirs(REFS_MODEL_DIR)

MIT_SEQS = os.path.join(NGSENV_BASEDIR, config.get("MIT_SEQS", "refs/mit_all.fasta"))  # to replace
PRIMER_SET = config.get("primer_set", {
    "set1": {
        "forward": "refs/MIT_species1_fwd.fasta",
        "reverse": "refs/MIT_species1_rv.fasta",
    },
    "set2": {
        "forward": "refs/MIT_species2_fwd.fasta",
        "reverse": "refs/MIT_species2_rv.fasta",
    }
})
PRIMER_TRIM_SET = config.get("primer_trim_set", {
    "forward": "refs/microhaps_pr_fwd.fasta",
    "reverse": "refs/microhaps_pr_rv.fasta",
})

rule all:
    input:
        f"{REFS_DIR}/mit_cons.fasta",
        BEDFILE,
        *[os.path.join(NGSENV_BASEDIR, primer_path) for primer_path in PRIMER_TRIM_SET.values()],
        *sum([
            [f"{REFS_MODEL_DIR}/{primer_set}.nc.pickle",
            f"{REFS_MODEL_DIR}/{primer_set}.cnb.pickle",
            f"{REFS_MODEL_DIR}/{primer_set}.pa.pickle",
            f"{REFS_MODEL_DIR}/{primer_set}.ensemble.pickle"] for primer_set in PRIMER_SET.keys()
        ],[]),
        PRIMER_TRIM_FILES

rule build_models:
    input:
        f"{REFS_MODEL_DIR}/{{primer_set}}_mit_spec_all.fasta",
        f"{REFS_MODEL_DIR}/{{primer_set}}_ordered_spec.txt"
    output:
        f"{REFS_MODEL_DIR}/{{primer_set}}.nc.pickle",
        f"{REFS_MODEL_DIR}/{{primer_set}}.cnb.pickle",
        f"{REFS_MODEL_DIR}/{{primer_set}}.pa.pickle",
        f"{REFS_MODEL_DIR}/{{primer_set}}.ensemble.pickle"
    shell:
        "ngs-pl create-model --type nc -o {output[0]} --cascade -m {input[1]} {input[0]} && "
        "ngs-pl create-model --type cnb -o {output[1]} --cascade -m {input[1]} {input[0]} && "
        "ngs-pl create-model --type pa -o {output[2]} --cascade -m {input[1]} {input[0]} && "
        "ngs-pl create-model --type ensemble -o {output[3]} --cascade -m {input[1]} {input[0]}"

rule concat_primer_for_trimming:
    localrule: True
    input:
        fwd = [os.path.join(NGSENV_BASEDIR, primer_fasta["forward"]) for primer_fasta in PRIMER_SET.values()],
        rv = [os.path.join(NGSENV_BASEDIR, primer_fasta["reverse"]) for primer_fasta in PRIMER_SET.values()],
    output:
        fwd = os.path.join(NGSENV_BASEDIR, PRIMER_TRIM_SET["forward"]),
        rv = os.path.join(NGSENV_BASEDIR, PRIMER_TRIM_SET["reverse"]),
    shell:
        "cat {input.fwd} > {output.fwd} && "
        "cat {input.rv} > {output.rv}"

rule combine_bed:
    localrule: True
    input:
        expand(f"{REFS_MODEL_DIR}/{{primer_set}}_mit_cons.bed", primer_set = PRIMER_SET.keys())
    output:
        BEDFILE
    shell:
        "cat {input} > {output}"

rule extract_all_hypothetical_products:
    input:
        cons = f"{REFS_DIR}/mit_cons.fasta",
        amplicon_product = f"{REFS_DIR}/{{primer_set}}_mit_cons_single_primersim.fasta",
        primer_out = f"{REFS_DIR}/{{primer_set}}_mit_cons_single_primersim.tsv",
        aligned_mit = f"{REFS_DIR}/mit_all.muscle.fasta",
        metadata = f"{REFS_DIR}/mit_metadata.tsv"
    output:
        temp(f"{REFS_MODEL_DIR}/{{primer_set}}_mit_spec_all.fasta"),
        temp(f"{REFS_MODEL_DIR}/{{primer_set}}_ordered_spec.txt"),
        temp(f"{REFS_MODEL_DIR}/{{primer_set}}_mit_cons.bed"),
    run:
        from seqpy.core.bioio import load, save, multisequence, biosequence
        import pandas as pd
        seqs = load(input.cons)
        pr_out = pd.read_csv(input.primer_out, sep="\t")
        
        start, end = (pr_out["end_upstream"][0], pr_out["begin_downstream"][0] -1)
        with open(output[2], "w+") as f:
            f.write(f"{seqs.seqs[0].label}\t{start}\t{end}\t{wildcards.primer_set}\n")
        cons_prod = seqs.seqs[0].seq[start:end]
        assert load(input.amplicon_product).seqs[0].seq.decode("utf-8") == cons_prod.decode("utf-8").upper(), "Amplicon product does not match consensus product"
        
        align_seqs = load(input.aligned_mit)
        all_prods = multisequence()
        label_order = []
        for seq in align_seqs:
            all_prods.append(
                biosequence(
                    seq.label,
                    seq.seq[start:end]
                )
            )
            label_order.append(seq.label)
        save(all_prods, output[0])

        meta = pd.read_table(input.metadata)
        assert set([a in meta["fasta_name"].tolist() for a in label_order]) == {True}, "Not all sequences in metadata"
        with open(output[1], "w+") as f:
            for label in label_order:
                f.write(meta[meta["fasta_name"] == label]["species"].values[0] + "\n")


rule primersim:
    input:
        cons = f"{REFS_DIR}/mit_cons.fasta",
        primers = f"{REFS_DIR}/concat_primers_{{primer_set}}.fasta",
    output:
        amplicon_product = temp(f"{REFS_DIR}/{{primer_set}}_mit_cons_single_primersim.fasta"),
        primer_out = temp(f"{REFS_DIR}/{{primer_set}}_mit_cons_single_primersim.tsv"),
    shell:
        "spcli $OPT/primersim.py --proof-read --strip-primers --template {input.cons} -o {output.primer_out} --outfragment {output.amplicon_product} {input.primers}"

def get_primer_set(wildcards):
    return PRIMER_SET[wildcards.primer_set].values()

rule concat_primers:
    localrule: True
    input:
        get_primer_set
    output:
        temp(f"{REFS_DIR}/concat_primers_{{primer_set}}.fasta")
    shell:
        "cat {input} > {output}""

rule extract_single_consensus:
    input:
        f"{REFS_DIR}/multi-consensus.fasta"
    output:
        f"{REFS_DIR}/mit_cons.fasta"
    run:
        from seqpy.core.bioio import load, save, multisequence
        seqs = load(input[0])
        b=multisequence()
        seqs.seqs[0].label = "MIT_CONS"
        b.append(seqs.seqs[0])
        save(b, output[0])

rule consensus_all_mit:
    input:
        f"{REFS_DIR}/mit_all.muscle.fasta"
    output:
        f"{REFS_DIR}/multi-consensus.fasta"
    shell:
        "spcli consensus --synthetic -o {output} {input}"

rule muscle_align_mit:
    input:
        f"{MIT_SEQS}"
    output:
        f"{REFS_DIR}/mit_all.muscle.fasta"
    shell:
        "muscle -super5 {input} -output {output}"