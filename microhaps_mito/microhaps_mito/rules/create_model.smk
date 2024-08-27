import os

NGSENV_BASEDIR = os.environ['NGSENV_BASEDIR']
REFS_DIR = os.path.join(NGSENV_BASEDIR, 'refs')
REFS_MODEL_DIR = os.path.join(REFS_DIR, 'model')
BEDFILE = config.get("bedfile", f"{REFS_MODEL_DIR}/mit_cons.bed")
if not os.path.exists(REFS_MODEL_DIR):
    os.makedirs(REFS_MODEL_DIR)

MIT_SEQS = config.get("MIT_SEQS", f"{REFS_DIR}/mit_all.fasta") # to replace
print(MIT_SEQS)
print(BEDFILE)


rule all:
    input:
        # *[f"{REFS_MODEL_DIR}/{model}.pkl" for model in [
        #     'nc', 'cnb', 'pa', 'ensemble'
        # ]],
        f"{REFS_DIR}/mit_cons.fasta",
        BEDFILE,
        f"{REFS_MODEL_DIR}/mit_spec_all.fasta",
        f"{REFS_MODEL_DIR}/nc.pickle",
        f"{REFS_MODEL_DIR}/cnb.pickle",
        f"{REFS_MODEL_DIR}/pa.pickle",
        f"{REFS_MODEL_DIR}/ensemble.pickle"

rule build_models:
    input:
        f"{REFS_MODEL_DIR}/mit_spec_all.fasta",
        f"{REFS_MODEL_DIR}/ordered_spec.txt"
    output:
        f"{REFS_MODEL_DIR}/nc.pickle",
        f"{REFS_MODEL_DIR}/cnb.pickle",
        f"{REFS_MODEL_DIR}/pa.pickle",
        f"{REFS_MODEL_DIR}/ensemble.pickle"
    shell:
        "ngs-pl create-model --type nc -o {output[0]} --cascade -m {input[1]} {input[0]} && "
        "ngs-pl create-model --type cnb -o {output[1]} --cascade -m {input[1]} {input[0]} && "
        "ngs-pl create-model --type pa -o {output[2]} --cascade -m {input[1]} {input[0]} && "
        "ngs-pl create-model --type ensemble -o {output[3]} --cascade -m {input[1]} {input[0]}"

rule extract_all_hypothetical_products:
    input:
        cons = f"{REFS_DIR}/mit_cons.fasta",
        amplicon_product = f"{REFS_DIR}/mit_cons_single_primersim.fasta",
        primer_out = f"{REFS_DIR}/mit_cons_single_primersim.tsv",
        aligned_mit = f"{REFS_DIR}/mit_all.muscle.fasta",
        metadata = f"{REFS_DIR}/mit_metadata.tsv"
    output:
        f"{REFS_MODEL_DIR}/mit_spec_all.fasta",
        f"{REFS_MODEL_DIR}/ordered_spec.txt",
        BEDFILE,
    run:
        from seqpy.core.bioio import load, save, multisequence, biosequence
        import pandas as pd
        seqs = load(input.cons)
        pr_out = pd.read_csv(input.primer_out, sep="\t")
        
        start, end = (pr_out["end_upstream"][0], pr_out["begin_downstream"][0] -1)
        with open(output[2], "w+") as f:
            f.write(f"{seqs.seqs[0].label}\t{start}\t{end}\n")
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
        primers = f"{REFS_DIR}/concat_primers.fasta",
    output:
        amplicon_product = temp(f"{REFS_DIR}/mit_cons_single_primersim.fasta"),
        primer_out = temp(f"{REFS_DIR}/mit_cons_single_primersim.tsv"),
    shell:
        "spcli $OPT/primersim.py --proof-read --strip-primers --template {input.cons} -o {output.primer_out} --outfragment {output.amplicon_product} {input.primers}"

rule concat_primers:
    input:
        f"{REFS_DIR}/microhap_pr_fwd.min_overlap.fasta",
        f"{REFS_DIR}/microhap_pr_rv.min_overlap.fasta"
    output:
        temp(f"{REFS_DIR}/concat_primers.fasta")
    shell:
        "cat {input} > {output}"

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
        temp(f"{REFS_DIR}/multi-consensus.fasta")
    shell:
        "spcli consensus --synthetic -o {output} {input}"

rule muscle_align_mit:
    input:
        f"{MIT_SEQS}"
    output:
        f"{REFS_DIR}/mit_all.muscle.fasta"
    shell:
        "muscle -super5 {input} -output {output}"