import streamlit as st
import subprocess
import os
import concurrent.futures
import pandas as pd
import io

st.set_page_config(page_title="COmapper Variant Calling", layout="centered")
st.title("COmapper: Local Variant Calling Test")

# --- Detect input files ---
fnas = sorted([f for f in os.listdir() if f.endswith(".fna")])
fqs = sorted([f for f in os.listdir() if f.endswith(".fq")])

# --- Require specific files ---
required_refs = {"A_lyrata.fna", "A_thaliana.fna"}
required_fqs = {"sim_lyrata.fq", "sim_thaliana.fq"}

if not required_refs.issubset(fnas) or not required_fqs.issubset(fqs):
    st.error("Missing one or more required files: A_lyrata.fna, A_thaliana.fna, sim_lyrata.fq, sim_thaliana.fq")
    st.stop()

# --- Input combinations ---
inputs = [
    ("A_lyrata.fna", "sim_lyrata.fq", "sim_lyrata.vcf"),
    ("A_thaliana.fna", "sim_thaliana.fq", "sim_thaliana.vcf")
]

# --- Display detected files ---
for ref, fq, _ in inputs:
    st.success(f"Found reference genome: {ref}")
    st.success(f"Found FASTQ file: {fq}")

# --- Variant calling pipeline function ---
def run_pipeline(ref, fq, vcf_out):
    sam = fq.replace(".fq", ".sam")
    bam = fq.replace(".fq", ".bam")
    sorted_bam = fq.replace(".fq", ".sort.bam")

    subprocess.run(["minimap2", "-a", ref, fq], stdout=open(sam, "w"), check=True)
    subprocess.run(["samtools", "view", "-bS", sam], stdout=open(bam, "wb"), check=True)
    subprocess.run(["samtools", "sort", "-o", sorted_bam, bam], check=True)
    subprocess.run(["samtools", "index", sorted_bam], check=True)

    with open(vcf_out, "w") as vcf:
        mpileup = subprocess.Popen(["bcftools", "mpileup", "-f", ref, sorted_bam], stdout=subprocess.PIPE)
        subprocess.run(["bcftools", "call", "-mv", "-Ov"], stdin=mpileup.stdout, stdout=vcf, check=True)

    return vcf_out

# --- Run pipelines when button is clicked ---
if st.button("Run both pipelines"):
    st.info("Running alignments and variant calling in parallel...")

    vcf_files = []
    with concurrent.futures.ThreadPoolExecutor() as executor:
        futures = [executor.submit(run_pipeline, *args) for args in inputs]
        for future in concurrent.futures.as_completed(futures):
            vcf_file = future.result()
            vcf_files.append(vcf_file)
            if os.path.exists(vcf_file):
                st.success(f" Completed: {vcf_file}")
                with open(vcf_file) as f:
                    st.download_button(f"Download {vcf_file}", f.read(), file_name=vcf_file)
            else:
                st.error(f" Failed to generate {vcf_file}")

# --- Optional: Display most recent VCF ---
for vcf in ["sim_lyrata.vcf", "sim_thaliana.vcf"]:
    if os.path.exists(vcf):
        st.subheader(f"Variant Table: {vcf}")
        with open(vcf) as f:
            lines = [line for line in f if not line.startswith("##")]
        df = pd.read_csv(io.StringIO("".join(lines)), sep='\t')
        st.dataframe(df)








    





















