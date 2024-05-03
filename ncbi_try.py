from Bio.Blast import NCBIWWW, NCBIXML


def run_ncbi_blast(sequence):
    # Perform BLAST search on NCBI
    result_handle = NCBIWWW.qblast("blastn", "nt", sequence)
    blast_records = NCBIXML.parse(result_handle)

    top_matches = []
    for record in blast_records:
        # Iterate through each BLAST record
        for alignment in record.alignments[:1]:  # Limit to top result
            # Extract information from alignment
            gene_name = alignment.title.split("|")[4].split()[0]  # Extract gene name
            accession = alignment.accession  # Extract accession number
            scientific_name = alignment.hit_def.split("[")[0].strip()  # Extract scientific name
            e_score = alignment.hsps[0].expect  # Extract E score
            alignment_score = alignment.hsps[0].score  # Extract alignment score
            identity = alignment.hsps[0].identities / alignment.hsps[0].align_length * 100  # Calculate identity (%)
            fasta_sequence = alignment.hsps[0].sbjct  # Extract fasta sequence

            # Append data to top_matches list
            top_matches.append(
                (gene_name, accession, scientific_name, e_score, alignment_score, identity, fasta_sequence))

    return top_matches


def get_genomic_sequence_from_user():
    sequence = input("Enter the genomic sequence: ")
    return sequence


# Get genomic sequence from user
sequence_A = get_genomic_sequence_from_user()

# Run NCBI BLAST search and retrieve top result
top_result = run_ncbi_blast(sequence_A)[0]

# Extract data from top result
gene_name, accession, scientific_name, e_score, alignment_score, identity, fasta_sequence = top_result

# Write data to a text file
with open("CLICKHERE/flask_file/blast_result.txt", "w") as file:
    file.write("Gene Name: " + gene_name + "\n")
    file.write("Accession Number: " + accession + "\n")
    file.write("Scientific Name: " + scientific_name + "\n")
    file.write("E Score: " + str(e_score) + "\n")
    file.write("Alignment Score: " + str(alignment_score) + "\n")
    file.write("Identity (%): " + str(identity) + "\n")
    file.write("Fasta Sequence: " + fasta_sequence + "\n")

print("Data written to blast_result.txt.")
