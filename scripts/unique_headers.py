# unique_headers.py

from Bio import SeqIO

def generate_header_mapping(fasta_file, mapping_file):
    headers = set()
    updated_headers = {}
    with open(mapping_file, 'w') as mapfile:
        mapfile.write("Old_Header\tNew_Header\n")
        for record in SeqIO.parse(fasta_file, "fasta"):
            header = record.id
            if header in headers:
                count = 1
                new_header = f"{header}_{count}"
                while new_header in headers:
                    count += 1
                    new_header = f"{header}_{count}"
                updated_headers[header] = new_header
                mapfile.write(f"{header}\t{new_header}\n")
            else:
                updated_headers[header] = header
                mapfile.write(f"{header}\t{header}\n")
            headers.add(updated_headers[header])
    return updated_headers
