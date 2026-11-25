from Bio.Seq import Seq

dna_seq = "ATGAAGTAA"
bases = ["A", "T", "G", "C"]
orig_protein = str(Seq(dna_seq).translate(to_stop=True))

def analyze_mutations(seq):
    results = []
    
    for i in range(len(seq)):
        original_base = seq[i]
        
        for new_base in bases:
            if new_base != original_base:
                new_seq = seq[:i] + new_base + seq[i+1:]
                mut_protein = str(Seq(new_seq).translate(to_stop=True))
                
                if new_seq[:3] != "ATG":
                    mut_type = "시작소실(Start-loss)"
                    mut_protein = "-"
                elif mut_protein == orig_protein:
                    mut_type = "침묵(Silent)"
                elif len(mut_protein) < len(orig_protein):
                    mut_type = "종결(Nonsense)"
                elif len(mut_protein) > len(orig_protein):
                    mut_type = "종결소실(Stop-loss)"
                else:
                    mut_type = "과오(Missense)"
                
                if mut_protein == "-":
                    p_len = 0
                else:
                    p_len = len(mut_protein)
                
                results.append((i + 1, original_base, new_base, mut_type, mut_protein, p_len))
                
    return results

mutation_results = analyze_mutations(dna_seq)

print(f"원래 DNA: {dna_seq} -> 원래 단백질: {orig_protein} (길이: {len(orig_protein)})\n")
for pos, orig, new, m_type, m_prot, p_len in mutation_results:
    print(f"{pos}번째 ({orig}->{new}) | {m_type:16s} | 단백질: {m_prot:4s} (길이: {p_len})")

