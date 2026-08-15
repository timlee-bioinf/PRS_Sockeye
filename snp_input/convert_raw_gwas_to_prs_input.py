from pathlib import Path
import csv
import re

src = Path('snp_input/gwas-association-downloaded_2026-08-15-pubmedId_36914875.tsv')
out = Path('snp_input/PRS_SNP_input_2023.tsv')

rows = []
seen = set()
with src.open('r', newline='', encoding='utf-8') as f:
    reader = csv.DictReader(f, delimiter='\t')
    for r in reader:
        snp = (r.get('SNPS') or r.get('SNP_ID_CURRENT') or r.get('STRONGEST SNP-RISK ALLELE') or '').strip()
        if not snp:
            continue

        # Accept rsIDs from either SNPS or SNP_ID_CURRENT; strip off allele suffixes like rs123-A
        rs = snp.split('-')[0].strip()
        if not rs.lower().startswith('rs'):
            continue

        strong = (r.get('STRONGEST SNP-RISK ALLELE') or '').strip()
        allele = ''
        if '-' in strong:
            allele = strong.rsplit('-', 1)[-1].strip()
        if allele not in {'A', 'C', 'G', 'T'}:
            continue

        weight_raw = (r.get('OR or BETA') or '').strip()
        m = re.search(r'([-+]?\d*\.?\d+(?:[Ee][-+]?\d+)?)', weight_raw)
        if not m:
            continue
        weight = float(m.group(1))

        key = (rs, allele)
        if key in seen:
            continue
        seen.add(key)
        rows.append({'rsID': rs, 'Risk allele': allele, 'Weight': weight})

with out.open('w', newline='', encoding='utf-8') as f:
    writer = csv.DictWriter(f, fieldnames=['rsID', 'Risk allele', 'Weight'], delimiter='\t')
    writer.writeheader()
    writer.writerows(rows)

print(f'Wrote {len(rows)} rows to {out}')
print('Preview:')
for row in rows[:3]:
    print(row)
