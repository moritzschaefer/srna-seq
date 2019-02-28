import pandas as pd
from pyliftover import LiftOver

lo = LiftOver('mm9', 'mm10')

df = pd.read_excel(snakemake.input[0]).set_index('Locus')

# the previously unknown miRNAs start with uc
df = df.loc[df.index.str.startswith('uc')]


def row_to_gffs(index, row):
    chr, start_end = row['hairpin_coordinates'].split(':')
    start, end = (lo.convert_coordinate(chr, int(v))[0][1]
                  for v in start_end.split('-'))
    if row.strand == '+':
        start_5p = start + 1
        end_5p = start + len(row.X5p_arm)
        start_3p = end - (len(row.X3p_arm) - 1)
        end_3p = end
    else:
        start_5p = end - (len(row.X5p_arm) - 1)
        end_5p = end
        start_3p = start
        end_3p = start + len(row.X3p_arm) - 1

    return [
        pd.Series({
            'chr': chr,
            'source': 'wen15',
            'type': 'exon',
            'start': start_5p,
            'end': end_5p,
            'score': 1000,
            'strand': row.strand,
            'phase': '.',
            'attributes': f'gene_id={index}-5p'
        }),
        pd.Series({
            'chr': chr,
            'source': 'wen15',
            'type': 'exon',
            'start': start_3p,
            'end': end_3p,
            'score': 1000,
            'strand': row.strand,
            'phase': '.',
            'attributes': f'gene_id={index}-3p'
        })
    ]


mirs = []

for index, row in df.iterrows():
    mirs.extend(row_to_gffs(index, row))

# chr,src,type,start,end,score,strand,phase,attributes...

pd.DataFrame(mirs).to_csv(
    snakemake.output[0], sep='\t', header=False, index=False)
