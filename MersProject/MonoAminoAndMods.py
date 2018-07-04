# MASS OF H20 is added for mass calculations
H20_MASS = 18.010565

# All possible modifications
modTable = {

    '4-hydroxynonenal (HNE)': ['C', 'H', 'K', 156.11504],
    'Acetylation (K)': ['K', 42.010567],
    'Beta-methylthiolation': ['C', 45.98772],
    'Carbamidomethylation': ['C', 57.021465],
    'Carboxylation (E)': ['E', 43.98983],
    'Carboxymethyl': ['C', 58.005478],
    'Citrullination': ['R', 0.984016],
    'Deamidation (NQ)': ['N', 'Q', 0.984016],
    'Dimethylation(KR)': ['K', 'R', 28.0313],
    'Dioxidation (M)': ['M', 31.989828],
    'FAD': ['C', 'H', 'Y', 783.1415],
    'Farnesylation': ['C', 204.1878],
    'Geranyl-geranyl': ['C', 272.2504],
    'Guanidination': ['K', 42.021797],
    'HexNAcylation (N)': ['N', 203.07938],
    'Hexose (NSY)': ['N', 'S', 'Y', 162.0528],
    'Lipoyl': ['K', 188.03296],
    'Methylation(KR)': ['K', 'R', 14.01565],
    'Methylation(others)': ['T', 'S', 'C', 'H', 'D', 'E', 14.01565],
    'Oxidation (HW)': ['H', 'W', 15.994915],
    'Oxidation (M)': ['M', 15.994915],
    'Palmitoylation': ['C', 'S', 'T', 'K', 238.22966],
    'Phosphopantetheine': ['S', 340.0858],
    'Phosphorylation (HCDR)': ['H', 'C', 'D', 'R', 79.96633],
    'Phosphorylation (STY)': ['S', 'T', 'Y', 79.96633],
    'Propionamide': ['C', 71.03712],
    'Pyridoxal phosphate': ['K', 229.014],
    'S-pyridylethylation': ['C', 105.057846],
    'Sulfation': ['Y', 'S', 'T', 79.95682],
    'Sulphone': ['M', 31.989828],
    'Ubiquitin': ['T', 'S', 'C', 'K', 114.04293],
    'Ubiquitination': ['K', 383.2281],


}

# Mono-isotopic mass
monoAminoMass = {
    'A': 71.03711,
    'R': 156.10111,
    'N': 114.0493,
    'D': 115.02694,
    'C': 103.00919,
    'E': 129.04259,
    'Q': 128.05858,
    'G': 57.02146,
    'H': 137.05891,
    'I': 113.08406,
    'L': 113.08406,
    'K': 128.09496,
    'M': 131.04049,
    'F': 147.06841,
    'P': 97.05276,
    'S': 87.03203,
    'T': 101.04768,
    'W': 186.07931,
    'Y': 163.06333,
    'V': 99.06841,

}


# Move to mers
modificationList = ['4-hydroxynonenal (HNE)', 'Acetylation (K)', 'Beta-methylthiolation', 'Carbamidomethylation',
                    'Carboxylation (E)', 'Carboxymethyl', 'Citrullination', 'Deamidation (NQ)', 'Dimethylation(KR)',
                    'Dioxidation (M)', 'FAD', 'Farnesylation', 'Geranyl-geranyl', 'Guanidination', 'HexNAcylation (N)',
                    'Hexose (NSY)', 'Lipoyl', 'Methylation(KR)', 'Methylation(others)', 'Oxidation (HW)',
                    'Oxidation (M)', 'Palmitoylation', 'Phosphopantetheine', 'Phosphorylation (HCDR)',
                    'Phosphorylation (STY)', 'Propionamide', 'Pyridoxal phosphate', 'S-pyridylethylation',
                    'Sulfation', 'Sulphone', 'Ubiquitin', 'Ubiquitination']

