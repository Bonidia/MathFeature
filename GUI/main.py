import os
import warnings
from gooey import Gooey, GooeyParser
warnings.filterwarnings('ignore')


path = os.path.dirname(__file__)
if path != '':
    path = path + '/'


descrip = {'Preprocessing (To eliminate noise from sequences) - DNA/RNA': 0,
           'Numerical Mapping - DNA/RNA (Not available yet)': 1,
           'Numerical Mapping - Protein (Not available yet)': 2,
           'Chaos Game Representation - DNA/RNA (Not available yet )': 3,
           'Numerical Mapping with Fourier Transform - DNA/RNA': 4,
           'Numerical Mapping with Fourier Transform - Protein': 5,
           'Shannon and Tsallis Entropy - DNA/RNA/Protein': 6,
           'Tsallis Entropy (entropic parameter) - DNA/RNA/Protein': 7,
           'Complex Networks - DNA/RNA/Protein': 8,
           'Complex Networks (Faster and without threshold) - DNA/RNA/Protein': 9,
           'Xmer k-Spaced Ymer Composition Frequency (kGap) - DNA/RNA/Protein': 10,
           'Basic k-mer - DNA/RNA': 11,
           'Nucleic acid composition (NAC) - DNA/RNA': 12,
           'Di-nucleotide composition (DNC) - DNA/RNA': 13,
           'Tri-nucleotide composition (TNC) - DNA/RNA': 14,
           'ORF Features or Coding Features - DNA/RNA': 15,
           'Fickett score - DNA/RNA': 16,
           'Pseudo K-tuple nucleotide composition - DNA/RNA': 17,
           'Basic k-mer - Protein': 18,
           'Amino acid composition (AAC) - Protein': 19,
           'Dipeptide composition (DPC) - Protein': 20,
           'Tripeptide composition (TPC) - Protein': 21,
           'AAindex Table - Protein': 22}


@Gooey(program_name='MathFeature',
       default_size=(800, 600),
       language='english',
       tabbed_groups=True,
       image_dir=path + 'img',
       menu=[{
           'name': 'File',
           'items': [{
               'type': 'AboutDialog',
               'menuTitle': 'About',
               'name': 'MathFeature',
               'description': 'Feature Extraction Package for Biological Sequences Based on Mathematical Descriptors',
               'version': '1.2.1',
               'copyright': '2021',
               'website': 'https://bonidia.github.io/MathFeature/',
               'developer': 'https://bonidia.github.io/website/'},
               {
                   'type': 'Link',
                   'menuTitle': 'Visit Our Site',
                   'url': 'https://bonidia.github.io/MathFeature/'
               }]},
           {
               'name': 'Help',
               'items': [{
                   'type': 'Link',
                   'menuTitle': 'Documentation',
                   'url': 'https://bonidia.github.io/MathFeature/'
               }]
           }])
def main():
    parser = GooeyParser(description='Feature Extraction Package for Biological Sequences')
    screen = parser.add_argument_group('Main screen - Choose descriptors')
    screen.add_argument('-t',
                        '--descriptor',
                        metavar='Descriptors',
                        required=True,
                        choices=descrip,
                        nargs='+')

    return parser.parse_args()


######################################################
######################################################
if __name__ == '__main__':
    args = main()
    descriptor = args.descriptor
    descriptor = ''.join(descriptor)
    if descrip[descriptor] == 0:
        os.system('python ' + path + 'preprocessing-screen.py')
    elif descrip[descriptor] == 1:
        print('Not available yet')
    elif descrip[descriptor] == 2:
        print('Not available yet')
    elif descrip[descriptor] == 3:
        print('Not available yet')
    elif descrip[descriptor] == 4:
        os.system('python ' + path + 'FourierClass-screen.py')
    elif descrip[descriptor] == 5:
        os.system('python ' + path + 'Mappings-Protein.py')
    elif descrip[descriptor] == 6:
        os.system('python ' + path + 'EntropyClass-screen.py')
    elif descrip[descriptor] == 7:
        os.system('python ' + path + 'TsallisEntropy-screen.py')
    elif descrip[descriptor] == 8:
        os.system('python ' + path + 'ComplexNetworksClass-screen.py')
    elif descrip[descriptor] == 9:
        os.system('python ' + path + 'ComplexNetworksClass-v2-screen.py')
    elif descrip[descriptor] == 10:
        os.system('python ' + path + 'Kgap-screen.py')
    elif descrip[descriptor] == 11:
        os.system('python ' + path + 'ExtractionTechniques-screen.py')
    elif descrip[descriptor] == 12:
        os.system('python ' + path + 'ExtractionTechniques-screen.py')
    elif descrip[descriptor] == 13:
        os.system('python ' + path + 'ExtractionTechniques-screen.py')
    elif descrip[descriptor] == 14:
        os.system('python ' + path + 'ExtractionTechniques-screen.py')
    elif descrip[descriptor] == 15:
        os.system('python ' + path + 'CodingClass-screen.py')
    elif descrip[descriptor] == 16:
        os.system('python ' + path + 'FickettScore-screen.py')
    elif descrip[descriptor] == 17:
        os.system('python ' + path + 'PseKNC-screen.py')
    elif descrip[descriptor] == 18:
        os.system('python ' + path + 'ExtractionTechniques-Protein-screen.py')
    elif descrip[descriptor] == 19:
        os.system('python ' + path + 'ExtractionTechniques-Protein-screen.py')
    elif descrip[descriptor] == 20:
        os.system('python ' + path + 'ExtractionTechniques-Protein-screen.py')
    elif descrip[descriptor] == 21:
        os.system('python ' + path + 'ExtractionTechniques-Protein-screen.py')
    elif descrip[descriptor] == 22:
        os.system('python ' + path + 'AAindexClass-screen.py')
######################################################
######################################################
