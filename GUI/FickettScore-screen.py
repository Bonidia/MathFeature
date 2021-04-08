import warnings
import os
import time
from Bio import SeqIO
from gooey import Gooey, GooeyParser
warnings.filterwarnings('ignore')
path = os.path.dirname(__file__)
if path != '':
	path = path + '/'


# Fickett TESTCODE data
# NAR 10(17) 5303-531

"""
fickett_value_orf: modified from source code of lncScore
downloaded from https://github.com/WGLab/lncScore

fickett_value_full_sequence: modified from source code of CPC2
downloaded from https://academic.oup.com/nar/article/45/W1/W12/3831091
"""

def look_up_position_prob(value, base, position_para, position_prob, position_weight):

	"""look up positional probability by base and value"""

	if float(value) < 0:
		return None
	for idx, val in enumerate(position_para):
		if float(value) >= val:
			return float(position_prob[base][idx]) * float(position_weight[base])


def look_up_content_prob(value, base, content_para, content_prob, content_weight):

	"""look up content probability by base and value"""

	if float(value) < 0:
		return None
	for idx, val in enumerate(content_para):
		if float(value) >= val:
			return float(content_prob[base][idx]) * float(content_weight[base])


def fickett_value_orf(seq, type_seq):

	"""calculate Fickett value. Input is DNA sequence"""

	position_prob = {
		'A': [0.94, 0.68, 0.84, 0.93, 0.58, 0.68, 0.45, 0.34, 0.20, 0.22],
		'C': [0.80, 0.70, 0.70, 0.81, 0.66, 0.48, 0.51, 0.33, 0.30, 0.23],
		'G': [0.90, 0.88, 0.74, 0.64, 0.53, 0.48, 0.27, 0.16, 0.08, 0.08],
		'T': [0.97, 0.97, 0.91, 0.68, 0.69, 0.44, 0.54, 0.20, 0.09, 0.09]}

	position_weight = {'A': 0.26, 'C': 0.18, 'G': 0.31, 'T': 0.33}
	position_para = [1.9, 1.8, 1.7, 1.6, 1.5, 1.4, 1.3, 1.2, 1.1, 0.0]

	content_prob = {
		'A': [0.28, 0.49, 0.44, 0.55, 0.62, 0.49, 0.67, 0.65, 0.81, 0.21],
		'C': [0.82, 0.64, 0.51, 0.64, 0.59, 0.59, 0.43, 0.44, 0.39, 0.31],
		'G': [0.40, 0.54, 0.47, 0.64, 0.64, 0.73, 0.41, 0.41, 0.33, 0.29],
		'T': [0.28, 0.24, 0.39, 0.40, 0.55, 0.75, 0.56, 0.69, 0.51, 0.58]}

	content_weight = {'A': 0.11, 'C': 0.12, 'G': 0.15, 'T': 0.14}
	content_para = [0.33, 0.31, 0.29, 0.27, 0.25, 0.23, 0.21, 0.17, 0]

	if len(seq) < 2:
		return 0
	fickett_score = 0
	seq = seq.upper()
	total_base = len(seq)
	A_content = float(seq.count('A')) / total_base
	C_content = float(seq.count('C')) / total_base
	G_content = float(seq.count('G')) / total_base
	if type_seq == 1:
		T_content = float(seq.count('T')) / total_base
	else:
		T_content = float(seq.count('U')) / total_base

	phase_0 = [seq[i] for i in range(0, len(seq)) if i % 3 == 0]
	phase_1 = [seq[i] for i in range(0, len(seq)) if i % 3 == 1]
	phase_2 = [seq[i] for i in range(0, len(seq)) if i % 3 == 2]
	
	A_position = max(phase_0.count('A'), phase_1.count('A'), phase_2.count('A')) / (min(phase_0.count('A'), phase_1.count('A'), phase_2.count('A')) + 1.0)

	C_position = max(phase_0.count('C'), phase_1.count('C'), phase_2.count('C')) / (min(phase_0.count('C'), phase_1.count('C'), phase_2.count('C')) + 1.0)

	G_position = max(phase_0.count('G'), phase_1.count('G'), phase_2.count('G')) / (min(phase_0.count('G'), phase_1.count('G'), phase_2.count('G')) + 1.0)

	if type_seq == 1:
		T_position = max(phase_0.count('T'), phase_1.count('T'), phase_2.count('T')) / (min(phase_0.count('T'), phase_1.count('T'), phase_2.count('T')) + 1.0)
	else:
		T_position = max(phase_0.count('U'), phase_1.count('U'), phase_2.count('U')) / (min(phase_0.count('U'), phase_1.count('U'), phase_2.count('U')) + 1.0)

	fickett_score += look_up_content_prob(A_content, 'A', content_para, content_prob, content_weight)
	fickett_score += look_up_content_prob(C_content, 'C', content_para, content_prob, content_weight)
	fickett_score += look_up_content_prob(G_content, 'G', content_para, content_prob, content_weight)
	fickett_score += look_up_content_prob(T_content, 'T', content_para, content_prob, content_weight)
	
	fickett_score += look_up_position_prob(A_position, 'A', position_para, position_prob, position_weight)
	fickett_score += look_up_position_prob(C_position, 'C', position_para, position_prob, position_weight)
	fickett_score += look_up_position_prob(G_position, 'G', position_para, position_prob, position_weight)
	fickett_score += look_up_position_prob(T_position, 'T', position_para, position_prob, position_weight)
	
	return fickett_score


def fickett_value_full_sequence(seq, type_seq):

	"""calculate Fickett from full sequence - CPC2"""

	position_para = [1.9, 1.8, 1.7, 1.6, 1.5, 1.4, 1.3, 1.2, 1.1, 0.0]
	content_para = [0.33, 0.31, 0.29, 0.27, 0.25, 0.23, 0.21, 0.19, 0.17, 0]

	position_prob = {
		'A': [0.51, 0.55, 0.57, 0.52, 0.48, 0.58, 0.57, 0.54, 0.50, 0.36],
		'C': [0.29, 0.44, 0.55, 0.49, 0.52, 0.60, 0.60, 0.56, 0.51, 0.38],
		'G': [0.62, 0.67, 0.74, 0.65, 0.61, 0.62, 0.52, 0.41, 0.31, 0.17],
		'T': [0.51, 0.60, 0.69, 0.64, 0.62, 0.67, 0.58, 0.48, 0.39, 0.24]}
	
	position_weight = {'A': 0.062, 'C': 0.093, 'G': 0.205, 'T': 0.154}
	content_weight = {'A': 0.084, 'C': 0.076, 'G': 0.081, 'T': 0.055}

	content_prob = {
		'A': [0.40, 0.55, 0.58, 0.58, 0.52, 0.48, 0.45, 0.45, 0.38, 0.19],
		'C': [0.50, 0.63, 0.59, 0.50, 0.46, 0.45, 0.47, 0.56, 0.59, 0.33],
		'G': [0.21, 0.40, 0.47, 0.50, 0.52, 0.56, 0.57, 0.52, 0.44, 0.23],
		'T': [0.30, 0.49, 0.56, 0.53, 0.48, 0.48, 0.52, 0.57, 0.60, 0.51]}

	if len(seq) < 2:
		return 0

	fickett_score = 0
	seq = seq.upper()
	total_base = len(seq)

	phase_0 = seq[::3]
	phase_1 = seq[1::3]
	phase_2 = seq[2::3]

	phase_0_A = phase_0.count('A')
	phase_1_A = phase_1.count('A')
	phase_2_A = phase_2.count('A')
	phase_0_C = phase_0.count('C')
	phase_1_C = phase_1.count('C')
	phase_2_C = phase_2.count('C')
	phase_0_G = phase_0.count('G')
	phase_1_G = phase_1.count('G')
	phase_2_G = phase_2.count('G')
	if type_seq == 1:
		phase_0_T = phase_0.count('T')
		phase_1_T = phase_1.count('T')
		phase_2_T = phase_2.count('T')
	else:
		phase_0_T = phase_0.count('U')
		phase_1_T = phase_1.count('U')
		phase_2_T = phase_2.count('U')

	A_content = float(phase_0_A + phase_1_A + phase_2_A) / total_base
	C_content = float(phase_0_C + phase_1_C + phase_2_C) / total_base
	G_content = float(phase_0_G + phase_1_G + phase_2_G) / total_base
	T_content = float(phase_0_T + phase_1_T + phase_2_T) / total_base
	A_position = max([phase_0_A, phase_1_A, phase_2_A]) / (min([phase_0_A, phase_1_A, phase_2_A]) + 1.0)
	C_position = max([phase_0_C, phase_1_C, phase_2_C]) / (min([phase_0_C, phase_1_C, phase_2_C]) + 1.0)
	G_position = max([phase_0_G, phase_1_G, phase_2_G]) / (min([phase_0_G, phase_1_G, phase_2_G]) + 1.0)
	T_position = max([phase_0_T, phase_1_T, phase_2_T]) / (min([phase_0_T, phase_1_T, phase_2_T]) + 1.0)

	fickett_score += look_up_content_prob(A_content, 'A', content_para, content_prob, content_weight)
	fickett_score += look_up_content_prob(C_content, 'C', content_para, content_prob, content_weight)
	fickett_score += look_up_content_prob(G_content, 'G', content_para, content_prob, content_weight)
	fickett_score += look_up_content_prob(T_content, 'T', content_para, content_prob, content_weight)

	fickett_score += look_up_position_prob(A_position, 'A', position_para, position_prob, position_weight)
	fickett_score += look_up_position_prob(C_position, 'C', position_para, position_prob, position_weight)
	fickett_score += look_up_position_prob(G_position, 'G', position_para, position_prob, position_weight)
	fickett_score += look_up_position_prob(T_position, 'T', position_para, position_prob, position_weight)

	return fickett_score


def header_fickett_score(foutput):
	file = open(foutput, 'a')
	file.write('%s,' % 'nameseq')
	file.write('%s,' % 'fickett_score-ORF')
	file.write('%s,' % 'fickett_score-full-sequence')
	file.write('label')
	file.write('\n')
	return


def calculate_sequences(finput, foutput, label_dataset, type_seq):
	header_fickett_score(foutput)
	for seq_record in SeqIO.parse(finput, 'fasta'):
		seq = seq_record.seq
		seq = seq.upper()
		name_seq = seq_record.name
		measure_orf = fickett_value_orf(seq, type_seq)
		measure_full = fickett_value_full_sequence(seq, type_seq)
		file = open(foutput, 'a')
		file.write('%s,' % name_seq)
		file.write('%s,' % str(measure_orf))
		file.write('%s,' % str(measure_full))
		file.write('%s' % str(label_dataset))
		file.write('\n')
		print('Recorded Sequence: %s' % name_seq)
	return


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
	parser = GooeyParser(description="Feature Extraction Package for Biological Sequences")
	screen = parser.add_argument_group('Fickett Score Descriptor - DNA/RNA')
	screen.add_argument('-i',
						'--input_file',
						metavar='Input File',
						help='Fasta format file, e.g., test.fasta',
						required=True,
						widget='FileChooser')
	screen.add_argument('-o',
						'--output_file',
						metavar='Output File',
						help='CSV format file, e.g., test.csv',
						required=True,
						widget='FileSaver')
	screen.add_argument('-l',
						'--label',
						metavar='Dataset Label',
						help='Dataset Label, e.g., 0, 1, lncRNA, mRNA, sncRNA.',
						required=True,
						action='store')
	screen.add_argument('-s',
						'--seq',
						metavar='Type of sequence',
						help='Type of sequence, e.g., DNA or RNA',
						required=True,
						choices=['DNA', 'RNA'],
						nargs='+')
	return parser.parse_args()


#############################################################################
if __name__ == '__main__':
	args = main()
	finput = str(args.input_file)
	foutput = str(args.output_file)
	label_dataset = str(args.label)
	seq = args.seq
	seq = ''.join(seq)
	start_time = time.time()
	calculate_sequences(finput, foutput, label_dataset, seq)
	print('Computation time %s seconds' % (time.time() - start_time))
#############################################################################
