
B = ['../3ANALIZA/TRANSPOSONS3.gff',
     'COMPARISON_1vs2_A.gff', 'COMPARISON_1vs2vs3_B1.gff',
     'COMPARISON_1vs2_B.gff', 'COMPARISON_1vs2vs3_B2.gff']

C = ['../2ANALIZA/TRANSPOSONS2.gff',
     'COMPARISON_1vs3_A.gff', 'COMPARISON_1vs2vs3_C1.gff',
     'COMPARISON_1vs3_B.gff', 'COMPARISON_1vs2vs3_C2.gff']

D = ['../1ANALIZA/TRANSPOSONS1.gff',
     'COMPARISON_2vs3_A.gff', 'COMPARISON_1vs2vs3_D1.gff',
     'COMPARISON_2vs3_B.gff', 'COMPARISON_1vs2vs3_D2.gff']

files_list = [B, C, D]


class Transposon:
    def __init__(self, chromosome, start, end, strand, line):
        self.chromosome = chromosome
        self.start = int(start)
        self.end = int(end)
        self.strand = strand
        self.line = line


def create_dictionary_with_TRANSPOSONS(TR_file):
    TRANSPOSONS = {}
    with open(TR_file, 'r') as TR_gff_file:
        list_for_chromosome = []
        for line in TR_gff_file:
            features = line.split()
            chr = features[0]
            start = features[3]
            end = features[4]
            strand = features[6]
            if chr not in TRANSPOSONS:
                list_for_chromosome = []
                TRANSPOSONS[chr] = list_for_chromosome
            transposon = Transposon(chr, start, end, strand, line)
            list_for_chromosome.append(transposon)
    return TRANSPOSONS


def if_intersect(s0, e0, s1, e1):
    if (s0 >= s1 and s0 <= e1) or (s1 >= s0 and s1 <= e0):
        return True
    else:
        return False


def find_intersections(TRANSPOSONS, file_in, file_out):
    with open(file_in, 'r') as input_file:
        with open(file_out, 'w') as output_file:
            for line in input_file:
                no_intersection = True
                features = line.split()
                chromosome = features[0]
                strand = features[6]
                start1 = int(features[3])
                end1 = int(features[4])
                start2 = int(features[11])
                end2 = int(features[12])
                for transposon in TRANSPOSONS[chromosome]:
                    if transposon.strand == strand:
                        start3 = transposon.start
                        end3 = transposon.end
                        if (if_intersect(start1, end1, start3, end3) or if_intersect(start2, end2, start3, end3)):
                            no_intersection = False
                            break
                if no_intersection:
                    output_file.write(line)


if __name__ == '__main__':
    for files in files_list:
        TRANSPOSONS = create_dictionary_with_TRANSPOSONS(files[0])
        find_intersections(TRANSPOSONS, files[1], files[2])
        find_intersections(TRANSPOSONS, files[3], files[4])