import os

TRANSPOSONS3 = {}


class Transposon3:
    def __init__(self, chromosome, name, start, end, strand, line):
        self.chromosome = chromosome
        self.name = name
        self.start = int(start)
        self.end = int(end)
        self.strand = strand
        self.line = line


def create_dictionary_with_TRANSPOSONS3():
    with open('../3ANALIZA/TRANSPOSONS3.gff', 'r') as TR3_gff_file:
        list_for_chromosome = []
        for line in TR3_gff_file:
            features = line.split()
            chr = features[0]
            name = features[2]
            start = features[3]
            end = features[4]
            strand = features[6]
            if chr not in TRANSPOSONS3:
                list_for_chromosome = []
                TRANSPOSONS3[chr] = list_for_chromosome
            transposon = Transposon3(chr, name, start, end, strand, line)
            list_for_chromosome.append(transposon)


def if_intersect(s0, e0, s1, e1):
    if (s0 >= s1 and s0 <= e1) or (s1 >= s0 and s1 <= e0):
        return True
    else:
        return False


def if_more_than_50percent(s0, e0, s1, e1):
    start_overlap = max(s0, s1)
    end_overlap = min(e0, e1)
    overlap = float(end_overlap - start_overlap + 1)
    length0 = float(e0 - s0 + 1)
    length1 = float(e1 - s1 + 1)
    if (((overlap/length0)*100) >= 50) or (((overlap/length1)*100) >= 50):
        return True
    else:
        return False


def find_intersections():
    with open('COMPARISON_1vs2_A.gff', 'r') as COMPARISON_1vs2_A_file:
        with open('COMPARISON_1vs2vs3_A1.gff', 'w') as output_A1_file:
            with open('COMPARISON_1vs2vs3_A2_unsorted.gff', 'w') as output_A2_file:
                for line_COM in COMPARISON_1vs2_A_file.readlines():
                    fields_COM = line_COM.split()
                    chr_COM = fields_COM[0]
                    strand_COM = fields_COM[6]
                    start1 = int(fields_COM[3])
                    end1 = int(fields_COM[4])
                    start2 = int(fields_COM[11])
                    end2 = int(fields_COM[12])
                    for transposon in TRANSPOSONS3[chr_COM]:
                        if transposon.strand == strand_COM:
                            start3 = transposon.start
                            end3 = transposon.end
                            if (if_intersect(start1, end1, start3, end3) and if_intersect(start2, end2, start3, end3)):
                                if (if_more_than_50percent(start1, end1, start3, end3) and if_more_than_50percent(start2, end2, start3, end3)):
                                    output_A1_file.write(line_COM[:-1] + '\t' + transposon.line)
                                else:
                                    output_A2_file.write(line_COM[:-1] + '\t' + transposon.line)
    with open('COMPARISON_1vs2_B.gff', 'r') as COMPARISON_1vs2_B_file:
        with open('COMPARISON_1vs2vs3_A2_unsorted.gff', 'a') as output_A2_file:
            for line_COM in COMPARISON_1vs2_B_file.readlines():
                fields_COM = line_COM.split()
                chr_COM = fields_COM[0]
                strand_COM = fields_COM[6]
                start1 = int(fields_COM[3])
                end1 = int(fields_COM[4])
                start2 = int(fields_COM[11])
                end2 = int(fields_COM[12])
                for transposon in TRANSPOSONS3[chr_COM]:
                    if transposon.strand == strand_COM:
                        start3 = transposon.start
                        end3 = transposon.end
                        if (if_intersect(start1, end1, start3, end3) and if_intersect(start2, end2, start3, end3)):
                            output_A2_file.write(line_COM[:-1] + '\t' + transposon.line)

if __name__ == '__main__':
    create_dictionary_with_TRANSPOSONS3()
    find_intersections()