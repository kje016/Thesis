
def get_realiability_sequence():
    with open('Reliability_Sequence.txt') as f:
        reliability_sequence = f.read()
    reliability_sequence = reliability_sequence.split(",")
    reliability_sequence = [int(x) for x in reliability_sequence]
    return reliability_sequence




def table_to_dict():
    dict = {}
    with open('BG.txt') as f:
        file = f.read().split('\n')
    for i, elem in enumerate(file):
        row = elem.split(' ')
        dict.update({int(row[0]): [int(a) for a in row[1:]]})
    print(f'{dict},')


table_to_dict()