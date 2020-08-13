CM_FORMAT = """#CMVIEW GRAPH FILE ver: 1.0
#SEQUENCE: {}
#PDB CHAIN CODE: A
#CHAIN: A
#MODEL: {}
#CT: Cb
#CUTOFF: 8.0
"""


def rr_to_cm(input):
    assert input[-3:] == '.rr'
    output = input[:-2] + "cm"

    with open(output, "w") as op:
        with open(input) as ip:
            for _ in range(4):
                ip.readline()
            model = ip.readline()[6:-1]
            seq = ip.readline()[:-1]
            op.write(CM_FORMAT.format(seq, model))
            for line in ip.readlines():
                if line[0].isnumeric():
                    a, b, c, d, e = line.split()
                    op.write(f'{a:s} {b:s} {e:s}\n')
