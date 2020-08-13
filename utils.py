
def rr_to_cm(input):
    assert input[:-3] == '.rr'
    output = input[:-2] + "cm"

    with open(output, "w") as f:
        for line in open(input):
            if not line[0].isnumeric():
                pass
            else:
                a, b, c, d, e = line.split()
                f.write(f'{a:s} {b:s} {e:s}\n')
