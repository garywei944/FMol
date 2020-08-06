def clean_data(file):
    """
    Generate a new file with only SMILES without dot
    Keep in a delete file the SMILES not accepted

    :param file: file containing SMILES, one per line.
    :return: None
    """
    with open(file, 'r') as em:
        data = em.readlines()
        for d in data:
            try:
                if '.' in d:
                    raise Exception(". in the SMILES ")
                with open(file + "_clean", "a") as clean:
                    clean.write(d)
            except Exception as e:
                with open(file + "_deleted", "a") as bad_out:
                    bad_out.write(repr(e))
                    bad_out.write(d)


if __name__ == '__main__':
    clean_data(file="emolecules")
