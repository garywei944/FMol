import os

command = "python alphafold.py -i ${} -o {} -m {} -r {} -t {} | cat &"


def launch_alphafold(target, target_file, model_dir="model", gpu_more_than_8g_memory=False):
    output_dir = target + "_out"

    print("Saving output to", )

    for replica in range(4):
        mode = ["D", "B"]
        if replica == 0:
            mode.append("T")
        for m in mode:
            print("Launching model:", m, replica)
            if not os.fork():
                os.system(command.format(target_file, output_dir, model_dir, replica, m))
        if not gpu_more_than_8g_memory:
            os.wait()

    print("All models running, waiting for them to complete")
    os.wait()

    print("Ensembling all replica outputs & Pasting contact maps")
    os.system("python alphafold.py -i {} -o {} -e".format(target_file, output_dir))
