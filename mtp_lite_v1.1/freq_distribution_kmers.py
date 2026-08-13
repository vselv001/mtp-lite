import sys

import numpy as np
import xlsxwriter


def write_histo_to_xcel_get_lower_upper(input_path, output_path):
    workbook = xlsxwriter.Workbook(output_path)
    worksheet = workbook.add_worksheet()

    freqs, no_of_kmers = [], []

    with open(input_path, "r") as histo:
        for row, line in enumerate(histo):
            args = line.strip().split()
            worksheet.write(row, 0, args[0])
            worksheet.write(row, 1, args[1])
            if args[0] == "0":
                continue
            freqs.append(int(args[0]))
            no_of_kmers.append(int(args[1]))

    workbook.close()

    start_index = 0
    while no_of_kmers[start_index] / no_of_kmers[start_index + 1] > 1.25:
        start_index += 1
    truncated_no_of_kmers = no_of_kmers[start_index:]
    max_index = no_of_kmers.index(max(truncated_no_of_kmers))
    end_index = start_index + 2 * (max_index - start_index)

    frequencies = np.array(no_of_kmers[start_index:end_index + 1])
    values = np.array(freqs[start_index:end_index + 1])

    mean = np.sum(values * frequencies) / np.sum(frequencies)
    variance = np.sum(frequencies * ((values - mean) ** 2)) / np.sum(frequencies)
    stdev = np.sqrt(variance)

    print(f"mean: {mean} std: {stdev}")

    lower = max(int(np.floor(mean - 3 * stdev)), 5)
    upper = int(np.ceil(mean + 3 * stdev))

    print(f"lower: {lower}")
    print(f"upper: {upper}")


if __name__ == "__main__":
    write_histo_to_xcel_get_lower_upper(sys.argv[1], sys.argv[2])
