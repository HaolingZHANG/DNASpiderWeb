from collections import Counter
from numpy import load, array, zeros, linspace, median, sum, max, log2, where
# noinspection PyPackageRequirements
from openpyxl import Workbook
from os.path import exists
from scipy.optimize import curve_fit

from dsw import obtain_vertices


def task1():
    if not exists(path="./raw/encode.xlsx"):
        record = load(file="./data/performance_evaluation_1.pkl", allow_pickle=True)
        saved_data = [[], [], [], []]
        for value in record:
            if value[4] > 0:
                saved_data[int(value[0] - 1)].append([int(value[1]), int(value[2]), int(value[3]), value[4]])
            else:
                saved_data[int(value[0] - 1)].append([int(value[1]), int(value[2]), int(value[3]), "N/A"])

        book = Workbook()

        for name, saved_values in zip(["DNA Fountain", "Yin-Yang Code", "HEDGES", "SPIDER-WEB"], saved_data):
            sheet = book.create_sheet(title=name)
            sheet.append(["constraint set index", "parameter index", "data index", "code rate"])
            for saved_value in saved_values:
                sheet.append(saved_value)

        book.save("./raw/encode.xlsx")


def task2():
    if not exists(path="./raw/correct.xlsx"):
        records = load(file="./data/correction_evaluation_1.pkl", allow_pickle=True)

        book = Workbook()

        sheet = book.create_sheet(title="correction rate")
        sheet.append(["error times", "correction flag", "runtime"])
        data_group = records["correction rate"]
        for error_times in [1, 2, 3, 4]:
            for data in data_group[error_times]:
                sheet.append([error_times, data[-2], data[-1]])

        sheet = book.create_sheet(title="minimum reads")
        data_group = records["minimum reads"]
        sheet.append(["reads for recovery", "error = 1", "errors = 2", "errors = 3", "errors = 4"])
        values = [[], [], [], []]
        for error_times in [1, 2, 3, 4]:
            values = array(sorted(Counter(data_group[error_times]).items(), key=lambda kv: (kv[0], kv[1])))
            counts = zeros(shape=100,)
            for location, count in array(values):
                counts[location - 1] = count
            for index in range(max(values[:, 0])):
                values[error_times - 1].append([counts[index]])
        for index in range(100):
            record = [index + 1]
            for error_times in [1, 2, 3, 4]:
                if index < len(values[error_times - 1]):
                    record += [values[error_times - 1][index]]
                else:
                    record += [""]
            sheet.append(record)

        sheet = book.create_sheet(title="frequency-based recovery")
        data_group = records["frequency-based recovery"]
        sheet.append(["order", "reads (error = 1)", "flag1", "reads (errors = 2)", "flag2",
                      "reads (errors = 3)", "flag3", "reads (errors = 4)", "flag4"])

        values, maximum_count = [[], [], [], []], 0
        for error_times in [1, 2, 3, 4]:
            right_set, obtained_set = data_group[error_times]
            right_set = set(right_set)
            counts = sorted(Counter(obtained_set).items(), key=lambda kv: (kv[1], kv[0]), reverse=True)
            maximum_count = max([maximum_count, len(counts)])
            for dna_string, count in counts:
                if dna_string in right_set:
                    values[error_times - 1].append([count, True])
                else:
                    values[error_times - 1].append([count, False])

        for index in range(maximum_count):
            record = [index + 1]
            for error_times in [1, 2, 3, 4]:
                if index < len(values[error_times - 1]):
                    record += [*values[error_times - 1][index]]
                else:
                    record += ["", ""]
            sheet.append(record)

        book.save("./raw/correct.xlsx")


def task3():
    def estimated_equation(x, a):
        return a ** x

    if not exists(path="./raw/encrypt.xlsx"):
        book = Workbook()

        reconstructions = load(file="./data/reconstruction_evaluation.npy")
        capacities, _ = load(file="./data/capacity_reference.pkl", allow_pickle=True)
        accessors = load(file="./data/coding_graphs.pkl", allow_pickle=True)

        location_group, percentages = [], linspace(0, 1, 101)
        for index, data in enumerate(reconstructions):
            medians = [0]
            for iteration_data in data.T[1:]:
                used_data = iteration_data[iteration_data > 0]
                if len(used_data) == 0:
                    break
                medians.append(median(used_data))

            number = len(obtain_vertices(accessors[str(index + 1).zfill(2)]))
            medians = array(medians) / number

            numbers = array(list(range(len(medians))))
            parameter = curve_fit(estimated_equation, xdata=medians[1:], ydata=numbers[1:], p0=(2,))[0][0]
            used_capacities = parameter ** percentages * 100 * capacities[index] / 8
            location_group.append(used_capacities)
        location_group = array(location_group).T

        sheet = book.create_sheet(title="reconstruction")
        sheet.append(["percentage"] + ["constraint " + str(index).zfill(2) for index in range(1, 13)])
        for percentage, values in zip(percentages, location_group):
            sheet.append([percentage, *values])

        rates = []
        for filter_index in ["01", "02", "03", "04", "05", "06", "07", "08", "09", "10", "11", "12"]:
            accessor = accessors[filter_index]
            out_degrees, counts = array(list(Counter(where(accessor >= 0)[0]).values())), zeros(shape=(3,), dtype=float)
            for out_degree in [2, 3, 4]:
                counts[out_degree - 2] = len(where(out_degrees == out_degree)[0])
            counts /= sum(counts)
            rates.append(counts[0] * log2(2.0) + counts[1] * log2(6.0) + counts[2] * log2(24.0))
        rates = array(rates)

        sheet = book.create_sheet(title="key strength")
        sheet.append(["string length"] + ["constraint " + str(index).zfill(2) for index in range(1, 13)])
        for length in linspace(0, 320, 321):
            sheet.append([length] + (length * rates).tolist())

        book.save("./raw/encrypt.xlsx")


if __name__ == "__main__":
    # task1()
    task2()
    # task3()
