import gzip


def get_part_length(parts):
    return len(parts[1].strip())


def gc_content(line):
    numerator = line.count('G') + line.count('C')
    denominator = line.count('G') + line.count('C') + line.count('A') + line.count('T') + line.count('N')
    return (numerator / denominator) * 100


def average_gc_content(parts):
    res = 0
    for part in parts:
        res += gc_content(part[1])
    return round(res / len(parts), 2)


def clean_repeats(parts):
    seen = set()
    nb_of_repeats = 0
    cleaned_parts = []
    nb_of_reads_with_n = 0
    average_percentage_of_n = 0
    for part in parts:
        if part[1].strip() not in seen:
            cleaned_parts.append(part)
            seen.add(part[1].strip())
            if 'N' in set(part[1].strip()):
                nb_of_reads_with_n += 1
                nb_of_n = part[1].strip().count('N')
                average_percentage_of_n += (nb_of_n / len(part[1].strip())) * 100
        else:
            nb_of_repeats += 1
    return cleaned_parts, nb_of_repeats, nb_of_reads_with_n, round(average_percentage_of_n / len(parts), 2)


def get_average_length(parts):
    lengths = sorted([get_part_length(part) for part in parts])
    lengths_set = sorted(list(set(lengths)))
    lengths_distribution = dict.fromkeys(lengths_set, 0)
    for length in lengths:
        lengths_distribution[length] += 1
    average_length = sum(lengths) / len(lengths)
    return average_length


def get_reads_from_file(filename):
    lines = []
    with open(filename, 'r') as file:
        for line in file:
            lines.append(line)
    return lines


def get_reads_from_tar_gz(filename):
    lines = []
    with gzip.open(filename, 'rb') as file:
        for line in file:
            lines.append(line.decode(encoding='utf8'))
    return lines


def get_data_quality_indicators(filename):
    lines = get_reads_from_tar_gz(filename)
    old_parts = [lines[i: i + 3] for i in range(0, len(lines), 4)]
    cleaned_parts, nb_of_repeats, nb_of_reads_with_n, average_n_per_read_sequence = clean_repeats(old_parts)
    average_length = get_average_length(cleaned_parts)
    return {
        'nb_of_reads': len(old_parts),
        'avg_length': average_length,
        'nb_of_repeats': nb_of_repeats,
        'nb_of_reads_with_n': nb_of_reads_with_n,
        'avg_gc_content': average_gc_content(old_parts),
        'avg_n_per_read_sequence': average_n_per_read_sequence
    }


def print_indicators(indicators):
    print(f'Reads in the file = {indicators["nb_of_reads"]}:')
    print(f'Reads sequence average length = {indicators["avg_length"]}')
    print()
    print(f'Repeats = {indicators["nb_of_repeats"]}')
    print(f'Reads with Ns = {indicators["nb_of_reads_with_n"]}')
    print()
    print(f'GC content average = {indicators["avg_gc_content"]}%')
    print(f"Ns per read sequence = {indicators['avg_n_per_read_sequence']}%")


def get_criteria_value(filename):
    indicators = get_data_quality_indicators(filename)
    return indicators["nb_of_repeats"], indicators["nb_of_reads_with_n"]


filenames = []
for _ in range(3):
    filename = input()
    filenames.append(filename)
filenames.sort(key=lambda filename: get_criteria_value(filename))
indicators = get_data_quality_indicators(filename)
print_indicators(get_data_quality_indicators(filenames[0]))
