from datetime import datetime


class Monitor(object):

    def __init__(self):
        """
        Initialize the monitor to identify the task progress.

        Example
            >>> from dsw import Monitor
            >>> monitor = Monitor()
            >>> monitor.output(current_state=1, total_state=10)
            \r|███                 | 10% ( 1/10) wait 0000:00:00.
            >>> monitor.output(current_state=5, total_state=10)
            \r|███████████         | 50% ( 5/10) wait 0000:00:00.
            >>> monitor.output(current_state=10, total_state=10)
            \r|████████████████████|100% (10/10) used 0000:00:00.
        """
        self.last_time = None

    def output(self, current_state, total_state, extra=None):
        """
        Output the current state of process.

        :param current_state: current state of process.
        :type current_state: int

        :param total_state: total state of process.
        :type total_state: int

        :param extra: extra vision information if required.
        :type extra: dict
        """
        if self.last_time is None:
            self.last_time = datetime.now()

        if current_state == 0:
            return

        position = int(current_state / total_state * 100)

        string = "|"

        for index in range(0, 100, 5):
            if position >= index:
                string += "█"
            else:
                string += " "

        string += "|"

        pass_time = (datetime.now() - self.last_time).total_seconds()
        wait_time = int(pass_time * (total_state - current_state) / current_state)

        string += " " * (3 - len(str(position))) + str(position) + "% ("

        string += " " * (len(str(total_state)) - len(str(current_state))) + str(current_state) + "/" + str(total_state)

        if current_state < total_state:
            minute, second = divmod(wait_time, 60)
            hour, minute = divmod(minute, 60)
            string += ") wait " + "%04d:%02d:%02d" % (hour, minute, second)
        else:
            minute, second = divmod(pass_time, 60)
            hour, minute = divmod(minute, 60)
            string += ") used " + "%04d:%02d:%02d" % (hour, minute, second)

        if extra is not None:
            string += " " + str(extra).replace("\'", "").replace("{", "(").replace("}", ")") + "."
        else:
            string += "."

        print("\r" + string, end="", flush=True)

        if current_state >= total_state:
            self.last_time = None
            print()


def calculus_addition(number, base):
    """
    Do huge number addition calculus with a small base value, as number + base.

    :param number: huge number.
    :type number: str

    :param base: base value:
    :type base: str

    :return: number + base.
    :rtype: str

    .. note::
        The integer of parameter "base" must less be than 10 in decimal system.

    Example
        >>> from dsw import calculus_addition
        >>> calculus_addition(number="99999999999999999999999999999999999999999999999999", base="2")
        '100000000000000000000000000000000000000000000000001'
    """
    number, base = list(number), list(base.zfill(len(number)))

    result = [0 for _ in range(len(number) + 1)]
    for index in range(len(number) - 1, -1, -1):
        sum_value = int(number[index]) + int(base[index]) + int(result[index + 1])
        if sum_value < 10:
            result[index + 1] = sum_value
        else:
            flag = 0
            while sum_value > 0:
                result[index + 1 - flag] = sum_value % 10
                sum_value //= 10
                flag += 1

    result = "".join(list(map(str, result)))

    return result if result[0] != "0" else result[1:]


def calculus_subtraction(number, base):
    """
    Do huge number subtraction calculus with a small base value, as number - base.

    :param number: huge number.
    :type number: str

    :param base: base value:
    :type base: str

    :return: number - base.
    :rtype: str

    .. note::
        The integer of parameter "base" must less be than 10 in decimal system.

    Example
        >>> from dsw import calculus_subtraction
        >>> calculus_subtraction(number="10000000000000000000000000000000000000000000000001", base="2")
        '9999999999999999999999999999999999999999999999999'
    """
    number, base = [int(item) for item in number], [int(item) for item in base]

    residue = ""
    for index in range(len(base)):
        flag_a, flag_b = len(number) - 1 - index, len(base) - 1 - index
        if int(number[flag_a]) >= int(base[flag_b]):
            residue = str(int(number[flag_a]) - int(base[flag_b])) + residue
        else:
            residue = str(10 + int(number[flag_a]) - int(base[flag_b])) + residue
            while number[flag_a - 1] == 0:
                number[flag_a - 1] = 9
                flag_a -= 1
            number[flag_a - 1] -= 1

        for flag in range(len(number) - 1 - index - 1, -1, -1):
            residue = str(number[flag]) + residue

    for index in range(len(residue)):
        if residue[index] != "0":
            return residue[index:]

    return "0"


def calculus_multiplication(number, base):
    """
    Do huge number multiplication calculus with a small base value, as number * base.

    :param number: huge number.
    :type number: str

    :param base: base value:
    :type base: str

    :return: number * base.
    :rtype: str

    .. note::
        The integer of parameter "base" must less be than 10 in decimal system.

    Example
        >>> from dsw import calculus_multiplication
        >>> calculus_multiplication(number="9999999999999999999999999999999999999999999999999", base="2")
        '19999999999999999999999999999999999999999999999998'
    """
    if base == "0":
        return "0"

    if base == "1":
        return number

    number = [int(item) for item in number]

    remainder = 0
    for index in range(len(number))[::-1]:
        current = number[index] * int(base) + remainder
        if current >= 10:
            number[index] = current % 10
            remainder = current // 10
        else:
            number[index] = current
            remainder = 0

    while remainder > 0:
        number.insert(0, remainder % 10)
        remainder //= 10

    result = "".join(list(map(str, number)))

    return result


def calculus_division(number, base):
    """
    Do huge number division calculus with a small base value, as number / base and number % base.

    :param number: huge number.
    :type number: str

    :param base: base value:
    :type base: str

    :return: number // base and number % base.
    :rtype: (str, str)

    Example
        >>> from dsw import calculus_division
        >>> calculus_division(number="9999999999999999999999999999999999999999999999999", base="2")
        ('4999999999999999999999999999999999999999999999999', '1')

    .. note::
        The integer of parameter "base" must less be than 10 in decimal system.

        When the divisor (base) is "0", the program will return "0", "0".
    """
    if base == "0":
        return "0", "0"

    if base == "1":
        return number, "0"

    if len(number) == 1 and number[0] < base:
        return "0", number[0]

    number, new_number, remainder = [int(item) for item in number], [], 0
    for index, quotient in enumerate(number):
        current = quotient + remainder * 10
        if current >= int(base):
            new_number.append(current // int(base))
            remainder = current - new_number[-1] * int(base)
        else:
            new_number.append(0)
            remainder = current

    quotient = "".join(list(map(str, new_number)))

    for index in range(len(quotient)):
        if quotient[index] != "0":
            return quotient[index:], str(remainder)

    return "0", str(remainder)


def bit_to_number(bit_array, is_string=True):
    """
    Transform a bit array to the equivalent decimal number.

    :param bit_array: bit array.
    :type bit_array: list

    :param is_string: type of equivalent decimal number is str.
    :type: is_string: bool

    :return: equivalent decimal number (may huge) of the inputted bit array.
    :rtype: str or int

    Example
        >>> from dsw import bit_to_number
        >>> bit_to_number(bit_array=[1, 1, 1, 1, 1, 0, 0, 1, 1, 1], is_string=True)
        '999'
        >>> bit_to_number(bit_array=[1, 1, 1, 1, 1, 0, 0, 1, 1, 1], is_string=False)
        999
    """
    if is_string:
        decimal_number = "0"

        for a_bit in bit_array:
            # multiply by 2
            decimal_number = calculus_multiplication(number=decimal_number, base="2")
            # add current bit
            decimal_number = calculus_addition(number=decimal_number, base=str(a_bit))
    else:
        decimal_number = 0

        for a_bit in bit_array:
            decimal_number = decimal_number * 2 + a_bit

    return decimal_number


def number_to_bit(decimal_number, bit_length):
    """
    Transform a decimal number to the equivalent bit array with specific length.

    :param decimal_number: decimal number (may huge) of the bit array.
    :type decimal_number: str or int

    :param bit_length: default length of the bit array.
    :type bit_length: int

    :return: bit array.
    :rtype: list

    Example
        >>> from dsw import number_to_bit
        >>> number_to_bit(decimal_number="999", bit_length=10)
        [1, 1, 1, 1, 1, 0, 0, 1, 1, 1]
        >>> number_to_bit(decimal_number=999, bit_length=10)
        [1, 1, 1, 1, 1, 0, 0, 1, 1, 1]
    """
    one_array = []
    if type(decimal_number) == str:
        while decimal_number != "0":  # decimal number > 0
            decimal_number, remainder = calculus_division(number=decimal_number, base="2")
            one_array.insert(0, int(remainder))
    elif type(decimal_number) == int:
        while decimal_number > 0:
            decimal_number, remainder = divmod(decimal_number, 2)
            one_array.insert(0, remainder)
    else:
        raise ValueError("No such type of decimal number (" + str(type(decimal_number)) + ")!")

    return [0] * (bit_length - len(one_array)) + one_array


def dna_to_number(dna_string, nucleotides=None, is_string=True):
    """
    Transform a DNA string to the equivalent decimal number.

    :param dna_string: required DNA string.
    :type dna_string: str

    :param nucleotides: usage of nucleotides.
    :type nucleotides: list

    :param is_string: type of equivalent decimal number is str.
    :type: is_string: bool

    :return: equivalent decimal number (may huge) of the inputted DNA string.
    :rtype: str or int

    Example
        >>> from dsw import dna_to_number
        >>> dna_to_number(dna_string="ACGTACGT", is_string=True)
        '6939'
        >>> dna_to_number(dna_string="ACGTACGT", is_string=False)
        6939
    """
    if nucleotides is None:
        nucleotides = ["A", "C", "G", "T"]

    nucleotide_values = list(map(nucleotides.index, dna_string))

    if is_string:
        decimal_number = "0"
        for nucleotide_value in nucleotide_values:
            # multiply by length of usage of nucleotides.
            decimal_number = calculus_multiplication(number=decimal_number, base=str(len(nucleotides)))
            # add current nucleotide value.
            decimal_number = calculus_addition(number=decimal_number, base=str(nucleotide_value))
    else:
        decimal_number = 0
        for nucleotide_value in nucleotide_values:
            decimal_number = decimal_number * 4 + nucleotide_value

    return decimal_number


def number_to_dna(decimal_number, dna_length, nucleotides=None):
    """
    Transform a decimal number to the equivalent DNA string with specific length.

    :param decimal_number: decimal number (may huge) of the DNA string.
    :type decimal_number: str or int

    :param dna_length: default length of the DNA string.
    :type dna_length: int

    :param nucleotides: usage of nucleotides.
    :type nucleotides: list

    :return: equivalent DNA string of the decimal number.
    :rtype: str

    Example
        >>> from dsw import number_to_dna
        >>> number_to_dna(decimal_number=6939, dna_length=8)
        'ACGTACGT'
        >>> number_to_dna(decimal_number="6939", dna_length=8)
        'ACGTACGT'
    """
    if nucleotides is None:
        nucleotides = ["A", "C", "G", "T"]

    one_array = []

    if type(decimal_number) == str:
        while decimal_number != "0":  # decimal number > 0
            decimal_number, remainder = calculus_division(number=decimal_number, base=str(len(nucleotides)))
            one_array.insert(0, nucleotides[int(remainder)])
    elif type(decimal_number) == int:
        while decimal_number > 0:
            decimal_number, remainder = divmod(decimal_number, len(nucleotides))
            one_array.insert(0, nucleotides[remainder])
    else:
        raise ValueError("No such type of decimal number (" + str(type(decimal_number)) + ")!")

    one_array = "".join(one_array)

    return nucleotides[0] * (dna_length - len(one_array)) + one_array
