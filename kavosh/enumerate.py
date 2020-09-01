import sys



def GenerateCompositions(N):
    # We generate all compositions using the following algorithm:
    # Write 1 ? 1 ? 1 ? ... ? 1 where there are N - 1 question marks.
    # Each question mark is either a plus sign, or a comma
    # Summing all 1s inbetween commas produces a composition
    # Iterating over all combinations of commas and plus signs gives the
    # complete set of unique combinations. We iterate over all 2 ** (N - 1)
    # numbers, write those numbers in binary, with 1s representing commas
    # and 0s representing plus signs.

    # there are 2^(N - 1) total compositions
    total_compositions = 2 ** (N - 1)

    for comp in range(total_compositions):
        # start with empty components in the composition
        components = []

        # the current value is one
        current_value = 1
        for i in range(N - 1):
            # want to start with low order bits first
            bit = (comp >> N - 2 - i) & 1

            # if the bit is one, treat as a comma and reset the current value
            if bit:
                components.append(current_value)
                current_value = 1
            # if the bit is zero, treat as a plus sign and add the two values together
            else:
                current_value += 1

        # append the remaining current value which has no comma at its end
        components.append(current_value)



def EnumerateSubgraphs(graph):
    pass
